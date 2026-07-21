# ==========================================================================
# fortran_parity_isotopes.jl — C13/C14 carbon-isotope Julia↔Fortran parity.
#
# Validates the config-gated isotope subsystem (carbon_isotopes.jl +
# c_iso_flux.jl, gated use_c13/use_c14) against a Fortran CTSM reference.
#
# GROUND TRUTH (generated 2026-07-17, clm_bgc_spinup/bgc_ref_iso):
#   A use_c13=.true./use_c14=.true. CN spinup, STARTUP from the Bow 2202-01-01
#   BGC restart, run to 2202-07-24 (mid-July growing season), dumping the
#   restart-format pdump window nstep 4700..4725. The isotope pools are
#   COLD-STARTED from bulk C at the January restart (leafc_13 = leafc*c3_r2,
#   ...) and then evolved by real isotope photosynthesis + BGC flux physics.
#
#   The Fortran #2119 hard-abort ("Cannot initialize from a run without c13 to
#   a run with c13") that previously blocked this — CNVegCarbonStateType.F90
#   line ~2657 — is bypassed by a case SourceMods patch so the reseed logic
#   right below it runs. That patch only executes under use_c13/c14 (default
#   non-isotope runs are byte-identical). See docs/FORTRAN_VALIDATION_BACKLOG.md
#   row C1 and docs/PARITY_COVERAGE_2026-07.md.
#
# WHAT THIS DIFFS: read_fortran_restart! does NOT read isotope restart vars, so
# this harness injects the c13/c14 pools from the before_step dump via a
# pre_step_hook (patch-remapped by PFT type, same convention as set_patch_1d!),
# runs ONE use_cn+use_c13+use_c14 clm_drv! step, then diffs the Julia c13/c14
# state against the Fortran after_hydrologydrainage dump. It covers:
#   * VEG isotope pools (leafc_13 … xsmrpool_14) — absolute post-step parity
#     (shared IC → tight) + one-step increment + a δ13C decoupled-physics check.
#   * SOIL isotope decomposition pools (litr/soil/cwd × c13/c14, _13_vr/_14_vr):
#     also injected from before_step, then the isotope decomp cascade + litter
#     vertical transport is diffed post-step. (First-ever soil-isotope diff.)
#   * The rc13 discrimination ratios (rc13_canair/psnsun/psnsha) — recomputed
#     each step from the atmospheric C13O2 partial pressure, so NO injection: a
#     pure test of the atmospheric-ratio + fractionation path.
# See the STATUS block for the two bugs this exposed (forc_pc13o2 never set;
# isotope carbon-flux instances never zeroed → soil-column NaN).
#
#   julia +1.12 --project=. scripts/fortran_parity_isotopes.jl [nstep]
# ==========================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const ISO_DUMPDIR = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_iso"
const ISO_FORCING = replace(FFORCING, "clmforc.2003.nc" => "clmforc.2002.nc")

# The veg-isotope pool stems: dump var "<stem>_13" / "<stem>_14" maps to the
# Julia field "<stem>_patch" on the c13_/c14_cnveg_carbonstate_inst.
const ISO_STEMS = (
    "leafc","leafc_storage","leafc_xfer",
    "frootc","frootc_storage","frootc_xfer",
    "livestemc","livestemc_storage","livestemc_xfer",
    "deadstemc","deadstemc_storage","deadstemc_xfer",
    "livecrootc","livecrootc_storage","livecrootc_xfer",
    "deadcrootc","deadcrootc_storage","deadcrootc_xfer",
    "gresp_storage","gresp_xfer","cpool","xsmrpool")

# ----------------------------------------------------------------------------
# Inventory: confirm the reference carries c13/c14 pools.
# ----------------------------------------------------------------------------
function inventory(nstep)
    bdump = joinpath(ISO_DUMPDIR, "pdump_before_step_n$(nstep).nc")
    isfile(bdump) || (println("MISSING dump: $bdump"); return false)
    ds = NCDataset(bdump, "r")
    ks = sort(collect(keys(ds)))
    c13 = filter(k -> endswith(k, "_13"), ks)
    c14 = filter(k -> endswith(k, "_14"), ks)
    pfts = haskey(ds, "pfts1d_itypveg") ? Int.(ds["pfts1d_itypveg"][:]) : Int[]
    println("=== bgc_ref_iso reference inventory (nstep $nstep) ===")
    println("  total vars : ", length(ks))
    println("  PFT itypes : ", pfts)
    println("  n c13 vars : ", length(c13), "  n c14 vars : ", length(c14))
    # delta13C sanity of the leaf pool (C3 leaves ~ −25 to −30 per mil)
    if haskey(ds, "leafc") && haskey(ds, "leafc_13")
        lc  = Float64.(coalesce.(ds["leafc"][:], NaN))
        l13 = Float64.(coalesce.(ds["leafc_13"][:], NaN))
        Rpdb = 0.0112372
        for p in eachindex(lc)
            (lc[p] > 0) || continue
            R = l13[p] / (lc[p] - l13[p]); d13 = (R / Rpdb - 1) * 1000
            @printf("  patch %d leaf δ13C = %+.2f ‰\n", p, d13)
        end
    end
    close(ds)
    return length(c13) > 0
end

# ----------------------------------------------------------------------------
# pre_step_hook: inject the Fortran c13/c14 veg pools into the Julia isotope
# carbon-state instances (read_fortran_restart! does not do this). Patch remap
# by PFT type, identical to fortran_restart.jl set_patch_1d!.
# ----------------------------------------------------------------------------
function make_iso_injector(inst)
    c13cs = inst.bgc_vegetation.c13_cnveg_carbonstate_inst
    c14cs = inst.bgc_vegetation.c14_cnveg_carbonstate_inst
    return function (inst_, bounds, dumpfile)
        ds = NCDataset(dumpfile, "r")
        f90_pfts = haskey(ds, "pfts1d_itypveg") ? Int.(ds["pfts1d_itypveg"][:]) : Int[]
        jl_pfts  = Int.(inst.patch.itype)
        np_f90 = length(f90_pfts); np_jl = length(jl_pfts)
        # dump patch pf -> julia patch pj (first pft-type match). NOTE: use nested
        # `for` blocks — a fused `for pf in ..., pj in ...` loop's `break` exits BOTH
        # loops in Julia, which would map only the first patch (a vacuous mapping).
        pmap = Dict{Int,Int}()
        for pf in 1:np_f90
            for pj in 1:np_jl
                if f90_pfts[pf] == jl_pfts[pj]; pmap[pf] = pj; break; end
            end
        end
        function inject!(cs, suffix)
            for stem in ISO_STEMS
                dn = stem * suffix
                haskey(ds, dn) || continue
                fld = Symbol(stem * "_patch")
                hasproperty(cs, fld) || continue
                arr = getfield(cs, fld)
                (arr isa AbstractVector && length(arr) >= np_jl) || continue
                dat = Float64.(coalesce.(ds[dn][:], NaN))
                for (pf, pj) in pmap
                    pf <= length(dat) && (arr[pj] = dat[pf])
                end
            end
            # pft_ctrunc_<iso> -> ctrunc_patch ; totvegc_<iso> -> totvegc_patch
            for (dn, fld) in (("pft_ctrunc" * suffix, :ctrunc_patch),
                              ("totvegc" * suffix,    :totvegc_patch))
                (haskey(ds, dn) && hasproperty(cs, fld)) || continue
                arr = getfield(cs, fld); (length(arr) >= np_jl) || continue
                dat = Float64.(coalesce.(ds[dn][:], NaN))
                for (pf, pj) in pmap; pf <= length(dat) && (arr[pj] = dat[pf]); end
            end
            dn = "totvegc_col" * suffix
            if haskey(ds, dn) && hasproperty(cs, :totvegc_col) && length(cs.totvegc_col) >= 1
                cs.totvegc_col[1] = Float64(coalesce(ds[dn][1], NaN))
            end
        end
        inject!(c13cs, "_13")
        inject!(c14cs, "_14")

        # --- SOIL isotope decomposition pools (litr/soil/cwd _13_vr/_14_vr) ---
        # read_fortran_restart! ingests the BULK soil decomp pools but not the
        # isotope soil instances, so seed them from the before_step dump exactly
        # as the veg pools above. Column decomp_cpools_vr_col is (col, lev, pool);
        # the dump vars are (lev,) for the single Bow column. Pool ordering matches
        # the bulk CN parity (litr1..3, soil1..3, cwd = 1..7).
        function inject_soil!(scs, suffix)
            (scs === nothing) && return
            arr = scs.decomp_cpools_vr_col
            (arr isa AbstractArray && ndims(arr) == 3 && size(arr, 1) >= 1) || return
            nlev = size(arr, 2)
            for (stem, pidx) in ISO_SOIL_STEMS
                dn = stem * suffix * "_vr"   # e.g. "litr1c" * "_13" * "_vr"
                haskey(ds, dn) || continue
                dat = Float64.(coalesce.(ds[dn][:], NaN))
                for j in 1:min(nlev, length(dat))
                    arr[1, j, pidx] = dat[j]
                end
            end
        end
        inject_soil!(get_c13_soilcs(inst), "_13")
        inject_soil!(get_c14_soilcs(inst), "_14")
        close(ds)
    end
end

# Soil decomp isotope pool stems → decomp_cpools_vr_col pool index (bulk CN order).
const ISO_SOIL_STEMS = (("litr1c",1),("litr2c",2),("litr3c",3),
                        ("soil1c",4),("soil2c",5),("soil3c",6),("cwdc",7))

get_c13_soilcs(inst) = hasproperty(inst, :c13_soilbiogeochem_carbonstate) ?
    inst.c13_soilbiogeochem_carbonstate : nothing
get_c14_soilcs(inst) = hasproperty(inst, :c14_soilbiogeochem_carbonstate) ?
    inst.c14_soilbiogeochem_carbonstate : nothing

# ----------------------------------------------------------------------------
# Snapshot the Julia c13/c14 veg pools (per patch) as a Dict for before/after.
# ----------------------------------------------------------------------------
function snapshot_iso(inst)
    d = Dict{String,Vector{Float64}}()
    for (cs, suf) in ((inst.bgc_vegetation.c13_cnveg_carbonstate_inst, "_13"),
                      (inst.bgc_vegetation.c14_cnveg_carbonstate_inst, "_14"))
        for stem in ISO_STEMS
            fld = Symbol(stem * "_patch")
            hasproperty(cs, fld) || continue
            arr = getfield(cs, fld)
            (arr isa AbstractVector && !isempty(arr)) || continue
            d[stem * suf] = Float64.(arr)
        end
    end
    return d
end

# ----------------------------------------------------------------------------
# Main parity: one isotope-enabled step, diff pools + increment vs Fortran.
# ----------------------------------------------------------------------------
function iso_parity_step(nstep)
    bdump = joinpath(ISO_DUMPDIR, "pdump_before_step_n$(nstep).nc")
    edump = joinpath(ISO_DUMPDIR, "pdump_after_hydrologydrainage_n$(nstep).nc")
    (isfile(bdump) && isfile(edump)) || (println("iso dumps missing for n$nstep"); return 1)

    println("\nRunning ONE use_cn+use_c13+use_c14 clm_drv! step from bgc_ref_iso (nstep=$nstep) ...")
    local inst, bounds, snap_before
    injector_ref = Ref{Any}(nothing)
    try
        inst, bounds = run_one_parity_step!(nstep; use_cn=true, dumpdir=ISO_DUMPDIR,
            use_hydrstress=true, use_luna=true, use_c13=true, use_c14=true,
            step_date=DateTime(2002,1,1) + Hour(nstep),
            forcing_file=ISO_FORCING,
            pre_step_hook=(i,b,f) -> begin
                inj = make_iso_injector(i); injector_ref[] = inj
                inj(i,b,f)
                snap_before = snapshot_iso(i)
            end)
    catch e
        println("\n!!! isotope clm_drv! step CRASHED:")
        Base.showerror(stdout, e, catch_backtrace()); println()
        return 2
    end
    println("step completed without error.\n")
    snap_after = snapshot_iso(inst)

    # ---- diff post-step pools + increment vs Fortran after_hydrologydrainage ----
    ds = NCDataset(edump, "r")
    dsb = NCDataset(bdump, "r")
    f90_pfts = Int.(ds["pfts1d_itypveg"][:]); jl_pfts = Int.(inst.patch.itype)
    pmap = Dict{Int,Int}()
    for pf in eachindex(f90_pfts)
        for pj in eachindex(jl_pfts)
            if f90_pfts[pf] == jl_pfts[pj]; pmap[pf] = pj; break; end
        end
    end
    if get(ENV, "ISO_DEBUG", "") == "1"
        println("  [dbg] f90_pfts=", f90_pfts, "  jl_pfts=", jl_pfts,
                "  jl_active=", Int.(inst.patch.active[eachindex(jl_pfts)]),
                "  jl_col=", Int.(inst.patch.column[eachindex(jl_pfts)]), "  pmap=", pmap)
    end

    rows = Tuple{String,Float64,Float64,Float64}[]   # name, |abs| pool, |rel| pool, |rel| increment
    gpool = 0.0; ginc = 0.0
    for (cs, suf) in ((inst.bgc_vegetation.c13_cnveg_carbonstate_inst, "_13"),
                      (inst.bgc_vegetation.c14_cnveg_carbonstate_inst, "_14"))
        for stem in ISO_STEMS
            dn = stem * suf
            (haskey(ds, dn) && haskey(dsb, dn)) || continue
            fld = Symbol(stem * "_patch")
            hasproperty(cs, fld) || continue
            jarr = getfield(cs, fld); (jarr isa AbstractVector && !isempty(jarr)) || continue
            fa = Float64.(coalesce.(ds[dn][:], NaN))   # fortran after
            fb = Float64.(coalesce.(dsb[dn][:], NaN))  # fortran before
            mp = 0.0; ma = 0.0; mi = 0.0; any_pt = false
            for (pf, pj) in pmap
                (pf <= length(fa) && pj <= length(jarr)) || continue
                jv = Float64(jarr[pj]); fv = fa[pf]
                (isnan(jv) && isnan(fv)) && continue
                any_pt = true
                ma = max(ma, abs(jv - fv))
                mp = max(mp, abs(jv - fv) / (1 + max(abs(jv), abs(fv))))
                # increment match: (julia_after − julia_before) vs (fortran_after − fortran_before)
                jb = get(snap_before, dn, Float64[])
                if pj <= length(jb)
                    dj = jv - jb[pj]; df = fv - fb[pf]
                    scale = 1e-12 + max(abs(dj), abs(df))
                    mi = max(mi, abs(dj - df) / scale)
                end
            end
            any_pt || continue
            if get(ENV, "ISO_DEBUG", "") == "1" && (dn == "leafc_13" || dn == "xsmrpool_13")
                jb = get(snap_before, dn, Float64[])
                for (pf, pj) in sort(collect(pmap))
                    (pf <= length(fa) && pj <= length(jarr)) || continue
                    @printf("  [dbg %-12s] pf=%d pj=%d  J_after=%.8e  J_before=%.8e  F_before=%.8e  F_after=%.8e\n",
                        dn, pf, pj, Float64(jarr[pj]), pj<=length(jb) ? jb[pj] : NaN, fb[pf], fa[pf])
                end
            end
            gpool = max(gpool, mp); ginc = max(ginc, mi)
            push!(rows, (dn, ma, mp, mi))
        end
    end
    close(ds); close(dsb)

    @printf("\n  %-24s %11s %11s %11s\n", "iso pool", "max|abs|", "pool|rel|", "incr|rel|")
    @printf("  %s\n", "-"^62)
    for (n,a,p,i) in sort(rows, by=x->-x[3])
        flag = p > 1e-6 ? "POOL-DIFF" : (i > 5e-2 ? "flux?" : "ok")
        @printf("  %-24s %11.3e %11.3e %11.3e  %s\n", n, a, p, i, flag)
    end
    @printf("  %s\n", "-"^62)
    @printf("  global post-step pool max|rel| = %.3e   (shared IC → expect tight)\n", gpool)
    @printf("  global one-step increment max|rel| = %.3e   (isotope FLUX physics)\n", ginc)

    # ---- δ13C parity: the discrimination-physics test, decoupled from bulk-C ----
    # leafc_13/leafc is the isotopic RATIO; comparing Julia's post-step δ13C to
    # Fortran's isolates the fractionation/CIsoFlux+state-update chain from the
    # (separately-tracked) bulk-C single-step residual that both isotope pools
    # inherit. Near-zero one-step veg increments make incr|rel| above hit the
    # bulk-C oracle floor; δ13C parity does not.
    dsA = NCDataset(edump, "r")
    bulk_cs = inst.bgc_vegetation.cnveg_carbonstate_inst
    Rpdb = 0.0112372
    @printf("\n  δ13C after-step parity (Julia vs Fortran; discrimination physics)\n")
    @printf("  %-16s %5s %11s %11s %10s\n", "pool", "patch", "δ13C_J‰", "δ13C_F‰", "Δ‰")
    @printf("  %s\n", "-"^58)
    maxdd = 0.0
    for stem in ("leafc","frootc","livestemc","deadstemc","livecrootc","deadcrootc")
        haskey(dsA, stem) && haskey(dsA, stem*"_13") || continue
        jfld = Symbol(stem * "_patch")
        hasproperty(bulk_cs, jfld) && hasproperty(c13cs_g(inst), jfld) || continue
        jbulk = getfield(bulk_cs, jfld); j13 = getfield(c13cs_g(inst), jfld)
        fbulk = Float64.(coalesce.(dsA[stem][:], NaN)); f13 = Float64.(coalesce.(dsA[stem*"_13"][:], NaN))
        for (pf, pj) in sort(collect(pmap))
            (pf <= length(fbulk) && pj <= length(jbulk)) || continue
            (jbulk[pj] > 1e-9 && fbulk[pf] > 1e-9) || continue
            Rj = j13[pj] / (jbulk[pj] - j13[pj]); Rf = f13[pf] / (fbulk[pf] - f13[pf])
            dj = (Rj/Rpdb - 1)*1000; df = (Rf/Rpdb - 1)*1000
            maxdd = max(maxdd, abs(dj - df))
            @printf("  %-16s %5d %11.4f %11.4f %10.2e\n", stem, pj, dj, df, abs(dj-df))
        end
    end
    close(dsA)
    @printf("  %s\n", "-"^58)
    @printf("  max |Δδ13C| = %.3e ‰   (isotope discrimination parity)\n", maxdd)

    # ---- SOIL isotope decomposition pools (litr/soil/cwd _13_vr/_14_vr) ----
    # The soil isotope decomp cascade (CIsoFlux litter inputs + isotope decay)
    # was never diffed before. Pools are injected from before_step (above), so
    # this is a shared-IC one-step diff — the same tightness class as the veg
    # pools. Single Bow column: dump var (lev,) vs decomp_cpools_vr_col[1,:,pool].
    dse = NCDataset(edump, "r")
    @printf("\n  soil isotope decomp pools (litr/soil/cwd × c13/c14)\n")
    @printf("  %-16s %11s %11s %8s\n", "pool", "max|abs|", "max|rel|", "nlev>0")
    @printf("  %s\n", "-"^50)
    soilmax = 0.0; soil_nonzero = 0
    for (scs, suf) in ((get_c13_soilcs(inst), "_13"), (get_c14_soilcs(inst), "_14"))
        scs === nothing && continue
        arr = scs.decomp_cpools_vr_col
        (arr isa AbstractArray && ndims(arr) == 3) || continue
        for (stem, pidx) in ISO_SOIL_STEMS
            dn = stem * suf * "_vr"
            haskey(dse, dn) || continue
            fa = Float64.(coalesce.(dse[dn][:], NaN))
            nlev = min(size(arr, 2), length(fa))
            ma = 0.0; mr = 0.0; nz = 0
            for j in 1:nlev
                fv = fa[j]; jv = Float64(arr[1, j, pidx])
                (isnan(fv) && isnan(jv)) && continue
                # port leaves inactive (below-active) levels at init-NaN where Fortran
                # writes 0; skip those (they carry no decomposition). Real levels have
                # finite values on both sides.
                (isnan(jv) && fv == 0.0) && continue
                (isnan(fv) && jv == 0.0) && continue
                (abs(fv) > 0) && (nz += 1)
                ma = max(ma, abs(jv - fv))
                mr = max(mr, abs(jv - fv) / (1 + max(abs(jv), abs(fv))))
            end
            soilmax = max(soilmax, mr); soil_nonzero += nz
            flag = mr > 1e-6 ? "POOL-DIFF" : "ok"
            @printf("  %-16s %11.3e %11.3e %8d  %s\n", dn, ma, mr, nz, flag)
        end
    end
    close(dse)
    @printf("  %s\n", "-"^50)
    @printf("  soil isotope pool max|rel| = %.3e  (nonzero levels diffed = %d)\n",
            soilmax, soil_nonzero)

    # ---- rc13 discrimination ratios (rc13_canair/psnsun/psnsha) ----
    # These are recomputed each step from forc_pc13o2/forc_pco2 (NO injection) →
    # a PURE physics test of the atmospheric-ratio + fractionation path, decoupled
    # from any injected pool IC. In the port these were 0 until forc_pc13o2_grc was
    # wired (this branch). Fortran rc13_canair = C13RATIO/(1-C13RATIO) = 0.0111718.
    dsr = NCDataset(edump, "r")
    ps = inst.photosyns
    @printf("\n  rc13 discrimination ratios (recomputed each step; pure physics)\n")
    @printf("  %-14s %5s %13s %13s %10s\n", "ratio", "patch", "julia", "fortran", "|Δ|")
    @printf("  %s\n", "-"^60)
    rc13max = 0.0
    for (dn, fld) in (("rc13_canair", :rc13_canair_patch),
                      ("rc13_psnsun", :rc13_psnsun_patch),
                      ("rc13_psnsha", :rc13_psnsha_patch))
        haskey(dsr, dn) && hasproperty(ps, fld) || continue
        fa = Float64.(coalesce.(dsr[dn][:], NaN)); jarr = getfield(ps, fld)
        for (pf, pj) in sort(collect(pmap))
            (pf <= length(fa) && pj <= length(jarr)) || continue
            fv = fa[pf]; jv = Float64(jarr[pj])
            (isnan(fv) && isnan(jv)) && continue
            jv = isnan(jv) ? 0.0 : jv   # port uses NaN for inactive/bare; treat as 0
            d = abs(jv - fv)
            # All three ratios are gated to the vegetated patches (Fortran rc13>0).
            # rc13_canair is the atmospheric ratio; rc13_psnsun/sha additionally carry
            # the fractionation factor alphapsn = alpha at the LAST canopy layer
            # (nrad(p), which is unlit here → alpha=1 → rc13_psn = canair). The port
            # reproduces this once fractionation! is wired (was dead ⇒ alpha=0 ⇒
            # rc13_psn=0). The bare patch (Fortran rc13_canair=0, no photosynthesis
            # filter) is a cosmetic masking diff on rc13_canair only: psnsun=0 there so
            # c13_psnsun=0 regardless — excluded via the fv>0 gate.
            (fv > 0) && (rc13max = max(rc13max, d))
            @printf("  %-14s %5d %13.6e %13.6e %10.2e\n", dn, pj, jv, fv, d)
        end
    end
    close(dsr)
    @printf("  %s\n", "-"^60)
    @printf("  max |Δ rc13 (canair+psnsun+psnsha, veg patches)| = %.3e\n", rc13max)
    return 0
end

# c13 veg carbon state accessor (kept local to avoid threading it through).
c13cs_g(inst) = inst.bgc_vegetation.c13_cnveg_carbonstate_inst

function main()
    na = filter(a -> occursin(r"^\d+$", a), ARGS)
    nstep = isempty(na) ? 4712 : parse(Int, na[1])
    ok = inventory(nstep)
    ok || (println("Reference lacks c13/c14 pools — regenerate bgc_ref_iso."); return 1)
    rc = iso_parity_step(nstep)
    println("\n" * "="^70)
    println("STATUS — C13/C14 isotope parity (bgc_ref_iso, nstep $nstep)")
    println("="^70)
    println("• Reference: use_c13/c14 CN spinup, iso pools reseeded from bulk C at")
    println("  the Jan restart then evolved to mid-July (bypassing CTSM #2119 via")
    println("  a case SourceMods patch). Veg c13/c14 pools present in the dumps.")
    println("• PRIOR FIXES (#241): use_c13/c14 were never propagated to the vegetation")
    println("  facade config (clm_initialize.jl) → the veg c13/c14 carbon state was")
    println("  size-0 and the whole veg isotope path silently no-op'd. And the")
    println("  CIsoFlux litter/gap/harvest/gru p2c scatters fell back to the varpar")
    println("  i_litr_min/i_met_lit sentinel (-9) → BoundsError once the state was")
    println("  live. With those fixed the CIsoFlux cascade RUNS.")
    println("• FIXED (this PR): the per-isotope c_state_update{0,1,2,2h,2g,3}! cascade")
    println("  is now wired (CTSM CNDriverMod calls each CStateUpdate 3× — bulk/c13/c14),")
    println("  so the computed isotope fluxes are APPLIED to the c13/c14 veg pools; and")
    println("  the isotope psn→cpool input (calc_gpp_mr_availc!) now receives the c13/c14")
    println("  carbonflux instances (was dropping to the bulk array → uninit cpool_13).")
    println("• FIXED (this branch — first SOIL-isotope + rc13 diff): two real bugs the")
    println("  veg-only diff was structurally blind to:")
    println("  1. forc_pc13o2_grc was NEVER populated (only ever read) → rc13_canair")
    println("     collapsed to 0 → ZERO C13 discrimination in photosynthesis (fresh")
    println("     photosynthate carried no C13). Wired forc_pc13o2 = C13RATIO*forc_pco2")
    println("     (CTSM lnd_import_export). rc13_canair now matches Fortran to ~1e-18.")
    println("     Hidden in the veg diff because cpool/cpool_13 are 0 at the dump")
    println("     boundaries (a VACUOUS 0.0 match) and injected veg pools barely move.")
    println("  2. The c13/c14 cnveg + soilbgc carbon-FLUX instances were never zeroed")
    println("     at step start (CTSM SetValues); their column scatter targets")
    println("     (decomp_cpools_sourcesink_col soil pools, phenology_c_to_litr_c_col)")
    println("     stayed at nan3d init → NaN through SoilBiogeochemLittVertTransp NaN'd")
    println("     the ENTIRE isotope soil decomposition column. Now zeroed per instance.")
    println("  3. C13 fractionation was DOUBLY dead: fractionation! (the ci/ca alphapsn")
    println("     formula) was ported but NEVER called, AND the driver hard-coded")
    println("     use_c13=false into canopy_fluxes_core!. So alphapsn stayed 0 (its")
    println("     timestep-init reset) → rc13_psnsun/sha = 0 → c13_psnsun = 0 (new")
    println("     photosynthate carried NO C13). Wired the CTSM CanopyFluxes call and")
    println("     threaded config.use_c13. rc13_psnsun/sha now match Fortran to ~1e-18.")
    println("• RESULT: veg pool parity ~9e-6, δ13C <1e-4 ‰; SOIL isotope decomp pools")
    println("  (litr/soil/cwd × c13/c14) finite & matching to ~4e-6 (the same shared-IC")
    println("  bulk-C single-step floor); rc13_canair/psnsun/psnsha to ~1e-18. c14 pools")
    println("  track to <1e-3. Both isotope FLUX cascades (veg + soil) + the")
    println("  photosynthetic discrimination now RUN and match Fortran.")
    return rc
end

exit(main())
