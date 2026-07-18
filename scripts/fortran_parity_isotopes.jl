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
# this harness injects the c13/c14 veg pools from the before_step dump via a
# pre_step_hook (patch-remapped by PFT type, same convention as
# set_patch_1d!), runs ONE use_cn+use_c13+use_c14 clm_drv! step, then diffs the
# Julia c13/c14 pools against the Fortran after_hydrologydrainage dump. It
# reports BOTH the absolute post-step pool parity (IC is shared → this is tight)
# AND the one-step increment (after−before) match, which is the sensitive test
# of the isotope FLUX physics (photosynthesis discrimination + CIsoFlux
# cascade).
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
        close(ds)
    end
end

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
    return 0
end

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
    println("• FIXED (this PR): use_c13/c14 were never propagated to the vegetation")
    println("  facade config (clm_initialize.jl) → the veg c13/c14 carbon state was")
    println("  size-0 and the whole veg isotope path silently no-op'd. And the")
    println("  CIsoFlux litter/gap/harvest/gru p2c scatters fell back to the varpar")
    println("  i_litr_min/i_met_lit sentinel (-9) → BoundsError once the state was")
    println("  live. With those fixed the CIsoFlux cascade now RUNS.")
    println("• REMAINING DIVERGENCE (documented, NOT fixed here): cn_driver.jl never")
    println("  calls c_state_update{0,1,2,2h,2g,3}! on the c13/c14 carbon state — CTSM")
    println("  CNDriverMod calls each CStateUpdate 3× (bulk, c13, c14). So the isotope")
    println("  FLUXES are computed but never applied to the veg POOLS: the one-step")
    println("  c13 pool increment diverges ~100% (Julia pools static; c14 changes only")
    println("  via c14_decay!). Wiring the per-isotope state-update cascade is the next")
    println("  isotope task — this harness validates it against bgc_ref_iso.")
    return rc
end

exit(main())
