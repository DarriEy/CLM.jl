# ==========================================================================
# fortran_parity_isotopes.jl — Tier-1 isotope / irrigation parity probe.
#
# Validates the Julia port against the Fortran reference run produced with
# irrigate=.true. (and an attempt at use_c13/use_c14=.true.) in the fresh dir
# tier1_iso/. See docs/FULL_PLATFORM_COVERAGE_PLAN.md Sprint 2.
#
# FINDINGS BAKED IN (see the report printed at the end):
#   * use_c13/use_c14 = .true. CRASH (SIGTRAP, rc=133) at restart read: the
#     branch restart 2202-07-16-57600 has NO c13_*/c14_* fields, and the
#     branch path (unlike init_interp's check_interp_non_ciso_to_ciso, which
#     is gated F here) aborts hard on the missing isotope restart variables.
#     => isotopes require a c13/c14-enabled spinup. This is a Tier-2 blocker.
#   * irrigate = .true. RUNS CLEAN (126 dumps, SUCCESSFUL TERMINATION) but is
#     a structural no-op at the Bow column: PFTs are {0 bare, 1 NET tree,
#     12 C3 grass} — NO crop CFTs, and IrrigationMod gates on
#     pftcon%irrigated(m)==1 (crop only). So n_irrig_steps_left==0 and
#     irrig_rate/irrig_rate_demand are NaN/inactive for every step.
#
# Because irrigation is inactive here, the only honest "parity" check is that
# the irrigate=.true. run reproduces the baseline (non-irrigation) CN /
# hydrology state through one Julia clm_drv! step. That is what this probe
# does, plus it dumps the irrigation/isotope inventory for the record.
#
#   julia +1.12 --project=. scripts/fortran_parity_isotopes.jl
# ==========================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const ISO_DUMPDIR = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/tier1_iso"
const NSTEP = 1757873   # first step of the 14-step tier1_iso run

# ----------------------------------------------------------------------------
# 1. Inventory: irrigation + isotope fields in the new reference dump.
# ----------------------------------------------------------------------------
function inventory()
    bdump = joinpath(ISO_DUMPDIR, "pdump_before_step_n$(NSTEP).nc")
    isfile(bdump) || (println("MISSING dump: $bdump"); return false)
    ds = NCDataset(bdump, "r")
    ks = sort(collect(keys(ds)))
    irr = filter(k -> occursin("irrig", lowercase(k)), ks)
    c13 = filter(k -> occursin("c13", lowercase(k)), ks)
    c14 = filter(k -> occursin("c14", lowercase(k)), ks)
    pfts = haskey(ds, "pfts1d_itypveg") ? Int.(ds["pfts1d_itypveg"][:]) : Int[]

    println("=== tier1_iso reference dump inventory (nstep $NSTEP) ===")
    println("  total vars         : ", length(ks))
    println("  PFT itypes         : ", pfts, "   (crop CFTs would be >=15)")
    println("  irrigation fields  : ", irr)
    println("  c13 fields         : ", isempty(c13) ? "NONE" : c13)
    println("  c14 fields         : ", isempty(c14) ? "NONE" : c14)

    # Irrigation state across all 14 steps
    println("\n  irrigation state per step (n_irrig_steps_left / irrig_rate):")
    any_active = false
    for n in NSTEP:(NSTEP+13)
        f = joinpath(ISO_DUMPDIR, "pdump_before_step_n$n.nc")
        isfile(f) || continue
        d = NCDataset(f, "r")
        nis = Int.(d["n_irrig_steps_left"][:])
        ir  = _dumpvar(d, "irrig_rate")
        active = any(nis .> 0) || (ir !== nothing && any(.!isnan.(ir) .& (ir .!= 0)))
        any_active |= active
        close(d)
    end
    println("    => any nonzero irrigation over the 14 steps? ", any_active,
            any_active ? "" : "  (inactive: no crop CFT at this column)")
    close(ds)
    return true
end

# ----------------------------------------------------------------------------
# 2. Run ONE Julia clm_drv! step from the irrigate dump and diff the standard
#    CN pools against the Fortran after_hydrologydrainage dump (same machinery
#    as fortran_parity_cn_summer.jl). Confirms irrigate=.true. == baseline.
# ----------------------------------------------------------------------------
function cn_parity_step()
    bdump = joinpath(ISO_DUMPDIR, "pdump_before_step_n$(NSTEP).nc")
    edump = joinpath(ISO_DUMPDIR, "pdump_after_hydrologydrainage_n$(NSTEP).nc")
    (isfile(bdump) && isfile(edump)) || (println("CN dumps missing"); return 1)

    # nstep 1757873 from the 2202-07-16-57600 restart. BGC run cycles datm
    # 2002-2009 (year_align=2002): model 2202 -> forcing 2002. Same mapping as
    # fortran_parity_cn_summer.jl (which used nstep 1757852).
    println("\nRunning ONE irrigate=true clm_drv! step from tier1_iso IC (nstep=$NSTEP) ...")
    local inst, bounds
    try
        inst, bounds = run_one_parity_step!(NSTEP; use_cn=true, dumpdir=ISO_DUMPDIR,
                            use_hydrstress=true, use_luna=true,
                            step_date=DateTime(2002,1,1)+Hour(NSTEP-1753153),
                            forcing_file=replace(FFORCING, "clmforc.2003.nc"=>"clmforc.2002.nc"))
    catch e
        println("\n!!! clm_drv! step CRASHED:")
        Base.showerror(stdout, e, catch_backtrace()); println()
        return 2
    end
    println("step completed without error.\n")

    ccs = inst.bgc_vegetation.cnveg_carbonstate_inst
    cns = inst.bgc_vegetation.cnveg_nitrogenstate_inst
    scs = inst.soilbiogeochem_carbonstate
    sns = inst.soilbiogeochem_nitrogenstate

    ds = NCDataset(edump, "r")
    jl_pfts = Int.(inst.patch.itype); f_pfts = Int.(ds["pfts1d_itypveg"][:])
    jl2f = Dict{Int,Int}()
    for pj in 1:length(jl_pfts), pf in 1:length(f_pfts)
        if f_pfts[pf] == jl_pfts[pj]; jl2f[pj] = pf; break; end
    end

    gmax = 0.0; rows = Tuple{String,Float64,Float64}[]
    function patchcmp(name, jlarr)
        d = _dumpvar(ds, name); d === nothing && return
        mabs=0.0; mrel=0.0
        for pj in 1:min(length(jlarr), length(f_pfts))
            haskey(jl2f,pj) || continue
            fv=d[jl2f[pj]]; jv=Float64(jlarr[pj])
            (isnan(fv)&&isnan(jv)) && continue
            a=abs(fv-jv); mabs=max(mabs,a); mrel=max(mrel, a/(1+max(abs(fv),abs(jv))))
        end
        gmax=max(gmax,mrel); push!(rows,(name,mabs,mrel))
    end
    function poolcmp(name, arr3, pidx)
        d=_dumpvar(ds,name); d===nothing && return
        nlev=min(length(d),size(arr3,2)); mabs=0.0; mrel=0.0
        for j in 1:nlev
            fv=d[j]; jv=Float64(arr3[1,j,pidx]); (isnan(fv)&&isnan(jv)) && continue
            a=abs(fv-jv); mabs=max(mabs,a); mrel=max(mrel,a/(1+max(abs(fv),abs(jv))))
        end
        gmax=max(gmax,mrel); push!(rows,(name,mabs,mrel))
    end
    function col2cmp(name, arr)
        d=_dumpvar(ds,name); d===nothing && return
        nlev=min(length(d),size(arr,2)); mabs=0.0; mrel=0.0
        for j in 1:nlev
            fv=d[j]; jv=Float64(arr[1,j]); (isnan(fv)&&isnan(jv)) && continue
            a=abs(fv-jv); mabs=max(mabs,a); mrel=max(mrel,a/(1+max(abs(fv),abs(jv))))
        end
        gmax=max(gmax,mrel); push!(rows,(name,mabs,mrel))
    end

    for s in ("leafc","frootc","livestemc","deadstemc","livecrootc","deadcrootc",
              "cpool","xsmrpool","gresp_storage")
        f=Symbol(s*"_patch"); hasproperty(ccs,f) && patchcmp(s, getfield(ccs,f))
    end
    for s in ("leafn","frootn","deadstemn","retransn","npool")
        f=Symbol(s*"_patch"); hasproperty(cns,f) && patchcmp(s, getfield(cns,f))
    end
    for (vn,p) in ("litr1c_vr"=>1,"litr2c_vr"=>2,"litr3c_vr"=>3,"soil1c_vr"=>4,
                   "soil2c_vr"=>5,"soil3c_vr"=>6,"cwdc_vr"=>7)
        poolcmp(vn, scs.decomp_cpools_vr_col, p)
    end
    for (vn,p) in ("litr1n_vr"=>1,"soil1n_vr"=>4,"cwdn_vr"=>7)
        poolcmp(vn, sns.decomp_npools_vr_col, p)
    end
    col2cmp("sminn_vr", sns.sminn_vr_col)
    col2cmp("smin_no3_vr", sns.smin_no3_vr_col)
    col2cmp("smin_nh4_vr", sns.smin_nh4_vr_col)
    close(ds)

    @printf("\n  %-16s %12s %12s\n", "CN field (irrigate=T)", "max|abs|", "max|rel|")
    @printf("  %s\n", "-"^46)
    for (n,a,r) in sort(rows, by=x->-x[3])
        @printf("  %-16s %12.3e %12.3e %s\n", n, a, r, r>1e-6 ? "DIFF" : "ok")
    end
    @printf("  %s\n  global CN max|rel| (irrigate path) = %.3e\n", "-"^46, gmax)
    return 0
end

function main()
    inventory() || return 1
    rc = cn_parity_step()
    println("\n" * "="^70)
    println("STATUS — Tier-1 isotope/irrigation parity")
    println("="^70)
    println("• use_c13/use_c14: CRASH at branch restart read (no c13/c14 on the")
    println("  2202-07-16 restart). Needs a c13/c14-enabled spinup => Tier 2.")
    println("• irrigate=.true.: clean run, but inactive at Bow (no crop CFT).")
    println("  Julia driver irrigation is a stub (clm_driver.jl ~752,~982); the")
    println("  IrrigationData/calc_irrigation_* port exists but is never called.")
    println("• CN/hydrology parity through one irrigate=true step is reported")
    println("  above (irrigation being a no-op here, this == the baseline).")
    return rc
end

exit(main())
