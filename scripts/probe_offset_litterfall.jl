# =============================================================================
# probe_offset_litterfall.jl — focused probe of the seasonal-deciduous leaf
# OFFSET litterfall ramp (CNPhenologyMod::CNOffsetLitterfall) against the
# Fortran autumn reference.
#
# Fortran ramp (CNPhenologyMod.F90, non-final-step branch):
#     t1 = dt * 2 / offset_counter^2
#     leafc_to_litter  = prev_leafc_to_litter  + t1*(leafc  - prev_leafc_to_litter *offset_counter)
#     frootc_to_litter = prev_frootc_to_litter + t1*(frootc - prev_frootc_to_litter*offset_counter)
# final step (|offset_counter - dt| <= dt/2): dump the whole pool.
#
# For each ramp step we report, for the grass patch (pft 12):
#   - the INJECTED (before_step) leafc / offset_counter / prev_leafc_to_litter
#   - Fortran's post-step leafc and prev_leafc_to_litter (= the flux it used)
#   - Julia's   post-step leafc and prev_leafc_to_litter
#   - the flux Fortran's own formula predicts from the injected state
# so we can see whether Julia's ramp flux, its state update, or neither, is off.
#
#   julia +1.12 --project=. scripts/probe_offset_litterfall.jl
# =============================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const AUT_DIR  = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_autumn"
const AUT_BASE = 1753153
const DT       = 3600.0
const GF       = 3          # Fortran grass patch index (pft 12)

steps = [1759931, 1759945, 1759993, 1760041, 1760113, 1760185, 1760233, 1760281, 1760290]

println("="^118)
println("OFFSET LITTERFALL RAMP PROBE — grass (pft 12), Bow, autumn 2202")
println("="^118)
@printf("%-9s %-12s | %9s %9s %11s | %9s %11s | %9s %11s | %11s\n",
        "nstep","date","inj leafc","offCnt_d","inj prevL2L","F leafc","F prevL2L","J leafc","J prevL2L","F formula")
println("-"^118)

for n in steps
    b = joinpath(AUT_DIR, "pdump_before_step_n$(n).nc")
    e = joinpath(AUT_DIR, "pdump_after_hydrologydrainage_n$(n).nc")
    (isfile(b) && isfile(e)) || continue

    dsb = NCDataset(b, "r")
    inj_leafc  = Float64(dsb["leafc"][GF])
    inj_frootc = Float64(dsb["frootc"][GF])
    inj_oc     = Float64(dsb["offset_counter"][GF])
    inj_of     = Float64(dsb["offset_flag"][GF])
    inj_prev   = Float64(dsb["prev_leafc_to_litter"][GF])
    close(dsb)

    dse = NCDataset(e, "r")
    f_leafc = Float64(dse["leafc"][GF])
    f_prev  = Float64(dse["prev_leafc_to_litter"][GF])
    f_oc    = Float64(dse["offset_counter"][GF])
    close(dse)

    inst, bounds = run_one_parity_step!(n; use_cn=true, dumpdir=AUT_DIR,
        use_hydrstress=true, use_luna=true,
        step_date=DateTime(2002,1,1) + Hour(n - AUT_BASE),
        forcing_file=replace(FFORCING, "clmforc.2003.nc" => "clmforc.2002.nc"))
    ccs = inst.bgc_vegetation.cnveg_carbonstate_inst
    cfx = inst.bgc_vegetation.cnveg_carbonflux_inst
    gj  = findfirst(==(12), Int.(inst.patch.itype))
    j_leafc = Float64(ccs.leafc_patch[gj])
    j_prev  = hasproperty(cfx, :prev_leafc_to_litter_patch) ?
              Float64(cfx.prev_leafc_to_litter_patch[gj]) : NaN
    j_l2l   = hasproperty(cfx, :leafc_to_litter_patch) ?
              Float64(cfx.leafc_to_litter_patch[gj]) : NaN

    # Fortran's own ramp formula evaluated on the INJECTED state.
    # NOTE offset_counter in the dump is the value BEFORE this step's decrement;
    # CNSeasonDecidPhenology decrements it, THEN CNOffsetLitterfall uses it.
    oc_used = inj_oc - DT
    f_formula = if abs(oc_used - DT) <= DT/2
        inj_leafc / DT                       # final step: dump the pool
    elseif oc_used > 0
        t1 = DT * 2.0 / (oc_used * oc_used)
        inj_prev + t1 * (inj_leafc - inj_prev * oc_used)
    else
        NaN
    end

    dt_s = DateTime(2202,1,1) + Hour(n - AUT_BASE)
    @printf("%-9d %-12s | %9.4f %9.4f %11.3e | %9.4f %11.3e | %9.4f %11.3e | %11.3e\n",
            n, Dates.format(dt_s, "mm-dd HH:MM"), inj_leafc, inj_oc/86400, inj_prev,
            f_leafc, f_prev, j_leafc, j_prev, f_formula)
    @printf("%-22s | Julia leafc_to_litter = %.4e   Fortran (from prevL2L) = %.4e   ratio J/F = %s\n",
            "", j_l2l, f_prev, f_prev != 0 ? @sprintf("%.4f", j_l2l/f_prev) : "n/a")
end
println("-"^118)
