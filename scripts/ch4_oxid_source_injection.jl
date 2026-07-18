# ==========================================================================
# CH4 OXIDATION SOURCE-TERM ISOLATION TEST
#
# Purpose: isolate the CH4 SOURCE recompute (oxidation) from transport and from
# the warm-CN production confound. Oxidation (`meth_oxid!` / ch4_oxid in
# ch4Mod.F90) is a pure function of INJECTABLE before-step state:
#   conc_ch4, conc_o2, watsat, h2osoi_vol, smp_l, t_soisno, jwt.
# Nothing from the CN driver enters it, so this is a fully clean recompute.
#
# IMPORTANT — the dumped CH4_OXID is POST-competition. ch4_tran (ch4Mod.F90:3609,
# 3623) overwrites ch4_oxid_depth(c,j) *= ch4stress(c,j) (or *= o2stress(c,j))
# and then += aereoxid*ch4_aere_depth. So the raw kinetic rate from ch4_oxid alone
# is NOT what the dump holds — at Bow's sat column O2STRESS_SAT≈0, scaling the raw
# rate down ~1000x. To make a VALID kernel check we apply Fortran's own dumped
# min(ch4stress,o2stress) to Julia's raw rate, isolating the oxidation kernel from
# the competition/transport solve.
#
# Method (oxidation is the ONLY thing that runs):
#   1. Build a Bow instance; inject the Fortran before-step biogeophysics
#      (T_SOISNO / H2OSOI / SMP via the 16-field oracle + read_fortran_restart!)
#      and the ch4 prognostics (conc_ch4_sat/unsat, conc_o2_sat/unsat) from the
#      pdump_before_step file — EXACTLY Fortran's start-of-ch4() state.
#   2. Compute jwt the way the driver does: get_jwt! for the unsat pass, zeros for
#      the sat pass (matching ch4Mod's `jwt(c)=0` for sat==1).
#   3. Call CLM.ch4_oxid! directly for sat=1 and sat=0 (RAW rate) — nothing else runs.
#   4. Apply Fortran's dumped competition factor min(ch4stress,o2stress) and diff
#      against Fortran's CH4_OXID_{SAT,UNSAT} (post-competition), per layer. A
#      residual here is a pure oxidation-kernel / Henry's-law / MM-kinetics defect;
#      the RAW/dump ratio (printed) is the competition scaling, a transport artifact.
#
# Usage: julia --project=. scripts/ch4_oxid_source_injection.jl [nstep]
# ==========================================================================

include(joinpath(@__DIR__, "fortran_parity_common.jl"))
using Printf

const DUMPDIR_CH4 = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_firech4"
const N0    = 4720
const DATE0 = DateTime(2003, 7, 16, 16)
const FNDEP = "/Users/darri.eythorsson/projects/cesm-inputdata/lnd/clm2/ndepdata/" *
              "fndep_clm_UNIFORM1e-8_0.9x1.25_yr2000_CLMjl-parity.nc"

_dv(ds, n) = haskey(ds, n) ? Array(ds[n]) : nothing

# --- inject the ch4 prognostic conc from the before-step pdump ---
function inject_ch4_before!(inst, bounds, dumpfile::String)
    ch4 = inst.ch4
    nc = bounds.endc
    nlev = CLM.varpar.nlevsoi
    NCDataset(dumpfile, "r") do ds
        col2d!(dst, name) = begin
            a = _dv(ds, name); a === nothing && return
            for c in 1:min(nc, size(a, 2)), j in 1:min(nlev, size(a, 1))
                v = a[j, c]; ismissing(v) || (dst[c, j] = Float64(v))
            end
        end
        col2d!(ch4.conc_ch4_sat_col,   "CONC_CH4_SAT")
        col2d!(ch4.conc_ch4_unsat_col, "CONC_CH4_UNSAT")
        col2d!(ch4.conc_o2_sat_col,    "CONC_O2_SAT")
        col2d!(ch4.conc_o2_unsat_col,  "CONC_O2_UNSAT")
    end
    return inst
end

function run_isolation(nstep::Int)
    step_date = DATE0 + Hour(nstep - N0)
    (inst, bounds, filt, tm) = build_bow_inst(; dtime=3600, start_date=step_date,
        use_cn=true, use_luna=true, use_lch4=true,
        fndep=FNDEP, cnfire_method=:nofire)
    if isempty(inst.photosyns.vcmx25_z_patch)
        inst.photosyns.vcmx25_z_patch = fill(30.0, bounds.endp, CLM.NLEVCAN)
        inst.photosyns.jmx25_z_patch  = fill(60.0, bounds.endp, CLM.NLEVCAN)
    end
    beforefile = joinpath(DUMPDIR_CH4, "pdump_before_step_n$(nstep).nc")
    afterfile  = joinpath(DUMPDIR_CH4, "bgcdump_after_ch4_n$(nstep).nc")
    inject_dump!(inst, bounds, beforefile)       # 16-field biogeophys oracle (T, H2OSOI, SMP)
    inject_ch4_before!(inst, bounds, beforefile) # ch4 conc prognostics (before-step)

    ch4  = inst.ch4
    ss   = inst.soilstate
    temp = inst.temperature
    wsb  = inst.water.waterstatebulk_inst
    nc = bounds.endc
    nlev    = CLM.varpar.nlevsoi
    nlevsno = CLM.varpar.nlevsno
    _sl = (nlevsno + 1):(nlevsno + nlev)   # snow-offset soil slice for t_soisno
    _dl = 1:nlev                            # soil slice

    # h2osoi_vol consistent with the injected h2osoi_liq/ice (as the driver keeps it)
    watsat = Matrix(ss.watsat_col[:, _dl])
    hv     = Matrix(wsb.ws.h2osoi_vol_col[:, _dl])
    smp    = Matrix(ss.smp_l_col[:, _dl])
    tsoi   = Matrix(temp.t_soisno_col[:, _sl])
    mask_soil = collect(Bool, filt.soilc)

    # jwt: get_jwt! for unsat, zeros for sat (matches ch4Mod driver)
    jwt_unsat = zeros(Int, nc)
    CLM.get_jwt!(jwt_unsat, mask_soil, watsat, hv, tsoi, nlev, inst.ch4_params)
    jwt_sat = zeros(Int, nc)

    # Push the sliced biogeophys into the inst arrays ch4_oxid! reads through the
    # struct-typed call (it re-slices watsat/h2osoi_vol/smp_l/t_soisno from args).
    CLM.ch4_oxid!(ch4, inst.ch4_params, mask_soil, watsat, hv, smp, tsoi,
                  jwt_sat, 1, false, nlev, 3600.0)
    CLM.ch4_oxid!(ch4, inst.ch4_params, mask_soil, watsat, hv, smp, tsoi,
                  jwt_unsat, 0, false, nlev, 3600.0)

    oxid_sat_J   = ch4.ch4_oxid_depth_sat_col[1, 1:nlev]
    oxid_unsat_J = ch4.ch4_oxid_depth_unsat_col[1, 1:nlev]

    (oxid_sat_F, oxid_unsat_F, ch4_sat, ch4_uns,
     o2s_sat, ch4s_sat, o2s_uns, ch4s_uns) =
        NCDataset(afterfile, "r") do ds
            (Float64.(_dv(ds, "CH4_OXID_SAT")[1:nlev, 1]),
             Float64.(_dv(ds, "CH4_OXID_UNSAT")[1:nlev, 1]),
             Float64.(_dv(ds, "CONC_CH4_SAT")[1:nlev, 1]),
             Float64.(_dv(ds, "CONC_CH4_UNSAT")[1:nlev, 1]),
             Float64.(_dv(ds, "O2STRESS_SAT")[1:nlev, 1]),
             Float64.(_dv(ds, "CH4STRESS_SAT")[1:nlev, 1]),
             Float64.(_dv(ds, "O2STRESS_UNSAT")[1:nlev, 1]),
             Float64.(_dv(ds, "CH4STRESS_UNSAT")[1:nlev, 1]))
        end

    # Apply Fortran's dumped competition factor min(ch4stress,o2stress) to Julia's
    # RAW rate (ch4_tran multiplies by ch4stress if ch4-limited else o2stress; the
    # min is the effective scaling and matches within the aereoxid add-back, small
    # at a dry sat column). This isolates the kernel from the transport competition.
    comp_sat = min.(ch4s_sat, o2s_sat)
    comp_uns = min.(ch4s_uns, o2s_uns)
    oxid_sat_Jc   = oxid_sat_J   .* comp_sat
    oxid_unsat_Jc = oxid_unsat_J .* comp_uns

    prof_relerr(jj, ff) = begin
        s = maximum(abs.(ff)); s <= 0 && return (0.0, 0.0)
        w = 0.0
        for i in eachindex(jj)
            (isfinite(jj[i]) && isfinite(ff[i])) || continue
            w = max(w, abs(jj[i] - ff[i]) / s)
        end
        (w, s)
    end

    println("="^90)
    println("CH4 OXIDATION SOURCE-TERM ISOLATION — nstep $nstep (Bow)  [only ch4_oxid! runs]")
    println("="^90)
    println("jwt_unsat = ", jwt_unsat[1], "   (sat pass uses jwt=0)")
    println("(J*comp = raw Julia rate * Fortran's dumped min(ch4stress,o2stress) competition)")
    println()
    @printf("%4s %12s %12s %12s | %12s %12s %12s\n",
            "lev", "OXID_SAT J*c", "OXID_SAT F", "comp_sat", "OXID_UNS J*c", "OXID_UNS F", "comp_uns")
    for j in 1:nlev
        @printf("%4d %12.4e %12.4e %12.4e | %12.4e %12.4e %12.4e\n",
                j, oxid_sat_Jc[j], oxid_sat_F[j], comp_sat[j],
                oxid_unsat_Jc[j], oxid_unsat_F[j], comp_uns[j])
    end
    println("-"^90)
    (wsr, _) = prof_relerr(oxid_sat_J, oxid_sat_F)     # raw (uncompensated) — shows the competition offset
    (wur, _) = prof_relerr(oxid_unsat_J, oxid_unsat_F)
    (ws, ss_) = prof_relerr(oxid_sat_Jc, oxid_sat_F)   # competition-compensated = kernel-only residual
    (wu, su)  = prof_relerr(oxid_unsat_Jc, oxid_unsat_F)
    @printf("RAW (no competition): OXID_SAT worst rel.err = %.3e   OXID_UNSAT = %.3e   (= the O2/CH4-stress scaling, a transport artifact)\n", wsr, wur)
    @printf("KERNEL-ONLY (competition-compensated): CH4_OXID_SAT   worst rel.err = %.4e   (|F| scale = %.4e)\n", ws, ss_)
    @printf("KERNEL-ONLY (competition-compensated): CH4_OXID_UNSAT worst rel.err = %.4e   (|F| scale = %.4e)\n", wu, su)
    return (ws, wu)
end

nstep = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 4721
run_isolation(nstep)
