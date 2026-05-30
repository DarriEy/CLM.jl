# Validate the multi-column harness: a batched run of N identical tiled gridcells
# must reproduce, tile-for-tile, a single-gridcell run (they are independent and
# identical). This is the data-parallel layout the GPU port batches over.
#
#   julia --project=. scripts/multicolumn_validate.jl

using CLM
using Printf

const FS = get(ENV, "CLM_FSURDAT",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PF = get(ENV, "CLM_PARAMFILE",
    "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")

function setup_forcing!(a2l, T, ng)
    for g in 1:ng
        a2l.forc_t_not_downscaled_grc[g] = T
        a2l.forc_pbot_not_downscaled_grc[g] = 85000.0
        a2l.forc_th_not_downscaled_grc[g] = T * (100000.0 / 85000.0)^(CLM.RAIR / CLM.CPAIR)
        a2l.forc_rho_not_downscaled_grc[g] = 85000.0 / (CLM.RAIR * T)
        a2l.forc_lwrad_not_downscaled_grc[g] = 300.0
        a2l.forc_vp_grc[g] = 800.0
        a2l.forc_hgt_grc[g] = 30.0
        a2l.forc_hgt_u_grc[g] = 30.0
        a2l.forc_hgt_t_grc[g] = 30.0
        a2l.forc_hgt_q_grc[g] = 30.0
        a2l.forc_topo_grc[g] = 0.0
        a2l.forc_wind_grc[g] = 3.0
        a2l.forc_u_grc[g] = 3.0
        a2l.forc_v_grc[g] = 0.0
        for b in 1:CLM.NUMRAD
            a2l.forc_solad_not_downscaled_grc[g, b] = 200.0
            a2l.forc_solai_grc[g, b] = 80.0
        end
        a2l.forc_solar_not_downscaled_grc[g] = 560.0
        a2l.forc_rain_not_downscaled_grc[g] = 0.0001
        a2l.forc_snow_not_downscaled_grc[g] = 0.0
    end
end

function run(ncopies, ksteps)
    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=FS, paramfile=PF, ncopies=ncopies)
    ng = bounds.endg
    config = CLM.CLMDriverConfig()
    fia = CLM.clump_filter_inactive_and_active
    (declin, _) = CLM.compute_orbital(120.0)
    nsw = 120.0 + 1800.0 / CLM.SECSPDAY
    setup_forcing!(inst.atm2lnd, 285.0, ng)
    CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    for n in 1:ksteps
        CLM.clm_drv!(config, inst, filt, fia, bounds, true, nsw, declin, declin,
            CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
            nstep=n, is_first_step=(n == 1), is_beg_curr_day=(n == 1),
            dtime=1800.0, mon=1, day=1, photosyns=inst.photosyns)
    end
    return inst, bounds
end

ksteps = 5
@printf("Single-gridcell run (ncopies=1)...\n")
inst1, b1 = run(1, ksteps)
nc1 = b1.endc
np1 = b1.endp
ref_tg = copy(inst1.temperature.t_grnd_col)
ref_lh = copy(inst1.energyflux.eflx_lh_tot_patch)
@printf("  nc=%d np=%d  t_grnd[1]=%.6f  LH[1]=%.4f\n", nc1, np1, ref_tg[1], ref_lh[1])

N = 3
@printf("Batched run (ncopies=%d)...\n", N)
instN, bN = run(N, ksteps)
@printf("  nc=%d (=%d x %d?)  np=%d (=%d x %d?)  ng=%d\n",
    bN.endc, N, nc1, bN.endp, N, np1, bN.endg)

# Columns/patches are ordered landunit-major (all gridcells' first landunit, then
# all gridcells' next landunit, ...), NOT grouped by gridcell. So group by
# gridcell and compare each gridcell's elements (in within-gridcell order) to the
# single run. NaN-aware: a NaN in the single run must coincide at the same slot.
function compare(arrN, ref, owner_grc, ngrid)
    err = 0.0
    for g in 1:ngrid
        idx = findall(==(g), owner_grc)
        length(idx) == length(ref) || return Inf
        for i in eachindex(ref)
            a = arrN[idx[i]]; b = ref[i]
            if isnan(a) || isnan(b)
                isequal(a, b) || (err = Inf)
            else
                err = max(err, abs(a - b))
            end
        end
    end
    return err
end
@printf("  col.gridcell    = %s\n", string(instN.column.gridcell[1:bN.endc]))
@printf("  col.landunit    = %s\n", string(instN.column.landunit[1:bN.endc]))
@printf("  col.itype       = %s\n", string(instN.column.itype[1:bN.endc]))
@printf("  nolakec filter  = %s\n", string(Int.(CLM.clump_filter.nolakec[1:bN.endc])))
@printf("  patch.column    = %s\n", string(instN.patch.column[1:bN.endp]))
@printf("  patch.gridcell  = %s\n", string(instN.patch.gridcell[1:bN.endp]))
@printf("  exposedvegp flt = %s\n", string(Int.(CLM.clump_filter.exposedvegp[1:bN.endp])))
tg_err = compare(instN.temperature.t_grnd_col, ref_tg, instN.column.gridcell[1:bN.endc], N)
lh_err = compare(instN.energyflux.eflx_lh_tot_patch, ref_lh, instN.patch.gridcell[1:bN.endp], N)
@printf("  max |t_grnd tile-vs-single| = %.3e\n", tg_err)
@printf("  max |LH tile-vs-single|     = %.3e\n", lh_err)

ok = bN.endc == N * nc1 && bN.endp == N * np1 && bN.endg == N &&
     tg_err < 1e-10 && lh_err < 1e-10
println(ok ? "\nMULTI-COLUMN HARNESS PASSED ✓" : "\nFAILED ✗")
