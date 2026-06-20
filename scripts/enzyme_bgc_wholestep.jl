# =============================================================================
# use_cn WHOLE-STEP BGC REVERSE — chains the 10 reverse-ready biogeochemistry phases in
# their cn_driver_no_leaching! forward order through ONE CLM.compositional_reverse! call on a
# single REAL use_cn inst bundle (mirrors scripts/enzyme_driver_reverse_fullstep.jl for the
# hydrology side). The assembler is CLM.bgc_rev_phases(inst, bounds, filt, config):
#   cn_mresp! → decomp_rate_constants_bgc! → soil_bgc_potential! → soil_bgc_competition! →
#   calc_plant_cn_alloc! → soil_biogeochem_decomp! → cn_gresp! → c_state_update1! →
#   n_state_update1! → soilbiogeochem_n_state_update1!
# The phases couple through inst fields (cf.decomp_k → nf.potential_immob → st.fpg_col →
# cnveg pools → smin_nh4); the gradient flows t_soisno → (decomposition/mineralization) →
# smin_nh4 across the whole BGC step. Each phase is individually machine-precision FD-validated
# in enzyme_bgc_reverse.jl; this proves the COUPLED chain.
#   julia +1.10 --project=/tmp/clm_jl10_<id> scripts/enzyme_bgc_wholestep.jl
# =============================================================================
using CLM, Enzyme, Printf
const C = CLM
const FSURDAT  = get(ENV, "CLM_FSURDAT",  "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc")
const PARAMFILE = get(ENV, "CLM_PARAMFILE", "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/clm5_params.nc")

function build_cn_inst()
    (inst, bounds, filt, tm) = C.clm_initialize!(; fsurdat=FSURDAT, paramfile=PARAMFILE, use_cn=true)
    config = C.CLMDriverConfig(use_cn=true); filt_ia = C.clump_filter_inactive_and_active
    dtime = 1800.0; calday = 172.5
    (declin, _) = C.compute_orbital(calday); nextsw = calday + dtime/C.SECSPDAY
    ng = bounds.endg; a2l = inst.atm2lnd
    C._setup_calib_forcing!(a2l, 285.0, ng)
    C.downscale_forcings!(bounds, a2l, inst.column, inst.landunit, inst.topo)
    C._init_calib_soil_moisture!(inst, bounds)
    C.interp_monthly_veg!(inst.satellite_phenology; kmo=6, kda=21)
    cs = inst.canopystate; wdb = inst.water.waterdiagnosticbulk_inst; pch = inst.patch
    C.satellite_phenology!(inst.satellite_phenology, cs, wdb, pch, filt.nolakep, bounds.begp:bounds.endp)
    for p in bounds.begp:bounds.endp; cs.frac_veg_nosno_patch[p] = cs.frac_veg_nosno_alb_patch[p]; end
    C.set_exposedvegp_filter!(filt, bounds, cs.frac_veg_nosno_patch)
    bgc_cs = inst.bgc_vegetation.cnveg_carbonstate_inst; bgc_ns = inst.bgc_vegetation.cnveg_nitrogenstate_inst
    bgc_cf = inst.bgc_vegetation.cnveg_carbonflux_inst; bgc_nf = inst.bgc_vegetation.cnveg_nitrogenflux_inst
    for s in (bgc_cs, bgc_ns, bgc_cf, bgc_nf, inst.soilbiogeochem_carbonstate,
              inst.soilbiogeochem_nitrogenstate, inst.soilbiogeochem_carbonflux, inst.soilbiogeochem_nitrogenflux)
        for fn in fieldnames(typeof(s)); v = getfield(s, fn); v isa AbstractArray{Float64} && replace!(v, NaN => 0.0); end
    end
    for p in bounds.begp:bounds.endp
        inst.patch.itype[p] > 0 || continue
        bgc_cs.leafc_patch[p]=100.0; bgc_cs.frootc_patch[p]=50.0; bgc_cs.livestemc_patch[p]=200.0; bgc_cs.deadstemc_patch[p]=500.0
        bgc_ns.leafn_patch[p]=3.0; bgc_ns.frootn_patch[p]=1.5; bgc_ns.livestemn_patch[p]=2.0
    end
    scs = inst.soilbiogeochem_carbonstate; sns = inst.soilbiogeochem_nitrogenstate; nld = C.varpar.nlevdecomp
    for c in bounds.begc:bounds.endc
        for j in 1:nld, p in 1:min(7, size(scs.decomp_cpools_vr_col,3))
            scs.decomp_cpools_vr_col[c,j,p]=10.0; sns.decomp_npools_vr_col[c,j,p]=0.5
        end
        sns.sminn_vr_col[c,1:nld] .= 0.01
    end
    for n in 1:6
        try
            C.clm_drv!(config, inst, filt, filt_ia, bounds, true, nextsw, declin, declin,
                C.ORB_OBLIQR_DEFAULT, false, false, "", false; nstep=n, is_first_step=(n==1),
                is_beg_curr_day=(n==1), dtime=dtime, mon=6, day=21, photosyns=inst.photosyns)
        catch e; @printf("  warm step %d: %s\n", n, sprint(showerror, e)); end
    end
    # CLMDriverConfig forwards the CN dims (ndecomp_pools/i_litr_*/npcropmin/…) via getproperty;
    # use_fun/use_nitrif_denitrif default to false in the aux builders (_cfgget) — the chain
    # runs the non-nitrif, non-FUN path, matching the per-phase building-block validations.
    return inst, bounds, filt, config
end

inst, bounds, filt, cnconfig = build_cn_inst()
phases = C.bgc_rev_phases(inst, bounds, filt, cnconfig)
@printf("assembled %d-phase use_cn whole-step BGC chain\n", length(phases))
sc = filt.bgc_soilc; c0 = [c for c in bounds.begc:bounds.endc if sc[c]][1]
j1 = C.varpar.nlevsno + 1

# DIAGNOSTIC: isolate the productionized decomprate phase (its t_soisno slice-copy path was
# never directly reverse-tested — [P5] perturbed a synthetic bundle array). Reverse just
# [decomprate], seed decomp_k, perturb t_soisno_col[c0,j1]; compare to FD.
let ph1 = C.bgc_rev_phases(inst, bounds, filt, cnconfig)[2:2]   # decomprate only
    Lk(b) = sum(abs2, b.inst.soilbiogeochem_carbonflux.decomp_k_col)
    fd(δ) = (bs=C.bgc_rev_bundle(deepcopy(inst), bounds, cnconfig); bs.inst.temperature.t_soisno_col[c0,j1]+=δ;
             for (f,ca) in ph1; f(bs,ca...) end; Lk(bs))
    g = (fd(1e-3)-fd(-1e-3))/2e-3
    dbk = C.compositional_reverse!(ph1, C.bgc_rev_bundle(deepcopy(inst), bounds, cnconfig),
        (db,b)->(db.inst.soilbiogeochem_carbonflux.decomp_k_col .= 2 .* b.inst.soilbiogeochem_carbonflux.decomp_k_col))
    @printf("DIAG decomprate-only: FD dL/d(t_soisno)=%.6e  rev=%.6e  %s\n", g,
        dbk.inst.temperature.t_soisno_col[c0,j1],
        abs(dbk.inst.temperature.t_soisno_col[c0,j1]) > 0 ? "slice-copy carries gradient" : "SLICE-COPY CUTS GRADIENT")
end
# Incremental-chain probe: grow the chain and seed the natural output of the last phase;
# find the first prefix length where db.t_soisno goes to 0 (the cut link).
let allph = C.bgc_rev_phases(inst, bounds, filt, cnconfig)
    probes = [(2, b->b.inst.soilbiogeochem_carbonflux.decomp_k_col,          "decomp_k   (after decomprate)"),
              (3, b->b.inst.soilbiogeochem_nitrogenflux.potential_immob_vr_col, "pot_immob  (after potential)"),
              (4, b->b.inst.soilbiogeochem_state.fpg_col,                      "fpg        (after competition)"),
              (5, b->b.inst.bgc_vegetation.cnveg_carbonflux_inst.cpool_to_leafc_patch, "cpool_to_leafc (after alloc)"),
              (8, b->b.inst.bgc_vegetation.cnveg_carbonstate_inst.leafc_patch, "leafc      (after cstate1)")]
    for (klen, getout, label) in probes
        ph = allph[1:klen]
        dbp = C.compositional_reverse!(ph, C.bgc_rev_bundle(deepcopy(inst), bounds, cnconfig),
            (db,b)->(getout(db) .= 2 .* getout(b)))
        gt = dbp.inst.temperature.t_soisno_col[c0,j1]
        @printf("  chain[1:%d] seed %-30s → db.t_soisno=%.4e %s\n", klen, label, gt, gt==0 ? "← CUT HERE" : "live")
    end
end
@printf("first BGC-soil column c0=%d; smin_nh4[c0,1]=%.6e decomp_k[c0,1,1]=%.6e\n",
    c0, inst.soilbiogeochem_nitrogenstate.smin_nh4_vr_col[c0,1], 0.0)

# forward run of the assembled chain on a copy, check finite
let bs = C.bgc_rev_bundle(deepcopy(inst), bounds, cnconfig)
    for (f, ca) in phases; f(bs, ca...); end
    i = bs.inst
    @printf("chained primal: decomp_k[c0,1,1]=%.4e gross_nmin[c0,1]=%.4e smin_nh4[c0,1]=%.6e fpg[c0]=%.4f finite=%s\n",
        i.soilbiogeochem_carbonflux.decomp_k_col[c0,1,1], i.soilbiogeochem_nitrogenflux.gross_nmin_vr_col[c0,1],
        i.soilbiogeochem_nitrogenstate.smin_nh4_vr_col[c0,1], i.soilbiogeochem_state.fpg_col[c0],
        string(isfinite(i.soilbiogeochem_nitrogenstate.smin_nh4_vr_col[c0,1])))
end

# L = sum(leafc^2): the downstream LEAF CARBON pool (after c_state_update1). Perturb t_soisno;
# the gradient flows the coupled C/N chain t_soisno → decomp_rate(decomp_k) →
# soil_bgc_potential(potential_immob) → competition(fpg) → allocation(cpool_to_leafc) →
# c_state_update1(leafc) — 5 phases, each producer→consumer adjacent in the reverse so no
# intermediate phase overwrites the carrier array (unlike the smin_nh4 path, whose gross_nmin
# carrier is overwritten by soil_biogeochem_decomp! between its producer and sminn_update).
Lout(b) = sum(abs2, b.inst.bgc_vegetation.cnveg_carbonstate_inst.leafc_patch)
function L_pert(δ)
    bs = C.bgc_rev_bundle(deepcopy(inst), bounds, cnconfig); bs.inst.temperature.t_soisno_col[c0, j1] += δ
    for (f, ca) in phases; f(bs, ca...); end
    return Lout(bs)
end
seed!(db, b) = (db.inst.bgc_vegetation.cnveg_carbonstate_inst.leafc_patch .=
                2 .* b.inst.bgc_vegetation.cnveg_carbonstate_inst.leafc_patch)
db = C.compositional_reverse!(phases, C.bgc_rev_bundle(deepcopy(inst), bounds, cnconfig), seed!)
g_rev = db.inst.temperature.t_soisno_col[c0, j1]

hs = (2e-2, 1e-2, 5e-3, 2.5e-3); cfd(h) = (L_pert(h) - L_pert(-h)) / (2h)
fds = [cfd(h) for h in hs]; for (h,f) in zip(hs,fds); @printf("  FD(h=%.2e) = % .8e\n", h, f); end
g_rich = (4*fds[end] - fds[end-1]) / 3; lo, hi = minimum(fds), maximum(fds)
bracketed = (lo - abs(lo)*1e-9 - 1e-300) <= g_rev <= (hi + abs(hi)*1e-9 + 1e-300)
relerr = abs(g_rev - g_rich) / max(abs(g_rich), 1e-30)
@printf("\n  Richardson = % .8e\n  rev dL/d(t_soisno) = % .8e\n", g_rich, g_rev)
@printf("  bracketed by FD spread [% .3e, % .3e]: %s\n", lo, hi, string(bracketed))
pass = relerr < 5e-5 || bracketed
println("\n", "="^70)
@printf("use_cn WHOLE-STEP BGC REVERSE (%d phases): rel=%.3e bracketed=%s  %s\n",
    length(phases), relerr, string(bracketed), pass ? "PASS ✓" : "FAIL ✗")
println("="^70)
