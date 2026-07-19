# gpu_validate_fun_e2e.jl — whole-cnfun!-on-GPU parity.
#
# Builds a use_fun FUN scenario (mirrors test/test_fun.jl make_fun_test_data),
# runs the WHOLE cnfun! (phase-1 kernels + setup + pack → substep → unpack + p2c)
# on the CPU, then adapts every struct + mask to GPU and runs cnfun!
# on-device, comparing the cnveg N/C flux outputs.
#
# Run: julia +1.12 --project=scripts scripts/gpu_validate_fun_e2e.jl

using CLM, Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))
include(joinpath(@__DIR__, "gpu_adapt.jl"))   # mf / mfs Float32 down-adaptors
# mf: arrays + ::FT scalars -> Float32 (keeps concrete-Float64 structs);
# mfs: also converts loose Float64 scalars (for ::FT-scalar structs).
ad(x)   = mf(device_array_type(), x)
ads(x)  = mfs(device_array_type(), x)
to(x)   = device_array_type()(x)
dmask(m) = to(collect(Bool, m))

# ---- Build the FUN scenario (np=1, nc=1, nlevdecomp=2) --------------------
function build()
    np, nc, nlevdecomp, nl, ng = 1, 1, 2, 1, 1
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    pftcon = CLM.PftConFUN(
        leafcn=[25.0], season_decid=[0.0], stress_decid=[0.0], a_fix=[-3.62],
        b_fix=[0.27], c_fix=[25.15], s_fix=[-6.0], akc_active=[1.0], akn_active=[1.0],
        ekc_active=[1.5], ekn_active=[1.5], kc_nonmyc=[2.0], kn_nonmyc=[2.0],
        perecm=[0.5], grperc=[0.3], fun_cn_flex_a=[5.0], fun_cn_flex_b=[200.0],
        fun_cn_flex_c=[80.0], FUN_fracfixers=[0.3], c3psn=[1.0])
    fun_params = CLM.FUNParams(ndays_off=21.0)
    patch = CLM.PatchData(); CLM.patch_init!(patch, np)
    patch.itype[1]=0; patch.column[1]=1; patch.wtcol[1]=1.0
    waterstate = CLM.WaterStateData(); CLM.waterstate_init!(waterstate, nc, np, nl, ng)
    for j in 1:nlevdecomp; waterstate.h2osoi_liq_col[1,j]=10.0; end
    temperature = CLM.TemperatureData(); CLM.temperature_init!(temperature, np, nc, nl, ng)
    for j in 1:nlevdecomp; temperature.t_soisno_col[1,j]=290.0; end
    soilstate = CLM.SoilStateData(); CLM.soilstate_init!(soilstate, np, nc)
    cr = [0.7, 0.3]; for j in 1:nlevdecomp; soilstate.crootfr_patch[1,j]=cr[j]; end
    canopystate = CLM.CanopyStateData(); CLM.canopystate_init!(canopystate, np)
    canopystate.tlai_patch[1]=3.0
    cnveg_state = CLM.CNVegStateData(); CLM.cnveg_state_init!(cnveg_state, np, nc)
    cnveg_state.leafcn_offset_patch[1]=25.0; cnveg_state.plantCN_patch[1]=30.0
    cnveg_state.onset_flag_patch[1]=0.0; cnveg_state.offset_flag_patch[1]=0.0
    cnveg_state.c_allometry_patch[1]=1.0; cnveg_state.n_allometry_patch[1]=0.033
    cnveg_cs = CLM.CNVegCarbonStateData(); CLM.cnveg_carbon_state_init!(cnveg_cs, np, nc, ng)
    cnveg_cs.leafc_patch[1]=100.0; cnveg_cs.frootc_patch[1]=50.0
    cnveg_cs.livestemc_patch[1]=200.0; cnveg_cs.livecrootc_patch[1]=100.0
    cnveg_cs.leafc_storage_patch[1]=10.0; cnveg_cs.leafc_storage_xfer_acc_patch[1]=0.0
    cnveg_cs.storage_cdemand_patch[1]=0.0
    cnveg_cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cnveg_cf, np, nc, ng)
    cnveg_cf.availc_patch[1]=0.001; cnveg_cf.leafc_to_litter_fun_patch[1]=0.001
    cnveg_cf.leafc_storage_to_xfer_patch[1]=0.0
    cnveg_ns = CLM.CNVegNitrogenStateData(); CLM.cnveg_nitrogen_state_init!(cnveg_ns, np, nc, ng)
    cnveg_ns.leafn_patch[1]=4.0; cnveg_ns.frootn_patch[1]=2.0; cnveg_ns.livestemn_patch[1]=1.0
    cnveg_ns.livecrootn_patch[1]=0.5; cnveg_ns.retransn_patch[1]=1.0
    cnveg_ns.leafn_storage_patch[1]=0.5; cnveg_ns.leafn_storage_xfer_acc_patch[1]=0.0
    cnveg_ns.storage_ndemand_patch[1]=0.0
    cnveg_nf = CLM.CNVegNitrogenFluxData()
    CLM.cnveg_nitrogen_flux_init!(cnveg_nf, np, nc, ng; nlevdecomp_full=nlevdecomp)
    cnveg_nf.plant_ndemand_patch[1]=0.0001; cnveg_nf.leafn_storage_to_xfer_patch[1]=0.0
    cnveg_nf.plant_ndemand_retrans_patch[1]=0.0
    soilbgc_nf = CLM.SoilBiogeochemNitrogenFluxData()
    CLM.soil_bgc_nitrogen_flux_init!(soilbgc_nf, nc, nlevdecomp, 1, 1)
    for j in 1:nlevdecomp
        soilbgc_nf.smin_no3_to_plant_vr_col[1,j]=0.5; soilbgc_nf.smin_nh4_to_plant_vr_col[1,j]=0.3
    end
    soilbgc_cf = CLM.SoilBiogeochemCarbonFluxData()
    CLM.soil_bgc_carbon_flux_init!(soilbgc_cf, nc, nlevdecomp, 1, 1)
    soilbgc_ns = CLM.SoilBiogeochemNitrogenStateData()
    CLM.soil_bgc_nitrogen_state_init!(soilbgc_ns, nc, ng, nlevdecomp, 1)
    dzsoi = fill(0.1, nlevdecomp)
    return (; pftcon, fun_params, patch, waterstate, temperature, soilstate, canopystate,
              cnveg_state, cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf, soilbgc_nf, soilbgc_cf,
              soilbgc_ns, dzsoi, np, nc, nlevdecomp)
end

run_fun!(d, msp, msc, dz, pft) = CLM.cnfun!(msp, msc, 1:d.np, 1:d.nc, d.fun_params, pft,
    d.patch, d.waterstate, d.temperature, d.soilstate, d.cnveg_state, d.cnveg_cs, d.cnveg_cf,
    d.cnveg_ns, d.cnveg_nf, d.soilbgc_nf, d.soilbgc_cf, d.canopystate, d.soilbgc_ns;
    dt=1800.0, nlevdecomp=d.nlevdecomp, dzsoi_decomp_vals=dz, use_flexiblecn=false)

# ---- CPU reference --------------------------------------------------------
cpu = build()
run_fun!(cpu, BitVector([true]), BitVector([true]), cpu.dzsoi, cpu.pftcon)
outs = [:Npassive_patch,:Nfix_patch,:Nuptake_patch,:Nactive_patch,:sminn_to_plant_fun_patch,
        :Necm_patch,:Nam_patch,:Nnonmyc_patch,:retransn_to_npool_patch,:free_retransn_to_npool_patch]
outc = [:npp_Nuptake_patch,:npp_growth_patch,:npp_Nfix_patch,:npp_Nretrans_patch,:soilc_change_patch]
cpu_nf = Dict(f => copy(getfield(cpu.cnveg_nf, f)) for f in outs)
cpu_cf = Dict(f => copy(getfield(cpu.cnveg_cf, f)) for f in outc)

# ---- Metal device run -----------------------------------------------------
if !gpu_functional()
    println("No GPU backend detected — skipping device leg."); exit(0)
end
dd = build()
pft_d = CLM.PftConFUN(; (k => to(Float32.(getfield(dd.pftcon,k))) for k in fieldnames(CLM.PftConFUN))...)
dev = (; patch=ad(dd.patch), waterstate=ads(dd.waterstate), temperature=ads(dd.temperature),
         soilstate=ad(dd.soilstate), canopystate=ad(dd.canopystate), cnveg_state=ads(dd.cnveg_state),
         cnveg_cs=ad(dd.cnveg_cs), cnveg_cf=ad(dd.cnveg_cf), cnveg_ns=ad(dd.cnveg_ns),
         cnveg_nf=ad(dd.cnveg_nf), soilbgc_nf=ad(dd.soilbgc_nf), soilbgc_cf=ad(dd.soilbgc_cf),
         soilbgc_ns=ad(dd.soilbgc_ns), fun_params=dd.fun_params, np=dd.np, nc=dd.nc, nlevdecomp=dd.nlevdecomp)
if !(dev.soilstate.crootfr_patch isa device_array_type())
    println("  BLOCKED: a struct did not move to the device under adapt."); exit(1)
end
run_fun!(dev, dmask(BitVector([true])), dmask(BitVector([true])), to(Float32.(dd.dzsoi)), pft_d)

# ---- Compare --------------------------------------------------------------
function compare()
    println("  field                          CPU            Metal          |Δ|")
    g = 0.0
    for f in outs
        a = cpu_nf[f]; b = Array(getfield(dev.cnveg_nf, f))
        dfd = maximum(abs.(Float64.(a) .- Float64.(b)); init=0.0); g = max(g, dfd)
        @printf("  nf.%-26s %13.6e %13.6e %10.2e\n", f, a[1], b[1], dfd)
    end
    for f in outc
        a = cpu_cf[f]; b = Array(getfield(dev.cnveg_cf, f))
        dfd = maximum(abs.(Float64.(a) .- Float64.(b)); init=0.0); g = max(g, dfd)
        @printf("  cf.%-26s %13.6e %13.6e %10.2e\n", f, a[1], b[1], dfd)
    end
    return g
end
gmax = compare()
@printf("\n  global max|CPU-Metal| = %.3e\n", gmax)
println(gmax < 1e-4 ? "★★ WHOLE cnfun! RUNS ON METAL — MATCHES CPU" : "RAN on Metal but Δ=$gmax")
