# gpu_validate_dynharvest_e2e.jl — whole cn_harvest! + cn_harvest_pft_to_column!
# on GPU parity.
#
# Mirrors the conservation scenario in test/test_dyn_harvest.jl (2 tree patches +
# 1 bare, on one column), runs the harvest mortality + p2c scatter on the CPU,
# then adapts every struct + mask + harvest state to GPU and reruns
# on-device, comparing the patch fluxes and the column litter/CWD/product pools.
#
# Run: julia +1.12 --project=scripts scripts/gpu_validate_dynharvest_e2e.jl

using CLM, Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))
include(joinpath(@__DIR__, "gpu_adapt.jl"))   # mf / mfs Float32 down-adaptors
ad(x)    = mf(device_array_type(), x)
to(x)    = device_array_type()(x)
dmask(m) = to(collect(Bool, m))

const NP, NC, NG, NLEV, NPOOL = 3, 1, 1, 1, 3
const ILMIN, ILMAX, IMET = 1, 3, 1
const DT = 1800.0

function build()
    patch = CLM.PatchData{Float64}()
    patch.itype = [1, 7, 0]; patch.column = [1, 1, 1]
    patch.gridcell = [1, 1, 1]; patch.wtcol = [0.4, 0.5, 0.1]

    lf_f = zeros(8, NPOOL); fr_f = zeros(8, NPOOL)
    for r in (2, 8); lf_f[r, :] .= [0.5, 0.3, 0.2]; fr_f[r, :] .= [0.4, 0.4, 0.2]; end
    pftcon = (lf_f = lf_f, fr_f = fr_f)

    sbs = CLM.SoilBiogeochemStateData{Float64}()
    sbs.leaf_prof_patch = ones(NP, NLEV); sbs.froot_prof_patch = ones(NP, NLEV)
    sbs.croot_prof_patch = ones(NP, NLEV); sbs.stem_prof_patch = ones(NP, NLEV)

    cs = CLM.CNVegCarbonStateData{Float64}(); ns = CLM.CNVegNitrogenStateData{Float64}()
    cpools_c = (:leafc_patch, :frootc_patch, :livestemc_patch, :deadstemc_patch,
        :livecrootc_patch, :deadcrootc_patch, :xsmrpool_patch, :leafc_storage_patch,
        :frootc_storage_patch, :livestemc_storage_patch, :deadstemc_storage_patch,
        :livecrootc_storage_patch, :deadcrootc_storage_patch, :gresp_storage_patch,
        :leafc_xfer_patch, :frootc_xfer_patch, :livestemc_xfer_patch, :deadstemc_xfer_patch,
        :livecrootc_xfer_patch, :deadcrootc_xfer_patch, :gresp_xfer_patch)
    for (k, f) in enumerate(cpools_c); setfield!(cs, f, Float64[10.0*k + p for p in 1:NP]); end
    npools_n = (:leafn_patch, :frootn_patch, :livestemn_patch, :deadstemn_patch,
        :livecrootn_patch, :deadcrootn_patch, :retransn_patch, :leafn_storage_patch,
        :frootn_storage_patch, :livestemn_storage_patch, :deadstemn_storage_patch,
        :livecrootn_storage_patch, :deadcrootn_storage_patch, :leafn_xfer_patch,
        :frootn_xfer_patch, :livestemn_xfer_patch, :deadstemn_xfer_patch,
        :livecrootn_xfer_patch, :deadcrootn_xfer_patch)
    for (k, f) in enumerate(npools_n); setfield!(ns, f, Float64[1.0*k + 0.1*p for p in 1:NP]); end

    cf = CLM.CNVegCarbonFluxData{Float64}(); nf = CLM.CNVegNitrogenFluxData{Float64}()
    cflux_patch = (:hrv_leafc_to_litter_patch, :hrv_frootc_to_litter_patch,
        :hrv_livestemc_to_litter_patch, :wood_harvestc_patch, :hrv_livecrootc_to_litter_patch,
        :hrv_deadcrootc_to_litter_patch, :hrv_xsmrpool_to_atm_patch,
        :hrv_leafc_storage_to_litter_patch, :hrv_frootc_storage_to_litter_patch,
        :hrv_livestemc_storage_to_litter_patch, :hrv_deadstemc_storage_to_litter_patch,
        :hrv_livecrootc_storage_to_litter_patch, :hrv_deadcrootc_storage_to_litter_patch,
        :hrv_gresp_storage_to_litter_patch, :hrv_leafc_xfer_to_litter_patch,
        :hrv_frootc_xfer_to_litter_patch, :hrv_livestemc_xfer_to_litter_patch,
        :hrv_deadstemc_xfer_to_litter_patch, :hrv_livecrootc_xfer_to_litter_patch,
        :hrv_deadcrootc_xfer_to_litter_patch, :hrv_gresp_xfer_to_litter_patch)
    for f in cflux_patch; setfield!(cf, f, zeros(NP)); end
    nflux_patch = (:hrv_leafn_to_litter_patch, :hrv_frootn_to_litter_patch,
        :hrv_livestemn_to_litter_patch, :wood_harvestn_patch, :hrv_livecrootn_to_litter_patch,
        :hrv_deadcrootn_to_litter_patch, :hrv_retransn_to_litter_patch,
        :hrv_leafn_storage_to_litter_patch, :hrv_frootn_storage_to_litter_patch,
        :hrv_livestemn_storage_to_litter_patch, :hrv_deadstemn_storage_to_litter_patch,
        :hrv_livecrootn_storage_to_litter_patch, :hrv_deadcrootn_storage_to_litter_patch,
        :hrv_leafn_xfer_to_litter_patch, :hrv_frootn_xfer_to_litter_patch,
        :hrv_livestemn_xfer_to_litter_patch, :hrv_deadstemn_xfer_to_litter_patch,
        :hrv_livecrootn_xfer_to_litter_patch, :hrv_deadcrootn_xfer_to_litter_patch)
    for f in nflux_patch; setfield!(nf, f, zeros(NP)); end
    cf.harvest_c_to_litr_c_col = zeros(NC, NLEV, NPOOL); cf.harvest_c_to_cwdc_col = zeros(NC, NLEV)
    cf.wood_harvestc_col = zeros(NC)
    nf.harvest_n_to_litr_n_col = zeros(NC, NLEV, NPOOL); nf.harvest_n_to_cwdn_col = zeros(NC, NLEV)
    nf.wood_harvestn_col = zeros(NC)
    return (; patch, pftcon, sbs, cs, ns, cf, nf)
end

run!(state, mask, d, pft) = CLM.cn_harvest!(state, mask, d.patch, pft, d.sbs, d.cs, d.ns,
    d.cf, d.nf; dt = DT, is_beg_curr_year = true, nlevdecomp = NLEV,
    i_litr_min = ILMIN, i_litr_max = ILMAX, i_met_lit = IMET)

# ---- CPU reference --------------------------------------------------------
cpu = build()
cpu_state = CLM.DynHarvestState(harvest = [0.25], do_harvest = true,
                                harvest_units = CLM.harvest_unitless_units)
run!(cpu_state, [true, true, true], cpu, cpu.pftcon)
outs_c = [:hrv_leafc_to_litter_patch, :wood_harvestc_patch, :harvest_c_to_litr_c_col,
          :harvest_c_to_cwdc_col, :wood_harvestc_col]
outs_n = [:hrv_retransn_to_litter_patch, :harvest_n_to_litr_n_col, :harvest_n_to_cwdn_col,
          :wood_harvestn_col]
cpu_c = Dict(f => copy(getfield(cpu.cf, f)) for f in outs_c)
cpu_n = Dict(f => copy(getfield(cpu.nf, f)) for f in outs_n)

# ---- Metal device run -----------------------------------------------------
if !gpu_functional()
    println("No GPU backend detected — skipping device leg."); exit(0)
end
dd = build()
dev = (; patch = ad(dd.patch), sbs = ad(dd.sbs), cs = ad(dd.cs), ns = ad(dd.ns),
         cf = ad(dd.cf), nf = ad(dd.nf))
pft_d = (lf_f = to(Float32.(dd.pftcon.lf_f)), fr_f = to(Float32.(dd.pftcon.fr_f)))
dev_state = CLM.DynHarvestState(harvest = to(Float32[0.25]), do_harvest = true,
                                harvest_units = CLM.harvest_unitless_units)
if !(dev.cs.leafc_patch isa device_array_type())
    println("  BLOCKED: a struct did not move to the device under adapt."); exit(1)
end
run!(dev_state, dmask([true, true, true]), dev, pft_d)

# ---- Compare --------------------------------------------------------------
function compare()
    println("  field                            max|Δ|")
    g = 0.0
    for f in outs_c
        a = cpu_c[f]; b = Array(getfield(dev.cf, f))
        dfd = maximum(abs.(Float64.(vec(a)) .- Float64.(vec(b))); init = 0.0); g = max(g, dfd)
        @printf("  cf.%-28s %10.2e\n", f, dfd)
    end
    for f in outs_n
        a = cpu_n[f]; b = Array(getfield(dev.nf, f))
        dfd = maximum(abs.(Float64.(vec(a)) .- Float64.(vec(b))); init = 0.0); g = max(g, dfd)
        @printf("  nf.%-28s %10.2e\n", f, dfd)
    end
    return g
end
gmax = compare()
@printf("\n  global max|CPU-Metal| = %.3e\n", gmax)
println(gmax < 1e-4 ? "★★ cn_harvest! + p2c scatter RUN ON METAL — MATCH CPU" :
                      "RAN on Metal but Δ=$gmax")
