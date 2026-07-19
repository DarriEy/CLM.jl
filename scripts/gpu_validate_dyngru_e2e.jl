# gpu_validate_dyngru_e2e.jl — whole cn_gross_unrep! + p2c scatter on GPU parity.
#
# Two natural-PFT patches + one bare (gated out) on one column. Runs the gross
# unrepresented-landcover mortality + p2c scatter on the CPU, then adapts every
# struct + mask + gross-unrep state to GPU and reruns on-device,
# comparing patch fluxes and the column litter/CWD/product accumulators.
#
# Run: julia +1.12 --project=scripts scripts/gpu_validate_dyngru_e2e.jl

using CLM, Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))
include(joinpath(@__DIR__, "gpu_adapt.jl"))
ad(x)    = mf(device_array_type(), x)
to(x)    = device_array_type()(x)
dmask(m) = to(collect(Bool, m))

const NP, NC, NG, NLEV, NLITR = 3, 1, 1, 1, 3
const ILMIN, ILMAX, IMET = 1, 3, 1
const DT, NOVEG, NC4 = 1800.0, 0, 14
const NATPFT = 15

function build()
    patch = CLM.PatchData{Float64}()
    patch.itype = [2, 7, 0]; patch.gridcell = [1, 1, 1]
    patch.column = [1, 1, 1]; patch.wtcol = [0.4, 0.5, 0.1]

    pconv = zeros(NATPFT); pconv[2+1] = 0.6; pconv[7+1] = 0.25
    pftcon = CLM.PftConGrossUnrep(pconv = pconv,
        lf_f = fill(1.0/NLITR, NATPFT, NLITR), fr_f = fill(1.0/NLITR, NATPFT, NLITR))

    sbs = CLM.SoilBiogeochemStateData{Float64}()
    sbs.leaf_prof_patch = ones(NP, NLEV); sbs.froot_prof_patch = ones(NP, NLEV)
    sbs.croot_prof_patch = ones(NP, NLEV); sbs.stem_prof_patch = ones(NP, NLEV)

    cs = CLM.CNVegCarbonStateData{Float64}(); ns = CLM.CNVegNitrogenStateData{Float64}()
    cpools = (:leafc_patch, :frootc_patch, :livestemc_patch, :deadstemc_patch,
        :livecrootc_patch, :deadcrootc_patch, :xsmrpool_patch, :leafc_storage_patch,
        :frootc_storage_patch, :livestemc_storage_patch, :deadstemc_storage_patch,
        :livecrootc_storage_patch, :deadcrootc_storage_patch, :gresp_storage_patch,
        :leafc_xfer_patch, :frootc_xfer_patch, :livestemc_xfer_patch, :deadstemc_xfer_patch,
        :livecrootc_xfer_patch, :deadcrootc_xfer_patch, :gresp_xfer_patch)
    for (k, f) in enumerate(cpools); setfield!(cs, f, Float64[10.0*k + p for p in 1:NP]); end
    npools = (:leafn_patch, :frootn_patch, :livestemn_patch, :deadstemn_patch,
        :livecrootn_patch, :deadcrootn_patch, :retransn_patch, :leafn_storage_patch,
        :frootn_storage_patch, :livestemn_storage_patch, :deadstemn_storage_patch,
        :livecrootn_storage_patch, :deadcrootn_storage_patch, :leafn_xfer_patch,
        :frootn_xfer_patch, :livestemn_xfer_patch, :deadstemn_xfer_patch,
        :livecrootn_xfer_patch, :deadcrootn_xfer_patch)
    for (k, f) in enumerate(npools); setfield!(ns, f, Float64[1.0*k + 0.1*p for p in 1:NP]); end

    cf = CLM.CNVegCarbonFluxData{Float64}(); nf = CLM.CNVegNitrogenFluxData{Float64}()
    cfl = (:gru_leafc_to_litter_patch, :gru_frootc_to_litter_patch, :gru_livestemc_to_atm_patch,
        :gru_deadstemc_to_atm_patch, :gru_wood_productc_gain_patch, :gru_livecrootc_to_litter_patch,
        :gru_deadcrootc_to_litter_patch, :gru_xsmrpool_to_atm_patch, :gru_leafc_storage_to_atm_patch,
        :gru_frootc_storage_to_atm_patch, :gru_livestemc_storage_to_atm_patch,
        :gru_deadstemc_storage_to_atm_patch, :gru_livecrootc_storage_to_atm_patch,
        :gru_deadcrootc_storage_to_atm_patch, :gru_gresp_storage_to_atm_patch,
        :gru_leafc_xfer_to_atm_patch, :gru_frootc_xfer_to_atm_patch, :gru_livestemc_xfer_to_atm_patch,
        :gru_deadstemc_xfer_to_atm_patch, :gru_livecrootc_xfer_to_atm_patch,
        :gru_deadcrootc_xfer_to_atm_patch, :gru_gresp_xfer_to_atm_patch)
    for f in cfl; setfield!(cf, f, zeros(NP)); end
    nfl = (:gru_leafn_to_litter_patch, :gru_frootn_to_litter_patch, :gru_livestemn_to_atm_patch,
        :gru_deadstemn_to_atm_patch, :gru_wood_productn_gain_patch, :gru_livecrootn_to_litter_patch,
        :gru_deadcrootn_to_litter_patch, :gru_retransn_to_litter_patch, :gru_leafn_storage_to_atm_patch,
        :gru_frootn_storage_to_atm_patch, :gru_livestemn_storage_to_atm_patch,
        :gru_deadstemn_storage_to_atm_patch, :gru_livecrootn_storage_to_atm_patch,
        :gru_deadcrootn_storage_to_atm_patch, :gru_leafn_xfer_to_atm_patch,
        :gru_frootn_xfer_to_atm_patch, :gru_livestemn_xfer_to_atm_patch, :gru_deadstemn_xfer_to_atm_patch,
        :gru_livecrootn_xfer_to_atm_patch, :gru_deadcrootn_xfer_to_atm_patch)
    for f in nfl; setfield!(nf, f, zeros(NP)); end
    cf.gru_c_to_litr_c_col = zeros(NC, NLEV, NLITR); cf.gru_c_to_cwdc_col = zeros(NC, NLEV)
    cf.gru_wood_productc_gain_col = zeros(NC)
    nf.gru_n_to_litr_n_col = zeros(NC, NLEV, NLITR); nf.gru_n_to_cwdn_col = zeros(NC, NLEV)
    nf.gru_wood_productn_gain_col = zeros(NC)
    return (; patch, pftcon, sbs, cs, ns, cf, nf)
end

grf() = (g = zeros(NG, NATPFT); g[1, 2+1] = 0.10; g[1, 7+1] = 0.20; g)

run!(state, mask, d, pft) = CLM.cn_gross_unrep!(mask, state, pft, d.patch, d.cs, d.ns, d.cf,
    d.nf, d.sbs; dt = DT, noveg = NOVEG, nc4_grass = NC4, is_beg_curr_year = true,
    nlevdecomp = NLEV, i_litr_min = ILMIN, i_litr_max = ILMAX, i_met_lit = IMET)

# ---- CPU reference --------------------------------------------------------
cpu = build()
cpu_state = CLM.DynGrossUnrepState(grossunrepfrac = grf(), do_grossunrep = true)
run!(cpu_state, [true, true, true], cpu, cpu.pftcon)
outs_c = [:gru_deadstemc_to_atm_patch, :gru_wood_productc_gain_patch, :gru_leafc_to_litter_patch,
          :gru_c_to_litr_c_col, :gru_c_to_cwdc_col, :gru_wood_productc_gain_col]
outs_n = [:gru_retransn_to_litter_patch, :gru_n_to_litr_n_col, :gru_n_to_cwdn_col,
          :gru_wood_productn_gain_col]
cpu_c = Dict(f => copy(getfield(cpu.cf, f)) for f in outs_c)
cpu_n = Dict(f => copy(getfield(cpu.nf, f)) for f in outs_n)

# ---- Metal device run -----------------------------------------------------
if !gpu_functional()
    println("No GPU backend detected — skipping device leg."); exit(0)
end
dd = build()
dev = (; patch = ad(dd.patch), sbs = ad(dd.sbs), cs = ad(dd.cs), ns = ad(dd.ns),
         cf = ad(dd.cf), nf = ad(dd.nf))
pft_d = CLM.PftConGrossUnrep(pconv = to(Float32.(dd.pftcon.pconv)),
    lf_f = to(Float32.(dd.pftcon.lf_f)), fr_f = to(Float32.(dd.pftcon.fr_f)))
dev_state = CLM.DynGrossUnrepState(grossunrepfrac = to(Float32.(grf())), do_grossunrep = true)
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
println(gmax < 1e-4 ? "★★ cn_gross_unrep! + p2c scatter RUN ON METAL — MATCH CPU" :
                      "RAN on Metal but Δ=$gmax")
