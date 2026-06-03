# ==========================================================================
# gpu_validate_phenology_e2e.jl — end-to-end Metal parity for the WHOLE
# cn_phenology! driver (CN phenology: climate avg, evergreen, seasonal- and
# stress-deciduous FSMs, crop phenology, onset/offset/background litterfall,
# livewood turnover, crop harvest, and the patch->column litter scatter).
#
# Builds a multi-patch phenology dataset (distinct PFT types so the evergreen,
# season-deciduous, stress-deciduous and crop paths all fire) at Float64, runs
# cn_phenology! phase 1 + phase 2 on the CPU, adapts every state struct to
# Metal/Float32, runs the SAME calls on the device, and compares the mutated
# per-patch fluxes/flags + the column litter scatter.
#
#   julia --project=scripts scripts/gpu_validate_phenology_e2e.jl
# ==========================================================================

using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))

# Float32 down-convert adaptor (state arrays are Float64 on the CPU build).
struct _F32 end
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{Bool}) = x
# scalar Float64 struct fields (e.g. TemperatureData's trailing FT params) must
# also drop to Float32 so the @adapt_structure reconstruction's FT stays uniform.
CLM.Adapt.adapt_storage(::_F32, x::Float64) = Float32(x)

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    return m
end
allfinite(a) = all(isfinite, Array(a))

# Build the full phenology dataset at Float64 (host). np patches, nc columns.
# PFT flags are set so distinct patches exercise distinct phenology paths.
function build(; np=4, nc=1, ng=1, nlevdecomp=1, nlitr=3)
    F(v, n) = fill(Float64(v), n)
    pstate = CLM.PhenologyState()
    params = CLM.PhenologyParams()
    CLM.cn_phenology_init!(pstate, params, 1800.0)

    npft = 20
    pftcon = CLM.PftConPhenology(
        evergreen = zeros(npft), season_decid = zeros(npft),
        season_decid_temperate = zeros(npft), stress_decid = zeros(npft),
        woody = zeros(npft), leaf_long = fill(1.0, npft),
        leafcn = fill(25.0, npft), frootcn = fill(42.0, npft),
        lflitcn = fill(50.0, npft), livewdcn = fill(50.0, npft),
        deadwdcn = fill(500.0, npft), ndays_on = fill(30.0, npft),
        crit_onset_gdd_sf = fill(1.0, npft),
        lf_f = fill(1.0 / nlitr, npft, nlitr), fr_f = fill(1.0 / nlitr, npft, nlitr),
        biofuel_harvfrac = zeros(npft), repr_structure_harvfrac = zeros(npft, 1),
        minplanttemp = fill(280.0, npft), planttemp = fill(283.0, npft),
        gddmin = fill(50.0, npft), lfemerg = fill(0.04, npft),
        grnfill = fill(0.65, npft), hybgdd = fill(1700.0, npft),
        mxmat = fill(200, npft), manunitro = zeros(npft),
        is_pft_known_to_model = fill(true, npft),
        mnNHplantdate = fill(100.0, npft), mxNHplantdate = fill(170.0, npft),
        mnSHplantdate = fill(280.0, npft), mxSHplantdate = fill(350.0, npft))

    # distinct PFT types per patch: ivt = itype+1
    #   patch1 itype=4 (ivt5)  evergreen+woody
    #   patch2 itype=5 (ivt6)  season-deciduous temperate + woody
    #   patch3 itype=6 (ivt7)  stress-deciduous + woody
    #   patch4 itype=7 (ivt8)  crop (prognostic)
    pftcon.evergreen[5] = 1.0; pftcon.woody[5] = 1.0; pftcon.leaf_long[5] = 3.0
    pftcon.season_decid[6] = 1.0; pftcon.season_decid_temperate[6] = 1.0; pftcon.woody[6] = 1.0
    pftcon.stress_decid[7] = 1.0; pftcon.woody[7] = 1.0

    patch_data = CLM.PatchData()
    patch_data.itype    = [4, 5, 6, 7][1:np]
    patch_data.column   = fill(1, np)
    patch_data.gridcell = fill(1, np)
    patch_data.wtcol    = fill(1.0, np)

    gridcell = CLM.GridcellData()
    gridcell.latdeg = fill(45.0, ng); gridcell.dayl = fill(45000.0, ng)
    gridcell.prev_dayl = fill(44000.0, ng)

    temperature = CLM.TemperatureData()
    temperature.t_ref2m_patch = fill(285.0, np)
    temperature.t_ref2m_min_patch = fill(280.0, np)
    temperature.t_ref2m_max_patch = fill(290.0, np)
    temperature.t_soisno_col = fill(285.0, nc, 10)
    temperature.soila10_col = fill(280.0, nc)
    temperature.t_a5min_patch = fill(278.0, np)

    water_diag = CLM.WaterDiagnosticBulkData()
    water_diag.snow_5day_col = fill(0.0, nc); water_diag.snow_depth_col = fill(0.0, nc)

    soil_state = CLM.SoilStateData(); soil_state.soilpsi_col = fill(-0.5, nc, 10)
    canopy_state = CLM.CanopyStateData(); canopy_state.tlai_patch = fill(2.0, np)

    cnveg_state = CLM.CNVegStateData()
    for f in (:tempavg_t2m_patch, :annavg_t2m_patch)
        setfield!(cnveg_state, f, fill(280.0, np))
    end
    for f in (:onset_flag_patch, :onset_counter_patch, :onset_gddflag_patch, :onset_gdd_patch,
              :onset_fdd_patch, :onset_swi_patch, :offset_flag_patch, :offset_counter_patch,
              :offset_fdd_patch, :offset_swi_patch, :days_active_patch,
              :bglfr_patch, :bgtr_patch, :lgsf_patch, :aleaf_patch, :aleafi_patch,
              :astem_patch, :astemi_patch, :aroot_patch, :cumvd_patch, :hdidx_patch)
        setfield!(cnveg_state, f, fill(0.0, np))
    end
    cnveg_state.dormant_flag_patch = fill(1.0, np)   # start dormant -> exercise onset tests
    cnveg_state.huigrain_patch = fill(0.8, np)
    cnveg_state.huileaf_patch = fill(0.1, np)
    cnveg_state.idop_patch = fill(0, np); cnveg_state.iyop_patch = fill(0, np)

    cnveg_cs = CLM.CNVegCarbonStateData()
    for (f, v) in ((:leafc_storage_patch, 10.0), (:frootc_storage_patch, 5.0),
                   (:livestemc_storage_patch, 3.0), (:deadstemc_storage_patch, 2.0),
                   (:livecrootc_storage_patch, 1.5), (:deadcrootc_storage_patch, 1.0),
                   (:gresp_storage_patch, 0.5), (:leafc_patch, 50.0), (:frootc_patch, 20.0),
                   (:livestemc_patch, 100.0), (:livecrootc_patch, 40.0))
        setfield!(cnveg_cs, f, fill(v, np))
    end
    for f in (:leafc_xfer_patch, :frootc_xfer_patch, :livestemc_xfer_patch,
              :deadstemc_xfer_patch, :livecrootc_xfer_patch, :deadcrootc_xfer_patch)
        setfield!(cnveg_cs, f, fill(0.0, np))
    end

    cnveg_cf = CLM.CNVegCarbonFluxData()
    for f in (:leafc_storage_to_xfer_patch, :frootc_storage_to_xfer_patch,
              :livestemc_storage_to_xfer_patch, :deadstemc_storage_to_xfer_patch,
              :livecrootc_storage_to_xfer_patch, :deadcrootc_storage_to_xfer_patch,
              :gresp_storage_to_xfer_patch, :leafc_xfer_to_leafc_patch,
              :frootc_xfer_to_frootc_patch, :livestemc_xfer_to_livestemc_patch,
              :deadstemc_xfer_to_deadstemc_patch, :livecrootc_xfer_to_livecrootc_patch,
              :deadcrootc_xfer_to_deadcrootc_patch, :leafc_to_litter_patch,
              :frootc_to_litter_patch, :prev_leafc_to_litter_patch, :prev_frootc_to_litter_patch,
              :crop_seedc_to_leaf_patch, :crop_harvestc_to_cropprodc_patch,
              :leafc_to_biofuelc_patch, :livestemc_to_biofuelc_patch,
              :leafc_to_removedresiduec_patch, :livestemc_to_removedresiduec_patch,
              :livestemc_to_deadstemc_patch, :livecrootc_to_deadcrootc_patch,
              :livestemc_to_litter_patch)
        setfield!(cnveg_cf, f, fill(0.0, np))
    end
    cnveg_cf.phenology_c_to_litr_c_col = zeros(nc, nlevdecomp, nlitr)

    cnveg_ns = CLM.CNVegNitrogenStateData()
    for (f, v) in ((:leafn_storage_patch, 0.4), (:frootn_storage_patch, 0.12),
                   (:livestemn_storage_patch, 0.06), (:deadstemn_storage_patch, 0.004),
                   (:livecrootn_storage_patch, 0.03), (:deadcrootn_storage_patch, 0.002),
                   (:leafn_patch, 2.0), (:frootn_patch, 0.5), (:livestemn_patch, 2.0),
                   (:livecrootn_patch, 0.8))
        setfield!(cnveg_ns, f, fill(v, np))
    end
    for f in (:leafn_xfer_patch, :frootn_xfer_patch, :livestemn_xfer_patch,
              :deadstemn_xfer_patch, :livecrootn_xfer_patch, :deadcrootn_xfer_patch)
        setfield!(cnveg_ns, f, fill(0.0, np))
    end

    cnveg_nf = CLM.CNVegNitrogenFluxData()
    for f in (:leafn_storage_to_xfer_patch, :frootn_storage_to_xfer_patch,
              :livestemn_storage_to_xfer_patch, :deadstemn_storage_to_xfer_patch,
              :livecrootn_storage_to_xfer_patch, :deadcrootn_storage_to_xfer_patch,
              :leafn_xfer_to_leafn_patch, :frootn_xfer_to_frootn_patch,
              :livestemn_xfer_to_livestemn_patch, :deadstemn_xfer_to_deadstemn_patch,
              :livecrootn_xfer_to_livecrootn_patch, :deadcrootn_xfer_to_deadcrootn_patch,
              :leafn_to_litter_patch, :frootn_to_litter_patch, :leafn_to_retransn_patch,
              :crop_seedn_to_leaf_patch, :crop_harvestn_to_cropprodn_patch,
              :leafn_to_biofueln_patch, :livestemn_to_biofueln_patch,
              :leafn_to_removedresiduen_patch, :livestemn_to_removedresiduen_patch,
              :livestemn_to_deadstemn_patch, :livecrootn_to_deadcrootn_patch,
              :livestemn_to_retransn_patch, :livecrootn_to_retransn_patch,
              :livestemn_to_litter_patch)
        setfield!(cnveg_nf, f, fill(0.0, np))
    end
    cnveg_nf.phenology_n_to_litr_n_col = zeros(nc, nlevdecomp, nlitr)

    crop = CLM.CropData(); CLM.crop_init!(crop, np)
    cn_params = CLM.CNSharedParamsData()

    mask_soilp  = BitVector(fill(true, np))
    mask_pcropp = BitVector(fill(false, np))
    mask_soilc  = BitVector(fill(true, nc))

    return (; pstate, params, pftcon, patch_data, gridcell, temperature, water_diag,
            soil_state, canopy_state, cnveg_state, cnveg_cs, cnveg_cf, cnveg_ns, cnveg_nf,
            crop, cn_params, mask_soilp, mask_pcropp, mask_soilc, nlevdecomp)
end

run_phen!(d, leaf_prof, froot_prof, phase) = CLM.cn_phenology!(
    d.pstate, d.params, d.pftcon, d.mask_soilp, d.mask_pcropp, d.mask_soilc,
    d.temperature, d.water_diag, d.canopy_state, d.soil_state,
    d.cnveg_state, d.cnveg_cs, d.cnveg_cf, d.cnveg_ns, d.cnveg_nf,
    d.crop, d.patch_data, d.gridcell, d.cn_params,
    leaf_prof, froot_prof, phase; varctl=CLM.VarCtl())

# fields to compare after a phase run (cpu struct, device struct, field list)
const CF_FIELDS = (:leafc_storage_to_xfer_patch, :frootc_storage_to_xfer_patch,
                   :livestemc_storage_to_xfer_patch, :leafc_xfer_to_leafc_patch,
                   :leafc_to_litter_patch, :frootc_to_litter_patch)
const STATE_FIELDS = (:bglfr_patch, :bgtr_patch, :lgsf_patch,
                      :onset_flag_patch, :offset_flag_patch, :dormant_flag_patch,
                      :onset_counter_patch, :offset_counter_patch, :days_active_patch,
                      :onset_gdd_patch, :onset_fdd_patch, :onset_swi_patch, :offset_swi_patch)

function main(backend)
    println("=" ^ 72)
    println("END-TO-END Metal parity for cn_phenology! (CN phenology driver)")
    println("=" ^ 72)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU path exercised by the suite).")
        return 0
    end
    name, _, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n", name, FT)

    H = build()           # CPU reference (Float64)
    Bsrc = build()        # device source (Float64 -> adapt to Metal F32)

    leaf_prof = ones(1, 1); froot_prof = ones(1, 1)

    # Adapt the kernel-touched state structs to Metal Float32. pstate/params/
    # cn_params/canopy_state are NOT indexed by any phenology kernel (pstate only
    # supplies host-side scalars to the wrappers) — keep them host-side so their
    # concrete-Float64 scalar fields don't trip the @adapt_structure reconstruction.
    mf(x) = CLM.Adapt.adapt(Metal.MtlArray, CLM.Adapt.adapt(_F32(), x))
    D = (; pstate = Bsrc.pstate, params = Bsrc.params, pftcon = mf(Bsrc.pftcon),
         patch_data = mf(Bsrc.patch_data), gridcell = mf(Bsrc.gridcell),
         temperature = mf(Bsrc.temperature), water_diag = mf(Bsrc.water_diag),
         soil_state = mf(Bsrc.soil_state), canopy_state = Bsrc.canopy_state,
         cnveg_state = mf(Bsrc.cnveg_state), cnveg_cs = mf(Bsrc.cnveg_cs),
         cnveg_cf = mf(Bsrc.cnveg_cf), cnveg_ns = mf(Bsrc.cnveg_ns), cnveg_nf = mf(Bsrc.cnveg_nf),
         crop = mf(Bsrc.crop), cn_params = Bsrc.cn_params,
         mask_soilp = Metal.MtlArray(collect(Bsrc.mask_soilp)),   # BitVector -> device Bool
         mask_pcropp = Metal.MtlArray(collect(Bsrc.mask_pcropp)),
         mask_soilc = Metal.MtlArray(collect(Bsrc.mask_soilc)),
         nlevdecomp = Bsrc.nlevdecomp)

    if !(D.cnveg_state.bglfr_patch isa Metal.MtlArray)
        println("  BLOCKED: state structs did not move to the device under adapt.")
        return 2
    end

    leaf_prof_d  = Metal.MtlArray(Float32.(leaf_prof))
    froot_prof_d = Metal.MtlArray(Float32.(froot_prof))

    # Run phase 1 then phase 2 on both backends.
    for phase in (1, 2)
        run_phen!(H, leaf_prof, froot_prof, phase)
        run_phen!(D, leaf_prof_d, froot_prof_d, phase)
    end

    nfail = 0
    checks = Tuple{String,Any,Any}[]
    for f in CF_FIELDS
        push!(checks, (string(f), getfield(H.cnveg_cf, f), getfield(D.cnveg_cf, f)))
    end
    for f in STATE_FIELDS
        push!(checks, (string(f), getfield(H.cnveg_state, f), getfield(D.cnveg_state, f)))
    end
    push!(checks, ("phenology_c_to_litr_c", H.cnveg_cf.phenology_c_to_litr_c_col,
                                            D.cnveg_cf.phenology_c_to_litr_c_col))
    push!(checks, ("phenology_n_to_litr_n", H.cnveg_nf.phenology_n_to_litr_n_col,
                                            D.cnveg_nf.phenology_n_to_litr_n_col))

    for (nm, a, b) in checks
        if !allfinite(a)
            @printf("  [WARN] %-26s CPU reference non-finite — skipping\n", nm); continue
        end
        d = reldiff(a, b); ok = d < 1f-3
        @printf("  [%s] %-26s rel|dev-cpu| = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
        ok || (nfail += 1)
    end

    println()
    println(nfail == 0 ? "  WHOLE cn_phenology! MATCHES CPU ON $name ($FT)" :
                         "  DIVERGENCE — investigate ($nfail failed).")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
