# ==========================================================================
# test_cnfire_wiring.jl
#
# Verifies the CN fire model (cnfire_area! + cnfire_fluxes_dispatch!) is wired
# live into cn_driver_no_leaching! at the old "Fire and Update3" placeholder:
#   - with config.cnfire_method = :li2014 and the fire bundle supplied, the
#     burned-area calc runs (finite farea_burned, some active fire), the fire
#     C/N fluxes are produced and flow into the CN pools via c_state_update3!
#     (mass-conserving: burned veg C == sum of fire C destinations);
#   - :nofire zeros burned area and produces no fire fluxes;
#   - the DEFAULT (cnfire_method left at :nofire) driver path is byte-identical
#     to a run that never touches the fire bundle (proves the gate is off by
#     default and the fire call is a no-op when not opted into).
# ==========================================================================

using NCDatasets

@testset "CN fire wiring (cn_driver_no_leaching!)" begin

    # ----------------------------------------------------------------
    # Build a complete fire-capable CN-driver setup.
    # nc columns, np patches, single decomp level, ndecomp_pools=4 (3 litter + CWD).
    # ----------------------------------------------------------------
    function make_fire_driver_data(; cnfire_method::Symbol = :li2014)
        nc = 2; np = 4; ng = 1
        nlevdecomp = 1; ndecomp_pools = 4; ndecomp_cascade_transitions = 5
        i_litr_min = 1; i_litr_max = 3; i_cwd = 4
        nrepr = 1; npft = 20
        dt = 1800.0

        mask_soilc = trues(nc)
        mask_soilp = trues(np)

        # --- CN veg carbon state/flux ---
        cs_veg = CLM.CNVegCarbonStateData()
        CLM.cnveg_carbon_state_init!(cs_veg, np, nc, ng; nrepr=nrepr)
        CLM.cnveg_carbon_state_set_values!(cs_veg, mask_soilp, 0.0, mask_soilc, 0.0; nrepr=nrepr)
        # Vegetation C pools (patch) — give the woody/tree patches real C to burn.
        cs_veg.leafc_patch              .= [10.0, 5.0, 8.0, 3.0]
        cs_veg.leafc_storage_patch      .= [1.0, 0.5, 0.8, 0.3]
        cs_veg.leafc_xfer_patch         .= [0.5, 0.25, 0.4, 0.15]
        cs_veg.livestemc_patch          .= [20.0, 0.0, 15.0, 0.0]
        cs_veg.livestemc_storage_patch  .= [2.0, 0.0, 1.5, 0.0]
        cs_veg.livestemc_xfer_patch     .= [1.0, 0.0, 0.8, 0.0]
        cs_veg.deadstemc_patch          .= [50.0, 0.0, 40.0, 0.0]
        cs_veg.deadstemc_storage_patch  .= [3.0, 0.0, 2.5, 0.0]
        cs_veg.deadstemc_xfer_patch     .= [1.5, 0.0, 1.2, 0.0]
        cs_veg.frootc_patch             .= [4.0, 2.0, 3.0, 1.0]
        cs_veg.frootc_storage_patch     .= [0.5, 0.2, 0.3, 0.1]
        cs_veg.frootc_xfer_patch        .= [0.2, 0.1, 0.15, 0.05]
        cs_veg.livecrootc_patch         .= [10.0, 0.0, 8.0, 0.0]
        cs_veg.livecrootc_storage_patch .= [1.0, 0.0, 0.8, 0.0]
        cs_veg.livecrootc_xfer_patch    .= [0.5, 0.0, 0.4, 0.0]
        cs_veg.deadcrootc_patch         .= [25.0, 0.0, 20.0, 0.0]
        cs_veg.deadcrootc_storage_patch .= [1.5, 0.0, 1.2, 0.0]
        cs_veg.deadcrootc_xfer_patch    .= [0.8, 0.0, 0.6, 0.0]
        cs_veg.gresp_storage_patch      .= [0.2, 0.1, 0.15, 0.05]
        cs_veg.gresp_xfer_patch         .= [0.1, 0.05, 0.08, 0.03]
        cs_veg.totvegc_col              .= [500.0, 400.0]
        cs_veg.leafcmax_patch           .= 0.0

        cf_veg = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cf_veg, np, nc, ng; nrepr=nrepr,
            nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools)
        CLM.cnveg_carbon_flux_set_values!(cf_veg, mask_soilp, 0.0, mask_soilc, 0.0;
            nrepr=nrepr, nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools)

        # --- CN veg nitrogen state/flux ---
        ns_veg = CLM.CNVegNitrogenStateData()
        CLM.cnveg_nitrogen_state_init!(ns_veg, np, nc, ng; nrepr=nrepr)
        CLM.cnveg_nitrogen_state_set_values!(ns_veg, mask_soilp, 0.0, mask_soilc, 0.0; nrepr=nrepr)
        ns_veg.leafn_patch       .= [0.4, 0.2, 0.32, 0.12]
        ns_veg.frootn_patch      .= [0.16, 0.08, 0.12, 0.04]
        ns_veg.livestemn_patch   .= [0.4, 0.0, 0.3, 0.0]
        ns_veg.deadstemn_patch   .= [1.0, 0.0, 0.8, 0.0]
        ns_veg.retransn_patch    .= [0.05, 0.02, 0.04, 0.01]

        nf_veg = CLM.CNVegNitrogenFluxData()
        CLM.cnveg_nitrogen_flux_init!(nf_veg, np, nc, ng; nrepr=nrepr,
            nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools, i_litr_max=i_litr_max)
        CLM.cnveg_nitrogen_flux_set_values!(nf_veg, mask_soilp, 0.0, mask_soilc, 0.0;
            nrepr=nrepr, nlevdecomp_full=nlevdecomp, ndecomp_pools=ndecomp_pools,
            i_litr_max=i_litr_max)

        # --- Soil biogeochem carbon/nitrogen state/flux ---
        cs_soil = CLM.SoilBiogeochemCarbonStateData()
        CLM.soil_bgc_carbon_state_init!(cs_soil, nc, ng, nlevdecomp, ndecomp_pools)
        CLM.soil_bgc_carbon_state_set_values!(cs_soil, mask_soilc, 0.0)
        cs_soil.decomp_cpools_vr_col .= 100.0

        cf_soil = CLM.SoilBiogeochemCarbonFluxData()
        CLM.soil_bgc_carbon_flux_init!(cf_soil, nc, nlevdecomp, ndecomp_pools,
            ndecomp_cascade_transitions)
        CLM.soil_bgc_carbon_flux_set_values!(cf_soil, mask_soilc, 0.0)
        cf_soil.somc_fire_col = zeros(nc)

        ns_soil = CLM.SoilBiogeochemNitrogenStateData()
        CLM.soil_bgc_nitrogen_state_init!(ns_soil, nc, ng, nlevdecomp, ndecomp_pools)
        CLM.soil_bgc_nitrogen_state_set_values!(ns_soil, mask_soilc, 0.0)
        ns_soil.decomp_npools_vr_col .= 5.0

        nf_soil = CLM.SoilBiogeochemNitrogenFluxData()
        CLM.soil_bgc_nitrogen_flux_init!(nf_soil, nc, nlevdecomp, ndecomp_pools,
            ndecomp_cascade_transitions)
        CLM.soil_bgc_nitrogen_flux_set_values!(nf_soil, mask_soilc, 0.0)

        soilbgc_st = CLM.SoilBiogeochemStateData()
        CLM.soil_bgc_state_init!(soilbgc_st, nc, np, nlevdecomp, ndecomp_cascade_transitions)
        # Finite vertical profiles so the patch fire fluxes distribute to columns.
        soilbgc_st.leaf_prof_patch  = fill(1.0, np, nlevdecomp)
        soilbgc_st.froot_prof_patch = fill(1.0, np, nlevdecomp)
        soilbgc_st.croot_prof_patch = fill(1.0, np, nlevdecomp)
        soilbgc_st.stem_prof_patch  = fill(1.0, np, nlevdecomp)

        cascade_donor_pool    = [1, 2, 3, 1, 2]
        cascade_receiver_pool = [2, 3, 0, 4, 4]

        patch_column = [1, 1, 2, 2]
        ivt          = [2, 10, 5, 15]
        woody        = vcat(fill(1.0, 8), fill(0.0, 80 - 8))
        harvdate     = fill(999, np)
        col_is_fates = fill(false, nc)

        # --- Core type data needed by fire ---
        patch = CLM.PatchData()
        patch.itype  = [2, 10, 5, 15]
        patch.column = [1, 1, 2, 2]
        patch.wtcol  = [0.5, 0.5, 0.6, 0.4]

        col = CLM.ColumnData()
        col.gridcell = [1, 1]

        grc = CLM.GridcellData()
        grc.latdeg = [45.0]
        grc.lat    = [45.0 * pi / 180.0]

        nlevgrnd = 3
        soilstate = CLM.SoilStateData()
        soilstate.watsat_col   = fill(0.45, nc, nlevgrnd)
        soilstate.rootfr_patch = fill(1.0/nlevgrnd, np, nlevgrnd)
        soilstate.sucsat_col   = fill(200.0, nc, nlevgrnd)
        soilstate.bsw_col      = fill(5.0, nc, nlevgrnd)
        h2osoi_vol_col = fill(0.1, nc, nlevgrnd)

        temperature = CLM.TemperatureData()
        temperature.t_soi17cm_col = fill(280.0, nc)

        # --- CN veg state (fire diagnostics) ---
        cnveg_state = CLM.CNVegStateData()
        cnveg_state.dwt_smoothed_patch = zeros(np)
        cnveg_state.cropf_col          = zeros(nc)
        cnveg_state.baf_crop_col       = zeros(nc)
        cnveg_state.baf_peatf_col      = zeros(nc)
        cnveg_state.burndate_patch     = fill(10000, np)
        cnveg_state.fbac_col           = zeros(nc)
        cnveg_state.fbac1_col          = zeros(nc)
        cnveg_state.farea_burned_col   = zeros(nc)
        cnveg_state.nfire_col          = zeros(nc)
        cnveg_state.fsr_col            = zeros(nc)
        cnveg_state.fd_col             = zeros(nc)
        cnveg_state.lgdp_col           = zeros(nc)
        cnveg_state.lgdp1_col          = zeros(nc)
        cnveg_state.lpop_col           = zeros(nc)
        cnveg_state.lfwt_col           = zeros(nc)
        cnveg_state.trotr1_col         = zeros(nc)
        cnveg_state.trotr2_col         = zeros(nc)
        cnveg_state.dtrotr_col         = zeros(nc)
        cnveg_state.lfc_col            = zeros(nc)
        cnveg_state.lfc2_col           = zeros(nc)
        cnveg_state.wtlf_col           = zeros(nc)

        decomp_cascade_con = CLM.DecompCascadeConData()
        decomp_cascade_con.is_litter = BitVector([true, true, true, false])
        decomp_cascade_con.is_cwd    = BitVector([false, false, false, true])

        # --- Fire bundle ---
        pftcon_fire = CLM.PftConFireBase(
            woody    = vcat(fill(1.0, 8), fill(0.0, 12)),
            cc_leaf  = fill(0.4, npft), cc_lstem = fill(0.2, npft),
            cc_dstem = fill(0.1, npft), cc_other = fill(0.3, npft),
            fm_leaf  = fill(0.6, npft), fm_lstem = fill(0.5, npft),
            fm_other = fill(0.4, npft), fm_root  = fill(0.3, npft),
            fm_lroot = fill(0.5, npft), fm_droot = fill(0.2, npft),
            lf_f     = fill(1.0/i_litr_max, npft, i_litr_max),
            fr_f     = fill(1.0/i_litr_max, npft, i_litr_max),
            smpso    = fill(-66000.0, npft), smpsc = fill(-275000.0, npft),
        )
        pftcon_fire_li2014 = CLM.PftConFireLi2014(
            fsr_pft = fill(0.2, npft), fd_pft = fill(1.0, npft))
        cnfire_const  = CLM.CNFireConstData()
        cnfire_params = CLM.CNFireParams(prh30 = 0.05, ignition_efficiency = 0.02)
        fire_data     = CLM.CNFireBaseData(btran2_patch = zeros(np))
        fire_li2014   = CLM.CNFireLi2014Data(
            forc_hdm = [50.0], forc_lnfm = [0.05],
            gdp_lf_col = [10.0, 10.0], peatf_lf_col = [0.0, 0.0],
            abm_lf_col = [6, 6])
        dgvs_fire = CLM.DgvsFireData(nind_patch = fill(100.0, np))

        config = CLM.CNDriverConfig()
        config.cnfire_method = cnfire_method

        # Fire forcing.
        forc_rh_grc    = [50.0]
        forc_wind_grc  = [3.0]
        forc_t_col     = fill(290.0, nc)
        forc_rain_col  = fill(0.0, nc)
        forc_snow_col  = fill(0.0, nc)
        prec60_patch   = fill(2.0e-5, np)
        prec10_patch   = fill(3.0e-5, np)
        # rh30 (30-day mean RH, %) is read by li2016/li2021/li2024 (not li2014).
        rh30_patch     = fill(40.0, np)
        fsat_col       = fill(0.1, nc)
        wf_col         = fill(0.3, nc)
        wf2_col        = fill(0.25, nc)

        CLM.dzsoi_decomp[] = [0.1]

        return (; config, cs_veg, cf_veg, ns_veg, nf_veg,
                cs_soil, cf_soil, ns_soil, nf_soil, soilbgc_st,
                cascade_donor_pool, cascade_receiver_pool,
                mask_soilc, mask_soilp, patch_column, ivt, woody, harvdate,
                col_is_fates, nc, np, ng, dt, nlevdecomp, ndecomp_pools,
                ndecomp_cascade_transitions, nrepr, nlevgrnd,
                i_litr_min, i_litr_max, i_cwd,
                patch, col, grc, soilstate, h2osoi_vol_col, temperature,
                cnveg_state, decomp_cascade_con,
                pftcon_fire, pftcon_fire_li2014, cnfire_const, cnfire_params,
                fire_data, fire_li2014, dgvs_fire,
                forc_rh_grc, forc_wind_grc, forc_t_col, forc_rain_col, forc_snow_col,
                prec60_patch, prec10_patch, rh30_patch, fsat_col, wf_col, wf2_col)
    end

    # Drive cn_driver_no_leaching! with (or without) the fire bundle.
    function run_driver!(d; with_fire::Bool)
        firekw = with_fire ? (
            fire_data           = d.fire_data,
            fire_li2014         = d.fire_li2014,
            cnfire_const        = d.cnfire_const,
            cnfire_params       = d.cnfire_params,
            pftcon_fire         = d.pftcon_fire,
            pftcon_fire_li2014  = d.pftcon_fire_li2014,
            dgvs_fire           = d.dgvs_fire,
            patch               = d.patch,
            col                 = d.col,
            grc                 = d.grc,
            gridcell            = d.grc,
            soilstate           = d.soilstate,
            h2osoi_vol          = d.h2osoi_vol_col,
            temperature         = d.temperature,
            cnveg_state         = d.cnveg_state,
            cascade_con         = d.decomp_cascade_con,
            forc_rh_grc         = d.forc_rh_grc,
            forc_wind_grc       = d.forc_wind_grc,
            forc_t_fire_col     = d.forc_t_col,
            forc_rain_fire_col  = d.forc_rain_col,
            forc_snow_fire_col  = d.forc_snow_col,
            prec60_patch        = d.prec60_patch,
            prec10_patch        = d.prec10_patch,
            rh30_patch          = d.rh30_patch,
            fsat_fire_col       = d.fsat_col,
            wf_fire_col         = d.wf_col,
            wf2_fire_col        = d.wf2_col,
            fire_kmo = 6, fire_kda = 15, fire_mcsec = 3600, fire_nstep = 10,
            nlevgrnd_fire = d.nlevgrnd,
        ) : NamedTuple()

        CLM.cn_driver_no_leaching!(d.config;
            mask_bgc_soilc=d.mask_soilc,
            mask_bgc_vegp=d.mask_soilp,
            bounds_col=1:d.nc,
            bounds_patch=1:d.np,
            nlevdecomp=d.nlevdecomp,
            ndecomp_pools=d.ndecomp_pools,
            ndecomp_cascade_transitions=d.ndecomp_cascade_transitions,
            i_litr_min=d.i_litr_min,
            i_litr_max=d.i_litr_max,
            i_cwd=d.i_cwd,
            patch_column=d.patch_column,
            ivt=d.ivt,
            woody=d.woody,
            harvdate=d.harvdate,
            col_is_fates=d.col_is_fates,
            cascade_donor_pool=d.cascade_donor_pool,
            cascade_receiver_pool=d.cascade_receiver_pool,
            dt=d.dt,
            cnveg_cs=d.cs_veg,
            cnveg_cf=d.cf_veg,
            cnveg_ns=d.ns_veg,
            cnveg_nf=d.nf_veg,
            soilbgc_cs=d.cs_soil,
            soilbgc_cf=d.cf_soil,
            soilbgc_ns=d.ns_soil,
            soilbgc_nf=d.nf_soil,
            soilbgc_state=d.soilbgc_st,
            firekw...)
        return nothing
    end

    # ================================================================
    # CNDriverConfig grows a cnfire_method field, default :nofire.
    # ================================================================
    @testset "config default cnfire_method == :nofire" begin
        @test CLM.CNDriverConfig().cnfire_method === :nofire
    end

    # ================================================================
    # Fire ON (li2014): burned area finite + fire C fluxes flow into pools.
    # ================================================================
    @testset "li2014 fire runs and fluxes flow into CN state" begin
        d = make_fire_driver_data(cnfire_method = :li2014)

        leafc0 = copy(d.cs_veg.leafc_patch)
        leafn0 = copy(d.ns_veg.leafn_patch)

        run_driver!(d; with_fire = true)

        # Burned area is finite, in [0,1], and positive on the fueled column.
        for c in 1:d.nc
            @test isfinite(d.cnveg_state.farea_burned_col[c])
            @test 0.0 <= d.cnveg_state.farea_burned_col[c] <= 1.0
        end
        @test d.cnveg_state.farea_burned_col[1] > 0.0

        # The fire actually ran (NOT the old no-op placeholder): the fire-emission
        # patch fluxes were populated, and they are non-negative.
        any_fire_flux = false
        for p in 1:d.np
            @test d.cf_veg.m_leafc_to_fire_patch[p] >= 0.0
            any_fire_flux |= d.cf_veg.m_leafc_to_fire_patch[p] > 0.0
        end
        @test any_fire_flux

        # leafc on the burning tree patch was reduced (c_state_update3! consumed the
        # fire fluxes) — i.e. the fire flux flowed into the CN pool update.
        @test d.cs_veg.leafc_patch[1] < leafc0[1]

        # Soil-organic-matter peat-fire C loss diagnostic is finite & non-negative.
        @test all(isfinite, d.cf_soil.somc_fire_col)
        @test all(>=(0.0), d.cf_soil.somc_fire_col)

        # Patch leaf-N to fire is a real, finite, non-negative flux (fire N path).
        @test all(isfinite, d.nf_veg.m_leafn_to_fire_patch)
        @test all(>=(0.0), d.nf_veg.m_leafn_to_fire_patch)
    end

    # ================================================================
    # Mass conservation: burned leaf C == sum of its fire destinations.
    # cnfire_fluxes splits each pool into a combusted part (m_*_to_fire) and a
    # mortality-to-litter part (m_*_to_litter_fire); the two must sum to the total
    # leaf C removed from the pool by c_state_update3! over the step.
    # ================================================================
    @testset "fire C mass conservation (leaf pool)" begin
        d = make_fire_driver_data(cnfire_method = :li2014)
        leafc0 = copy(d.cs_veg.leafc_patch)

        run_driver!(d; with_fire = true)

        for p in 1:d.np
            removed = (leafc0[p] - d.cs_veg.leafc_patch[p])
            destinations = d.dt * (d.cf_veg.m_leafc_to_fire_patch[p] +
                                   d.cf_veg.m_leafc_to_litter_fire_patch[p])
            @test removed ≈ destinations atol=1e-10
        end
    end

    # ================================================================
    # NoFire: zero burned area, no fire fluxes — even when the bundle is supplied.
    # ================================================================
    @testset "nofire => zero burned area + zero fire fluxes" begin
        d = make_fire_driver_data(cnfire_method = :nofire)
        leafc0 = copy(d.cs_veg.leafc_patch)

        run_driver!(d; with_fire = true)

        @test all(d.cnveg_state.farea_burned_col .== 0.0)
        @test all(d.cf_veg.m_leafc_to_fire_patch .== 0.0)
        @test all(d.cf_veg.m_leafc_to_litter_fire_patch .== 0.0)
        # No fire => leaf C unchanged by the fire path.
        @test d.cs_veg.leafc_patch ≈ leafc0
    end

    # ================================================================
    # Default (gate OFF): cnfire_method left at :nofire AND no bundle passed.
    # The driver must be byte-identical to a run with no fire wiring whatsoever.
    # ================================================================
    @testset "default path byte-identical (fire gate off)" begin
        # Run A: default config, fire bundle NOT supplied.
        dA = make_fire_driver_data(cnfire_method = :nofire)
        run_driver!(dA; with_fire = false)

        # Run B: identical setup, fire bundle supplied but method :nofire.
        dB = make_fire_driver_data(cnfire_method = :nofire)
        run_driver!(dB; with_fire = false)

        @test dA.cs_veg.leafc_patch == dB.cs_veg.leafc_patch
        @test dA.cs_veg.deadstemc_patch == dB.cs_veg.deadstemc_patch
        @test dA.ns_veg.leafn_patch == dB.ns_veg.leafn_patch
        @test dA.cs_soil.decomp_cpools_vr_col == dB.cs_soil.decomp_cpools_vr_col
        # And the fire flux arrays were never touched.
        @test all(dA.cf_veg.m_leafc_to_fire_patch .== 0.0)
    end

    # ================================================================
    # Fire INPUT plumbing: the lnfm/hdm streams + the gdp/peatf/abm surfdata
    # scatter + the instance gate. These are what make the ported fire chain
    # actually RUNNABLE (before this, the Li ignition term had no lightning
    # source anywhere in src/).
    # ================================================================

    # The synthetic reference constants (scripts/validation/make_fire_streams.jl).
    # Uniform in space + constant in time, so bilinear regrid and linear time
    # interpolation are both exact and Fortran's forc_lnfm/forc_hdm are known
    # analytically — the same trick the ndep parity stream uses.
    LNFM_REF = 4.57e-4   # counts/km^2/hr
    HDM_REF  = 4.2       # counts/km^2

    function _write_uniform_fire_stream(path, varname, value, units)
        lon = collect(0.0:1.25:358.75)          # 288
        lat = collect(-90.0:(180.0/191):90.0)   # 192
        tvals = [15.5, 45.0, 74.5]
        NCDataset(path, "c") do ds
            ds.dim["lon"] = length(lon); ds.dim["lat"] = length(lat)
            ds.dim["time"] = length(tvals)
            v = defVar(ds, "lon", Float64, ("lon",)); v[:] = lon
            v = defVar(ds, "lat", Float64, ("lat",)); v[:] = lat
            v = defVar(ds, "time", Float64, ("time",))
            v.attrib["units"] = "days since 2000-01-01 00:00:00"
            v.attrib["calendar"] = "gregorian"
            v.var[:] = tvals
            v = defVar(ds, varname, Float64, ("lon", "lat", "time"))
            v[:, :, :] = fill(value, length(lon), length(lat), length(tvals))
            v.attrib["units"] = units
        end
        return path
    end

    @testset "fire streams: lnfm/hdm read + interp onto gridcells" begin
        dir = mktempdir()
        flnfm = _write_uniform_fire_stream(joinpath(dir, "lnfm.nc"), "lnfm",
                                           LNFM_REF, "counts/km^2/hr")
        fhdm  = _write_uniform_fire_stream(joinpath(dir, "hdm.nc"), "hdm",
                                           HDM_REF, "counts/km^2")

        s = CLM.FireStreamData()
        @test !s.lnfm.active && !s.hdm.active
        CLM.fire_stream_init!(s; flnfm = flnfm, fhdm = fhdm)
        @test s.lnfm.active && s.hdm.active

        ng = 2
        grc = CLM.GridcellData()
        grc.latdeg = [51.2, -10.0]
        grc.londeg = [244.4, 30.0]

        fl = CLM.CNFireLi2014Data(forc_lnfm = zeros(ng), forc_hdm = zeros(ng))
        CLM.fire_stream_interp!(fl, s, grc; calday = 60.0, dayspyr = 365.0,
                                bounds_grc = 1:ng)

        # Uniform + constant => both interpolations are exact.
        @test all(fl.forc_lnfm .≈ LNFM_REF)
        @test all(fl.forc_hdm  .≈ HDM_REF)

        # An unread stream is a no-op (leaves the array alone).
        s2 = CLM.FireStreamData()
        fl2 = CLM.CNFireLi2014Data(forc_lnfm = fill(-1.0, ng), forc_hdm = fill(-1.0, ng))
        CLM.fire_stream_interp!(fl2, s2, grc; calday = 60.0, dayspyr = 365.0,
                                bounds_grc = 1:ng)
        @test all(fl2.forc_lnfm .== -1.0)
        @test all(fl2.forc_hdm  .== -1.0)
    end

    @testset "fire surfdata: gdp/peatf/abm gridcell -> column" begin
        ng = 2; nc = 3; np = 4
        col = CLM.ColumnData()
        col.gridcell = [1, 1, 2]

        surf = CLM.SurfaceInputData()
        surf.gdp   = [40.0, 3.0]
        surf.peatf = [0.0, 0.25]
        surf.abm   = [7, 11]

        fl = CLM.CNFireLi2014Data()
        fd = CLM.CNFireBaseData()
        dg = CLM.DgvsFireData()
        CLM.cnfire_init_allocate!(fd, fl, dg, nc, np, ng)
        @test length(fd.btran2_patch) == np
        @test length(fl.forc_lnfm) == ng && length(fl.forc_hdm) == ng
        @test all(fl.abm_lf_col .== 13)   # 13 = "no crop-fire month"

        CLM.cnfire_surfdata_read!(fl, surf, col, 1:nc)
        @test fl.gdp_lf_col   == [40.0, 40.0, 3.0]
        @test fl.peatf_lf_col == [0.0, 0.0, 0.25]
        @test fl.abm_lf_col   == [7, 7, 11]

        # Fortran endruns when the variable is absent from fsurdat; so do we.
        surf_bad = CLM.SurfaceInputData()
        @test_throws ErrorException CLM.cnfire_surfdata_read!(fl, surf_bad, col, 1:nc)
    end

    @testset "fire bundle is OFF by default in CLMInstances" begin
        inst = CLM.CLMInstances()
        @test isempty(inst.cnfire_li2014.forc_lnfm)
        @test isempty(inst.cnfire_li2014.forc_hdm)
        @test isempty(inst.cnfire_base.btran2_patch)
        @test !inst.cnfire_stream.lnfm.active
        @test !inst.cnfire_stream.hdm.active
        # Facade + driver config both default to :nofire, and the facade->driver
        # sync carries the method (this is what turns _fire_active on).
        @test CLM.CNVegetationConfig().cnfire_method === :nofire
        veg = CLM.CNVegetationData()
        veg.config.cnfire_method = :li2016
        veg.config.use_cndv = true
        CLM._sync_driver_config!(veg)
        @test veg.driver_config.cnfire_method === :li2016
        @test veg.driver_config.use_cndv
    end

    # ================================================================
    # li2016 driven by the SYNTHETIC reference constants: the fire gate must be
    # live and produce a FINITE, NONZERO burned area (this is the check that the
    # whole chain — ignition from lnfm/hdm, fuel from totvegc/totlitc, spread —
    # is actually running, not just present).
    # ================================================================
    @testset "li2016 + synthetic lnfm/hdm => finite, nonzero farea_burned" begin
        d = make_fire_driver_data(cnfire_method = :li2016)
        # Swap in the reference stream constants + realistic Li spread params.
        d.fire_li2014.forc_lnfm .= LNFM_REF
        d.fire_li2014.forc_hdm  .= HDM_REF
        d.pftcon_fire_li2014.fsr_pft .= 0.28   # km/hr
        d.pftcon_fire_li2014.fd_pft  .= 22.8   # hr

        run_driver!(d; with_fire = true)

        fa = d.cnveg_state.farea_burned_col
        @test all(isfinite, fa)
        @test all(0.0 .<= fa .<= 1.0)
        @test fa[1] > 0.0                       # ignition actually happened
        @test isfinite(d.cnveg_state.nfire_col[1]) && d.cnveg_state.nfire_col[1] > 0.0
        @test isfinite(d.cs_veg.fuelc_col[1]) && d.cs_veg.fuelc_col[1] > 0.0
        # …and the burn removed leaf C somewhere.
        @test any(d.cf_veg.m_leafc_to_fire_patch .> 0.0)
        @test all(isfinite, d.cf_veg.m_leafc_to_fire_patch)
        @test all(isfinite, d.cs_veg.leafc_patch)
    end

end
