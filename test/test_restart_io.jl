@testset "Restart I/O (restFileMod)" begin
    ng, nl, nc, np = 2, 3, 6, 12

    # Vertical level params must be set before allocating column/state arrays.
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)

    # Build a fully-allocated instance and bounds.
    function build_inst()
        inst = CLM.CLMInstances()
        CLM.clm_instInit!(inst; ng=ng, nl=nl, nc=nc, np=np)
        return inst
    end
    bounds = CLM.BoundsType(begc=1, endc=nc, begp=1, endp=np,
                            begg=1, endg=ng, begl=1, endl=nl)

    # Deterministic, distinct, finite fill so the round-trip has real values.
    fill_seq!(A, base) = (A .= reshape(base .+ (1:length(A)), size(A)); A)

    # Fill the prognostic state that the biogeophys registry covers.
    function seed_state!(inst)
        T = inst.temperature
        fill_seq!(T.t_soisno_col, 250.0)
        fill_seq!(T.t_grnd_col, 270.0)
        fill_seq!(T.t_lake_col, 275.0)
        fill_seq!(T.t_h2osfc_col, 274.0)
        fill_seq!(T.t_veg_patch, 290.0)
        fill_seq!(T.t_ref2m_patch, 288.0)

        ws = inst.water.waterstatebulk_inst.ws
        fill_seq!(ws.h2osoi_liq_col, 10.0)
        fill_seq!(ws.h2osoi_ice_col, 5.0)
        fill_seq!(ws.h2osno_no_layers_col, 1.0)
        fill_seq!(ws.h2osfc_col, 0.5)
        fill_seq!(ws.wa_col, 4000.0)
        fill_seq!(ws.liqcan_patch, 0.1)
        fill_seq!(ws.snocan_patch, 0.2)
        fill_seq!(inst.water.waterstatebulk_inst.int_snow_col, 2.0)

        wd = inst.water.waterdiagnosticbulk_inst
        fill_seq!(wd.snow_depth_col, 0.05)
        fill_seq!(wd.frac_sno_col, 0.3)
        fill_seq!(wd.frac_sno_eff_col, 0.25)

        fill_seq!(inst.soilhydrology.zwt_col, 3.0)
        fill_seq!(inst.soilhydrology.zwt_perched_col, 2.0)

        inst.column.snl .= -3  # integer snow-layer count
        fill_seq!(inst.column.dz, 0.01)
        fill_seq!(inst.column.z, 0.1)
        fill_seq!(inst.column.zi, 0.0)

        fill_seq!(inst.lakestate.lake_icefrac_col, 0.0)
        fill_seq!(inst.lakestate.savedtke1_col, 0.5)

        fill_seq!(inst.canopystate.tlai_patch, 1.0)
        fill_seq!(inst.canopystate.elai_patch, 0.9)
        fill_seq!(inst.canopystate.htop_patch, 5.0)
        return inst
    end

    mktempdir() do dir
        file = joinpath(dir, "clm_restart.nc")

        # --- biogeophys-only round-trip ---
        src = seed_state!(build_inst())
        CLM.write_restart(src, file; bounds=bounds)
        @test isfile(file)

        dst = build_inst()  # fresh, NaN-initialised
        CLM.read_restart!(dst, file; bounds=bounds)

        # Bit-for-bit comparison of every covered field.
        @test dst.temperature.t_soisno_col == src.temperature.t_soisno_col
        @test dst.temperature.t_grnd_col == src.temperature.t_grnd_col
        @test dst.temperature.t_lake_col == src.temperature.t_lake_col
        @test dst.temperature.t_h2osfc_col == src.temperature.t_h2osfc_col
        @test dst.temperature.t_veg_patch == src.temperature.t_veg_patch
        @test dst.temperature.t_ref2m_patch == src.temperature.t_ref2m_patch

        sws = src.water.waterstatebulk_inst.ws
        dws = dst.water.waterstatebulk_inst.ws
        @test dws.h2osoi_liq_col == sws.h2osoi_liq_col
        @test dws.h2osoi_ice_col == sws.h2osoi_ice_col
        @test dws.h2osno_no_layers_col == sws.h2osno_no_layers_col
        @test dws.h2osfc_col == sws.h2osfc_col
        @test dws.wa_col == sws.wa_col
        @test dws.liqcan_patch == sws.liqcan_patch
        @test dws.snocan_patch == sws.snocan_patch
        @test dst.water.waterstatebulk_inst.int_snow_col ==
              src.water.waterstatebulk_inst.int_snow_col

        @test dst.water.waterdiagnosticbulk_inst.snow_depth_col ==
              src.water.waterdiagnosticbulk_inst.snow_depth_col
        @test dst.water.waterdiagnosticbulk_inst.frac_sno_col ==
              src.water.waterdiagnosticbulk_inst.frac_sno_col

        @test dst.soilhydrology.zwt_col == src.soilhydrology.zwt_col
        @test dst.soilhydrology.zwt_perched_col == src.soilhydrology.zwt_perched_col

        # Integer snow-layer count survives the Float64 round-trip exactly.
        @test dst.column.snl == src.column.snl
        @test eltype(dst.column.snl) <: Integer
        @test dst.column.dz == src.column.dz
        @test dst.column.z == src.column.z
        @test dst.column.zi == src.column.zi

        @test dst.lakestate.lake_icefrac_col == src.lakestate.lake_icefrac_col
        @test dst.lakestate.savedtke1_col == src.lakestate.savedtke1_col

        @test dst.canopystate.tlai_patch == src.canopystate.tlai_patch
        @test dst.canopystate.elai_patch == src.canopystate.elai_patch
        @test dst.canopystate.htop_patch == src.canopystate.htop_patch
    end

    # --- NaN <-> fill-value round-trips as NaN (not the sentinel) ---
    mktempdir() do dir
        file = joinpath(dir, "clm_restart_nan.nc")
        src = build_inst()
        src.temperature.t_grnd_col .= NaN
        src.temperature.t_grnd_col[1] = 271.5
        CLM.write_restart(src, file; bounds=bounds)
        dst = build_inst()
        CLM.read_restart!(dst, file; bounds=bounds)
        @test dst.temperature.t_grnd_col[1] == 271.5
        @test all(isnan, dst.temperature.t_grnd_col[2:end])
    end

    # --- CN pools round-trip when use_cn=true ---
    mktempdir() do dir
        file = joinpath(dir, "clm_restart_cn.nc")
        src = build_inst()
        scs = src.soilbiogeochem_carbonstate
        sns = src.soilbiogeochem_nitrogenstate
        fill_seq!(scs.decomp_cpools_vr_col, 100.0)
        fill_seq!(scs.ctrunc_vr_col, 0.0)
        fill_seq!(sns.decomp_npools_vr_col, 10.0)
        fill_seq!(sns.sminn_vr_col, 1.0)
        fill_seq!(sns.smin_no3_vr_col, 0.5)
        fill_seq!(sns.smin_nh4_vr_col, 0.3)

        vcs = src.bgc_vegetation.cnveg_carbonstate_inst
        vns = src.bgc_vegetation.cnveg_nitrogenstate_inst
        fill_seq!(vcs.leafc_patch, 50.0)
        fill_seq!(vcs.frootc_patch, 40.0)
        fill_seq!(vcs.livestemc_patch, 30.0)
        fill_seq!(vcs.deadstemc_patch, 20.0)
        fill_seq!(vcs.cpool_patch, 1.0)
        fill_seq!(vns.leafn_patch, 5.0)
        fill_seq!(vns.frootn_patch, 4.0)
        fill_seq!(vns.npool_patch, 0.1)

        CLM.write_restart(src, file; bounds=bounds, use_cn=true)
        dst = build_inst()
        CLM.read_restart!(dst, file; bounds=bounds, use_cn=true)

        dcs = dst.soilbiogeochem_carbonstate
        dns = dst.soilbiogeochem_nitrogenstate
        @test dcs.decomp_cpools_vr_col == scs.decomp_cpools_vr_col
        @test dcs.ctrunc_vr_col == scs.ctrunc_vr_col
        @test dns.decomp_npools_vr_col == sns.decomp_npools_vr_col
        @test dns.sminn_vr_col == sns.sminn_vr_col
        @test dns.smin_no3_vr_col == sns.smin_no3_vr_col
        @test dns.smin_nh4_vr_col == sns.smin_nh4_vr_col

        dvcs = dst.bgc_vegetation.cnveg_carbonstate_inst
        dvns = dst.bgc_vegetation.cnveg_nitrogenstate_inst
        @test dvcs.leafc_patch == vcs.leafc_patch
        @test dvcs.frootc_patch == vcs.frootc_patch
        @test dvcs.livestemc_patch == vcs.livestemc_patch
        @test dvcs.deadstemc_patch == vcs.deadstemc_patch
        @test dvcs.cpool_patch == vcs.cpool_patch
        @test dvns.leafn_patch == vns.leafn_patch
        @test dvns.frootn_patch == vns.frootn_patch
        @test dvns.npool_patch == vns.npool_patch

        # use_cn=false reader must NOT touch CN fields (they stay NaN-init).
        dst2 = build_inst()
        CLM.read_restart!(dst2, file; bounds=bounds, use_cn=false)
        @test all(isnan, dst2.bgc_vegetation.cnveg_carbonstate_inst.leafc_patch)
    end

    # ----------------------------------------------------------------------
    # CN flux/accumulator + phenology + crop + isotope state round-trips.
    # These need the CN/crop/isotope arrays actually allocated, so we flip
    # the facade config flags and re-run cn_vegetation_init! to materialize
    # the (otherwise length-0) pools before seeding.
    # ----------------------------------------------------------------------
    function build_cn_inst(; use_cn=true, use_crop=true, use_c13=false, use_c14=false)
        inst = build_inst()
        v = inst.bgc_vegetation
        v.config.use_cn   = use_cn
        v.config.use_crop = use_crop
        v.config.use_c13  = use_c13
        v.config.use_c14  = use_c14
        CLM.cn_vegetation_init!(v, np, nc, ng;
                                 nlevdecomp=10, ndecomp_pools=7,
                                 ndecomp_cascade_transitions=5)
        return inst
    end

    # --- CN flux/accumulator + phenology + crop round-trip ---
    mktempdir() do dir
        file = joinpath(dir, "clm_restart_cnflux.nc")
        src = build_cn_inst()

        vcs = src.bgc_vegetation.cnveg_carbonstate_inst
        vns = src.bgc_vegetation.cnveg_nitrogenstate_inst
        vss = src.bgc_vegetation.cnveg_state_inst
        cr  = src.crop

        # Carbon flux/storage/transfer pools
        fill_seq!(vcs.leafc_storage_patch, 3.0)
        fill_seq!(vcs.leafc_xfer_patch, 4.0)
        fill_seq!(vcs.livecrootc_patch, 6.0)
        fill_seq!(vcs.gresp_storage_patch, 2.0)
        fill_seq!(vcs.gresp_xfer_patch, 1.5)
        fill_seq!(vcs.xsmrpool_patch, 1.0)
        fill_seq!(vcs.cropseedc_deficit_patch, -0.5)
        fill_seq!(vcs.reproductivec_patch, 7.0)           # 2-D (patch × nrepr)
        # Nitrogen flux/storage/transfer pools
        fill_seq!(vns.leafn_storage_patch, 0.3)
        fill_seq!(vns.livestemn_patch, 0.6)
        fill_seq!(vns.livestemn_storage_patch, 0.5)
        fill_seq!(vns.retransn_patch, 0.4)
        fill_seq!(vns.reproductiven_patch, 0.7)           # 2-D (patch × nrepr)
        # Phenology counters
        fill_seq!(vss.dormant_flag_patch, 0.0)
        fill_seq!(vss.days_active_patch, 30.0)
        fill_seq!(vss.onset_counter_patch, 10.0)
        fill_seq!(vss.offset_flag_patch, 0.0)
        fill_seq!(vss.tempsum_potential_gpp_patch, 5.0)
        fill_seq!(vss.annsum_counter_col, 100.0)
        fill_seq!(vss.annavg_t2m_col, 285.0)
        # Crop block (CNVegState side)
        fill_seq!(vss.gddmaturity_patch, 1500.0)
        fill_seq!(vss.cumvd_patch, 12.0)
        vss.idop_patch    .= fill(120, np)
        vss.peaklai_patch .= fill(1, np)
        fill_seq!(vss.gddmaturity_thisyr, 50.0)           # 2-D (patch × mxharvests)
        # CropType state
        cr.croplive_patch .= [isodd(i) for i in 1:np]
        cr.sown_in_this_window .= [iseven(i) for i in 1:np]
        cr.harvdate_patch .= fill(200, np)
        cr.nyrs_crop_active_patch .= collect(1:np)
        fill_seq!(cr.vf_patch, 0.5)
        fill_seq!(cr.cphase_patch, 2.0)
        fill_seq!(cr.gddaccum_thisyr_patch, 9.0)          # 2-D (patch × mxharvests)
        fill_seq!(cr.hui_thisyr_patch, 8.0)

        CLM.write_restart(src, file; bounds=bounds, use_cn=true, use_crop=true)
        dst = build_cn_inst()
        CLM.read_restart!(dst, file; bounds=bounds, use_cn=true, use_crop=true)

        dvcs = dst.bgc_vegetation.cnveg_carbonstate_inst
        dvns = dst.bgc_vegetation.cnveg_nitrogenstate_inst
        dvss = dst.bgc_vegetation.cnveg_state_inst
        dcr  = dst.crop

        # Carbon flux pools
        @test dvcs.leafc_storage_patch == vcs.leafc_storage_patch
        @test dvcs.leafc_xfer_patch == vcs.leafc_xfer_patch
        @test dvcs.livecrootc_patch == vcs.livecrootc_patch
        @test dvcs.gresp_storage_patch == vcs.gresp_storage_patch
        @test dvcs.gresp_xfer_patch == vcs.gresp_xfer_patch
        @test dvcs.xsmrpool_patch == vcs.xsmrpool_patch
        @test dvcs.cropseedc_deficit_patch == vcs.cropseedc_deficit_patch
        @test dvcs.reproductivec_patch == vcs.reproductivec_patch
        # Nitrogen flux pools
        @test dvns.leafn_storage_patch == vns.leafn_storage_patch
        @test dvns.livestemn_patch == vns.livestemn_patch
        @test dvns.livestemn_storage_patch == vns.livestemn_storage_patch
        @test dvns.retransn_patch == vns.retransn_patch
        @test dvns.reproductiven_patch == vns.reproductiven_patch
        # Phenology counters
        @test dvss.dormant_flag_patch == vss.dormant_flag_patch
        @test dvss.days_active_patch == vss.days_active_patch
        @test dvss.onset_counter_patch == vss.onset_counter_patch
        @test dvss.offset_flag_patch == vss.offset_flag_patch
        @test dvss.tempsum_potential_gpp_patch == vss.tempsum_potential_gpp_patch
        @test dvss.annsum_counter_col == vss.annsum_counter_col
        @test dvss.annavg_t2m_col == vss.annavg_t2m_col
        # Crop block (CNVegState)
        @test dvss.gddmaturity_patch == vss.gddmaturity_patch
        @test dvss.cumvd_patch == vss.cumvd_patch
        @test dvss.idop_patch == vss.idop_patch
        @test eltype(dvss.idop_patch) <: Integer
        @test dvss.peaklai_patch == vss.peaklai_patch
        @test dvss.gddmaturity_thisyr == vss.gddmaturity_thisyr
        # CropType — integer + boolean survive the Float64 round-trip exactly
        @test dcr.croplive_patch == cr.croplive_patch
        @test eltype(dcr.croplive_patch) <: Bool
        @test dcr.sown_in_this_window == cr.sown_in_this_window
        @test dcr.harvdate_patch == cr.harvdate_patch
        @test eltype(dcr.harvdate_patch) <: Integer
        @test dcr.nyrs_crop_active_patch == cr.nyrs_crop_active_patch
        @test dcr.vf_patch == cr.vf_patch
        @test dcr.cphase_patch == cr.cphase_patch
        @test dcr.gddaccum_thisyr_patch == cr.gddaccum_thisyr_patch
        @test dcr.hui_thisyr_patch == cr.hui_thisyr_patch

        # use_crop=false reader must NOT touch crop fields (stay at cold-start).
        dst3 = build_cn_inst()
        # Cold-start croplive is all-false; planting date sentinel is typemax(Int).
        CLM.read_restart!(dst3, file; bounds=bounds, use_cn=true, use_crop=false)
        @test !any(dst3.crop.croplive_patch)            # untouched (cold default)
        @test dst3.crop.harvdate_patch != cr.harvdate_patch
        # ...but the CN flux pools (use_cn block) DID round-trip.
        @test dst3.bgc_vegetation.cnveg_carbonstate_inst.xsmrpool_patch ==
              vcs.xsmrpool_patch
    end

    # --- C13/C14 isotope pools round-trip, gated on use_c13/use_c14 ---
    mktempdir() do dir
        file = joinpath(dir, "clm_restart_iso.nc")
        src = build_cn_inst(use_c13=true, use_c14=true)
        c13 = src.bgc_vegetation.c13_cnveg_carbonstate_inst
        c14 = src.bgc_vegetation.c14_cnveg_carbonstate_inst
        fill_seq!(c13.leafc_patch, 11.0)
        fill_seq!(c13.leafc_storage_patch, 11.5)
        fill_seq!(c13.xsmrpool_patch, 12.0)
        fill_seq!(c13.reproductivec_patch, 1.1)
        fill_seq!(c14.leafc_patch, 13.0)
        fill_seq!(c14.xsmrpool_patch, 14.0)

        CLM.write_restart(src, file; bounds=bounds,
                          use_cn=true, use_crop=true, use_c13=true, use_c14=true)
        dst = build_cn_inst(use_c13=true, use_c14=true)
        CLM.read_restart!(dst, file; bounds=bounds,
                          use_cn=true, use_crop=true, use_c13=true, use_c14=true)

        d13 = dst.bgc_vegetation.c13_cnveg_carbonstate_inst
        d14 = dst.bgc_vegetation.c14_cnveg_carbonstate_inst
        @test d13.leafc_patch == c13.leafc_patch
        @test d13.leafc_storage_patch == c13.leafc_storage_patch
        @test d13.xsmrpool_patch == c13.xsmrpool_patch
        @test d13.reproductivec_patch == c13.reproductivec_patch
        @test d14.leafc_patch == c14.leafc_patch
        @test d14.xsmrpool_patch == c14.xsmrpool_patch

        # Isotope flags OFF on read must leave the isotope pools untouched
        # even though the suffixed vars are present in the file.
        dst4 = build_cn_inst(use_c13=true, use_c14=true)
        CLM.read_restart!(dst4, file; bounds=bounds,
                          use_cn=true, use_crop=true, use_c13=false, use_c14=false)
        @test all(isnan, dst4.bgc_vegetation.c13_cnveg_carbonstate_inst.leafc_patch)
        @test all(isnan, dst4.bgc_vegetation.c14_cnveg_carbonstate_inst.leafc_patch)
    end
end
