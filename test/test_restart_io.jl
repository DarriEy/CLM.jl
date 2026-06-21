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
end
