@testset "init_interp (cross-grid restart start)" begin
    ng, nl, nc, np = 2, 3, 6, 12

    # Build a fully-allocated instance at the CURRENT varpar level structure.
    function build_inst()
        inst = CLM.CLMInstances()
        CLM.clm_instInit!(inst; ng=ng, nl=nl, nc=nc, np=np)
        return inst
    end
    bounds = CLM.BoundsType(begc=1, endc=nc, begp=1, endp=np,
                            begg=1, endg=ng, begl=1, endl=nl)

    fill_seq!(A, base) = (A .= reshape(base .+ (1:length(A)), size(A)); A)

    # Save the global soil-layer structure so we leave it untouched for the
    # rest of the suite (these tests deliberately mutate it).
    _orig_soilstruct = CLM.varctl.soil_layerstruct_predefined

    # Set a monotone-increasing z coordinate for every column over the soil
    # levels (the snow levels above stay NaN, treated as missing). Returns the
    # number of vertical levels in the z matrix.
    function set_monotone_z!(inst)
        nlev = size(inst.column.z, 2)
        for c in 1:size(inst.column.z, 1)
            for k in 1:nlev
                inst.column.z[c, k] = 0.05 * k          # 0.05, 0.10, ... (m)
            end
        end
        return nlev
    end

    # Seed the prognostic state with distinct finite values + a known
    # temperature profile that is a linear function of depth z (so vertical
    # interpolation has an analytic, bracketed answer).
    function seed_state!(inst)
        T = inst.temperature
        fill_seq!(T.t_grnd_col, 270.0)
        fill_seq!(T.t_h2osfc_col, 274.0)
        fill_seq!(T.t_veg_patch, 290.0)
        fill_seq!(T.t_ref2m_patch, 288.0)
        # t_soisno: linear in z with a per-column offset so the interpolated dest
        # equals an analytic line that DIFFERS per column (catches any collapse
        # of columns onto a single source profile).
        nlev = size(T.t_soisno_col, 2)
        for c in 1:size(T.t_soisno_col, 1)
            for k in 1:nlev
                z = inst.column.z[c, k]
                T.t_soisno_col[c, k] = 250.0 + 5.0 * c + 100.0 * z  # linear T(z), per-col offset
            end
        end

        ws = inst.water.waterstatebulk_inst.ws
        fill_seq!(ws.h2osno_no_layers_col, 1.0)
        fill_seq!(ws.h2osfc_col, 0.5)
        fill_seq!(ws.wa_col, 4000.0)
        fill_seq!(ws.liqcan_patch, 0.1)
        fill_seq!(ws.snocan_patch, 0.2)
        fill_seq!(inst.water.waterstatebulk_inst.int_snow_col, 2.0)
        fill_seq!(inst.soilhydrology.zwt_col, 3.0)
        fill_seq!(inst.canopystate.tlai_patch, 1.0)
        return inst
    end

    # --- (a) identical grid → must match read_restart! exactly --------------
    @testset "identical grid reduces to read_restart!" begin
        CLM.varctl.soil_layerstruct_predefined = "10SL_3.5m"
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
        mktempdir() do dir
            file = joinpath(dir, "finidat.nc")
            src = build_inst()
            set_monotone_z!(src)
            seed_state!(src)
            CLM.write_restart(src, file; bounds=bounds)

            dst_ii = build_inst()
            CLM.init_interp!(dst_ii, file; bounds=bounds)

            dst_rr = build_inst()
            CLM.read_restart!(dst_rr, file; bounds=bounds)

            # init_interp! on a matching grid is byte-identical to read_restart!
            @test isequal(dst_ii.temperature.t_soisno_col,
                          dst_rr.temperature.t_soisno_col)
            @test isequal(dst_ii.temperature.t_grnd_col, dst_rr.temperature.t_grnd_col)
            @test isequal(dst_ii.temperature.t_veg_patch, dst_rr.temperature.t_veg_patch)
            @test isequal(dst_ii.soilhydrology.zwt_col, dst_rr.soilhydrology.zwt_col)
            @test isequal(dst_ii.water.waterstatebulk_inst.ws.wa_col,
                          dst_rr.water.waterstatebulk_inst.ws.wa_col)
            @test isequal(dst_ii.canopystate.tlai_patch, dst_rr.canopystate.tlai_patch)
        end
    end

    # --- (b) different #soil levels → vertical interpolation ----------------
    @testset "different soil levels → vertical interp" begin
        # Source grid: 4 soil layers (nlevgrnd=4).
        CLM.varctl.soil_layerstruct_predefined = "4SL_2m"
        CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
        src = build_inst()
        nlev_s = set_monotone_z!(src)
        seed_state!(src)

        # Capture source soil profile + z for the analytic check before we
        # change the global level structure.
        src_t = copy(src.temperature.t_soisno_col)
        src_z = copy(src.column.z)

        mktempdir() do dir
            file = joinpath(dir, "finidat_coarse.nc")
            CLM.write_restart(src, file; bounds=bounds)

            # Destination grid: 10 soil layers (nlevgrnd=15) → finer vertical.
            CLM.varctl.soil_layerstruct_predefined = "10SL_3.5m"
            CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
            dst = build_inst()
            nlev_d = set_monotone_z!(dst)
            @test nlev_d != nlev_s   # grids genuinely differ vertically

            CLM.init_interp!(dst, file; bounds=bounds)

            dt = dst.temperature.t_soisno_col
            dz = dst.column.z

            # Finite everywhere a destination z level exists.
            @test all(isfinite, dt)

            # Source min/max over real (z-defined) levels, per column.
            for c in 1:nc
                # source valid levels (z finite — all here)
                smin = minimum(src_t[c, :])
                smax = maximum(src_t[c, :])
                for k in 1:nlev_d
                    # bracketed by source range (clamp at the ends)
                    @test dt[c, k] >= smin - 1e-9
                    @test dt[c, k] <= smax + 1e-9
                end
                # monotone-preserving: source T(z) increasing → dest increasing
                for k in 2:nlev_d
                    @test dt[c, k] >= dt[c, k-1] - 1e-9
                end
                # analytic: T was linear (250 + 100*z); within the source z
                # range the interpolant must reproduce the line exactly.
                zs_lo = minimum(src_z[c, :]); zs_hi = maximum(src_z[c, :])
                for k in 1:nlev_d
                    z = dz[c, k]
                    if zs_lo <= z <= zs_hi
                        @test isapprox(dt[c, k], 250.0 + 5.0 * c + 100.0 * z; atol=1e-6)
                    end
                end
            end
        end

    end

    # Restore the global soil-layer structure for the rest of the suite.
    CLM.varctl.soil_layerstruct_predefined = _orig_soilstruct
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
end
