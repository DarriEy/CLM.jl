# Water-table reverse phase: the CPU HOST phase list (watertable_rev_phase!), the device
# decomposition (watertable_dev_rev_phase!, run on the host here) and production
# clm_drv! must take the same branch — aquifer water_table! with use_aquifer_layer,
# PerchedWaterTable + ThetaBasedWaterTable without it, then the bedrock clip — so they
# give the same forward state and the same Enzyme gradients through zwt / qflx_drain.
# (The host phase used to call the aquifer water_table! unconditionally, which made
# host-list gradients through zwt / qflx_drain differ from production when
# use_aquifer_layer = false, the CLM5 default.)

using Test, Random
import CLM
include(joinpath(@__DIR__, "testdata.jl"))

@testset "water-table reverse phase: host == device == production" begin
    enzyme_ok = try
        @eval using Enzyme
        true
    catch e
        @info "Enzyme.jl not available — skipping water-table reverse test" exception=e
        false
    end
    if !enzyme_ok
        @test_skip true
        return
    end
    fsurdat, paramfile = bow_params()
    if !isfile(fsurdat) || !isfile(paramfile)
        testdata_missing("water-table reverse phase", fsurdat, paramfile) && return
    end

    (inst, bounds, filt, tm) = CLM.clm_initialize!(; fsurdat=fsurdat, paramfile=paramfile)
    for g in 1:bounds.endg
        a = inst.atm2lnd; a.forc_t_not_downscaled_grc[g] = 285.0
        a.forc_pbot_not_downscaled_grc[g] = 85000.0
        a.forc_th_not_downscaled_grc[g] = 285.0 * (100000.0/85000.0)^(CLM.RAIR/CLM.CPAIR)
        a.forc_rho_not_downscaled_grc[g] = 85000.0/(CLM.RAIR*285.0)
        a.forc_lwrad_not_downscaled_grc[g] = 300.0; a.forc_vp_grc[g] = 800.0
        a.forc_hgt_grc[g] = 30.0; a.forc_wind_grc[g] = 3.0
        a.forc_hgt_u_grc[g] = 30.0; a.forc_hgt_t_grc[g] = 30.0; a.forc_hgt_q_grc[g] = 30.0
        a.forc_topo_grc[g] = 0.0
        for b in 1:CLM.NUMRAD; a.forc_solad_not_downscaled_grc[g,b] = 200.0; a.forc_solai_grc[g,b] = 80.0; end
        a.forc_solar_not_downscaled_grc[g] = 560.0
        a.forc_rain_not_downscaled_grc[g] = 0.0001; a.forc_snow_not_downscaled_grc[g] = 0.0
    end
    CLM.downscale_forcings!(bounds, inst.atm2lnd, inst.column, inst.landunit, inst.topo)
    config0 = CLM.CLMDriverConfig(); filt_ia = CLM.clump_filter_inactive_and_active
    (declin, _) = CLM.compute_orbital(120.0); nextsw = 120.0 + 1800.0/CLM.SECSPDAY
    for n in 1:3
        CLM.clm_drv!(config0, inst, filt, filt_ia, bounds, true, nextsw, declin, declin,
            CLM.ORB_OBLIQR_DEFAULT, false, false, "", false;
            nstep=n, is_first_step=(n==1), is_beg_curr_day=(n==1), dtime=1800.0,
            mon=1, day=1, photosyns=inst.photosyns)
    end
    bc_col = bounds.begc:bounds.endc; nlevsoi = CLM.varpar.nlevsoi
    hc = filt.hydrologyc; c0 = first(c for c in bc_col if hc[c])

    outs(i) = (i.soilhydrology.zwt_col, i.soilhydrology.zwt_perched_col,
               i.soilhydrology.frost_table_col,
               i.water.waterfluxbulk_inst.wf.qflx_drain_col,
               i.water.waterfluxbulk_inst.wf.qflx_drain_perched_col,
               i.water.waterfluxbulk_inst.wf.qflx_rsub_sat_col,
               i.water.waterstatebulk_inst.ws.wa_col)

    # production clm_drv! water-table block (clm_driver.jl, HydrologyNoDrainage)
    function production!(i, aquifer)
        sh = i.soilhydrology; ss = i.soilstate; col = i.column
        if aquifer
            CLM.water_table!(sh, ss, i.temperature.t_soisno_col, i.water.waterstatebulk_inst.ws,
                i.water.waterfluxbulk_inst, col.dz, col.z, col.zi, hc, bc_col, nlevsoi, 1800.0;
                recompute_frost_table = true)
        else
            CLM.perched_water_table!(sh, ss, i.temperature.t_soisno_col,
                i.water.waterstatebulk_inst.ws, col.dz, col.z, col.zi, hc, bc_col, nlevsoi)
            CLM.theta_based_water_table!(sh, ss, i.water.waterstatebulk_inst.ws,
                col.dz, col.z, col.zi, col.nbedrock, hc, bc_col, nlevsoi)
        end
        CLM.clamp_zwt_to_bedrock!(sh.zwt_col, hc, col.nbedrock, col.zi, CLM.varpar.nlevsno)
    end

    # leaves of the gradient bundle that are inputs of the water-table phase
    grad_leaves(db) = (db.inst.soilhydrology, db.inst.soilstate,
        db.inst.water.waterstatebulk_inst.ws, db.inst.water.waterfluxbulk_inst.wf,
        db.inst.temperature, db.inst.column)
    function float_fields(x)
        out = Pair{String,Any}[]
        for f in fieldnames(typeof(x))
            v = getfield(x, f)
            v isa AbstractArray{<:AbstractFloat} && push!(out, string(f) => v)
        end
        out
    end

    Enzyme.API.strictAliasing!(false)
    revmode = Enzyme.set_runtime_activity(Enzyme.Reverse)
    rng = Random.MersenneTwister(11)

    for aquifer in (false, true)
        config = CLM.CLMDriverConfig(use_aquifer_layer = aquifer)
        aux = CLM.watertable_rev_aux(bounds, filt, config)

        # ---- forward: host phase == device phase == production
        bh = CLM.driver_rev_bundle(deepcopy(inst)); CLM.watertable_rev_phase!(bh, aux)
        bd = CLM.driver_rev_bundle(deepcopy(inst))
        CLM.watertable_dev_rev_phase!(CLM.watertable_rev_view(bd), aux)
        ip = deepcopy(inst); production!(ip, aquifer)
        for (h, d, p) in zip(outs(bh.inst), outs(bd.inst), outs(ip))
            @test isequal(h, p)
            @test isequal(d, p)
        end
        @test isfinite(bh.inst.soilhydrology.zwt_col[c0])

        # ---- reverse (Enzyme, CPU): host phase == device phase, seeds on zwt,
        # zwt_perched and qflx_drain
        seeds = Dict(k => 0.5 .+ rand(rng, size(a)...) for (k, a) in
            (("zwt", inst.soilhydrology.zwt_col), ("zwtp", inst.soilhydrology.zwt_perched_col),
             ("drain", inst.water.waterfluxbulk_inst.wf.qflx_drain_col)))
        function seed!(db)
            copyto!(db.inst.soilhydrology.zwt_col, seeds["zwt"])
            copyto!(db.inst.soilhydrology.zwt_perched_col, seeds["zwtp"])
            copyto!(db.inst.water.waterfluxbulk_inst.wf.qflx_drain_col, seeds["drain"])
        end
        b1 = CLM.driver_rev_bundle(deepcopy(inst)); d1 = Enzyme.make_zero(b1); seed!(d1)
        CLM._phase_autodiff!(revmode, CLM.watertable_rev_phase!, b1, d1, (aux,))
        b2 = CLM.driver_rev_bundle(deepcopy(inst)); d2 = Enzyme.make_zero(b2); seed!(d2)
        CLM._phase_autodiff!(revmode, CLM.watertable_dev_rev_phase!,
            CLM.watertable_rev_view(b2), CLM.watertable_rev_view(d2), (aux,))
        nnz = 0
        for (x1, x2) in zip(grad_leaves(d1), grad_leaves(d2))
            for ((f, g1), (_, g2)) in zip(float_fields(x1), float_fields(x2))
                m = isfinite.(g1) .& isfinite.(g2)
                @test isequal(isfinite.(g1), isfinite.(g2))
                sc = max(maximum(abs, g1[m]; init = 0.0), 1e-300)
                @test maximum(abs.(g1[m] .- g2[m]); init = 0.0) <= 1e-12 * sc
                nnz += count(!iszero, g1[m])
            end
        end
        @test nnz > 0

        # ---- host gradient vs central FD of L = seeds · (zwt, zwt_perched, qflx_drain)
        # w.r.t. the soil-liquid entry with the largest |gradient| in column c0
        gl = d1.inst.water.waterstatebulk_inst.ws.h2osoi_liq_col
        j = argmax(abs.(gl[c0, :])); g = gl[c0, j]
        x0 = inst.water.waterstatebulk_inst.ws.h2osoi_liq_col[c0, j]
        hstep = 1e-6 * max(abs(x0), 1.0)
        function lossd(δ)
            i = deepcopy(inst); i.water.waterstatebulk_inst.ws.h2osoi_liq_col[c0, j] += δ
            b = CLM.driver_rev_bundle(i); CLM.watertable_rev_phase!(b, aux)
            (b.inst.soilhydrology.zwt_col, b.inst.soilhydrology.zwt_perched_col,
             b.inst.water.waterfluxbulk_inst.wf.qflx_drain_col)
        end
        op = lossd(hstep); om = lossd(-hstep)
        fd = 0.0
        for (k, (a, b)) in zip(("zwt", "zwtp", "drain"), zip(op, om))
            d = a .- b; m = isfinite.(d); fd += sum(seeds[k][m] .* d[m])
        end
        fd /= 2hstep
        @test abs(g - fd) <= 1e-6 * max(abs(fd), 1e-8)
    end
end
