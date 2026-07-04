# ==========================================================================
# test_fortran_parity.jl — gated per-timestep Julia↔Fortran parity regression.
#
# Runs the single-step shared-IC harness for a daytime and a night step on the
# Bow column and asserts the post-step prognostic state matches the Fortran
# reference within per-field tolerances. The tolerances are set to PASS the
# current code yet FAIL if any of the localized parity fixes regress:
#   - canopy radiative-transfer injection (SABV would collapse to 0)
#   - canopy_fluxes_params a_coef/a_exp wiring (night T_VEG would be ~2 K off)
#   - frost_table use_aquifer_layer gate (ZWT_PERCH would be ~0.27 m off)
#   - T_STEM injection (would be ~20 K off)
#   - Bow daytime coupled canopy leaf temperature (machine-level T_VEG parity)
#   - the cold-start pedotransfer / forcing / canopy radiation fixes (T_GRND
#     over-heating, ~16 K).
#
# GATED: the Fortran reference dumps + Bow input files live outside the repo
# (machine-specific). When absent, the test is skipped (not failed) so CI and
# other machines stay green.
# ==========================================================================

using Test

@testset "Fortran parity (per-step, gated)" begin
    common = joinpath(@__DIR__, "..", "scripts", "fortran_parity_common.jl")
    dumpdir = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_parity_run"

    # Steps to check + the field tolerances (abs). Daytime (n13461, peak sun)
    # exercises the radiation/canopy path; night (n13470) exercises the
    # under-canopy resistance (a_coef) fix.
    steps = (13461, 13470)
    have_dumps = isdir(dumpdir) && all(
        isfile(joinpath(dumpdir, "pdump_before_step_n$(n).nc")) &&
        isfile(joinpath(dumpdir, "pdump_after_hydrologydrainage_n$(n).nc"))
        for n in steps)

    if !(isfile(common) && have_dumps)
        @info "Fortran parity test SKIPPED (reference dumps not present on this machine)"
        @test_skip have_dumps
    else
        Base.include(@__MODULE__, common)   # build_bow_inst / run_one_parity_step! / DUMPDIR

        # (dump_name, kind, abs_tol, inst->array)
        rd(v) = ismissing(v) ? NaN : Float64(v)
        fields(inst, nstep) = begin
            ws = inst.water.waterstatebulk_inst.ws
            wd = inst.water.waterdiagnosticbulk_inst
            te = inst.temperature; sh = inst.soilhydrology; sa = inst.solarabs
            ef = inst.energyflux
            [
                ("T_GRND",     :col1d, 0.20,  te.t_grnd_col),
                ("T_SOISNO",   :col2d, 0.20,  te.t_soisno_col),
                ("T_VEG",      :patch, nstep == 13461 ? 1e-8 : 1.20, te.t_veg_patch),
                ("T_STEM",     :patch, 0.50,  te.t_stem_patch),
                ("SABV_P",     :patch, 5.0,   sa.sabv_patch),
                ("SABG_P",     :patch, 5.0,   sa.sabg_patch),
                ("EFLX_GNET_P",:patch, 6.0,   ef.eflx_gnet_patch),
                ("ZWT",        :col1d, 0.02,  sh.zwt_col),
                ("ZWT_PERCH",  :col1d, 0.05,  sh.zwt_perched_col),
                ("H2OSOI_LIQ", :col2d, 0.05,  ws.h2osoi_liq_col),
                ("H2OSOI_ICE", :col2d, 0.05,  ws.h2osoi_ice_col),
                ("WA",         :col1d, 1.0,   ws.wa_col),
                ("H2OSFC",     :col1d, 1e-3,  ws.h2osfc_col),
                ("SNOW_DEPTH", :col1d, 1e-3,  wd.snow_depth_col),
                ("frac_sno",   :col1d, 1e-3,  wd.frac_sno_col),
            ]
        end

        maxabs(ds, name, kind, jl) = begin
            haskey(ds, name) || return NaN
            fv = Float64[]; jv = Float64[]
            if kind == :col1d
                push!(fv, rd(ds[name][1])); push!(jv, Float64(jl[1]))
            elseif kind == :col2d
                d = ds[name][:, :]; nlev = min(size(d, 1), size(jl, 2))
                for k in 1:nlev; push!(fv, rd(d[k, 1])); push!(jv, Float64(jl[1, k])); end
            elseif kind == :patch
                d = ds[name][:]; n = min(length(d), length(jl))
                for i in 1:n; push!(fv, rd(d[i])); push!(jv, Float64(jl[i])); end
            end
            keep = .!(isnan.(fv) .| isnan.(jv))
            fv = fv[keep]; jv = jv[keep]
            isempty(fv) ? NaN : maximum(abs.(fv .- jv))
        end

        for n in steps
            @testset "step n$n" begin
                inst, _ = run_one_parity_step!(n; use_hydrstress=true)
                ds = NCDataset(joinpath(dumpdir, "pdump_after_hydrologydrainage_n$(n).nc"), "r")
                for (name, kind, tol, jl) in fields(inst, n)
                    m = maxabs(ds, name, kind, jl)
                    isnan(m) && continue   # field absent in this dump
                    @test m <= tol
                end
                close(ds)
            end
        end
    end
end
