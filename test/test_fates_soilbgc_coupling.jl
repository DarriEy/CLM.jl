# test_fates_soilbgc_coupling.jl
# FATES -> CLM soil-biogeochemistry host coupling (use_fates_bgc).
#
#   A. fates_update_c_litter_fluxes! (UpdateCLitterFluxes, clmfates_interfaceMod.F90
#      :1367-1441): FATES litter fragmentation reaches decomp_cpools_sourcesink_col of
#      the metabolic / cellulose / lignin pools of FATES columns only, conserving mass,
#      and diagnoses fates_litter_flux; fates_update_n_litter_fluxes! zeroes the FATES
#      column's N source/sink (the ctsm5.3 UpdateNLitterFluxes).
#   B. Configuration: CLMDriverConfig(use_fates=true) carries soil BGC exactly like
#      CTSM controlMod (use_fates_bgc = .not. use_fates_sp) with the 6-pool /
#      8-transition FATES CENTURY cascade; the CN and SP configurations are unchanged.
#   C. fates_column_sites / fates_site_columns_from_flags (f2hmap hsites / fcolumn).
#   D. (data-gated) Live run on the fates_parity_1pt site: soil HR is positive, FATES
#      litter enters the litter pools after the first daily step, and the soil column
#      carbon closes step by step: d(TOTSOMC + TOTLITC) = (fates_litter_flux - HR) dt;
#      and HR / soil + litter C track the Fortran CTSM run of the same case
#      (scripts/validation/fates_fortran_parity/fates_soilbgc_fortran.csv, produced by
#      stage_soilbgc_run.sh; FATES_SOILBGC_FORTRAN_H0 may point at a fresh h0 .nc).

using Test
using CLM
using Dates

include("testdata.jl")

# Minimal PRT stub (FluxIntoLitterPools only reads PRT state on the MIMICS path).
mutable struct _SBCStubPRT <: CLM.AbstractPRTVartypes end
CLM.GetState(::_SBCStubPRT, organ_id::Integer, element_id::Integer,
             position_id::Union{Nothing,Integer}=nothing) = 0.0

@testset "FATES soil-BGC coupling" begin

    @testset "A. litter fluxes into the CLM decomposition pools" begin
        old_numpft   = CLM.numpft[]
        old_numel    = CLM.num_elements[]
        old_ellist   = copy(CLM.element_list)
        old_elpos    = copy(CLM.element_pos)
        old_decomp   = CLM.hlm_decomp[]
        old_edparams = CLM.EDParams[]
        old_edpft    = CLM.EDPftvarcon_inst[]
        try
            npft = 2; nlevsoil = 3; nlevdecomp = 3; ndp = 7
            CLM.numpft[] = npft
            CLM.num_elements[] = 1
            empty!(CLM.element_list); push!(CLM.element_list, CLM.carbon12_element)
            fill!(CLM.element_pos, 0); CLM.element_pos[CLM.carbon12_element] = 1
            CLM.hlm_decomp[] = "CENTURY"
            ep = CLM.ed_params_type(); ep.ED_val_cwd_fcel = 0.76; ep.ED_val_cwd_flig = 0.24
            CLM.EDParams[] = ep
            pc = CLM.EDPftvarcon_type()
            pc.lf_flab = [0.25, 0.25]; pc.lf_fcel = [0.50, 0.50]; pc.lf_flig = [0.25, 0.25]
            CLM.EDPftvarcon_inst[] = pc

            dz = [0.02, 0.04, 0.06]
            function mksite(scale)
                site = CLM.ed_site_type()
                patch = CLM.fates_patch_type()
                patch.area = CLM.area
                patch.nitr_repro_stoich = zeros(CLM.maxpft)
                litt = CLM.litter_type()
                CLM.InitAllocate!(litt, npft, nlevsoil, CLM.carbon12_element)
                CLM.ZeroFlux!(litt)
                CLM.InitConditions!(litt, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
                litt.ag_cwd_frag .= scale .* [0.01, 0.02, 0.0, 0.0]
                litt.bg_cwd_frag .= 0.0; litt.bg_cwd_frag[1, 2] = scale * 0.003
                litt.leaf_fines_frag .= scale .* [0.02, 0.03, 0.01]
                litt.root_fines_frag .= 0.0
                litt.root_fines_frag[CLM.ilabile, 1] = scale * 0.004
                litt.root_fines_frag[CLM.icellulose, 2] = scale * 0.006
                litt.root_fines_frag[CLM.ilignin, 3] = scale * 0.002
                litt.seed_decay .= scale .* [0.01, 0.0]
                litt.seed_germ_decay .= scale .* [0.0, 0.005]
                patch.litter = [litt]
                site.oldest_patch = patch
                bc_in = CLM.bc_in_type()
                CLM.allocate_bcin!(bc_in; npatches=1, nlevsoil=nlevsoil,
                                   nlevdecomp=nlevdecomp, max_comp=4)
                bc_in.nlevsoil = nlevsoil; bc_in.nlevdecomp = nlevdecomp
                bc_in.dz_sisl = copy(dz); bc_in.dz_decomp_sisl = copy(dz)
                bc_in.decomp_id = [1, 2, 3]
                bc_in.max_rooting_depth_index_col = nlevsoil
                bc_out = CLM.bc_out_type()
                CLM.allocate_bcout!(bc_out; npatches=1, nlevsoil=nlevsoil,
                                    nlevdecomp=nlevdecomp, max_comp=4)
                total = sum(litt.ag_cwd_frag) + sum(litt.bg_cwd_frag) +
                        sum(litt.leaf_fines_frag) + sum(litt.root_fines_frag) +
                        sum(litt.seed_decay) + sum(litt.seed_germ_decay)
                return site, bc_in, bc_out, total
            end
            s1, bi1, bo1, tot1 = mksite(1.0)
            s2, bi2, bo2, tot2 = mksite(2.0)
            fates = (nsites = 2, sites = [s1, s2], bc_in = [bi1, bi2], bc_out = [bo1, bo2])

            # 4 columns: 1 = non-FATES, 2 = FATES site 1, 3 = FATES site 2 (masked
            # out of the bgc filter this step), 4 = FATES but beyond nsites.
            nc = 4
            col_is_fates = [false, true, true, false]
            mask = BitVector([true, true, false, true])
            dt = 3600.0
            sentinel = 7.0
            ss = fill(sentinel, nc, nlevdecomp, ndp)
            cf = (decomp_cpools_sourcesink_col = ss, fates_litter_flux = fill(NaN, nc))
            # Fortran adds onto the SetValues-zeroed array; start the FATES column at 0
            # except for one cell that proves the update is additive.
            ss[2, :, :] .= 0.0
            ss[2, 1, 4] = 0.5

            CLM.fates_update_c_litter_fluxes!(fates, cf, col_is_fates, mask;
                dtime=dt, nlevdecomp=nlevdecomp, decomp_method=1, i_met_lit=1)

            # untouched: non-FATES column, masked-out FATES column
            @test all(ss[1, :, :] .== sentinel)
            @test all(ss[3, :, :] .== sentinel)
            @test all(ss[4, :, :] .== sentinel)
            @test isnan(cf.fates_litter_flux[1]) && isnan(cf.fates_litter_flux[3])

            # pool mapping: labile -> met (1), cellulose -> cel (2), lignin -> lig (3)
            for j in 1:nlevdecomp
                @test ss[2, j, 1] == bo1.litt_flux_lab_c_si[j] * dt
                @test ss[2, j, 2] == bo1.litt_flux_cel_c_si[j] * dt
                @test ss[2, j, 3] == bo1.litt_flux_lig_c_si[j] * dt
            end
            @test ss[2, 1, 4] == 0.5                       # additive, SOM untouched
            @test all(ss[2, :, 5:ndp] .== 0.0)
            # mass conservation: kg/m2/day of fragmentation == depth-integrated input
            conv = CLM.days_per_sec * CLM.g_per_kg       # kg/m2/day -> g/m2/s
            col_in = sum(ss[2, j, i] * dz[j] for j in 1:nlevdecomp, i in 1:3) / dt
            @test col_in ≈ tot1 * conv rtol = 1e-12
            @test cf.fates_litter_flux[2] ≈ tot1 * conv rtol = 1e-12
            # the masked-out site was never asked for its litter
            @test all(bo2.litt_flux_lab_c_si .== 0.0)

            # MIMICS: lignin shares the structural (cellulose) pool
            ss2 = zeros(nc, nlevdecomp, ndp)
            cf2 = (decomp_cpools_sourcesink_col = ss2, fates_litter_flux = zeros(nc))
            CLM.fates_update_c_litter_fluxes!(fates, cf2, col_is_fates, trues(nc);
                dtime=dt, nlevdecomp=nlevdecomp, decomp_method=CLM.MIMICS_DECOMP)
            for j in 1:nlevdecomp
                @test ss2[2, j, 2] ≈ (bo1.litt_flux_cel_c_si[j] + bo1.litt_flux_lig_c_si[j]) * dt
                @test ss2[2, j, 3] == 0.0
            end
            @test cf2.fates_litter_flux[3] ≈ tot2 * conv rtol = 1e-12

            # N: the ctsm5.3 UpdateNLitterFluxes zeroes the FATES column's N source/sink
            nss = fill(3.0, nc, nlevdecomp, ndp)
            nf = (decomp_npools_sourcesink_col = nss, fates_litter_flux = fill(NaN, nc))
            CLM.fates_update_n_litter_fluxes!(fates, nf, col_is_fates, mask)
            @test all(nss[2, :, :] .== 0.0)
            @test all(nss[1, :, :] .== 3.0) && all(nss[3, :, :] .== 3.0)
            @test all(nf.fates_litter_flux .== 0.0)

            # no FATES instance => no-op
            ss3 = fill(1.0, nc, nlevdecomp, ndp)
            CLM.fates_update_c_litter_fluxes!(nothing,
                (decomp_cpools_sourcesink_col = ss3, fates_litter_flux = zeros(nc)),
                col_is_fates, mask; dtime=dt, nlevdecomp=nlevdecomp, decomp_method=1)
            @test all(ss3 .== 1.0)
        finally
            CLM.numpft[] = old_numpft
            CLM.num_elements[] = old_numel
            empty!(CLM.element_list); append!(CLM.element_list, old_ellist)
            copyto!(CLM.element_pos, old_elpos)
            CLM.hlm_decomp[] = old_decomp
            CLM.EDParams[] = old_edparams
            CLM.EDPftvarcon_inst[] = old_edpft
        end
    end

    @testset "B. driver configuration" begin
        f = CLM.CLMDriverConfig(use_fates=true)
        @test f.use_fates && f.use_fates_bgc && !f.use_fates_sp
        @test f.ndecomp_pools == 6
        @test f.ndecomp_cascade_transitions == 8
        fsp = CLM.CLMDriverConfig(use_fates=true, use_fates_sp=true)
        @test fsp.use_fates_sp && !fsp.use_fates_bgc
        @test fsp.ndecomp_pools == 7 && fsp.ndecomp_cascade_transitions == 10
        fno = CLM.CLMDriverConfig(use_fates=true, use_fates_bgc=false)
        @test !fno.use_fates_bgc
        cn = CLM.CLMDriverConfig(use_cn=true)
        @test !cn.use_fates_bgc && cn.ndecomp_pools == 7 && cn.ndecomp_cascade_transitions == 10
        sp = CLM.CLMDriverConfig()
        @test !sp.use_fates && !sp.use_fates_bgc && sp.ndecomp_pools == 7
    end

    @testset "C. site/column maps" begin
        flags = [false, true, false, true, true]
        @test CLM.fates_column_sites(flags, 2) == [0, 1, 0, 2, 0]
        @test CLM.fates_site_columns_from_flags(flags, 2) == [2, 4]
        @test CLM.fates_site_columns_from_flags(flags, 4) == [2, 4, 5, 0]
    end

    @testset "D. live FATES-BGC run (data-gated)" begin
        sitedir = symfluence_path("domain_Bow_at_Banff_lumped_era5", "settings", "CLM",
                                  "parameters")
        fsurdat = joinpath(sitedir, "surfdata_clm.nc")
        fparam  = joinpath(sitedir, "clm5_params.nc")
        forcdir = symfluence_path("domain_Bow_at_Banff_lumped_era5", "data", "forcing",
                                  "CLM_input")
        forcing = joinpath(forcdir, "clmforc.2002.nc")
        if !(isfile(fsurdat) && isfile(fparam) && isfile(forcing))
            testdata_missing("FATES soil-BGC live run", fsurdat, fparam, forcing)
            return
        end
        harness = joinpath(@__DIR__, "..", "scripts", "validation", "fates_fortran_parity",
                           "fates_soilbgc_parity.jl")
        # The harness (and fates_fortran_parity.jl it includes) read the site from
        # SYMFLUENCE_DATA; point it at the same root the test data resolved to.
        ENV["SYMFLUENCE_DATA"] = symfluence_data_root()
        m = Module(:_FatesSoilBGCHarness)
        Core.eval(m, :(include(p::AbstractString) = Base.include($m, p)))
        Base.include(m, harness)
        nsteps = 36                      # crosses the first daily step (nstep 25)
        recs = Base.invokelatest(getglobal(m, :run_julia), nsteps)
        hr   = [r["HR"] for r in recs]
        flit = [r["FLIT"] for r in recs]
        soilc = [r["TOTSOMC"] + r["TOTLITC"] for r in recs]
        @test all(isfinite, hr) && all(>=(0), hr)
        @test maximum(hr) > 1e-8                              # decomposition is live
        @test recs[end]["FATES_HET_RESP"] ≈ hr[end] * 1e-3 rtol = 1e-12   # g -> kg
        @test recs[end]["TOTLITC"] > 0                         # FATES litter arrived
        @test flit[end] > 0
        @test recs[end]["LITTERC_HR"] > 0                      # and is decomposing
        # Soil column C closes every step (no leaching of SOM C at this site/window):
        #   d(TOTSOMC + TOTLITC) = (fates_litter_flux - HR) * dt
        dtime = 3600.0
        err = [abs((soilc[k] - soilc[k-1]) - (flit[k] - hr[k]) * dtime) for k in 3:length(recs)]
        @test maximum(err) < 1e-9

        # Fortran reference: the committed hourly extract of the fates_parity_1pt case
        # (override with an h0 .nc via FATES_SOILBGC_FORTRAN_H0).
        ref = get(ENV, "FATES_SOILBGC_FORTRAN_H0",
                  getglobal(m, :FORTRAN_CSV))
        if isfile(ref)
            ft = Base.invokelatest(getglobal(m, :read_fortran), ref)
            n = min(length(recs), length(ft["HR"]))
            # Post-thaw (first day boundary onward) the soil-BGC thread tracks Fortran
            # to within the host soil-temperature/moisture residual.
            for k in 26:n
                @test abs(hr[k] - ft["HR"][k]) <= 0.10 * ft["HR"][k]
                @test abs(soilc[k] - (ft["TOTSOMC"][k] + ft["TOTLITC"][k])) < 0.05
            end
            @test abs(recs[n]["TOTLITC"] - ft["TOTLITC"][n]) <= 0.05 * ft["TOTLITC"][n]
        end
    end
end
