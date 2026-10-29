# test_fates_niche_modes.jl
# The FATES research modes (scripts/validation/fates_fortran_parity/README.md "RESEARCH MODES"): namelist validity checks, the
# switches reaching the FATES control globals, host seed-dispersal wiring, inventory
# initialization through clm_fates_init!, and the FATES restart file round-trip
# (write -> rebuild via create_patchcohort_structure! -> identical state).

using Test
using CLM

# clm_fates_init! rewrites the FATES module globals (control Refs, parameter tables,
# dimension counts). Snapshot every module-level Ref (+ element_list) and restore it.
function _niche_snapshot()
    snap = Dict{Symbol,Any}()
    for n in names(CLM; all = true)
        isdefined(CLM, n) || continue
        v = getfield(CLM, n)
        v isa Base.RefValue && (snap[n] = deepcopy(v[]))
    end
    snap[:__element_list] = copy(CLM.element_list)
    return snap
end
function _niche_restore!(snap)
    for (n, v) in snap
        n === :__element_list && continue
        getfield(CLM, n)[] = v
    end
    empty!(CLM.element_list); append!(CLM.element_list, snap[:__element_list])
    return nothing
end

function _niche_cohorts(site)
    out = NamedTuple[]
    p = site.oldest_patch
    while p !== nothing
        c = p.tallest
        while c !== nothing
            push!(out, (; patch_age = p.age, patch_area = p.area, pft = c.pft, n = c.n,
                          dbh = c.dbh, height = c.height, layer = c.canopy_layer,
                          coage = c.coage, crowndamage = c.crowndamage,
                          leafc = CLM.GetState(c.prt, CLM.leaf_organ, CLM.carbon12_element),
                          structc = CLM.GetState(c.prt, CLM.struct_organ, CLM.carbon12_element),
                          storec = CLM.GetState(c.prt, CLM.store_organ, CLM.carbon12_element)))
            c = c.shorter
        end
        p = p.younger
    end
    return out
end

@testset "FATES research modes" begin
    snap = _niche_snapshot()
    try
        @testset "host-side validity checks (CTSM build-namelist + check_allset)" begin
            chk = CLM.fates_check_niche_modes
            @test chk() === nothing
            @test chk(; use_fates_tree_damage = true, use_fates_cohort_age_tracking = true,
                        fates_seeddisp_cadence = 3) === nothing
            @test_throws ErrorException chk(; fates_seeddisp_cadence = 4)
            @test_throws ErrorException chk(; fates_seeddisp_cadence = -1)
            @test_throws ErrorException chk(; use_fates_inventory_init = true)
            @test_throws ErrorException chk(; use_fates_inventory_init = true,
                                              fates_inventory_ctrl_filename = "/no/such/file")
            mktemp() do path, io
                close(io)
                @test chk(; use_fates_inventory_init = true,
                            fates_inventory_ctrl_filename = path) === nothing
                @test_throws ErrorException chk(; use_fates_inventory_init = true,
                                                  fates_inventory_ctrl_filename = path,
                                                  use_fates_cohort_age_tracking = true)
            end
            @test_throws ErrorException chk(; use_fates_ed_st3 = true,
                                              use_fates_ed_prescribed_phys = true)
            @test_throws ErrorException chk(; use_fates_tree_damage = true,
                                              fates_parteh_mode = CLM.prt_cnp_flex_allom_hyp)
        end

        @testset "switches reach the FATES control globals" begin
            CLM.clm_fates_init!(CLM.CLMInstances(); nsites = 1, nlevsoil = 3, nlevdecomp = 1,
                                use_tree_damage = true, use_ed_st3 = true,
                                use_cohort_age_tracking = true, seeddisp_cadence = 2)
            @test CLM.hlm_use_tree_damage[] == CLM.itrue
            @test CLM.hlm_use_ed_st3[] == CLM.itrue
            @test CLM.hlm_use_ed_prescribed_phys[] == CLM.ifalse
            @test CLM.hlm_use_cohort_age_tracking[] == CLM.itrue
            @test CLM.hlm_seeddisp_cadence[] == 2
            @test CLM.hlm_use_inventory_init[] == CLM.ifalse
            # FATES's own check_allset still rejects ST3 + prescribed physiology
            @test_throws Exception CLM.clm_fates_init!(CLM.CLMInstances(); nsites = 1,
                nlevsoil = 3, nlevdecomp = 1, use_ed_st3 = true, use_ed_prescribed_phys = true)
            # defaults: every mode off
            CLM.clm_fates_init!(CLM.CLMInstances(); nsites = 1, nlevsoil = 3, nlevdecomp = 1)
            for r in (CLM.hlm_use_tree_damage, CLM.hlm_use_ed_st3, CLM.hlm_use_ed_prescribed_phys,
                      CLM.hlm_use_cohort_age_tracking, CLM.hlm_use_inventory_init)
                @test r[] == CLM.ifalse
            end
            @test CLM.hlm_seeddisp_cadence[] == 0
            @test CLM.CLMInstances().fates === nothing
        end

        @testset "host seed dispersal (WrapUpdateFatesSeedInOut / WrapGlobalSeedDispersal)" begin
            inst = CLM.CLMInstances()
            fates = CLM.clm_fates_init!(inst; nsites = 1, nlevsoil = 3, nlevdecomp = 1,
                                        seeddisp_cadence = 1)
            @test fates.seeds !== nothing
            @test size(fates.seeds.incoming_global) == (CLM.numpft[], 1)
            @test CLM.dispersal_date[] == 0          # the Fortran's unset-clock init date
            site = fates.sites[1]
            np = CLM.numpft[]
            # single gridcell -> no neighbours -> nothing ever comes back in
            @test isempty(CLM.lneighbors[][1].neighbor_indices)

            CLM.fates_get_and_set_time!(2002, 7, 1, 0)
            site.seed_out .= 1.0; site.seed_in .= 5.0
            CLM.fates_wrap_update_seed_inout!(fates)    # day 1 /= dispersal_date 0 -> reset
            @test all(iszero, site.seed_out) && all(iszero, site.seed_in)
            site.seed_out .= 2.0                        # today's SeedUpdate export
            CLM.fates_wrap_seed_outgoing!(fates)
            @test all(==(2.0), fates.seeds.outgoing_local[:, 1])
            CLM.fates_wrap_global_seed_dispersal!(fates) # sets the flag + the date
            @test CLM.dispersal_date[] == 1
            @test all(iszero, fates.seeds.outgoing_local)
            @test all(iszero, fates.seeds.incoming_global)

            # next day: the flag carried from yesterday's dispersal resets seed_out again
            CLM.fates_get_and_set_time!(2002, 7, 2, 0)
            site.seed_out .= 3.0
            CLM.fates_wrap_update_seed_inout!(fates)
            @test all(iszero, site.seed_out)

            # monthly cadence: within the month seed_out accumulates (no reset)
            CLM.hlm_seeddisp_cadence[] = CLM.fates_dispersal_cadence_monthly
            CLM.dispersal_date[] = 7; CLM.dispersal_flag[] = false
            site.seed_out .= 4.0
            CLM.fates_wrap_update_seed_inout!(fates)
            @test all(==(4.0), site.seed_out)
            @test all(iszero, site.seed_in)
            @test length(site.seed_out) == np
        end

        @testset "inventory init through clm_fates_init!" begin
            mktempdir() do d
                pss = joinpath(d, "site.pss"); css = joinpath(d, "site.css")
                write(pss, "time patch trk age area\n2000.0 p1 2 10.0 0.4\n2000.0 p2 2 60.0 0.6\n")
                write(css, "time patch dbh height pft nplant\n" *
                           "2000.0 p1 8.0 -1.0 2 0.05\n2000.0 p1 3.0 -1.0 6 0.10\n" *
                           "2000.0 p2 30.0 -1.0 2 0.02\n2000.0 p2 15.0 -1.0 2 0.04\n")
                ctl = joinpath(d, "ctl.txt")
                # file longitude negative, host longitude negative too: the Fortran
                # aborts here (no proximal site); the port wraps both into [0, 360).
                write(ctl, "type latitude longitude pss-name css-name\n1 51.36 -116.01 $pss $css\n")
                inst = CLM.CLMInstances()
                push!(inst.gridcell.lat, deg2rad(51.36)); push!(inst.gridcell.lon, deg2rad(-116.01))
                fates = CLM.clm_fates_init!(inst; nsites = 1, nlevsoil = 3, nlevdecomp = 1,
                                            use_inventory_init = true, inventory_ctrl_file = ctl)
                coh = _niche_cohorts(fates.sites[1])
                A = CLM.area
                @test sort(unique(round.([c.patch_area / A for c in coh], digits = 6))) == [0.4, 0.6]
                @test sort(unique([c.pft for c in coh])) == [2, 6]
                # density: css n is per m2 of the patch -> n * patch area
                @test sum(c.n for c in coh if c.pft == 6) ≈ 0.10 * 0.4 * A
                @test sum(c.n for c in coh if c.pft == 2) ≈ (0.05 * 0.4 + (0.02 + 0.04) * 0.6) * A
                @test maximum(c.dbh for c in coh) ≈ 30.0
                # coage 0, not NaN: a NaN coage blocks fuse_cohorts for these cohorts
                @test all(c.coage == 0.0 for c in coh)
                # too far away -> FATES endrun
                inst2 = CLM.CLMInstances()
                push!(inst2.gridcell.lat, deg2rad(10.0)); push!(inst2.gridcell.lon, deg2rad(10.0))
                @test_throws Exception CLM.clm_fates_init!(inst2; nsites = 1, nlevsoil = 3,
                    nlevdecomp = 1, use_inventory_init = true, inventory_ctrl_file = ctl)
            end
        end

        @testset "FATES restart file round-trip" begin
            inst = CLM.CLMInstances()
            fates = CLM.clm_fates_init!(inst; nsites = 1, nlevsoil = 3, nlevdecomp = 1,
                                        use_cohort_age_tracking = true)
            site = fates.sites[1]
            # perturb the state so the round trip carries non-cold-start values
            let c = site.oldest_patch.tallest
                k = 0
                while c !== nothing
                    k += 1; c.coage = 0.5 * k; c = c.shorter
                end
            end
            site.oldest_patch.age = 3.25
            ref = _niche_cohorts(site)
            ref_site = (site.cstatus, copy(site.dstatus), site.grow_deg_days,
                        copy(site.vegtemp_memory), site.spread)
            mktempdir() do d
                path = joinpath(d, "fates_rst.nc")
                CLM.write_fates_restart(path, fates)
                @test isfile(path)
                old_site = fates.sites[1]
                CLM.read_fates_restart!(path, inst)
                new_site = fates.sites[1]
                @test new_site !== old_site            # rebuilt from the file
                got = _niche_cohorts(new_site)
                @test length(got) == length(ref)
                @test all(isequal(a, b) for (a, b) in zip(got, ref))
                @test (new_site.cstatus, new_site.dstatus, new_site.grow_deg_days,
                       new_site.vegtemp_memory, new_site.spread) == ref_site
                # restarted cohorts get Fortran's NanValues + ZeroValues: the fields it
                # zeroes (fire mortality, ...) are 0 and the N/P demand sentinels are -9,
                # not the NaN defaults that broke the first post-restart mass balance
                let c = new_site.oldest_patch.tallest
                    while c !== nothing
                        @test c.fire_mort == 0.0 && c.crownfire_mort == 0.0 &&
                              c.daily_n_demand == -9.0
                        c = c.shorter
                    end
                end
                # init_site_vars!/zero_site! ran on the rebuilt site (the closed TODO)
                @test length(new_site.seed_out) == CLM.numpft[]
                @test length(new_site.dz_soil) == 3
            end
        end
    finally
        _niche_restore!(snap)
    end
end
