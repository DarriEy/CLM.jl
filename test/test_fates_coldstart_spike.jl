# test_fates_coldstart_spike.jl
# FATES live-driver wiring W1 + W2 SPIKE — the ownership-model proof.
#
# Proves that a finite, cold-started, carbon-only single FATES site can be attached
# to a CLMInstances via the `fates::Union{AbstractFatesInterface,Nothing}` field,
# and that this attached state survives the AD dual-copy (_make_dual_instances) and
# the GPU adapt UNTOUCHED (excluded, like surfdata) — the documented v1 boundary
# that FATES columns are CPU-only / non-differentiable.

using Test
using CLM
using Adapt
import ForwardDiff
using ForwardDiff: Dual

# Walk youngest -> oldest collecting patches.
function _spike_patch_list(site)
    ps = CLM.fates_patch_type[]
    p = site.youngest_patch
    while p !== nothing
        push!(ps, p)
        p = p.older
    end
    return ps
end

# Walk tallest -> shortest collecting cohorts of a patch.
function _spike_cohort_list(patch)
    cs = CLM.fates_cohort_type[]
    c = patch.tallest
    while c !== nothing
        push!(cs, c)
        c = c.shorter
    end
    return cs
end

@testset "FATES W1+W2 cold-start spike" begin
    # clm_fates_init! mutates the FATES module-global control Refs + parameter
    # tables. Save/restore so the surrounding default suite is unperturbed.
    old_global   = CLM.prt_global[]
    old_numel    = CLM.num_elements[]
    old_ellist   = copy(CLM.element_list)
    old_parteh   = CLM.hlm_parteh_mode[]
    old_use_sp   = CLM.hlm_use_sp[]
    old_hydro    = CLM.hlm_use_planthydro[]
    old_agetrk   = CLM.hlm_use_cohort_age_tracking[]
    old_nleafage = CLM.nleafage[]
    old_numpft   = CLM.numpft[]
    old_tod      = CLM.hlm_current_tod[]
    old_nlevsc   = CLM.nlevsclass[]
    old_nlevca   = CLM.nlevcoage[]
    old_nlevage  = CLM.nlevage[]
    old_nlevdam  = CLM.nlevdamage[]
    old_nlevhgt  = CLM.nlevheight[]
    old_sfp      = CLM.SFParams[]
    old_restart  = CLM.hlm_is_restart[]
    old_nocomp   = CLM.hlm_use_nocomp[]
    old_biogeog  = CLM.hlm_use_fixed_biogeog[]
    old_luh      = CLM.hlm_use_luh[]
    old_inv      = CLM.hlm_use_inventory_init[]
    old_damage   = CLM.hlm_use_tree_damage[]
    old_doy      = CLM.hlm_day_of_year[]
    old_dpy      = CLM.hlm_days_per_year[]
    old_freqday  = CLM.hlm_freq_day[]
    old_nharv    = CLM.hlm_num_lu_harvest_cats[]
    old_ch4      = CLM.hlm_use_ch4[]
    old_vert     = CLM.hlm_use_vertsoilc[]
    old_spit     = CLM.hlm_spitfire_mode[]
    old_paramfp  = CLM.fates_maxPatchesPerSite[]
    # read_fates_params! overwrites the FATES parameter globals; save/restore them
    # so the surrounding default suite is unperturbed (a real param table now).
    old_edparams = CLM.EDParams[]
    old_edpft    = CLM.EDPftvarcon_inst[]
    old_paramd   = CLM.ParamDerived[]

    try
        nlevsoil = 3
        nlevdecomp = 1     # non-vertsoilc => one decomp layer

        # =================================================================
        # W1 + W2: build + cold-start a single FATES site, attach to inst.
        # numpft now comes from the REAL FATES param file (read_fates_params!),
        # not a hand-set count — the carbon-only cold start must be finite on
        # the real per-PFT trait/allometry table.
        # =================================================================
        inst = CLM.CLMInstances()
        @test inst.fates === nothing   # default OFF — inert

        fates = CLM.clm_fates_init!(inst; nsites = 1,
                                    nlevsoil = nlevsoil, nlevdecomp = nlevdecomp)
        numpft = CLM.numpft[]          # resolved from the param file (14)
        @test numpft == 14

        # ---- attachment + container shape ----
        @test inst.fates !== nothing
        @test inst.fates === fates
        @test inst.fates isa CLM.fates_interface_type
        @test inst.fates isa CLM.AbstractFatesInterface
        @test inst.fates.nsites == 1
        @test length(inst.fates.sites) == 1
        @test length(inst.fates.bc_in) == 1
        @test length(inst.fates.bc_out) == 1

        site = inst.fates.sites[1]

        # ---- finite-site assertions: patches ----
        patches = _spike_patch_list(site)
        @test length(patches) >= 1
        # NBG default cold start => exactly one primaryland patch holding the area.
        @test length(patches) == 1
        p1 = patches[1]
        @test p1.land_use_label == CLM.primaryland

        # Sum of patch areas ≈ the FATES notional site area (~10000 m²).
        total_area = sum(p.area for p in patches)
        @test total_area ≈ CLM.area rtol = 1e-9
        @test CLM.area ≈ 10000.0

        # ---- finite-site assertions: cohorts ----
        cohorts = vcat((_spike_cohort_list(p) for p in patches)...)
        @test length(cohorts) >= 1
        # One cohort per PFT seeded at cold start.
        @test length(cohorts) == numpft
        for c in cohorts
            @test isfinite(c.dbh)    && c.dbh    > 0.0
            @test isfinite(c.n)      && c.n      > 0.0
            @test isfinite(c.height) && c.height > 0.0
        end
        @test sort([c.pft for c in cohorts]) == collect(1:numpft)

        # ---- no NaN in the site's key arrays ----
        @test all(isfinite, site.zi_soil)
        @test all(isfinite, site.dz_soil)
        @test all(isfinite, site.z_soil)
        @test all(isfinite, site.area_PFT)
        @test site.mass_balance[1].old_stock > 0.0
        @test isfinite(site.mass_balance[1].old_stock)
        # bc_out aggregates that ed_update_site fills must be finite (no NaN leaked).
        bc_out = inst.fates.bc_out[1]
        @test all(isfinite, bc_out.elai_pa)
        @test all(isfinite, bc_out.htop_pa)

        # =================================================================
        # OWNERSHIP-MODEL PROOF #1 — the AD dual-copy SKIPS .fates.
        # =================================================================
        D = Dual{Nothing, Float64, 1}
        inst_d = CLM._make_dual_instances(inst, D)
        @test inst_d isa CLM.CLMInstances
        # The FATES site is attached by reference, NOT Float64->Dual converted.
        @test inst_d.fates === inst.fates
        @test inst_d.fates !== nothing

        # =================================================================
        # OWNERSHIP-MODEL PROOF #2 — CPU identity adapt leaves .fates intact.
        # =================================================================
        inst_a = adapt(Array, inst)
        @test inst_a.fates === inst.fates
        @test inst_a.fates !== nothing
        @test inst_a.fates.nsites == 1

        # =================================================================
        # FIXED-BIOGEOG CLIMATE SCREEN — seed only climate-appropriate PFTs.
        # The default cold start above seeded all 14 PFTs (fixed-biogeog OFF).
        # The :drop_cold_deciduous screen restricts seeding + use_this_pft to
        # the non-cold-deciduous (evergreen + drought-deciduous) PFTs, which
        # retires the cold-deciduous-driven multi-year boom→bust die-back.
        # =================================================================
        season_decid = copy(CLM.prt_params.season_decid)   # from the real param file
        want = [ft for ft in 1:numpft if season_decid[ft] == CLM.ifalse]
        @test !isempty(want) && length(want) < numpft   # a real, partial screen

        inst2 = CLM.CLMInstances()
        CLM.clm_fates_init!(inst2; nsites = 1, nlevsoil = nlevsoil,
                            nlevdecomp = nlevdecomp,
                            fates_biogeog_screen = :drop_cold_deciduous)
        @test CLM.hlm_use_fixed_biogeog[] == CLM.itrue
        site2 = inst2.fates.sites[1]
        # use_this_pft flags exactly the non-cold-deciduous PFTs.
        for ft in 1:numpft
            @test site2.use_this_pft[ft] ==
                  (season_decid[ft] == CLM.ifalse ? CLM.itrue : CLM.ifalse)
        end
        # Cold-start cohorts cover exactly the screened PFT set — no cold-deciduous.
        coh2 = vcat((_spike_cohort_list(p) for p in _spike_patch_list(site2))...)
        @test sort(unique(c.pft for c in coh2)) == want
        @test all(season_decid[c.pft] == CLM.ifalse for c in coh2)
        @test length(coh2) == length(want)

        # Explicit per-FATES-PFT presence vector (general mechanism): seed PFT 1 only.
        inst3 = CLM.CLMInstances()
        areafrac = zeros(numpft); areafrac[1] = 1.0
        CLM.clm_fates_init!(inst3; nsites = 1, nlevsoil = nlevsoil,
                            nlevdecomp = nlevdecomp, fates_pft_areafrac = areafrac)
        @test CLM.hlm_use_fixed_biogeog[] == CLM.itrue
        site3 = inst3.fates.sites[1]
        @test site3.use_this_pft[1] == CLM.itrue
        @test all(site3.use_this_pft[ft] == CLM.ifalse for ft in 2:numpft)
        coh3 = vcat((_spike_cohort_list(p) for p in _spike_patch_list(site3))...)
        @test all(c.pft == 1 for c in coh3)

        # Guard: an all-zero (empty) presence vector is rejected.
        @test_throws ErrorException CLM.clm_fates_init!(CLM.CLMInstances(); nsites = 1,
            nlevsoil = nlevsoil, nlevdecomp = nlevdecomp,
            fates_pft_areafrac = zeros(numpft))

    finally
        CLM.prt_global[]                  = old_global
        CLM.num_elements[]                = old_numel
        empty!(CLM.element_list); append!(CLM.element_list, old_ellist)
        CLM.hlm_parteh_mode[]             = old_parteh
        CLM.hlm_use_sp[]                  = old_use_sp
        CLM.hlm_use_planthydro[]          = old_hydro
        CLM.hlm_use_cohort_age_tracking[] = old_agetrk
        CLM.nleafage[]                    = old_nleafage
        CLM.numpft[]                      = old_numpft
        CLM.hlm_current_tod[]             = old_tod
        CLM.nlevsclass[]                  = old_nlevsc
        CLM.nlevcoage[]                   = old_nlevca
        CLM.nlevage[]                     = old_nlevage
        CLM.nlevdamage[]                  = old_nlevdam
        CLM.nlevheight[]                  = old_nlevhgt
        CLM.SFParams[]                    = old_sfp
        CLM.hlm_is_restart[]              = old_restart
        CLM.hlm_use_nocomp[]              = old_nocomp
        CLM.hlm_use_fixed_biogeog[]       = old_biogeog
        CLM.hlm_use_luh[]                 = old_luh
        CLM.hlm_use_inventory_init[]      = old_inv
        CLM.hlm_use_tree_damage[]         = old_damage
        CLM.hlm_day_of_year[]             = old_doy
        CLM.hlm_days_per_year[]           = old_dpy
        CLM.hlm_freq_day[]                = old_freqday
        CLM.hlm_num_lu_harvest_cats[]     = old_nharv
        CLM.hlm_use_ch4[]                 = old_ch4
        CLM.hlm_use_vertsoilc[]           = old_vert
        CLM.hlm_spitfire_mode[]           = old_spit
        CLM.fates_maxPatchesPerSite[]     = old_paramfp
        CLM.EDParams[]                    = old_edparams
        CLM.EDPftvarcon_inst[]            = old_edpft
        CLM.ParamDerived[]                = old_paramd
    end
end
