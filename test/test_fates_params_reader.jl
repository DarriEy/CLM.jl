# test_fates_params_reader.jl
# Real FATES parameter-file reader (read_fates_params!) tests.
#
# Verifies that read_fates_params! parses the committed FATES default parameter
# file (data/fates/fates_params_default.cdl) and populates the FATES parameter
# globals with the physically-meaningful values (matching the CDL), and that a
# carbon-only FATES site still cold-starts FINITE on the real 14-PFT params.

using Test
using CLM

# Walk youngest -> oldest collecting patches.
function _rdr_patch_list(site)
    ps = CLM.fates_patch_type[]
    p = site.youngest_patch
    while p !== nothing
        push!(ps, p)
        p = p.older
    end
    return ps
end

# Walk tallest -> shortest collecting cohorts of a patch.
function _rdr_cohort_list(patch)
    cs = CLM.fates_cohort_type[]
    c = patch.tallest
    while c !== nothing
        push!(cs, c)
        c = c.shorter
    end
    return cs
end

@testset "FATES real param-file reader" begin
    # read_fates_params! mutates the FATES module-global parameter Refs + dim
    # Refs. Save/restore so the surrounding suite is unperturbed.
    old_parteh   = CLM.hlm_parteh_mode[]
    old_numpft   = CLM.numpft[]
    old_nleafage = CLM.nleafage[]
    old_nlevsc   = CLM.nlevsclass[]
    old_nlevage  = CLM.nlevage[]
    old_nlevca   = CLM.nlevcoage[]
    old_nlevdam  = CLM.nlevdamage[]
    old_nlevhgt  = CLM.nlevheight[]
    old_edparams = CLM.EDParams[]
    old_edpft    = CLM.EDPftvarcon_inst[]
    old_paramd   = CLM.ParamDerived[]
    old_sfp      = CLM.SFParams[]
    old_global   = CLM.prt_global[]
    old_numel    = CLM.num_elements[]
    old_ellist   = copy(CLM.element_list)

    try
        # parteh mode must be set before the read (PRTDerivedParams! + the
        # consistency checks branch on it). Carbon-only, like clm_fates_init!.
        CLM.hlm_parteh_mode[] = CLM.prt_carbon_allom_hyp

        # =================================================================
        # read_fates_params! runs and reports coverage.
        # =================================================================
        s = CLM.read_fates_params!()

        # The default file declares 14 FATES PFTs.
        @test s.npft == 14
        @test CLM.numpft[] == 14

        # Every Julia-consumer (registered) parameter must exist in the file.
        @test isempty(s.registered_missing)

        # Coverage: every registered param was filled from the file.
        @test s.n_populated == s.n_registered

        # The only CDL fates_* vars with no Julia consumer are the descriptive
        # char-string name arrays (no numeric field). Everything numeric maps.
        @test s.n_cdl_unmapped == 5
        @test Set(s.cdl_unmapped) == Set([
            "fates_pftname", "fates_hydro_organ_name", "fates_alloc_organ_name",
            "fates_landuseclass_name", "fates_litterclass_name"])

        # =================================================================
        # Sample known CDL values match the populated globals.
        # =================================================================
        p  = CLM.prt_params
        ev = CLM.EDPftvarcon_inst[]
        edp = CLM.EDParams[]
        sf  = CLM.SFParams[]

        # 1. PRT scalars-per-PFT (carbon-to-biomass, wood density, growth resp).
        @test all(==(2.0), p.c2b)                       # fates_c2b = 2 (all PFTs)
        @test p.wood_density[1] ≈ 0.548327              # fates_wood_density[1]
        @test p.grperc[1] ≈ 0.11                        # fates_grperc[1]

        # 2. EDPftvarcon per-PFT (recruit height min, C3/C4 pathway).
        @test ev.hgt_min[1] ≈ 1.3                       # fates_recruit_height_min[1]
        @test ev.hgt_min[7] ≈ 0.2                       # shrub PFT 7
        @test ev.c3psn[1] ≈ 1.0                         # PFT 1 is C3
        @test ev.c3psn[14] ≈ 0.0                        # PFT 14 (c4_grass) is C4

        # 3. EDPftvarcon PFT x leafage 2-D (vcmax25top) — dim order proof: the CDL
        #    declares (leafage, pft); the [pft, leafage] field must read PFT-first.
        @test ev.vcmax25top[1, 1] ≈ 50.0                # fates_leaf_vcmax25top pft1
        @test ev.vcmax25top[12, 1] ≈ 86.0               # pft12
        @test ev.vcmax25top[14, 1] ≈ 78.0               # pft14

        # 4. PRT 2-D (PFT x organ) — stoich_nitr; leaf organ (col 1) pft1 = 0.033.
        @test p.nitr_stoich_p1[1, 1] ≈ 0.033

        # 5. PRT organ id mapping (organ-dimensioned 1-D).
        @test p.organ_id == [1, 2, 3, 6]

        # 6. EDParams scalars (maxcohort, q10, regen/rad/stomatal model switches).
        @test edp.max_cohort_per_patch == 100           # fates_maxcohort
        @test edp.q10_mr ≈ 1.5                           # fates_q10_mr
        @test edp.regeneration_model == 1               # fates_regeneration_model
        @test edp.radiation_model == 1                  # fates_rad_model (Norman)
        @test edp.stomatal_model == 1                   # fates_leaf_stomatal_model

        # 7. SPITFIRE per-CWD fraction + scalar threshold.
        @test sf.SF_val_CWD_frac ≈ [0.045, 0.075, 0.21, 0.67]
        @test sf.SF_val_fire_threshold ≈ 50.0

        # 8. History-bin dim Refs sized from the file's bin-edge arrays.
        @test CLM.nlevsclass[] == 13   # 13 sizeclass edges
        @test CLM.nlevage[]    == 7    # 7 ageclass edges
        @test CLM.nlevcoage[]  == 2    # 2 coage edges
        @test CLM.nlevdamage[] == 2    # 2 damage edges
        @test CLM.nleafage[]   == 1    # leafage class

        # =================================================================
        # No NaN in the populated core arrays.
        # =================================================================
        @test !any(isnan, p.c2b)
        @test !any(isnan, p.wood_density)
        @test !any(isnan, p.slatop)
        @test !any(isnan, p.nitr_stoich_p1)
        @test !any(isnan, ev.vcmax25top)
        @test !any(isnan, ev.hgt_min)
        @test !any(isnan, ev.initd)
        @test !any(isnan, sf.SF_val_CWD_frac)
        @test !any(isnan, edp.dinc_vai)   # VAI bins computed from file params

        # Derived params (jmax25 = 1.67*vcmax25; branch_frac = sum CWD_frac[1:3]).
        pd = CLM.ParamDerived[]
        @test pd.jmax25top[1, 1] ≈ 1.67 * 50.0
        @test pd.branch_frac[1] ≈ sum(sf.SF_val_CWD_frac[1:3])
        @test !any(isnan, pd.jmax25top)
        @test !any(isnan, pd.branch_frac)

    finally
        CLM.hlm_parteh_mode[] = old_parteh
        CLM.numpft[]          = old_numpft
        CLM.nleafage[]        = old_nleafage
        CLM.nlevsclass[]      = old_nlevsc
        CLM.nlevage[]         = old_nlevage
        CLM.nlevcoage[]       = old_nlevca
        CLM.nlevdamage[]      = old_nlevdam
        CLM.nlevheight[]      = old_nlevhgt
        CLM.EDParams[]        = old_edparams
        CLM.EDPftvarcon_inst[] = old_edpft
        CLM.ParamDerived[]    = old_paramd
        CLM.SFParams[]        = old_sfp
        CLM.prt_global[]      = old_global
        CLM.num_elements[]    = old_numel
        empty!(CLM.element_list); append!(CLM.element_list, old_ellist)
    end
end

@testset "FATES cold-start finite on real params" begin
    # clm_fates_init! mutates many FATES module-global control + parameter Refs.
    # Save/restore the ones it touches so the surrounding suite is unperturbed.
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
    old_edparams = CLM.EDParams[]
    old_edpft    = CLM.EDPftvarcon_inst[]
    old_paramd   = CLM.ParamDerived[]
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

    try
        nlevsoil = 3
        nlevdecomp = 1

        inst = CLM.CLMInstances()
        @test inst.fates === nothing

        fates = CLM.clm_fates_init!(inst; nsites = 1,
                                    nlevsoil = nlevsoil, nlevdecomp = nlevdecomp)
        numpft = CLM.numpft[]
        @test numpft == 14   # PFT count from the real param file

        @test inst.fates === fates
        @test inst.fates.nsites == 1

        site = inst.fates.sites[1]

        # ---- finite site: one primaryland patch holding the whole area ----
        patches = _rdr_patch_list(site)
        @test length(patches) == 1
        @test patches[1].land_use_label == CLM.primaryland
        @test sum(p.area for p in patches) ≈ CLM.area rtol = 1e-9

        # ---- finite cohorts: one per PFT, positive dbh/n/height ----
        cohorts = vcat((_rdr_cohort_list(p) for p in patches)...)
        @test length(cohorts) == numpft
        for c in cohorts
            @test isfinite(c.dbh)    && c.dbh    > 0.0
            @test isfinite(c.n)      && c.n      > 0.0
            @test isfinite(c.height) && c.height > 0.0
        end
        @test sort([c.pft for c in cohorts]) == collect(1:numpft)

        # ---- no NaN in site/bc_out aggregates ----
        @test all(isfinite, site.zi_soil)
        @test all(isfinite, site.dz_soil)
        @test all(isfinite, site.z_soil)
        @test all(isfinite, site.area_PFT)
        @test site.mass_balance[1].old_stock > 0.0
        @test isfinite(site.mass_balance[1].old_stock)
        bc_out = inst.fates.bc_out[1]
        @test all(isfinite, bc_out.elai_pa)
        @test all(isfinite, bc_out.htop_pa)

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
        CLM.EDParams[]                    = old_edparams
        CLM.EDPftvarcon_inst[]            = old_edpft
        CLM.ParamDerived[]                = old_paramd
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
    end
end
