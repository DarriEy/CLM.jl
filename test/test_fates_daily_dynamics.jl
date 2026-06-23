# test_fates_daily_dynamics.jl
# FATES live-driver wiring W5 — the DAILY demographic step.
#
# Cold-starts a single carbon-only FATES site (via clm_fates_init!), packs a
# representative daily bc_in, and drives one (then a few) `ed_ecosystem_dynamics` +
# `ed_update_site` daily steps through the W5 helper `fates_daily_dynamics_step!`
# (the same routine the gated `clm_drv` driver hook calls at the day boundary).
#
# Asserts: the step completes without error; the final TotalBalanceCheck passes
# (carbon mass conserved within FATES' 1e-5 tolerance — it throws otherwise);
# cohort/patch state stays finite and physical (n >= 0, dbh > 0, patch areas still
# sum to AREA); and the unpacked canopy structure (elai/htop) is finite. A 3-day
# loop confirms stability (no NaN / blowup).
#
# clm_fates_init! mutates FATES module-global control Refs + parameter tables;
# save/restore them so the surrounding default suite is unperturbed.

using Test
using CLM

@testset "FATES W5 daily dynamics" begin
    # ---- save FATES globals (same set as the W3/W4 driver-hooks test) ----
    old_global   = CLM.prt_global[]
    old_numel    = CLM.num_elements[]
    old_ellist   = copy(CLM.element_list)
    old_parteh   = CLM.hlm_parteh_mode[]
    old_use_sp   = CLM.hlm_use_sp[]
    old_st3      = CLM.hlm_use_ed_st3[]
    old_hydro    = CLM.hlm_use_planthydro[]
    old_agetrk   = CLM.hlm_use_cohort_age_tracking[]
    old_nleafage = CLM.nleafage[]
    old_numpft   = CLM.numpft[]
    old_tod      = CLM.hlm_current_tod[]
    old_nlevsc   = CLM.nlevsclass[]
    old_nlevca   = CLM.nlevcoage[]
    old_nlevage  = CLM.nlevage[]
    old_nlevdam  = CLM.nlevdamage[]
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
    old_nlevsno  = CLM.varpar.nlevsno
    old_nlevsoi  = CLM.varpar.nlevsoi
    old_nlevgrnd = CLM.varpar.nlevgrnd
    old_nlevmaxu = CLM.varpar.nlevmaxurbgrnd

    try
        numpft   = 2
        nlevsoil = 3
        nlevdecomp = 1

        # ---- minimal subgrid: 1 gridcell, 1 landunit, 1 column, 2 patches ----
        # patch 1 = bare ground (col.patchi[c]+0), patch 2 = vegetated (+1).
        ng, nl, nc, np = 1, 1, 1, 2
        CLM.varpar.nlevsno = 5
        CLM.varpar.nlevsoi = nlevsoil
        CLM.varpar.nlevgrnd = nlevsoil
        CLM.varpar.nlevmaxurbgrnd = nlevsoil

        inst = CLM.CLMInstances()
        CLM.clm_instInit!(inst; ng=ng, nl=nl, nc=nc, np=np)

        col = inst.column
        col.gridcell[1] = 1
        col.patchi[1]   = 1          # bareground patch index 1
        col.npatches[1] = 2
        col.is_fates[1] = true
        inst.patch.column[1] = 1
        inst.patch.column[2] = 1
        inst.patch.is_fates[1] = true
        inst.patch.is_fates[2] = true

        cveg = 1                      # FATES column
        pveg = col.patchi[cveg] + 1   # vegetated patch = 2

        # ---- attach a cold-started carbon-only FATES site ----
        fates = CLM.clm_fates_init!(inst; nsites = 1, numpft_in = numpft,
                                    nlevsoil = nlevsoil, nlevdecomp = nlevdecomp)
        @test inst.fates !== nothing
        @test inst.fates.nsites == 1

        # ---- representative daily soil-column state ----
        ss   = inst.soilstate
        temp = inst.temperature
        wdb  = inst.water.waterdiagnosticbulk_inst
        cs   = inst.canopystate

        joff = CLM.varpar.nlevsno
        for j in 1:nlevsoil
            temp.t_soisno_col[cveg, joff + j] = 290.0
            wdb.h2osoi_liqvol_col[cveg, j]    = 0.30
            ss.watsat_col[cveg, j]            = 0.45
            ss.eff_porosity_col[cveg, j]      = 0.45
            ss.sucsat_col[cveg, j]            = 100.0
            ss.bsw_col[cveg, j]               = 6.0
        end

        site = inst.fates.sites[1]

        # ---- record pre-step carbon stock + cohort/patch census ----
        function _census(site)
            ncoh = 0; min_n = Inf; min_dbh = Inf; area_sum = 0.0
            pp = site.oldest_patch
            while pp !== nothing
                area_sum += pp.area
                cc = pp.tallest
                while cc !== nothing
                    ncoh += 1
                    min_n   = min(min_n, cc.n)
                    min_dbh = min(min_dbh, cc.dbh)
                    cc = cc.shorter
                end
                pp = pp.younger
            end
            return ncoh, min_n, min_dbh, area_sum
        end

        ncoh0, _, _, area0 = _census(site)
        @test ncoh0 >= 1                      # cold start seeded at least one cohort
        @test area0 ≈ CLM.area atol = 1e-6    # patch areas sum to AREA

        # ================================================================
        # ONE daily demographic step via the W5 helper (the driver hook).
        # TotalBalanceCheck(-1) inside throws on mass-imbalance > 1e-5, so a
        # clean return IS the conservation assertion.
        # ================================================================
        CLM.fates_daily_dynamics_step!(inst; nlevsoil=nlevsoil,
                                       nlevdecomp=nlevdecomp, use_fates_bgc=false)

        # The final mass-conservation audit passed (else it would have thrown);
        # confirm err_fates (net_flux - change_in_stock) is within tolerance.
        for el in 1:CLM.num_elements[]
            err = site.mass_balance[el].err_fates
            @test isfinite(err)
            @test abs(err) <= 1.0e-5          # FATES TotalBalanceCheck tolerance
        end

        # Cohort / patch state finite + physical.
        ncoh1, min_n1, min_dbh1, area1 = _census(site)
        @test ncoh1 >= 1
        @test isfinite(min_n1)   && min_n1   >= 0.0
        @test isfinite(min_dbh1) && min_dbh1 >  0.0
        @test area1 ≈ CLM.area atol = 1e-6    # patch areas still sum to AREA

        # Unpacked canopy structure finite.
        @test isfinite(cs.elai_patch[pveg])
        @test isfinite(cs.esai_patch[pveg])
        @test isfinite(cs.htop_patch[pveg]) && cs.htop_patch[pveg] >= 0.0
        @test isfinite(cs.hbot_patch[pveg])
        @test isfinite(cs.z0m_patch[pveg])

        # ================================================================
        # STABILITY: a few more daily steps — no NaN / blowup, mass conserved.
        # ================================================================
        for _step in 1:3
            CLM.fates_daily_dynamics_step!(inst; nlevsoil=nlevsoil,
                                           nlevdecomp=nlevdecomp, use_fates_bgc=false)
            for el in 1:CLM.num_elements[]
                @test abs(site.mass_balance[el].err_fates) <= 1.0e-5
            end
            _, min_n, min_dbh, area_s = _census(site)
            @test isfinite(min_n)   && min_n   >= 0.0
            @test isfinite(min_dbh) && min_dbh >  0.0
            @test area_s ≈ CLM.area atol = 1e-6
            @test isfinite(cs.elai_patch[pveg])
            @test isfinite(cs.htop_patch[pveg])
        end

    finally
        CLM.prt_global[]                  = old_global
        CLM.num_elements[]                = old_numel
        empty!(CLM.element_list); append!(CLM.element_list, old_ellist)
        CLM.hlm_parteh_mode[]             = old_parteh
        CLM.hlm_use_sp[]                  = old_use_sp
        CLM.hlm_use_ed_st3[]              = old_st3
        CLM.hlm_use_planthydro[]          = old_hydro
        CLM.hlm_use_cohort_age_tracking[] = old_agetrk
        CLM.nleafage[]                    = old_nleafage
        CLM.numpft[]                      = old_numpft
        CLM.hlm_current_tod[]             = old_tod
        CLM.nlevsclass[]                  = old_nlevsc
        CLM.nlevcoage[]                   = old_nlevca
        CLM.nlevage[]                     = old_nlevage
        CLM.nlevdamage[]                  = old_nlevdam
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
        CLM.varpar.nlevsno                = old_nlevsno
        CLM.varpar.nlevsoi                = old_nlevsoi
        CLM.varpar.nlevgrnd               = old_nlevgrnd
        CLM.varpar.nlevmaxurbgrnd         = old_nlevmaxu
    end
end
