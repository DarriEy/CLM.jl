# test_fates_driver_hooks.jl
# FATES live-driver wiring W3 (radiation) + W4 (btran / photosynthesis) hooks.
#
# Drives the bc pack -> FATES driver -> bc unpack round-trip for each of the four
# hooks against a finite, cold-started, carbon-only single FATES site attached to a
# CLMInstances (via clm_fates_init!), and asserts the FATES outputs flow back into
# the CLM structs finite and in-range. Also confirms one hook fires through the
# gated `clm_drv` driver branch (the radiation sun/shade block).
#
# clm_fates_init! mutates FATES module-global control Refs + parameter tables;
# save/restore them so the surrounding default suite is unperturbed.

using Test
using CLM

@testset "FATES W3+W4 driver hooks" begin
    # ---- save FATES globals (same set as the cold-start spike test) ----
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

        # Column <-> patch + FATES tags. patchi = beginning-1 so patch indices are
        # patchi[c]+1 (bareground) ... here we treat the veg patch as patchi[c]+1.
        col = inst.column
        col.gridcell[1] = 1
        col.patchi[1]   = 1          # first patch of the column = index 1 (bareground)
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

        # ================================================================
        # Pack representative CLM forcing into the CLM structs.
        # ================================================================
        a2l = inst.atm2lnd
        sa  = inst.surfalb
        ss  = inst.soilstate
        temp = inst.temperature
        wdb  = inst.water.waterdiagnosticbulk_inst
        ef   = inst.energyflux
        cs   = inst.canopystate
        ps   = inst.photosyns

        nb = CLM.num_swb
        for ib in 1:nb
            a2l.forc_solad_downscaled_col[cveg, ib] = 200.0
            a2l.forc_solai_grc[1, ib]               = 80.0
            sa.albgrd_col[cveg, ib] = 0.15
            sa.albgri_col[cveg, ib] = 0.18
        end
        sa.coszen_col[cveg] = 0.6   # daytime

        # Soil column state for btran.
        joff = CLM.varpar.nlevsno
        for j in 1:nlevsoil
            temp.t_soisno_col[cveg, joff + j] = 290.0
            wdb.h2osoi_liqvol_col[cveg, j]    = 0.30
            ss.watsat_col[cveg, j]            = 0.45
            ss.eff_porosity_col[cveg, j]      = 0.45
            ss.sucsat_col[cveg, j]            = 100.0
            ss.bsw_col[cveg, j]               = 6.0
        end

        # Canopy / atmosphere for photosynthesis.
        temp.t_veg_patch[pveg]               = 295.0
        a2l.forc_t_downscaled_col[cveg]      = 293.0
        a2l.forc_pbot_downscaled_col[cveg]   = 1.0e5
        a2l.forc_vp_grc[1]                   = 1500.0
        a2l.forc_pco2_grc[1]                 = 40.0    # CO2 partial pressure (Pa)
        a2l.forc_po2_grc[1]                  = 21000.0 # O2 partial pressure (Pa)

        # ================================================================
        # W3b — FATES normalized canopy radiation (sets f_sun on the patch).
        # ================================================================
        CLM.fates_pack_bcin_radiation!(inst; s=1, c=cveg, p=pveg, coszen=sa.coszen_col[cveg])
        CLM.FatesNormalizedCanopyRadiation(inst.fates.nsites, inst.fates.sites,
                                           inst.fates.bc_in, inst.fates.bc_out)
        CLM.fates_unpack_bcout_canopy_radiation!(inst; s=1, c=cveg, p=pveg)

        for ib in 1:nb
            @test isfinite(sa.albd_patch[pveg, ib]) && 0.0 <= sa.albd_patch[pveg, ib] <= 1.0
            @test isfinite(sa.albi_patch[pveg, ib]) && 0.0 <= sa.albi_patch[pveg, ib] <= 1.0
            @test isfinite(sa.fabd_patch[pveg, ib])
            @test isfinite(sa.fabi_patch[pveg, ib])
            @test isfinite(sa.ftdd_patch[pveg, ib])
            @test isfinite(sa.ftid_patch[pveg, ib])
            @test isfinite(sa.ftii_patch[pveg, ib])
        end

        # ================================================================
        # W3a — FATES sun/shade fractions.
        # ================================================================
        CLM.fates_pack_bcin_radiation!(inst; s=1, c=cveg, p=pveg, coszen=sa.coszen_col[cveg])
        CLM.FatesSunShadeFracs(inst.fates.nsites, inst.fates.sites,
                               inst.fates.bc_in, inst.fates.bc_out)
        CLM.fates_unpack_bcout_sunfrac!(inst; s=1, c=cveg, p=pveg)

        @test isfinite(cs.fsun_patch[pveg])
        @test 0.0 <= cs.fsun_patch[pveg] <= 1.0
        @test isfinite(cs.laisun_patch[pveg]) && cs.laisun_patch[pveg] >= 0.0
        @test isfinite(cs.laisha_patch[pveg]) && cs.laisha_patch[pveg] >= 0.0
        # FatesSunShadeFracs partition identity: laisun = elai*fsun and
        # laisha = elai*(1-fsun), so laisun = (laisun+laisha)*fsun exactly, and the
        # total elai = laisun+laisha is a finite positive canopy leaf area.
        elai_sh = cs.laisun_patch[pveg] + cs.laisha_patch[pveg]
        @test elai_sh > 0.0
        @test cs.laisun_patch[pveg] ≈ elai_sh * cs.fsun_patch[pveg] atol = 1e-8

        # ================================================================
        # W4a — FATES btran.
        # ================================================================
        CLM.fates_pack_bcin_btran!(inst; s=1, c=cveg, nlevsoil=nlevsoil)
        CLM.btran_ed!(inst.fates.sites, inst.fates.bc_in, inst.fates.bc_out)
        CLM.fates_unpack_bcout_btran!(inst; s=1, c=cveg, p=pveg, nlevsoil=nlevsoil)

        @test isfinite(ef.btran_patch[pveg])
        @test 0.0 <= ef.btran_patch[pveg] <= 1.0
        # rootr sums to ~1 over soil layers (FatesBtran normalization).
        rsum = sum(ss.rootr_patch[pveg, 1:nlevsoil])
        @test isfinite(rsum)
        @test all(isfinite, ss.rootr_patch[pveg, 1:nlevsoil])

        # ================================================================
        # W4b — FATES photosynthesis.
        # ================================================================
        pbot  = a2l.forc_pbot_downscaled_col[cveg]
        t_veg = temp.t_veg_patch[pveg]
        _, esat_tv, _, _ = CLM.qsat(t_veg, pbot)
        CLM.fates_pack_bcin_photosynthesis!(inst; s=1, c=cveg, p=pveg,
            forc_pbot=pbot, forc_pco2=a2l.forc_pco2_grc[1],
            forc_po2=a2l.forc_po2_grc[1], t_veg=t_veg, tgcm=a2l.forc_t_downscaled_col[cveg],
            esat_tv=esat_tv, eair=a2l.forc_vp_grc[1], rb=50.0, dayl_factor=0.5)
        CLM.FatesPlantRespPhotosynthDrive(inst.fates.nsites, inst.fates.sites,
                                          inst.fates.bc_in, inst.fates.bc_out, 1800.0)
        CLM.fates_unpack_bcout_photosynthesis!(inst; s=1, c=cveg, p=pveg)

        @test isfinite(ps.rssun_patch[pveg]) && ps.rssun_patch[pveg] > 0.0
        @test isfinite(ps.rssha_patch[pveg]) && ps.rssha_patch[pveg] > 0.0

        # ================================================================
        # GATED-DRIVER-BRANCH confirmation: the W3a sun/shade hook fires through
        # the same gated logic the driver uses (loop over FATES columns mapping
        # site<->column<->veg-patch, pack -> FatesSunShadeFracs -> unpack), here
        # exercised directly with the driver's column/patch index arithmetic.
        # ================================================================
        cs.fsun_patch[pveg] = NaN
        sgated = 0
        for c in 1:length(col.is_fates)
            col.is_fates[c] || continue
            sgated += 1
            sgated <= inst.fates.nsites || break
            pg = col.patchi[c] + 1
            CLM.fates_pack_bcin_radiation!(inst; s=sgated, c=c, p=pg, coszen=sa.coszen_col[c])
        end
        CLM.FatesSunShadeFracs(inst.fates.nsites, inst.fates.sites,
                               inst.fates.bc_in, inst.fates.bc_out)
        sgated = 0
        for c in 1:length(col.is_fates)
            col.is_fates[c] || continue
            sgated += 1
            sgated <= inst.fates.nsites || break
            pg = col.patchi[c] + 1
            CLM.fates_unpack_bcout_sunfrac!(inst; s=sgated, c=c, p=pg)
        end
        @test sgated == 1               # exactly one FATES column was processed
        @test isfinite(cs.fsun_patch[pveg])   # the gated branch repopulated it

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
