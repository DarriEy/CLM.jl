using Test
using CLM

# ============================================================================
# Crop tails vs Fortran ground truth generated on the replacement box
# (ctsm5.3.012, use_cn + use_crop + use_fertilizer, use_fun = .false.,
#  use_flexiblecn = .true.; docs/CROP_LIFECYCLE_MAP.md §3i-§3k,
#  docs/CROP_PARITY.md § "Crop tails (2026-10-28)").
#
# Every number below is Fortran output (or the Fortran input of the same call),
# extracted at full precision by
# scripts/validation/crop_parity/extract_crop_tails_fixture.jl. The full per-step
# comparison over all ~3 300 dumped steps is scripts/fortran_parity_crop_tails.jl.
#
#  1. grain-fill allocation (the recursive astem/aleaf decay) -- chained 24 steps
#     from the leaf-emergence -> grain-fill transition, Julia carrying its OWN
#     aleaf/astem, for corn / spring wheat / soybean / tropical corn
#  2. grain C flux: calc_allometry! + FlexibleCN calc_plant_cn_alloc_flexiblecn!
#  3. harvest: cn_offset_litterfall! crop branch (grain -> food / seed, C and N)
#     -- this branch was MISSING from the port until this change
#  4. n_soyfix! against the real CNSoyfix driven in-model over a branch sweep
#  5. vernalization! on every winter-wheat call, plus the call gate (the port's
#     `isfinite` gate admitted SPVAL and reset cumvd every step)
#  6. soil mineral-N update carries fert/soyfix N into NH4 under use_crop
#  7. itype 16 (c3_irrigated) is a GENERIC crop (npcropmin is 0-based)
# ============================================================================

include(joinpath(@__DIR__, "fixtures", "crop_tails_fortran_ref.jl"))

const _P = CT_PARAMS
_npft() = length(_P.arooti)

function _ct_pfta()
    CLM.PftConAllocation{Float64}(
        woody = copy(_P.woody), froot_leaf = copy(_P.froot_leaf), croot_stem = copy(_P.croot_stem),
        stem_leaf = copy(_P.stem_leaf), flivewd = copy(_P.flivewd), leafcn = copy(_P.leafcn),
        frootcn = copy(_P.frootcn), livewdcn = copy(_P.livewdcn), deadwdcn = copy(_P.deadwdcn),
        graincn = copy(_P.graincn), grperc = copy(_P.grperc),
        arooti = copy(_P.arooti), arootf = copy(_P.arootf), fleafi = copy(_P.fleafi),
        bfact = copy(_P.bfact), astemf = copy(_P.astemf), declfact = copy(_P.declfact),
        allconss = copy(_P.allconss), aleaff = copy(_P.aleaff), allconsl = copy(_P.allconsl))
end

function _ct_patch(itypes)
    pt = CLM.PatchData()
    pt.itype = collect(itypes); pt.column = collect(1:length(itypes))
    pt.gridcell = ones(Int, length(itypes)); pt.wtcol = ones(length(itypes))
    return pt
end

@testset "Crop tails vs ctsm5.3.012 (grain fill, harvest, soyfix, vernalization)" begin

    # ------------------------------------------------------------------ 1
    @testset "grain-fill allocation, chained ($(r.itype))" for r in CT_GRAINFILL
        pfta = _ct_pfta(); pt = _ct_patch([r.itype])
        cr = CLM.CropData(); CLM.crop_init!(cr, 1); cr.croplive_patch .= true
        vs = CLM.CNVegStateData(); CLM.cnveg_state_init!(vs, 1, 1)
        vs.aleaf_patch[1] = r.prior.aleaf; vs.astem_patch[1] = r.prior.astem
        vs.aleafi_patch[1] = r.prior.aleafi; vs.astemi_patch[1] = r.prior.astemi
        vs.aroot_patch[1] = r.prior.aroot; vs.arepr_patch[1, end] = r.prior.arepr
        ngf = 0
        for k in eachindex(r.nsteps)
            cr.hui_patch[1] = r.hui[k]; cr.gddtsoi_patch[1] = r.gddtsoi[k]
            vs.huileaf_patch[1] = r.huileaf[k]; vs.huigrain_patch[1] = r.huigrain[k]
            vs.gddmaturity_patch[1] = r.gddmaturity[k]; vs.peaklai_patch[1] = r.peaklai[k]
            r.hui[k] >= r.huigrain[k] && (ngf += 1)
            CLM.calc_crop_allocation_fractions!(trues(1), 1:1, pfta, pt, cr, vs)
            # Julia carries its own aleaf/astem into the next step (the recursion)
            @test vs.aleaf_patch[1] ≈ r.aleaf[k] rtol = 1e-12 atol = 1e-15
            @test vs.astem_patch[1] ≈ r.astem[k] rtol = 1e-12 atol = 1e-15
            @test vs.aroot_patch[1] ≈ r.aroot[k] rtol = 1e-12 atol = 1e-15
            @test vs.arepr_patch[1, end] ≈ r.arepr[k] rtol = 1e-12 atol = 1e-15
        end
        @test ngf == length(r.nsteps)          # every step really is in the grain-fill branch
        # astem never grows in grain fill and grain allocation never shrinks
        @test r.astem[end] <= r.astem[1]
        @test r.arepr[end] >= r.arepr[1] >= 0
    end

    @testset "grain-fill fixture is non-vacuous" begin
        # the recursion is live (strict astem decay) in some chains, and the
        # max(astemf, ...) floor is hit in others (spring wheat 10 days in)
        @test count(r -> r.astem[end] < r.astem[1], CT_GRAINFILL) >= 3
        @test any(r -> r.astem[end] == _P.astemf[r.itype + 1], CT_GRAINFILL)
        @test maximum(r -> r.arepr[end], CT_GRAINFILL) > 0.3
    end

    # ------------------------------------------------------------------ 2
    @testset "grain C flux (allometry + FlexibleCN allocation)" begin
        pfnc = CLM.PftConNutrientCompetition{Float64}(
            woody = copy(_P.woody), froot_leaf = copy(_P.froot_leaf), croot_stem = copy(_P.croot_stem),
            stem_leaf = copy(_P.stem_leaf), flivewd = copy(_P.flivewd), leafcn = copy(_P.leafcn),
            frootcn = copy(_P.frootcn), livewdcn = copy(_P.livewdcn), deadwdcn = copy(_P.deadwdcn),
            fcur = copy(_P.fcur), graincn = copy(_P.graincn), grperc = copy(_P.grperc),
            grpnow = copy(_P.grpnow), fleafcn = copy(_P.fleafcn), ffrootcn = copy(_P.ffrootcn),
            fstemcn = copy(_P.fstemcn), astemf = copy(_P.astemf), season_decid = copy(_P.season_decid),
            stress_decid = copy(_P.stress_decid), evergreen = copy(_P.evergreen))
        sp = CLM.CNSharedParamsData{Float64}(); sp.use_fun = false; sp.use_flexiblecn = true
        @test length(CT_CNALLOC) >= 10
        for r in CT_CNALLOC
            pt = _ct_patch([r.itype])
            cr = CLM.CropData(); CLM.crop_init!(cr, 1); cr.croplive_patch .= true
            vs = CLM.CNVegStateData(); CLM.cnveg_state_init!(vs, 1, 1)
            vs.aleaf_patch .= r.aleaf; vs.astem_patch .= r.astem; vs.aroot_patch .= r.aroot
            vs.arepr_patch[1, end] = r.arepr
            cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf, 1, 1, 1)
            cf.annsum_npp_patch .= 0.0
            CLM.calc_allometry!(trues(1), 1:1, _ct_pfta(), sp, pt, cf, vs)
            @test vs.c_allometry_patch[1] ≈ r.c_allometry rtol = 1e-13
            @test vs.n_allometry_patch[1] ≈ r.n_allometry rtol = 1e-13
            cs = CLM.CNVegCarbonStateData(); CLM.cnveg_carbon_state_init!(cs, 1, 1, 1)
            ns = CLM.CNVegNitrogenStateData(); CLM.cnveg_nitrogen_state_init!(ns, 1, 1, 1)
            nf = CLM.CNVegNitrogenFluxData(); CLM.cnveg_nitrogen_flux_init!(nf, 1, 1, 1)
            for f in (:leafc_storage_patch, :frootc_storage_patch, :livestemc_storage_patch, :livecrootc_storage_patch)
                getfield(cs, f) .= 0.0
            end
            for f in (:leafn_storage_patch, :frootn_storage_patch, :livestemn_storage_patch, :livecrootn_storage_patch)
                getfield(ns, f) .= 0.0
            end
            ns.npool_patch .= r.npool
            cf.availc_patch .= r.availc; cf.gpp_before_downreg_patch .= 1.0; cf.npp_growth_patch .= 0.0
            nf.plant_ndemand_patch .= 0.0; nf.retransn_to_npool_patch .= 0.0; nf.sminn_to_plant_fun_patch .= 0.0
            can = CLM.CanopyStateData{Float64}(); CLM.canopystate_init!(can, 1)
            can.laisun_patch .= 1.0; can.laisha_patch .= 1.0
            CLM.calc_plant_cn_alloc_flexiblecn!(trues(1), 1:1, pfnc, sp, pt, cr, can, vs, cs, cf, ns, nf;
                                                fpg_col = [1.0], dt = 3600.0, carbon_resp_opt = 0)
            @test cf.plant_calloc_patch[1] == r.plant_calloc          # = availc (non-FUN FlexibleCN)
            @test cf.cpool_to_reproductivec_patch[1, end] ≈ r.cpool_to_reprc rtol = 1e-13
            @test cf.cpool_to_reproductivec_storage_patch[1, end] ≈ r.cpool_to_reprc_stor atol = 1e-20
            @test cf.cpool_to_leafc_patch[1] ≈ r.cpool_to_leafc rtol = 1e-13
            @test cf.cpool_to_livestemc_patch[1] ≈ r.cpool_to_livestemc rtol = 1e-13
            @test cf.cpool_to_frootc_patch[1] ≈ r.cpool_to_frootc rtol = 1e-13
        end
    end

    # ------------------------------------------------------------------ 3
    @testset "harvest: grain to food / seed ($(r.itype) n=$(r.n))" for r in CT_HARVEST
        pt = _ct_patch([r.itype])
        pf = CLM.PftConPhenology{Float64}(lflitcn = copy(_P.lflitcn), leafcn = copy(_P.leafcn),
                 frootcn = copy(_P.frootcn), biofuel_harvfrac = copy(_P.biofuel_harvfrac),
                 repr_structure_harvfrac = zeros(_npft(), 1))
        ps = CLM.PhenologyState(); ps.dt = 3600.0
        vs = CLM.CNVegStateData(); CLM.cnveg_state_init!(vs, 1, 1)
        vs.offset_flag_patch .= 1.0; vs.offset_counter_patch .= r.offset_counter
        cr = CLM.CropData(); CLM.crop_init!(cr, 1); cr.harvest_count .= r.harvest_count
        cs = CLM.CNVegCarbonStateData(); CLM.cnveg_carbon_state_init!(cs, 1, 1, 1)
        ns = CLM.CNVegNitrogenStateData(); CLM.cnveg_nitrogen_state_init!(ns, 1, 1, 1)
        cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf, 1, 1, 1)
        nf = CLM.CNVegNitrogenFluxData(); CLM.cnveg_nitrogen_flux_init!(nf, 1, 1, 1)
        for a in (cf.leafc_to_litter_patch, cf.frootc_to_litter_patch, cf.livestemc_to_litter_patch,
                  cf.prev_leafc_to_litter_patch, cf.prev_frootc_to_litter_patch, cf.leafc_to_litter_fun_patch,
                  nf.leafn_to_litter_patch, nf.leafn_to_retransn_patch, nf.frootn_to_litter_patch,
                  nf.livestemn_to_litter_patch, cf.repr_grainc_to_food_thisyr_patch,
                  cf.repr_grainc_to_seed_thisyr_patch, nf.repr_grainn_to_food_thisyr_patch,
                  nf.repr_grainn_to_seed_thisyr_patch)
            a .= 0.0
        end
        cs.leafc_patch .= r.leafc; cs.frootc_patch .= r.frootc; cs.livestemc_patch .= r.livestemc
        cs.reproductivec_patch[1, end] = r.reprc; cs.cropseedc_deficit_patch .= r.cropseedc_deficit
        ns.leafn_patch .= r.leafn; ns.frootn_patch .= r.frootn; ns.livestemn_patch .= r.livestemn
        ns.reproductiven_patch[1, end] = r.reprn; ns.cropseedn_deficit_patch .= r.cropseedn_deficit
        cf.cpool_to_leafc_patch .= r.cpool_to_leafc; cf.cpool_to_livestemc_patch .= r.cpool_to_livestemc
        cf.cpool_to_frootc_patch .= r.cpool_to_frootc; cf.cpool_to_reproductivec_patch[1, end] = r.cpool_to_reprc
        nf.npool_to_reproductiven_patch[1, end] = r.npool_to_reprn
        CLM.cn_offset_litterfall!(ps, trues(1), pf, vs, cs, ns, cf, nf, cr, pt;
                                  mask_pcropp = trues(1), use_fun = false, CNratio_floating = true)
        @test r.repr_grainc_to_food > 0                       # a real harvest, not a vacuous zero
        @test cf.repr_grainc_to_food_patch[1, end] ≈ r.repr_grainc_to_food rtol = 1e-13
        @test cf.repr_grainc_to_seed_patch[1, end] ≈ r.repr_grainc_to_seed rtol = 1e-13
        @test nf.repr_grainn_to_food_patch[1, end] ≈ r.repr_grainn_to_food rtol = 1e-13
        @test nf.repr_grainn_to_seed_patch[1, end] ≈ r.repr_grainn_to_seed rtol = 1e-13
        @test cf.leafc_to_litter_patch[1] ≈ r.leafc_to_litter rtol = 1e-13
        @test cf.frootc_to_litter_patch[1] ≈ r.frootc_to_litter rtol = 1e-13
        @test cf.livestemc_to_litter_patch[1] ≈ r.livestemc_to_litter rtol = 1e-13
        @test nf.livestemn_to_litter_patch[1] ≈ r.livestemn_to_litter rtol = 1e-13
        @test nf.leafn_to_litter_patch[1] ≈ r.leafn_to_litter rtol = 1e-13
        @test nf.leafn_to_retransn_patch[1] ≈ r.leafn_to_retransn rtol = 1e-12 atol = 1e-20
        @test nf.frootn_to_litter_patch[1] ≈ r.frootn_to_litter rtol = 1e-13
    end

    # ------------------------------------------------------------------ 4
    @testset "n_soyfix! vs in-model CNSoyfix sweep" begin
        S = CT_SOYFIX
        np = length(S.itype); nc = np
        ncase = size(S.fpg, 2)
        nz = 0; nfx = 0
        for k in 1:ncase
            pt = _ct_patch(S.itype); pt.wtcol = copy(S.wtcol)
            cr = CLM.CropData(); CLM.crop_init!(cr, np); cr.croplive_patch .= true
            cr.hui_patch .= S.hui[:, k]
            vs = CLM.CNVegStateData(); CLM.cnveg_state_init!(vs, np, nc); vs.gddmaturity_patch .= S.gddmaturity
            nf = CLM.CNVegNitrogenFluxData(); CLM.cnveg_nitrogen_flux_init!(nf, np, nc, 1)
            nf.plant_ndemand_patch .= S.plant_ndemand
            sbs = CLM.SoilBiogeochemStateData(); CLM.soil_bgc_state_init!(sbs, nc, np, 1, 1); sbs.fpg_col .= S.fpg[:, k]
            sbns = CLM.SoilBiogeochemNitrogenStateData(); CLM.soil_bgc_nitrogen_state_init!(sbns, nc, 1, 1, 1)
            sbns.sminn_col .= S.sminn[:, k]
            sbnf = CLM.SoilBiogeochemNitrogenFluxData(); CLM.soil_bgc_nitrogen_flux_init!(sbnf, nc, 1, 1, 1)
            wd = CLM.WaterDiagnosticBulkData{Float64}(); wd.wf_col = copy(S.wf[:, k])
            CLM.n_soyfix!(sbnf, nf, sbs, sbns, vs, cr, wd;
                mask_soilc = trues(nc), bounds = 1:nc, mask_soilp = trues(np), bounds_p = 1:np,
                patch = pt, ntmp_soybean = CLM.ntmp_soybean, nirrig_tmp_soybean = CLM.nirrig_tmp_soybean,
                ntrp_soybean = CLM.ntrp_soybean, nirrig_trp_soybean = CLM.nirrig_trp_soybean)
            for p in 1:np
                @test nf.soyfixn_patch[p] ≈ S.soyfixn[p, k] rtol = 1e-14 atol = 1e-25
                S.soyfixn[p, k] > 0 && (nz += 1)
                S.soyfixn[p, k] > 0 && S.fpg[p, k] < 1 && (nfx += 1)
            end
            for c in 1:nc
                @test sbnf.soyfixn_to_sminn_col[c] ≈ S.soyfixn_to_sminn[c, k] rtol = 1e-14 atol = 1e-25
            end
        end
        @test nz >= 20      # the Fortran fixation is non-zero across many cases
        @test nfx == nz     # and only ever behind the fpg < 1 gate
    end

    # ------------------------------------------------------------------ 5
    @testset "vernalization! on every winter-wheat call" begin
        @test length(CT_VERN) >= 20
        @test any(r -> r.vf == 1.0, CT_VERN) || maximum(r -> r.vf, CT_VERN) > 0.95   # reaches (near) full vernalization
        @test any(r -> r.cumvd > r.cumvd0, CT_VERN)                                    # vernalization accumulates
        @test any(r -> r.t_ref2m_min <= CLM.TFRZ - 3.0, CT_VERN)                       # hardening branch (frozen) exercised
        @test any(r -> r.t_ref2m < CLM.TFRZ && r.snow_depth > 0, CT_VERN)              # snow-insulated crown temperature
        pstate = CLM.PhenologyState()
        for r in CT_VERN
            pt = _ct_patch([CLM.nwwheat])
            can = CLM.CanopyStateData{Float64}(); CLM.canopystate_init!(can, 1); can.tlai_patch .= r.tlai
            tm = CLM.TemperatureData{Float64}()
            tm.t_ref2m_patch = [r.t_ref2m]; tm.t_ref2m_min_patch = [r.t_ref2m_min]; tm.t_ref2m_max_patch = [r.t_ref2m_max]
            wd = CLM.WaterDiagnosticBulkData{Float64}(); wd.snow_depth_col = [r.snow_depth]
            vs = CLM.CNVegStateData(); CLM.cnveg_state_init!(vs, 1, 1)
            vs.cumvd_patch .= r.cumvd0; vs.hdidx_patch .= r.hdidx0
            cr = CLM.CropData(); CLM.crop_init!(cr, 1); cr.vf_patch .= r.vf0
            @test CLM._vernalization_gate(r.t_ref2m_min, r.vf0, CLM.nwwheat)
            CLM.vernalization!(1, pstate, can, tm, wd, vs, cr, pt)
            @test cr.vf_patch[1] ≈ r.vf rtol = 1e-14
            @test vs.cumvd_patch[1] ≈ r.cumvd rtol = 1e-14 atol = 1e-15
            @test vs.hdidx_patch[1] ≈ r.hdidx rtol = 1e-14 atol = 1e-15
        end
        # the gate: SPVAL daily minimum (every step but the end-of-day one) must NOT
        # call vernalization -- `isfinite(1e36)` is true, which is what was wrong
        @test !CLM._vernalization_gate(CLM.SPVAL, 0.9, CLM.nwwheat)
        @test !CLM._vernalization_gate(NaN, 0.9, CLM.nwwheat)
        @test !CLM._vernalization_gate(270.0, 1.0, CLM.nwwheat)
        @test !CLM._vernalization_gate(270.0, 0.9, CLM.ntmp_soybean)
        @test CLM._vernalization_gate(270.0, 0.9, CLM.nirrig_wwheat)
    end

    # ------------------------------------------------------------------ 6
    @testset "soil N update carries crop fert + soyfix N into NH4" begin
        nc = 1; nl = 2
        ns = CLM.SoilBiogeochemNitrogenStateData(); CLM.soil_bgc_nitrogen_state_init!(ns, nc, 1, nl, 1)
        nf = CLM.SoilBiogeochemNitrogenFluxData(); CLM.soil_bgc_nitrogen_flux_init!(nf, nc, nl, 1, 1)
        st = CLM.SoilBiogeochemStateData(); CLM.soil_bgc_state_init!(st, nc, 1, nl, 1)
        for a in (ns.smin_nh4_vr_col, ns.smin_no3_vr_col, ns.sminn_vr_col)
            a .= 0.0
        end
        for f in (:ndep_to_sminn_col, :nfix_to_sminn_col, :ffix_to_sminn_col)
            getfield(nf, f) .= 0.0
        end
        for f in (:gross_nmin_vr_col, :actual_immob_nh4_vr_col, :actual_immob_no3_vr_col,
                  :sminn_to_plant_fun_nh4_vr_col, :sminn_to_plant_fun_no3_vr_col,
                  :smin_nh4_to_plant_vr_col, :smin_no3_to_plant_vr_col, :f_nit_vr_col,
                  :f_denit_vr_col, :supplement_to_sminn_vr_col)
            getfield(nf, f) .= 0.0
        end
        nf.fert_to_sminn_col .= 2.0e-6; nf.soyfixn_to_sminn_col .= 3.0e-7
        st.ndep_prof_col .= [0.7 0.3]; st.nfixation_prof_col .= [0.4 0.6]
        dt = 3600.0
        CLM.soilbiogeochem_n_state_update1!(ns, nf, st; mask_bgc_soilc = trues(1), bounds_col = 1:1,
                                            nlevdecomp = nl, dt = dt, use_fun = false, use_crop = true)
        @test ns.smin_nh4_vr_col[1, 1] ≈ 2.0e-6 * dt * 0.7 + 3.0e-7 * dt * 0.4
        @test ns.smin_nh4_vr_col[1, 2] ≈ 2.0e-6 * dt * 0.3 + 3.0e-7 * dt * 0.6
        # and the default (no crops) path is unchanged
        ns.smin_nh4_vr_col .= 0.0
        CLM.soilbiogeochem_n_state_update1!(ns, nf, st; mask_bgc_soilc = trues(1), bounds_col = 1:1,
                                            nlevdecomp = nl, dt = dt, use_fun = false)
        @test all(ns.smin_nh4_vr_col .== 0.0)
    end

    # ------------------------------------------------------------------ 7
    @testset "itype 16 (c3_irrigated) is a generic crop in calc_allometry!" begin
        pfta = _ct_pfta(); sp = CLM.CNSharedParamsData{Float64}(); sp.use_fun = false
        pt = _ct_patch([16, 17])
        vs = CLM.CNVegStateData(); CLM.cnveg_state_init!(vs, 2, 2)
        vs.aleaf_patch .= 0.5; vs.astem_patch .= 0.3; vs.aroot_patch .= 0.2; vs.arepr_patch .= 0.0
        cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf, 2, 2, 1); cf.annsum_npp_patch .= 0.0
        CLM.calc_allometry!(trues(2), 1:2, pfta, sp, pt, cf, vs)
        f1 = _P.froot_leaf[17]          # 1-based row for itype 16
        g1 = _P.grperc[17]
        @test vs.c_allometry_patch[1] ≈ 1 + g1 + f1 + f1 * g1                   # Fortran generic branch
        @test vs.c_allometry_patch[2] ≈ 1.25 * (1 + 0.2 / 0.5 + (0.3 / 0.5) * (1 + _P.croot_stem[18]))
    end
end
