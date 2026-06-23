# test_fates_soilbgcflux.jl
# Tests for FATES Batch 11 (Tier F): FatesSoilBGCFluxMod — the FATES <-> host
# soil-biogeochem flux coupling.
#
# Strategy:
#   * Define a lightweight stub PRT object (subtype of CLM.AbstractPRTVartypes)
#     with known per-(organ,element) state and net-alloc, so we can exercise the
#     flux-coupling math without standing up the full PRT global descriptor.
#   * Build a synthetic single-site / single-patch / two-cohort system with known
#     litter fragmentation pools + cohort fine-root carbon and N/P demand drivers,
#     and a minimal bc_in / bc_out.
#   * Run:
#       - FluxIntoLitterPools  : assert litter -> soil flux MASS CONSERVATION
#         (every kg of frag flux sent into the patch must appear in the depth-
#         integrated cellulose+lignin+labile soil fluxes, after the g/m3/s ->
#         kg/m2/day back-conversion).
#       - PrepNutrientAquisitionBCs : assert the per-competitor fine-root carbon
#         (veg_rootc) aggregates correctly and num_plant_comps == #cohorts.
#       - UnPackNutrientAquisitionBCs : assert per-cohort daily N/P demand +
#         uptake unpack correctly from the host BCs and that BCs are zeroed.
#       - EffluxIntoLitterPools : assert root efflux conserves into the root-fines
#         labile pool.

using Test
using CLM

# ---------------------------------------------------------------------------
# A minimal stub PRT object. We only need GetState (organ,element)->mass and
# GetNetAlloc (organ,element)->net alloc.
# ---------------------------------------------------------------------------
mutable struct _StubPRT <: CLM.AbstractPRTVartypes
    state::Dict{Tuple{Int,Int},Float64}     # (organ, element) -> mass [kg]
    netalloc::Dict{Tuple{Int,Int},Float64}  # (organ, element) -> net alloc
end
_StubPRT() = _StubPRT(Dict{Tuple{Int,Int},Float64}(), Dict{Tuple{Int,Int},Float64}())

function CLM.GetState(this::_StubPRT, organ_id::Integer, element_id::Integer,
                      position_id::Union{Nothing,Integer}=nothing)
    return get(this.state, (Int(organ_id), Int(element_id)), 0.0)
end
function CLM.GetNetAlloc(this::_StubPRT, organ_id::Integer, element_id::Integer,
                         position_id::Union{Nothing,Integer}=nothing)
    return get(this.netalloc, (Int(organ_id), Int(element_id)), 0.0)
end

@testset "FATES Batch 11: FatesSoilBGCFluxMod" begin

    # --- Preserve and restore mutated module-global state --------------------
    old_numpft   = CLM.numpft[]
    old_numel    = CLM.num_elements[]
    old_ellist   = copy(CLM.element_list)
    old_elpos    = copy(CLM.element_pos)
    old_npscale  = CLM.fates_np_comp_scaling[]
    old_parteh   = CLM.hlm_parteh_mode[]
    old_use_ch4  = CLM.hlm_use_ch4[]
    old_nu_com   = CLM.hlm_nu_com[]
    old_decomp   = CLM.hlm_decomp[]
    old_edparams = CLM.EDParams[]
    old_edpft    = CLM.EDPftvarcon_inst[]

    try
        npft = 2
        nlevsoil = 3
        nlevdecomp = 3

        CLM.numpft[]       = npft
        CLM.num_elements[] = 1                  # carbon-only litter flux
        empty!(CLM.element_list); push!(CLM.element_list, CLM.carbon12_element)
        fill!(CLM.element_pos, 0)
        CLM.element_pos[CLM.carbon12_element] = 1

        CLM.hlm_parteh_mode[] = CLM.prt_cnp_flex_allom_hyp   # not c-only (so unpack runs)
        CLM.hlm_use_ch4[]     = CLM.ifalse
        CLM.hlm_nu_com[]      = "RD"
        CLM.hlm_decomp[]      = "CTC"                        # not MIMICS
        CLM.fates_np_comp_scaling[] = CLM.coupled_np_comp_scaling

        # --- EDParams: cwd cellulose/lignin fractions + uptake modes ---------
        ep = CLM.ed_params_type()
        ep.ED_val_cwd_fcel = 0.75
        ep.ED_val_cwd_flig = 0.25     # fcel + flig == 1 -> CWD fully partitioned
        ep.n_uptake_mode   = CLM.coupled_n_uptake
        ep.p_uptake_mode   = CLM.coupled_p_uptake
        CLM.EDParams[] = ep

        # --- EDPftvarcon: leaf-litter chem fracs + uptake kinetics -----------
        pft_con = CLM.EDPftvarcon_type()
        pft_con.lf_flab = [0.30, 0.30]
        pft_con.lf_fcel = [0.45, 0.45]
        pft_con.lf_flig = [0.25, 0.25]     # sum to 1 -> seeds fully partitioned
        pft_con.vmax_nh4 = [1.0e-7, 2.0e-7]
        pft_con.vmax_no3 = [0.5e-7, 1.0e-7]
        pft_con.vmax_p   = [3.0e-8, 4.0e-8]
        pft_con.prescribed_nuptake = [0.0, 0.0]
        pft_con.decompmicc = [100.0, 100.0]
        CLM.EDPftvarcon_inst[] = pft_con

        # --- prt_params: allocate + set woody / agb / fine-root profile ------
        pp = CLM.prt_params
        CLM.allocate_prt_params!(pp, npft, CLM.num_organ_types, 1)
        pp.woody          .= CLM.itrue
        pp.allom_agb_frac .= 0.6
        # Fine-root profile: 1-parameter exponential (type 2), scale a.
        pp.fnrt_prof_mode .= 2.0
        pp.fnrt_prof_a    .= 4.0
        pp.fnrt_prof_b    .= 0.0
        # organ_param_id maps each PRT organ -> param-file slot (use organ index).
        for io in 1:CLM.num_organ_types
            pp.organ_param_id[io] = io
        end
        # N stoichiometry (only used on the MIMICS C-only path; set sane values).
        pp.nitr_stoich_p1 .= 0.02

        # ---------------------------------------------------------------------
        # Build a site with soil layering.
        # ---------------------------------------------------------------------
        site = CLM.ed_site_type()
        site.nlevsoil = nlevsoil
        site.dz_soil  = [0.1, 0.2, 0.3]
        site.z_soil   = [0.05, 0.20, 0.45]
        # zi has a zero index (surface) -> length nlevsoil+1
        site.zi_soil  = [0.0, 0.1, 0.3, 0.6]
        site.rootfrac_scr = zeros(Float64, nlevsoil)
        site.ema_npp  = 0.0

        # ---------------------------------------------------------------------
        # Build one patch with a carbon litter pool that has KNOWN frag fluxes.
        # ---------------------------------------------------------------------
        patch = CLM.fates_patch_type()
        patch.area = CLM.area               # patch fills the whole notional area
        patch.nocomp_pft_label = 1
        patch.nitr_repro_stoich = zeros(Float64, CLM.maxpft)

        litt = CLM.litter_type()
        CLM.InitAllocate!(litt, npft, nlevsoil, CLM.carbon12_element)
        CLM.ZeroFlux!(litt)
        CLM.InitConditions!(litt, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        patch.litter = [litt]

        # Known fragmentation fluxes [kg/m2/day]:
        #   ag_cwd_frag: one value per cwd class (surface profile distributed)
        litt.ag_cwd_frag .= [0.01, 0.02, 0.0, 0.0]
        #   bg_cwd_frag: cwd x soil
        litt.bg_cwd_frag .= 0.0
        litt.bg_cwd_frag[1, 1] = 0.005
        litt.bg_cwd_frag[1, 2] = 0.003
        #   leaf fines (dcmpy): labile/cellulose/lignin (surface)
        litt.leaf_fines_frag .= [0.02, 0.03, 0.01]
        #   root fines (dcmpy x soil)
        litt.root_fines_frag .= 0.0
        litt.root_fines_frag[CLM.ilabile, 1]    = 0.004
        litt.root_fines_frag[CLM.icellulose, 2] = 0.006
        litt.root_fines_frag[CLM.ilignin, 3]    = 0.002
        #   decaying seeds (pft): partitioned by lf_flab/fcel/flig (surface)
        litt.seed_decay      .= [0.01, 0.0]
        litt.seed_germ_decay .= [0.0, 0.005]

        # Total known mass sent into the soil per day [kg/m2/day]:
        total_in = sum(litt.ag_cwd_frag) + sum(litt.bg_cwd_frag) +
                   sum(litt.leaf_fines_frag) + sum(litt.root_fines_frag) +
                   sum(litt.seed_decay) + sum(litt.seed_germ_decay)

        # ---------------------------------------------------------------------
        # Two cohorts (tallest -> shorter), each with a stub PRT.
        # ---------------------------------------------------------------------
        prt1 = _StubPRT()
        prt1.state[(CLM.fnrt_organ,  CLM.carbon12_element)] = 0.5  # kg fineroot C
        prt1.state[(CLM.leaf_organ,  CLM.carbon12_element)] = 1.0
        prt1.state[(CLM.sapw_organ,  CLM.carbon12_element)] = 2.0
        prt1.state[(CLM.struct_organ,CLM.carbon12_element)] = 3.0

        prt2 = _StubPRT()
        prt2.state[(CLM.fnrt_organ,  CLM.carbon12_element)] = 0.25

        coh1 = CLM.fates_cohort_type()
        coh1.pft = 1
        coh1.n   = 100.0
        coh1.prt = prt1
        coh1.isnew = false
        coh1.c_area = 5.0
        coh1.resp_acc_hold = 0.0

        coh2 = CLM.fates_cohort_type()
        coh2.pft = 2
        coh2.n   = 50.0
        coh2.prt = prt2
        coh2.isnew = false
        coh2.c_area = 2.0
        coh2.resp_acc_hold = 0.0

        coh1.shorter = coh2
        patch.tallest = coh1

        site.oldest_patch = patch

        # ---------------------------------------------------------------------
        # bc_in / bc_out
        # ---------------------------------------------------------------------
        bc_in  = CLM.bc_in_type()
        CLM.allocate_bcin!(bc_in; npatches=1, nlevsoil=nlevsoil,
                           nlevdecomp=nlevdecomp, max_comp=4)
        bc_in.dz_sisl        = [0.1, 0.2, 0.3]
        bc_in.dz_decomp_sisl = [0.1, 0.2, 0.3]
        bc_in.decomp_id      = [1, 2, 3]
        bc_in.max_rooting_depth_index_col = nlevsoil

        bc_out = CLM.bc_out_type()
        CLM.allocate_bcout!(bc_out; npatches=1, nlevsoil=nlevsoil,
                            nlevdecomp=nlevdecomp, max_comp=4)

        # =====================================================================
        # 1) FluxIntoLitterPools — litter -> soil MASS CONSERVATION
        # =====================================================================
        CLM.FluxIntoLitterPools(site, bc_in, bc_out)

        # Back-convert g/m3/s -> kg/m2/day and integrate over depth:
        #   kg/m2/day = (g/m3/s) / (days_per_sec * g_per_kg / dz)
        #            => multiply by dz / (days_per_sec * g_per_kg)
        recovered = 0.0
        for id in 1:nlevdecomp
            fcel = bc_out.litt_flux_cel_c_si[id]
            flig = bc_out.litt_flux_lig_c_si[id]
            flab = bc_out.litt_flux_lab_c_si[id]
            recovered += (fcel + flig + flab) * bc_in.dz_decomp_sisl[id] /
                         (CLM.days_per_sec * CLM.g_per_kg)
        end

        @test isfinite(recovered)
        @test recovered ≈ total_in rtol=1e-10
        # all three chemical fractions received some mass
        @test sum(bc_out.litt_flux_cel_c_si) > 0
        @test sum(bc_out.litt_flux_lig_c_si) > 0
        @test sum(bc_out.litt_flux_lab_c_si) > 0

        # =====================================================================
        # 2) PrepNutrientAquisitionBCs — per-competitor fine-root C aggregation
        # =====================================================================
        CLM.PrepNutrientAquisitionBCs(site, bc_in, bc_out)

        @test bc_out.num_plant_comps == 2
        @test bc_out.ft_index[1] == 1
        @test bc_out.ft_index[2] == 2

        # veg_rootc[icomp,id] summed over depth, back to [gC/m2] then [kgC]:
        #   veg_rootc = fnrt_c * n * rootfrac * area_inv * g_per_kg / dz
        #   => fnrt_c*n*g_per_kg*area_inv = sum_id( veg_rootc[icomp,id]*dz )
        # so total competitor fine-root C [kgC over notional area]:
        for (icomp, coh) in enumerate((coh1, coh2))
            fnrt_c = CLM.GetState(coh.prt, CLM.fnrt_organ, CLM.carbon12_element)
            integ = 0.0
            for id in 1:nlevdecomp
                integ += bc_out.veg_rootc[icomp, id] * site.dz_soil[id]
            end
            # integ is in gC/m2; convert to kgC over notional area then per plant:
            expected_density = fnrt_c * coh.n * CLM.area_inv * CLM.g_per_kg  # gC/m2 (root frac sums to 1)
            @test integ ≈ expected_density rtol=1e-8
        end

        # =====================================================================
        # 3) UnPackNutrientAquisitionBCs — per-cohort N/P demand + uptake unpack
        # =====================================================================
        # Host returns per-competitor uptake fluxes [g/m2/day] (icomp ordering =
        # tallest->shorter walk, matching PrepNutrientAquisitionBCs).
        bc_in.plant_nh4_uptake_flux[1, 1] = 0.02
        bc_in.plant_nh4_uptake_flux[2, 1] = 0.01
        bc_in.plant_no3_uptake_flux[1, 1] = 0.005
        bc_in.plant_no3_uptake_flux[2, 1] = 0.004
        bc_in.plant_p_uptake_flux[1, 1]   = 0.003
        bc_in.plant_p_uptake_flux[2, 1]   = 0.002

        sites  = [site]
        bcins  = [bc_in]
        CLM.UnPackNutrientAquisitionBCs(sites, bcins)

        edp = CLM.EDPftvarcon_inst[]
        for (icomp, coh) in enumerate((coh1, coh2))
            pft = coh.pft
            fnrt_c = CLM.GetState(coh.prt, CLM.fnrt_organ, CLM.carbon12_element)
            exp_ndem = fnrt_c * (edp.vmax_nh4[pft] + edp.vmax_no3[pft]) * CLM.sec_per_day
            @test coh.daily_n_demand ≈ exp_ndem rtol=1e-12
            @test coh.daily_p_demand ≈ fnrt_c * edp.vmax_p[pft] * CLM.sec_per_day rtol=1e-12

            exp_nh4 = bc_in.plant_nh4_uptake_flux[icomp, 1] * CLM.kg_per_g * CLM.area / coh.n
            # bc_in zeroed already -> recompute from the known inputs
        end
        # daily uptake unpack (use saved inputs since bc_in is now zeroed)
        @test coh1.daily_nh4_uptake ≈ 0.02 * CLM.kg_per_g * CLM.area / coh1.n rtol=1e-12
        @test coh2.daily_nh4_uptake ≈ 0.01 * CLM.kg_per_g * CLM.area / coh2.n rtol=1e-12
        @test coh1.daily_p_gain     ≈ 0.003 * CLM.kg_per_g * CLM.area / coh1.n rtol=1e-12

        # BCs zeroed after unpack
        @test all(bc_in.plant_nh4_uptake_flux .== 0.0)
        @test all(bc_in.plant_no3_uptake_flux .== 0.0)
        @test all(bc_in.plant_p_uptake_flux .== 0.0)

        # =====================================================================
        # 4) EffluxIntoLitterPools — root efflux conserves into root-fines labile
        # =====================================================================
        # Fresh litter so we can measure the increment cleanly.
        litt2 = CLM.litter_type()
        CLM.InitAllocate!(litt2, npft, nlevsoil, CLM.carbon12_element)
        CLM.ZeroFlux!(litt2)
        CLM.InitConditions!(litt2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        patch.litter = [litt2]

        coh1.daily_c_efflux = 1.0e-4  # kgC/plant/day
        before = sum(litt2.root_fines_frag)
        CLM.EffluxIntoLitterPools(site, patch, coh1, bc_in)
        after = sum(litt2.root_fines_frag)

        # All efflux (root frac sums to 1) -> labile root fines, area-normalized:
        expected_efflux = coh1.daily_c_efflux * coh1.n * CLM.area_inv
        @test (after - before) ≈ expected_efflux rtol=1e-10
        # only the labile dcmpy row received mass
        @test sum(litt2.root_fines_frag[CLM.ilabile, :]) ≈ expected_efflux rtol=1e-10
        @test sum(litt2.root_fines_frag[CLM.icellulose, :]) == 0.0
        @test sum(litt2.root_fines_frag[CLM.ilignin, :]) == 0.0

    finally
        # Restore module-global state.
        CLM.numpft[]       = old_numpft
        CLM.num_elements[] = old_numel
        empty!(CLM.element_list); append!(CLM.element_list, old_ellist)
        copyto!(CLM.element_pos, old_elpos)
        CLM.fates_np_comp_scaling[] = old_npscale
        CLM.hlm_parteh_mode[] = old_parteh
        CLM.hlm_use_ch4[]     = old_use_ch4
        CLM.hlm_nu_com[]      = old_nu_com
        CLM.hlm_decomp[]      = old_decomp
        CLM.EDParams[]        = old_edparams
        CLM.EDPftvarcon_inst[] = old_edpft
    end
end
