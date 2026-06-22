# test_fates_edpftvarcon.jl
# Tests for FATES EDPftvarcon (Tier F, Batch 3) — the FATES PFT-indexed
# parameter container: construct/allocate for a small synthetic nPFT, check the
# registered fates_* name set, the numrad 2-D assembly, derived relationships,
# and a FatesCheckParams consistency throw.

using Test
using CLM

@testset "FATES EDPftvarcon" begin

    npft = 2          # small synthetic PFT count
    nleafage = 1      # leaf-age classes
    norgan = 4        # hydraulic organs (leaf, stem, transp root, abs root)
    nhlmpft = 2       # HLM PFT count for hlm_pft_map

    # -----------------------------------------------------------------------
    @testset "lower-bound constants" begin
        @test CLM.lower_bound_pft == 1
        @test CLM.lower_bound_general == 1
    end

    # -----------------------------------------------------------------------
    @testset "Register! lists expected fates_* names" begin
        ev = CLM.EDPftvarcon_type()
        fp = CLM.fates_parameters_type()
        CLM.Register!(ev, fp)

        registered = Set(fp.parameters[i].name for i in 1:fp.num_parameters)

        # representative names across every registration group, verbatim
        expected = [
            "fates_mort_freezetol",
            "fates_recruit_height_min",
            "fates_leaf_c3psn",
            "fates_rad_leaf_xl",
            "fates_frag_leaf_flab",
            "fates_frag_fnrt_flig",
            "fates_phen_flush_fraction",
            "fates_cnp_eca_km_nh4",
            "fates_cnp_prescribed_nuptake",
            "fates_landuse_luc_pprod10",
            "fates_dev_arbitrary_pft",
            "fates_hlm_pft_map",             # 2-D PFT x HLM-PFT
            "fates_rad_leaf_rhovis",         # numrad (assembled into rhol)
            "fates_rad_stem_taunir",
            "fates_leaf_vcmax25top",         # 2-D PFT x leaf-age
            "fates_hydro_kmax_node",         # 2-D PFT x organ
            "fates_hydro_vg_alpha_node",
        ]
        for nm in expected
            @test nm in registered
        end

        # count = 107 1-D PFT names + 1 (hlm_pft_map) + 8 numrad + 12 hydr-organ
        # + 1 leaf-age
        @test fp.num_parameters ==
              length(CLM._edpft_pft_param_names) + 1 +
              length(CLM._edpft_numrad_param_names) +
              length(CLM._edpft_hydr_organ_param_names) + 1
    end

    # -----------------------------------------------------------------------
    # Build a fully-populated registry for a synthetic nPFT, run Receive!, and
    # check shapes + the numrad vis/nir assembly + a 2-D round-trip.
    @testset "Receive! shapes + numrad/2-D assembly" begin
        ev = CLM.EDPftvarcon_type()
        fp = CLM.fates_parameters_type()
        CLM.Register!(ev, fp)

        # populate every registered parameter with finite data of the right shape
        for i in 1:fp.num_parameters
            par = fp.parameters[i]
            if par.dimension_shape == CLM.dimension_shape_1d
                par.dimension_sizes[1] = npft
                CLM.SetData1D!(fp, i, fill(0.1, npft))
            elseif par.dimension_shape == CLM.dimension_shape_2d
                # second dim depends on which 2-D param this is
                ncol = if par.name == "fates_leaf_vcmax25top"
                    nleafage
                elseif par.name == "fates_hlm_pft_map"
                    nhlmpft
                else
                    norgan
                end
                CLM.SetData2D!(fp, i, fill(0.2, npft, ncol))
            end
        end

        # override the two numrad leaf-reflectance vectors so we can verify the
        # vis/nir column assembly distinctly
        let i = CLM.FindIndex(fp, "fates_rad_leaf_rhovis")
            CLM.SetData1D!(fp, i, Float64[0.11, 0.12])
        end
        let i = CLM.FindIndex(fp, "fates_rad_leaf_rhonir")
            CLM.SetData1D!(fp, i, Float64[0.41, 0.42])
        end

        CLM.Receive!(ev, fp)

        # 1-D arrays sized to npft
        @test length(ev.freezetol) == npft
        @test length(ev.c3psn) == npft
        @test length(ev.prescribed_nuptake) == npft
        @test length(ev.landusechange_pprod10) == npft

        # numrad assembled into PFT x num_swb matrices
        @test size(ev.rhol) == (npft, CLM.num_swb)
        @test size(ev.taus) == (npft, CLM.num_swb)
        @test ev.rhol[:, CLM.ivis] == [0.11, 0.12]
        @test ev.rhol[:, CLM.inir] == [0.41, 0.42]

        # 2-D leaf-age + organ + hlm-pft-map shapes
        @test size(ev.vcmax25top) == (npft, nleafage)
        @test size(ev.hydr_kmax_node) == (npft, norgan)
        @test size(ev.hlm_pft_map) == (npft, nhlmpft)
    end

    # -----------------------------------------------------------------------
    # GetDecompyFrac maps the litter pool indices to the correct leaf/fnrt
    # fractions, and throws on bad indices.
    @testset "GetDecompyFrac mapping + throws" begin
        ev = CLM.EDPftvarcon_type()
        ev.lf_flab = [0.25, 0.30]
        ev.lf_fcel = [0.45, 0.50]
        ev.lf_flig = [0.30, 0.20]
        ev.fr_flab = [0.15, 0.10]
        ev.fr_fcel = [0.55, 0.60]
        ev.fr_flig = [0.30, 0.30]
        CLM.EDPftvarcon_inst[] = ev

        @test CLM.GetDecompyFrac(1, CLM.leaf_organ, CLM.ilabile)   == 0.25
        @test CLM.GetDecompyFrac(2, CLM.leaf_organ, CLM.icellulose) == 0.50
        @test CLM.GetDecompyFrac(1, CLM.leaf_organ, CLM.ilignin)   == 0.30
        @test CLM.GetDecompyFrac(2, CLM.fnrt_organ, CLM.ilabile)   == 0.10
        @test CLM.GetDecompyFrac(1, CLM.fnrt_organ, CLM.icellulose) == 0.55

        # unknown decomposability pool index
        @test_throws ErrorException CLM.GetDecompyFrac(1, CLM.leaf_organ, 99)
        # unknown organ index
        @test_throws ErrorException CLM.GetDecompyFrac(1, CLM.sapw_organ, CLM.ilabile)
    end

    # -----------------------------------------------------------------------
    # FatesCheckParams: a deliberately bad parameter triggers an endrun (error).
    @testset "FatesCheckParams consistency throw" begin
        # set up the global EDParams so the early model-switch checks pass
        edp = CLM.ed_params()
        edp.radiation_model = CLM.norman_solver
        edp.regeneration_model = CLM.default_regeneration
        edp.dayl_switch = CLM.itrue
        edp.logging_mechanical_frac = 0.0
        edp.logging_collateral_frac = 0.0
        edp.logging_direct_frac = 0.0

        # carbon-only PARTEH mode skips the CNP nutrient checks
        CLM.hlm_parteh_mode[] = CLM.prt_carbon_allom_hyp
        CLM.hlm_use_fixed_biogeog[] = CLM.ifalse
        CLM.hlm_use_inventory_init[] = CLM.itrue   # so initd==0 does not abort
        CLM.hlm_use_nocomp[] = CLM.ifalse
        CLM.hlm_use_sp[] = CLM.ifalse

        # build a minimal, valid EDPftvarcon_inst then break ONE thing: xl out
        # of the allowed [-0.4, 0.6] range.
        ev = CLM.EDPftvarcon_type()
        np = 1
        ev.freezetol = [-10.0]
        ev.initd = [0.0]
        ev.xl = [0.9]                      # <-- invalid (> 0.6)
        ev.c3psn = [1.0]
        ev.allom_frbstor_repro = [0.5]
        ev.phenflush_fraction = [0.5]
        ev.seed_dispersal_pdf_scale = fill(CLM.fates_check_param_set + 1.0, np)
        ev.seed_dispersal_pdf_shape = fill(CLM.fates_check_param_set + 1.0, np)
        ev.seed_dispersal_max_dist  = fill(CLM.fates_check_param_set + 1.0, np)
        ev.seed_dispersal_fraction  = fill(CLM.fates_check_param_set + 1.0, np)
        ev.mort_ip_age_senescence  = fill(CLM.fates_check_param_set + 1.0, np)
        ev.mort_r_age_senescence   = fill(CLM.fates_check_param_set + 1.0, np)
        ev.mort_ip_size_senescence = [1.0]
        ev.mort_r_size_senescence  = [0.05]
        ev.prescribed_nuptake = [0.0]
        ev.prescribed_puptake = [0.0]
        ev.eca_km_nh4 = [1.0]
        ev.eca_km_no3 = [1.0]
        ev.eca_lambda_ptase = [0.0]
        ev.eca_alpha_ptase = [0.0]
        ev.hlm_pft_map = reshape([1.0], np, 1)
        ev.maintresp_leaf_atkin2017_baserate = [1.4e-6]
        CLM.EDPftvarcon_inst[] = ev

        # prt_params must be sized for at least 1 PFT (evergreen/woody/smode +
        # nitr_stoich_p1 / organ_param_id used by the per-PFT checks)
        CLM.allocate_prt_params!(CLM.prt_params, np, 6, 1)
        CLM.prt_params.evergreen[1] = CLM.itrue       # skip phenflush deciduous check
        CLM.prt_params.woody[1] = CLM.ifalse
        CLM.prt_params.allom_smode[1] = 1
        # organ_param_id[leaf_organ] must index a valid nitr_stoich_p1 column
        CLM.prt_params.organ_param_id[CLM.leaf_organ] = 1

        # is_master=false -> no-op (returns without checking)
        @test CLM.FatesCheckParams(false) === nothing

        # is_master=true -> the bad xl triggers endrun (error)
        @test_throws ErrorException CLM.FatesCheckParams(true)

        # fixing xl lets the check pass cleanly
        ev.xl = [0.1]
        @test CLM.FatesCheckParams(true) === nothing
    end
end
