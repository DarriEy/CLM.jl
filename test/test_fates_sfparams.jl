# test_fates_sfparams.jl
# Tests for FATES SFParamsMod (Tier F, Batch 3) — the SPITFIRE fire-model
# parameter holder, its registration, retrieval, init + consistency checks.

using Test
using CLM

@testset "FATES SFParamsMod" begin

    @testset "parameter names preserved verbatim" begin
        @test CLM.SF_name_fdi_alpha          == "fates_fire_fdi_alpha"
        @test CLM.SF_name_miner_total        == "fates_fire_miner_total"
        @test CLM.SF_name_fuel_energy        == "fates_fire_fuel_energy"
        @test CLM.SF_name_part_dens          == "fates_fire_part_dens"
        @test CLM.SF_name_miner_damp         == "fates_fire_miner_damp"
        @test CLM.SF_name_max_durat          == "fates_fire_max_durat"
        @test CLM.SF_name_durat_slope        == "fates_fire_durat_slope"
        @test CLM.SF_name_drying_ratio       == "fates_fire_drying_ratio"
        @test CLM.SF_name_fire_threshold     == "fates_fire_threshold"
        @test CLM.SF_name_CWD_frac           == "fates_frag_cwd_frac"
        @test CLM.SF_name_max_decomp         == "fates_frag_maxdecomp"
        @test CLM.SF_name_SAV                == "fates_fire_SAV"
        @test CLM.SF_name_FBD                == "fates_fire_FBD"
        @test CLM.SF_name_min_moisture       == "fates_fire_min_moisture"
        @test CLM.SF_name_mid_moisture       == "fates_fire_mid_moisture"
        @test CLM.SF_name_low_moisture_Coeff == "fates_fire_low_moisture_Coeff"
        @test CLM.SF_name_low_moisture_Slope == "fates_fire_low_moisture_Slope"
        @test CLM.SF_name_mid_moisture_Coeff == "fates_fire_mid_moisture_Coeff"
        @test CLM.SF_name_mid_moisture_Slope == "fates_fire_mid_moisture_Slope"
    end

    @testset "SpitFireParamsInit! defaults + array shapes" begin
        CLM.SpitFireParamsInit!()
        p = CLM.sf_params()

        # scalars reset to the unset sentinel
        @test p.SF_val_fdi_alpha      == CLM.fates_unset_r8
        @test p.SF_val_miner_total    == CLM.fates_unset_r8
        @test p.SF_val_fuel_energy    == CLM.fates_unset_r8
        @test p.SF_val_part_dens      == CLM.fates_unset_r8
        @test p.SF_val_miner_damp     == CLM.fates_unset_r8
        @test p.SF_val_max_durat      == CLM.fates_unset_r8
        @test p.SF_val_durat_slope    == CLM.fates_unset_r8
        @test p.SF_val_drying_ratio   == CLM.fates_unset_r8
        @test p.SF_val_fire_threshold == CLM.fates_unset_r8

        # per-CWD array sized ncwd
        @test length(p.SF_val_CWD_frac) == CLM.ncwd
        @test all(p.SF_val_CWD_frac .== CLM.fates_unset_r8)

        # per-fuel-class arrays sized num_fuel_classes
        nfsc_arrays = (
            p.SF_val_max_decomp, p.SF_val_SAV, p.SF_val_FBD,
            p.SF_val_min_moisture, p.SF_val_mid_moisture,
            p.SF_val_low_moisture_Coeff, p.SF_val_low_moisture_Slope,
            p.SF_val_mid_moisture_Coeff, p.SF_val_mid_moisture_Slope,
        )
        for a in nfsc_arrays
            @test length(a) == CLM.num_fuel_classes
            @test all(a .== CLM.fates_unset_r8)
        end
    end

    @testset "SpitFireRegisterParams! registers expected fates_fire_* names" begin
        fp = CLM.fates_parameters_type()
        CLM.SpitFireRegisterParams!(fp)

        # 9 scalar + 1 per-CWD + 9 per-fuel-class = 19 registrations
        @test fp.num_parameters == 19

        registered = [strip(fp.parameters[i].name) for i in 1:fp.num_parameters]
        expected = [
            CLM.SF_name_fdi_alpha, CLM.SF_name_miner_total, CLM.SF_name_fuel_energy,
            CLM.SF_name_part_dens, CLM.SF_name_miner_damp, CLM.SF_name_max_durat,
            CLM.SF_name_durat_slope, CLM.SF_name_drying_ratio, CLM.SF_name_fire_threshold,
            CLM.SF_name_CWD_frac,
            CLM.SF_name_SAV, CLM.SF_name_FBD, CLM.SF_name_min_moisture,
            CLM.SF_name_mid_moisture, CLM.SF_name_low_moisture_Coeff,
            CLM.SF_name_low_moisture_Slope, CLM.SF_name_mid_moisture_Coeff,
            CLM.SF_name_mid_moisture_Slope, CLM.SF_name_max_decomp,
        ]
        for nm in expected
            @test nm in registered
        end

        # The CWD frac registers under the NCWD dimension; the fuel-class arrays
        # under the litterclass (fsc) dimension.
        icwd = CLM.FindIndex(fp, CLM.SF_name_CWD_frac)
        @test fp.parameters[icwd].dimension_shape == CLM.dimension_shape_1d
        @test fp.parameters[icwd].dimension_names[1] == CLM.dimension_name_cwd
        isav = CLM.FindIndex(fp, CLM.SF_name_SAV)
        @test fp.parameters[isav].dimension_shape == CLM.dimension_shape_1d
        @test fp.parameters[isav].dimension_names[1] == CLM.dimension_name_fsc
        ifdi = CLM.FindIndex(fp, CLM.SF_name_fdi_alpha)
        @test fp.parameters[ifdi].dimension_shape == CLM.dimension_shape_scalar
    end

    @testset "register -> set -> receive round-trip" begin
        fp = CLM.fates_parameters_type()
        CLM.SpitFireRegisterParams!(fp)

        # Populate the storage the way the (stubbed) host file read would.
        # Scalars:
        scalar_vals = Dict(
            CLM.SF_name_fdi_alpha      => 0.301,
            CLM.SF_name_miner_total    => 0.055,
            CLM.SF_name_fuel_energy    => 18000.0,
            CLM.SF_name_part_dens      => 513.0,
            CLM.SF_name_miner_damp     => 0.41739,
            CLM.SF_name_max_durat      => 240.0,
            CLM.SF_name_durat_slope    => -11.06,
            CLM.SF_name_drying_ratio   => 66000.0,
            CLM.SF_name_fire_threshold => 50.0,
        )
        for (nm, v) in scalar_vals
            i = CLM.FindIndex(fp, nm)
            CLM.SetDataScalar!(fp, i, v)
        end

        # per-CWD (ncwd) — must sum to unity for the consistency check below.
        cwd_frac = [0.045, 0.075, 0.21, 0.67]
        @test length(cwd_frac) == CLM.ncwd
        icwd = CLM.FindIndex(fp, CLM.SF_name_CWD_frac)
        fp.parameters[icwd].dimension_sizes[1] = CLM.ncwd
        CLM.SetData1D!(fp, icwd, cwd_frac)

        # per-fuel-class (num_fuel_classes)
        nfsc = CLM.num_fuel_classes
        nfsc_vals = Dict(
            CLM.SF_name_SAV                => collect(range(13.0,  step=1.0,  length=nfsc)),
            CLM.SF_name_FBD                => collect(range(15.0,  step=2.0,  length=nfsc)),
            CLM.SF_name_min_moisture       => fill(0.18, nfsc),
            CLM.SF_name_mid_moisture       => fill(0.5,  nfsc),
            CLM.SF_name_low_moisture_Coeff => fill(1.18, nfsc),
            CLM.SF_name_low_moisture_Slope => fill(0.0096, nfsc),
            CLM.SF_name_mid_moisture_Coeff => fill(1.22, nfsc),
            CLM.SF_name_mid_moisture_Slope => fill(0.0027, nfsc),
            CLM.SF_name_max_decomp         => fill(0.5,  nfsc),
        )
        for (nm, v) in nfsc_vals
            i = CLM.FindIndex(fp, nm)
            fp.parameters[i].dimension_sizes[1] = nfsc
            CLM.SetData1D!(fp, i, v)
        end

        # Retrieve into the module-global SFParams instance.
        CLM.SpitFireReceiveParams!(fp)
        p = CLM.sf_params()

        @test p.SF_val_fdi_alpha      ≈ 0.301
        @test p.SF_val_fuel_energy    ≈ 18000.0
        @test p.SF_val_fire_threshold ≈ 50.0
        @test p.SF_val_CWD_frac       ≈ cwd_frac
        @test p.SF_val_SAV            ≈ nfsc_vals[CLM.SF_name_SAV]
        @test p.SF_val_max_decomp     ≈ nfsc_vals[CLM.SF_name_max_decomp]

        # Consistency checks: CWD already sums to unity, threshold valid,
        # decomp rates ≥ 0 → should pass and leave the (already-unity) CWD
        # fractions essentially unchanged.
        CLM.SpitFireCheckParams(true)
        @test sum(p.SF_val_CWD_frac) ≈ 1.0

        # is_master=false is a no-op.
        @test CLM.SpitFireCheckParams(false) === nothing
    end

    @testset "SpitFireCheckParams correction + failure modes" begin
        # Small (< 1e-5) CWD residual gets corrected into the largest pool.
        CLM.SpitFireParamsInit!()
        p = CLM.sf_params()
        p.SF_val_max_decomp .= 0.5
        p.SF_val_fire_threshold = 50.0
        p.SF_val_CWD_frac .= [0.045, 0.075, 0.21, 0.6699965]  # sums to 0.9999965
        CLM.SpitFireCheckParams(true)
        @test sum(p.SF_val_CWD_frac) ≈ 1.0
        @test p.SF_val_CWD_frac[4] ≈ 0.67  # correction lands in the largest pool

        # Negative decomposition rate -> abort.
        CLM.SpitFireParamsInit!()
        p = CLM.sf_params()
        p.SF_val_max_decomp .= 0.5
        p.SF_val_max_decomp[2] = -1.0
        p.SF_val_fire_threshold = 50.0
        p.SF_val_CWD_frac .= [0.045, 0.075, 0.21, 0.67]
        @test_throws Exception CLM.SpitFireCheckParams(true)

        # CWD fractions wildly off unity -> abort.
        CLM.SpitFireParamsInit!()
        p = CLM.sf_params()
        p.SF_val_max_decomp .= 0.5
        p.SF_val_fire_threshold = 50.0
        p.SF_val_CWD_frac .= [0.1, 0.1, 0.1, 0.1]  # sums to 0.4
        @test_throws Exception CLM.SpitFireCheckParams(true)

        # Fire threshold below the minimum -> abort.
        CLM.SpitFireParamsInit!()
        p = CLM.sf_params()
        p.SF_val_max_decomp .= 0.5
        p.SF_val_CWD_frac .= [0.045, 0.075, 0.21, 0.67]
        p.SF_val_fire_threshold = 0.0  # below min_fire_threshold
        @test_throws Exception CLM.SpitFireCheckParams(true)
    end

end
