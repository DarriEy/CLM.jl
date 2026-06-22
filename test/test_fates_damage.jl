# test_fates_damage.jl
# Tests for FATES Batch 5 (Tier F): DamageMainMod — the tree-damage module.
#   * IsItDamageTime!     — calendar-driven damage-event flag (event codes).
#   * GetDamageFrac       — damage-class transition fraction (param_derived).
#   * GetCrownReduction   — crown-loss fraction for a damage class.
#   * GetDamageMortality  — crown-loss-dependent logistic mortality.

using Test
using CLM

@testset "FATES Batch 5: DamageMainMod" begin

    numpft = 3
    nla    = 2     # nleafage
    nld    = 4     # nlevdamage

    CLM.nleafage[]   = nla
    CLM.nlevdamage[] = nld

    # --- Synthetic EDPftvarcon_inst. ---
    pft = CLM.EDPftvarcon_type()
    pft.vcmax25top     = [10.0 20.0; 30.0 40.0; 50.0 60.0]
    pft.damage_frac    = [0.0, 0.2, 0.5]
    # per-PFT damage mortality params (inflection p1, rate p2)
    pft.damage_mort_p1 = [0.5, 0.4, 0.3]
    pft.damage_mort_p2 = [2.0, 5.0, 10.0]
    CLM.EDPftvarcon_inst[] = pft

    # --- Synthetic SF CWD fractions (needed by param_derived Init!). ---
    sfp = CLM.sf_params_type()
    sfp.SF_val_CWD_frac = [0.1, 0.2, 0.3, 0.4]
    CLM.SFParams[] = sfp

    # --- Synthetic ed_params: damage-bin edges (percent crown lost per class). ---
    edp = CLM.ed_params_type()
    edp.ED_val_history_damage_bin_edges = [0.0, 20.0, 50.0, 80.0]
    CLM.EDParams[] = edp

    # --- Build the derived damage-transition matrices into the live param_derived. ---
    pd = CLM.param_derived_type()
    CLM.Init!(pd, numpft)
    CLM.ParamDerived[] = pd

    # ---------------------------------------------------------------------------
    # GetCrownReduction — fraction of crown lost = bin_edge / 100
    # ---------------------------------------------------------------------------
    @testset "GetCrownReduction" begin
        @test CLM.GetCrownReduction(1) ≈ 0.0
        @test CLM.GetCrownReduction(2) ≈ 0.20
        @test CLM.GetCrownReduction(3) ≈ 0.50
        @test CLM.GetCrownReduction(4) ≈ 0.80
        # bounds: always in [0, 1] given edges in [0, 100]
        for cd in 1:nld
            cr = CLM.GetCrownReduction(cd)
            @test 0.0 <= cr <= 1.0
        end
    end

    # ---------------------------------------------------------------------------
    # GetDamageFrac — transition fraction read from param_derived matrix
    # ---------------------------------------------------------------------------
    @testset "GetDamageFrac" begin
        # PFT 1 has damage_frac = 0 -> identity matrix.
        for i in 1:nld, j in 1:nld
            expected = (i == j) ? 1.0 : 0.0
            @test CLM.GetDamageFrac(i, j, 1) ≈ expected
        end

        # Matches the underlying param_derived matrix exactly for every entry.
        for ft in 1:numpft, i in 1:nld, j in 1:nld
            @test CLM.GetDamageFrac(i, j, ft) ==
                  CLM.param_derived().damage_transitions[i, j, ft]
        end

        # PFT 2 (damage_frac 0.2): diagonal stay = 1 - df for non-terminal rows.
        @test CLM.GetDamageFrac(1, 1, 2) ≈ 0.8
        # Each source row is a valid distribution.
        for ft in 1:numpft, i in 1:nld
            row = [CLM.GetDamageFrac(i, j, ft) for j in 1:nld]
            @test all(row .>= 0.0)
            @test sum(row) ≈ 1.0
        end
    end

    # ---------------------------------------------------------------------------
    # GetDamageMortality — logistic in crown loss; undamaged class -> 0
    # ---------------------------------------------------------------------------
    @testset "GetDamageMortality" begin
        # Undamaged class is exactly zero for every PFT.
        for ft in 1:numpft
            @test CLM.GetDamageMortality(1, ft) == 0.0
        end

        # Closed-form check for a damaged class.
        ft = 2
        cd = 3
        crown_loss = 50.0 / 100.0
        p1 = pft.damage_mort_p1[ft]
        p2 = pft.damage_mort_p2[ft]
        expected = 1.0 / (1.0 + exp(-1.0 * p2 * (crown_loss - p1)))
        @test CLM.GetDamageMortality(cd, ft) ≈ expected

        # Mortality is a probability rate in (0, 1) for all damaged classes,
        # and monotonically non-decreasing with crown damage class.
        for ft in 1:numpft
            prev = -1.0
            for cd in 2:nld
                m = CLM.GetDamageMortality(cd, ft)
                @test 0.0 < m < 1.0
                @test m >= prev
                prev = m
            end
        end
    end

    # ---------------------------------------------------------------------------
    # IsItDamageTime! — calendar-driven event flag for each event code
    # ---------------------------------------------------------------------------
    @testset "IsItDamageTime! calendar logic" begin
        is_master = CLM.ifalse   # suppress log output

        set_code(c) = (edp.damage_event_code = c; CLM.EDParams[] = edp)
        set_cal(; day=1, month=1, year=1, modelday=1, doy=1) = begin
            CLM.hlm_current_day[]   = day
            CLM.hlm_current_month[] = month
            CLM.hlm_current_year[]  = year
            CLM.hlm_model_day[]     = Float64(modelday)
            CLM.hlm_day_of_year[]   = doy
        end

        # code 1: damage off, never true.
        set_code(1)
        set_cal(; modelday=1, day=1)
        CLM.IsItDamageTime!(is_master)
        @test CLM.damage_time[] == false

        # code 2: only on first model day.
        set_code(2)
        set_cal(; modelday=1)
        CLM.IsItDamageTime!(is_master)
        @test CLM.damage_time[] == true
        set_cal(; modelday=2)
        CLM.IsItDamageTime!(is_master)
        @test CLM.damage_time[] == false

        # code 3: every day.
        set_code(3)
        set_cal(; modelday=42, day=17)
        CLM.IsItDamageTime!(is_master)
        @test CLM.damage_time[] == true

        # code 4: first day of month only.
        set_code(4)
        set_cal(; day=1)
        CLM.IsItDamageTime!(is_master)
        @test CLM.damage_time[] == true
        set_cal(; day=2)
        CLM.IsItDamageTime!(is_master)
        @test CLM.damage_time[] == false

        # negative code (> -366): specific day of year.
        set_code(-100)
        set_cal(; doy=100)
        CLM.IsItDamageTime!(is_master)
        @test CLM.damage_time[] == true
        set_cal(; doy=101)
        CLM.IsItDamageTime!(is_master)
        @test CLM.damage_time[] == false

        # YYYYMMDD specific event: 2025-06-22 -> 20250622.
        set_code(20250622)
        set_cal(; day=22, month=6, year=2025)
        CLM.IsItDamageTime!(is_master)
        @test CLM.damage_time[] == true
        set_cal(; day=23, month=6, year=2025)
        CLM.IsItDamageTime!(is_master)
        @test CLM.damage_time[] == false

        # invalid code raises (Fortran endrun).
        set_code(7)
        set_cal()
        @test_throws Exception CLM.IsItDamageTime!(is_master)

        # undamaged_class flag exported and == 1.
        @test CLM.undamaged_class == 1
    end

end
