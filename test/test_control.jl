@testset "Control Module" begin

    # ------------------------------------------------------------------
    # RUNTYP labels
    # ------------------------------------------------------------------
    @testset "RUNTYP labels" begin
        @test CLM.RUNTYP[CLM.NSRSTARTUP + 1]  == "initial"
        @test CLM.RUNTYP[CLM.NSRCONTINUE + 1] == "restart"
        @test CLM.RUNTYP[CLM.NSRBRANCH + 1]   == "branch"
    end

    # ------------------------------------------------------------------
    # control_init! — default VarCtl should pass with minimal overrides
    # ------------------------------------------------------------------
    @testset "control_init! defaults" begin
        ctl = CLM.VarCtl(nsrest=CLM.NSRSTARTUP)
        CLM.control_init!(ctl)

        # anoxia follows use_lch4 (default true → anoxia = true)
        @test ctl.anoxia == ctl.use_lch4

        # nfix_timeconst should be resolved from sentinel
        @test ctl.nfix_timeconst != -1.2345

        # use_nitrif_denitrif default is true → nfix_timeconst = 10
        @test ctl.nfix_timeconst == 10.0

        # nrevsn cleared for startup
        @test ctl.nrevsn == ""

        # FATES defaults are false, so fates_bgc and fates_sp stay false
        @test ctl.use_fates_bgc == false
        @test ctl.use_fates_sp  == false
    end

    # ------------------------------------------------------------------
    # control_init! — nfix_timeconst without nitrif_denitrif
    # ------------------------------------------------------------------
    @testset "nfix_timeconst without nitrif_denitrif" begin
        ctl = CLM.VarCtl(nsrest=CLM.NSRSTARTUP, use_nitrif_denitrif=false, use_lch4=false)
        CLM.control_init!(ctl)
        @test ctl.nfix_timeconst == 0.0
    end

    # ------------------------------------------------------------------
    # control_init! — FATES overrides
    # ------------------------------------------------------------------
    @testset "FATES overrides" begin
        ctl = CLM.VarCtl(
            nsrest=CLM.NSRSTARTUP,
            use_fates=true,
            use_fates_sp=false,
            use_lch4=false,
        )
        CLM.control_init!(ctl)

        @test ctl.spinup_matrixcn == false
        @test ctl.hist_wrt_matrixcn_diag == false
        @test ctl.use_fates_bgc == true   # !use_fates_sp
        @test ctl.anoxia == false          # use_lch4 = false
    end

    @testset "FATES SP mode" begin
        ctl = CLM.VarCtl(
            nsrest=CLM.NSRSTARTUP,
            use_fates=true,
            use_fates_sp=true,
            use_lch4=false,
        )
        CLM.control_init!(ctl)
        @test ctl.use_fates_bgc == false
    end

    # ------------------------------------------------------------------
    # control_init! — Continue run sets nrevsn
    # ------------------------------------------------------------------
    @testset "continue run sets nrevsn" begin
        ctl = CLM.VarCtl(nsrest=CLM.NSRCONTINUE)
        CLM.control_init!(ctl)
        @test ctl.nrevsn == "set by restart pointer file"
    end

    # ------------------------------------------------------------------
    # control_init! — error paths
    # ------------------------------------------------------------------
    @testset "error: invalid co2_type" begin
        ctl = CLM.VarCtl(nsrest=CLM.NSRSTARTUP, co2_type="bogus")
        @test_throws ErrorException CLM.control_init!(ctl)
    end

    @testset "error: nsrest unset" begin
        ctl = CLM.VarCtl()  # nsrest defaults to IUNDEF
        @test_throws ErrorException CLM.control_init!(ctl)
    end

    @testset "error: branch without nrevsn" begin
        ctl = CLM.VarCtl(nsrest=CLM.NSRBRANCH, nrevsn="")
        @test_throws ErrorException CLM.control_init!(ctl)
    end

    @testset "error: co2_ppmv out of range" begin
        ctl = CLM.VarCtl(nsrest=CLM.NSRSTARTUP, co2_ppmv=0.0)
        @test_throws ErrorException CLM.control_init!(ctl)

        ctl2 = CLM.VarCtl(nsrest=CLM.NSRSTARTUP, co2_ppmv=5000.0)
        @test_throws ErrorException CLM.control_init!(ctl2)
    end

    @testset "error: SNICAR dual internal mixing" begin
        ctl = CLM.VarCtl(nsrest=CLM.NSRSTARTUP, snicar_snobc_intmix=true, snicar_snodst_intmix=true)
        @test_throws ErrorException CLM.control_init!(ctl)
    end

    @testset "error: use_fates + use_cn" begin
        ctl = CLM.VarCtl(nsrest=CLM.NSRSTARTUP, use_fates=true, use_cn=true, use_lch4=false)
        @test_throws ErrorException CLM.control_init!(ctl)
    end

    @testset "error: use_fates + use_crop" begin
        ctl = CLM.VarCtl(nsrest=CLM.NSRSTARTUP, use_fates=true, use_crop=true, use_lch4=false)
        @test_throws ErrorException CLM.control_init!(ctl)
    end

    @testset "error: use_crop without create_crop_landunit" begin
        ctl = CLM.VarCtl(nsrest=CLM.NSRSTARTUP, use_crop=true, create_crop_landunit=false)
        @test_throws ErrorException CLM.control_init!(ctl)
    end

    @testset "error: single_column without lat/lon" begin
        ctl = CLM.VarCtl(nsrest=CLM.NSRSTARTUP, single_column=true)
        @test_throws ErrorException CLM.control_init!(ctl)
    end

    @testset "error: invalid spinup_state" begin
        ctl = CLM.VarCtl(nsrest=CLM.NSRSTARTUP, spinup_state=5)
        @test_throws ErrorException CLM.control_init!(ctl)
    end

    @testset "error: glc_do_dynglacier + collapse_urban" begin
        ctl = CLM.VarCtl(nsrest=CLM.NSRSTARTUP, glc_do_dynglacier=true, collapse_urban=true)
        @test_throws ErrorException CLM.control_init!(ctl)
    end

    # ------------------------------------------------------------------
    # control_print — smoke test (should not error)
    # ------------------------------------------------------------------
    @testset "control_print smoke test" begin
        ctl = CLM.VarCtl(nsrest=CLM.NSRSTARTUP)
        CLM.control_init!(ctl)
        buf = IOBuffer()
        CLM.control_print(ctl; io=buf)
        output = String(take!(buf))
        @test contains(output, "define run:")
        @test contains(output, "initial")
        @test contains(output, "co2_type")
    end

    @testset "control_print with FATES" begin
        ctl = CLM.VarCtl(nsrest=CLM.NSRSTARTUP, use_fates=true, use_lch4=false)
        CLM.control_init!(ctl)
        buf = IOBuffer()
        CLM.control_print(ctl; io=buf)
        output = String(take!(buf))
        @test contains(output, "use_fates = true")
        @test contains(output, "fates_spitfire_mode")
    end

end
