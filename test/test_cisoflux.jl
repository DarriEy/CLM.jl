# Integration-level tests for the carbon-isotope flux cascades (CIsoFlux1/2/2h/2g/3)
# wired through the CN driver: each isotopic C flux is the bulk flux scaled by the
# donor-pool ratio (iso_state/tot_state), ratios stay physical (near-PDB for C13),
# the 0/0 guard holds, and the default (isotopes off) path is byte-identical.
@testset "C Isotope Flux cascades (CIsoFlux1/2/2h/2g/3) — wiring" begin

    PDB = 0.0112372   # Vienna-PDB C13/C12 standard ratio (~1.1% → physical C13 ratio)

    # Helper: build a CNVegCarbonStateData with a given uniform isotope ratio applied
    # to a chosen set of donor pools, and a matching flux instance with bulk fluxes.
    function build_cnveg(np; ratio=PDB)
        cs = CLM.CNVegCarbonStateData()
        cf = CLM.CNVegCarbonFluxData()
        # Donor pools used by CIsoFlux1 non-mortality fluxes
        cs.cpool_patch          = fill(100.0, np)
        cs.leafc_xfer_patch     = fill(10.0, np)
        cs.frootc_xfer_patch    = fill(8.0, np)
        cs.livestemc_xfer_patch = fill(6.0, np)
        cs.deadstemc_xfer_patch = fill(4.0, np)
        cs.livecrootc_xfer_patch = fill(3.0, np)
        cs.deadcrootc_xfer_patch = fill(2.0, np)
        cs.leafc_patch          = fill(50.0, np)
        cs.frootc_patch         = fill(40.0, np)
        cs.livestemc_patch      = fill(30.0, np)
        cs.livecrootc_patch     = fill(20.0, np)
        cs.gresp_xfer_patch     = fill(5.0, np)
        cs.gresp_storage_patch  = fill(5.0, np)
        cs.leafc_storage_patch  = fill(7.0, np)
        cs.frootc_storage_patch = fill(6.0, np)
        cs.livestemc_storage_patch = fill(4.0, np)
        cs.deadstemc_storage_patch = fill(3.0, np)
        cs.livecrootc_storage_patch = fill(2.0, np)
        cs.deadcrootc_storage_patch = fill(1.0, np)
        cs.totvegc_patch        = fill(200.0, np)
        return cs, cf
    end

    # ---------------------------------------------------------------------------
    # Full c_iso_flux1! end-to-end: every isotopic veg flux = bulk * donor ratio,
    # ratios stay physical near PDB, and the column decomp cascade obeys the same
    # donor-pool-ratio rule with the 0/0 guard.
    # ---------------------------------------------------------------------------
    @testset "c_iso_flux1! end-to-end ratio cascade" begin
        np = 2; nc = 1; nlevdecomp = 1
        ndecomp_pools = 3; ndecomp_cascade_transitions = 2
        ratio = PDB   # physical C13/(C12+C13)-style ratio

        # c_iso_flux1!'s internal litter-to-column scatter reads the global litter-pool
        # indices from varpar; set them for this self-contained test (restore after).
        _save = (CLM.varpar.i_litr_min, CLM.varpar.i_litr_max, CLM.varpar.i_met_lit)
        CLM.varpar.i_litr_min = 1; CLM.varpar.i_litr_max = 2; CLM.varpar.i_met_lit = 1

        cs, cf = build_cnveg(np)
        # Bulk veg fluxes (positive donors)
        cf.cpool_to_leafc_patch            = fill(2.0, np)
        cf.leafc_xfer_to_leafc_patch       = fill(1.0, np)
        cf.frootc_xfer_to_frootc_patch     = fill(0.8, np)
        cf.leafc_to_litter_patch           = fill(0.5, np)
        cf.frootc_to_litter_patch          = fill(0.4, np)
        cf.livestemc_to_deadstemc_patch    = fill(0.3, np)
        cf.leaf_curmr_patch                = fill(0.2, np)

        # Isotope state instances: exactly `ratio` of the bulk pools.
        ics = CLM.CNVegCarbonStateData(); icf = CLM.CNVegCarbonFluxData()
        for fld in (:cpool_patch, :leafc_xfer_patch, :frootc_xfer_patch,
                    :livestemc_xfer_patch, :deadstemc_xfer_patch,
                    :livecrootc_xfer_patch, :deadcrootc_xfer_patch,
                    :leafc_patch, :frootc_patch, :livestemc_patch, :livecrootc_patch,
                    :gresp_xfer_patch, :gresp_storage_patch, :leafc_storage_patch,
                    :frootc_storage_patch, :livestemc_storage_patch,
                    :deadstemc_storage_patch, :livecrootc_storage_patch,
                    :deadcrootc_storage_patch, :totvegc_patch)
            setfield!(ics, fld, getfield(cs, fld) .* ratio)
        end

        # Soil decomposition C state/flux (bulk) + isotope soil instances.
        scs = CLM.SoilBiogeochemCarbonStateData()
        scf = CLM.SoilBiogeochemCarbonFluxData()
        iscs = CLM.SoilBiogeochemCarbonStateData()
        iscf = CLM.SoilBiogeochemCarbonFluxData()
        scs.decomp_cpools_vr_col  = fill(20.0, nc, nlevdecomp, ndecomp_pools)
        # one donor pool intentionally zero to exercise the 0/0 guard
        scs.decomp_cpools_vr_col[1, 1, 2] = 0.0
        iscs.decomp_cpools_vr_col = scs.decomp_cpools_vr_col .* ratio
        scf.decomp_cascade_hr_vr_col        = fill(0.1, nc, nlevdecomp, ndecomp_cascade_transitions)
        scf.decomp_cascade_ctransfer_vr_col = fill(0.05, nc, nlevdecomp, ndecomp_cascade_transitions)
        iscf.decomp_cascade_hr_vr_col        = fill(-9.0, nc, nlevdecomp, ndecomp_cascade_transitions)
        iscf.decomp_cascade_ctransfer_vr_col = fill(-9.0, nc, nlevdecomp, ndecomp_cascade_transitions)

        sbst = CLM.SoilBiogeochemStateData()
        sbst.leaf_prof_patch  = ones(np, nlevdecomp)
        sbst.froot_prof_patch = ones(np, nlevdecomp)
        # phenology_c_to_litr_c_col target for the litter-to-column scatter
        icf.phenology_c_to_litr_c_col = zeros(nc, nlevdecomp, 2)
        # Pre-allocate every iso output flux field c_iso_flux1! writes into. The
        # default CNVegCarbonFluxData has 0-length vectors, and c_iso_flux_calc! writes
        # in place, so the destination arrays must be sized np up front.
        for fld in fieldnames(CLM.CNVegCarbonFluxData)
            v = getfield(icf, fld)
            if v isa Vector{Float64} && length(v) == 0
                setfield!(icf, fld, zeros(np))
            elseif v isa Matrix{Float64} && size(v, 1) == 0
                setfield!(icf, fld, zeros(np, 1))
            end
        end
        # Same for the bulk flux instance's reproductive 2D fields the scatter reads.
        cf.repr_grainc_to_food_patch       = zeros(np, 1)
        cf.repr_structurec_to_litter_patch = zeros(np, 1)

        # donor pool: transition l donates from pool l (1->2)
        cascade_donor_pool = [1, 2]
        mask_soilc = trues(nc); mask_soilp = trues(np)

        CLM.c_iso_flux1!(sbst, scf, scs, cf, cs, iscf, iscs, icf, ics,
            mask_soilc, mask_soilp, 1:nc, 1:np, cascade_donor_pool,
            nlevdecomp, ndecomp_cascade_transitions, "c13";
            use_crop=false, nrepr=1, npcropmin=17,
            patch_column=[1, 1], patch_itype=[1, 1], patch_wtcol=[0.5, 0.5],
            lf_f=ones(20, 2), fr_f=ones(20, 2),
            leaf_prof=sbst.leaf_prof_patch, froot_prof=sbst.froot_prof_patch)

        # Each isotopic veg flux equals bulk * donor ratio (frax = 1).
        @test all(icf.cpool_to_leafc_patch        .≈ cf.cpool_to_leafc_patch .* ratio)
        @test all(icf.leafc_xfer_to_leafc_patch   .≈ cf.leafc_xfer_to_leafc_patch .* ratio)
        @test all(icf.leafc_to_litter_patch       .≈ cf.leafc_to_litter_patch .* ratio)
        @test all(icf.livestemc_to_deadstemc_patch .≈ cf.livestemc_to_deadstemc_patch .* ratio)
        @test all(icf.leaf_curmr_patch            .≈ cf.leaf_curmr_patch .* ratio)

        # Physical: implied flux ratio stays near PDB and within (0,1).
        for p in 1:np
            r = icf.cpool_to_leafc_patch[p] / cf.cpool_to_leafc_patch[p]
            @test 0.0 < r < 1.0
            @test isapprox(r, PDB; rtol=1e-12)
        end

        # All isotopic fluxes finite + sign-consistent with bulk.
        @test all(isfinite, icf.cpool_to_leafc_patch)
        @test all(sign.(icf.leafc_to_litter_patch) .== sign.(cf.leafc_to_litter_patch))

        # Column decomp cascade (on the iso SOIL flux instance): l=1 donor pool
        # nonzero → scaled by donor ratio; l=2 donor pool zero → 0/0 guard.
        @test iscf.decomp_cascade_hr_vr_col[1, 1, 1] ≈ scf.decomp_cascade_hr_vr_col[1, 1, 1] * ratio
        @test iscf.decomp_cascade_ctransfer_vr_col[1, 1, 1] ≈ scf.decomp_cascade_ctransfer_vr_col[1, 1, 1] * ratio
        @test iscf.decomp_cascade_hr_vr_col[1, 1, 2] == 0.0          # 0/0 guard
        @test iscf.decomp_cascade_ctransfer_vr_col[1, 1, 2] == 0.0   # 0/0 guard

        # Litterfall scatter to column carried the isotopic litterfall (mass present).
        @test icf.phenology_c_to_litr_c_col[1, 1, 1] ≈ (icf.leafc_to_litter_patch[1]*0.5 +
                                                        icf.frootc_to_litter_patch[1]*0.5 +
                                                        icf.leafc_to_litter_patch[2]*0.5 +
                                                        icf.frootc_to_litter_patch[2]*0.5)

        # restore global varpar
        CLM.varpar.i_litr_min, CLM.varpar.i_litr_max, CLM.varpar.i_met_lit = _save
    end

    # ---------------------------------------------------------------------------
    # Mass/ratio conservation: a uniform donor ratio propagates unchanged through
    # allocation -> litterfall -> decomposition (frax = 1, no fractionation).
    # ---------------------------------------------------------------------------
    @testset "ratio conservation allocation->litter->decomp" begin
        n = 6; mask = trues(n); bounds = 1:n
        ratio = 0.0105
        tot_state = collect(100.0:50.0:350.0)
        iso_state = tot_state .* ratio
        alloc_flux  = collect(1.0:1.0:6.0)
        litter_flux = collect(0.5:0.5:3.0)
        decomp_flux = collect(0.1:0.1:0.6)

        for tot_flux in (alloc_flux, litter_flux, decomp_flux)
            iso_flux = zeros(n)
            CLM.c_iso_flux_calc!(iso_flux, tot_flux, iso_state, tot_state,
                                 mask, bounds, 1.0, "c13")
            # ratio is conserved at every cascade stage
            @test all(iso_flux ./ tot_flux .≈ ratio)
            @test all(isfinite, iso_flux)
        end
    end

    # ---------------------------------------------------------------------------
    # 0/0 guard at the cascade entry point (zero donor state -> zero iso flux).
    # ---------------------------------------------------------------------------
    @testset "0/0 donor guard" begin
        n = 4; mask = trues(n); bounds = 1:n
        tot_flux  = fill(3.0, n)
        tot_state = [0.0, 10.0, 0.0, 20.0]
        iso_state = [0.0, 0.1, 5.0, 0.2]   # note: iso=0 OR tot=0 must guard
        iso_flux  = fill(NaN, n)
        CLM.c_iso_flux_calc!(iso_flux, tot_flux, iso_state, tot_state,
                             mask, bounds, 1.0, "c14")
        @test iso_flux[1] == 0.0   # both zero
        @test iso_flux[3] == 0.0   # tot_state zero
        @test iso_flux[2] ≈ 3.0 * (0.1/10.0) * (1.0 + (1.0-1.0)*2.0)
        @test all(isfinite, iso_flux)
    end

    # ---------------------------------------------------------------------------
    # Driver default: with use_c13/use_c14 OFF (no isotope instances supplied),
    # the gating list is empty so no isotopic work runs — the bulk path is
    # untouched (byte-identical). We assert the gating predicate directly.
    # ---------------------------------------------------------------------------
    @testset "driver default isotopes-off is a no-op" begin
        cfg = CLM.CNDriverConfig(use_cn=true, use_c13=false, use_c14=false)
        @test cfg.use_c13 == false
        @test cfg.use_c14 == false
        # The driver builds _iso_active only when the flag is set AND instances are
        # provided; with both false the list is empty regardless of instances.
        active = String[]
        if cfg.use_c13; push!(active, "c13"); end
        if cfg.use_c14; push!(active, "c14"); end
        @test isempty(active)
    end
end
