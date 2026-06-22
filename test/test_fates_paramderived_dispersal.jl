# test_fates_paramderived_dispersal.jl
# Tests for FATES Batch 4 (Tier F):
#   * FatesParameterDerivedMod — parameters statically derived at init from the
#     raw EDPftvarcon / EDParams / SFParams (jmax25top/tpu25top/kp25top from
#     vcmax25top, branch_frac from CWD fractions, damage-class transition
#     matrices).
#   * FatesDispersalMod — seed dispersal: kernel probability densities,
#     neighbor topology, and the calendar-driven IsItDispersalTime logic.

using Test
using CLM

@testset "FATES Batch 4: ParameterDerived + Dispersal" begin

    # ---------------------------------------------------------------------------
    # FatesParameterDerivedMod
    # ---------------------------------------------------------------------------
    @testset "FatesParameterDerivedMod" begin
        numpft = 3
        nla    = 2     # nleafage
        nld    = 4     # nlevdamage

        # Set the FATES module-global dimensions used by the allocators.
        CLM.nleafage[]   = nla
        CLM.nlevdamage[] = nld

        # --- Populate a small synthetic EDPftvarcon_inst. ---
        pft = CLM.EDPftvarcon_type()
        # vcmax25top is (numpft x nleafage)
        vcmax = [10.0  20.0;
                 30.0  40.0;
                 50.0  60.0]
        pft.vcmax25top  = vcmax
        # per-PFT annual damage fraction. Note: damage_frac == 1.0 would make the
        # TERMINAL damage class an undefined 0/0 row in the Fortran formula too
        # (diag = 1-df = 0 and no i+1 split -> NaN); we use < 1.0 so every row is
        # well-defined and faithfully comparable to hand-computed values.
        pft.damage_frac = [0.0, 0.2, 0.5]
        CLM.EDPftvarcon_inst[] = pft

        # --- Populate synthetic SF_val_CWD_frac (first 3 sum to branch_frac). ---
        sfp = CLM.sf_params_type()
        sfp.SF_val_CWD_frac = [0.1, 0.2, 0.3, 0.4]   # ncwd = 4
        CLM.SFParams[] = sfp

        # --- Populate synthetic damage-bin edges on ed_params. ---
        edp = CLM.ed_params_type()
        # nld = 4 edges; extended with 100 -> widths derived below
        edp.ED_val_history_damage_bin_edges = [0.0, 20.0, 50.0, 80.0]
        CLM.EDParams[] = edp

        # --- Run Init. ---
        pd = CLM.param_derived_type()
        CLM.Init!(pd, numpft)

        @testset "derived photosynthesis constants" begin
            @test size(pd.jmax25top) == (numpft, nla)
            @test size(pd.tpu25top)  == (numpft, nla)
            @test size(pd.kp25top)   == (numpft, nla)
            for ft in 1:numpft, iage in 1:nla
                v = vcmax[ft, iage]
                @test pd.jmax25top[ft, iage] ≈ 1.67    * v
                @test pd.tpu25top[ft, iage]  ≈ 0.167   * v
                @test pd.kp25top[ft, iage]   ≈ 20000.0 * v
            end
        end

        @testset "branch_frac = sum of first 3 CWD fractions" begin
            @test length(pd.branch_frac) == numpft
            expected = 0.1 + 0.2 + 0.3
            for ft in 1:numpft
                @test pd.branch_frac[ft] ≈ expected
            end
        end

        @testset "damage transition matrix: shape, rows sum to one, structure" begin
            @test size(pd.damage_transitions) == (nld, nld, numpft)

            # Every row sums to one (probabilities — they must go somewhere).
            for ft in 1:numpft, i in 1:nld
                @test sum(pd.damage_transitions[i, :, ft]) ≈ 1.0
            end

            # No probability mass into LESS-damaged classes (lower triangle = 0).
            for ft in 1:numpft, i in 1:nld, j in 1:(i - 1)
                @test pd.damage_transitions[i, j, ft] == 0.0
            end

            # PFT 1: damage_frac = 0 -> identity matrix (all mass on diagonal).
            for i in 1:nld
                @test pd.damage_transitions[i, i, 1] ≈ 1.0
            end

            # Hand-computed PFT 2, damage_frac = 0.2.
            # widths from edges [0,20,50,80,100] -> [20,30,30,20].
            widths = [20.0, 30.0, 30.0, 20.0]
            df = 0.2
            for i in 1:nld
                # Pre-normalization the diagonal is (1 - df); since the off-diagonal
                # mass also sums to df (when i < nld) the row already sums to one,
                # so normalization is a no-op and the diagonal stays (1 - df).
                if i < nld
                    @test pd.damage_transitions[i, i, 2] ≈ (1.0 - df)
                    denom = sum(widths[(i + 1):nld])
                    for j in (i + 1):nld
                        @test pd.damage_transitions[i, j, 2] ≈ df * widths[j] / denom
                    end
                else
                    # last class: only the diagonal, renormalized to 1.
                    @test pd.damage_transitions[i, i, 2] ≈ 1.0
                end
            end

            # PFT 3: damage_frac = 0.5 -> diagonal (stay) = 0.5 for every row that
            # has somewhere to send the damaged fraction (i < nld), and the last
            # class is renormalized to 1.
            df3 = 0.5
            for i in 1:nld
                if i < nld
                    @test pd.damage_transitions[i, i, 3] ≈ (1.0 - df3)
                    denom = sum(widths[(i + 1):nld])
                    for j in (i + 1):nld
                        @test pd.damage_transitions[i, j, 3] ≈ df3 * widths[j] / denom
                    end
                else
                    @test pd.damage_transitions[i, i, 3] ≈ 1.0
                end
            end
        end

        @testset "module singleton accessor wiring" begin
            CLM.ParamDerived[] = pd
            @test CLM.param_derived() === pd
        end
    end

    # ---------------------------------------------------------------------------
    # FatesDispersalMod
    # ---------------------------------------------------------------------------
    @testset "FatesDispersalMod" begin

        # --- Synthetic PFT dispersal parameters. ---
        numpft = 2
        pft = CLM.EDPftvarcon_type()
        # scale (param_a) and shape (param_b) per PFT
        pft.seed_dispersal_pdf_scale = [0.5, 1.0]
        pft.seed_dispersal_pdf_shape = [2.0, 1.5]
        CLM.EDPftvarcon_inst[] = pft

        @testset "exponential kernel" begin
            CLM.fates_dispersal_kernel_mode[] = CLM.fates_dispersal_kernel_exponential
            for ipft in 1:numpft
                a = pft.seed_dispersal_pdf_scale[ipft]
                for dist in (0.0, 1.0, 5.0)
                    @test CLM.ProbabilityDensity(ipft, dist) ≈ exp(-a * dist)
                end
                # at zero distance, density is exactly 1
                @test CLM.ProbabilityDensity(ipft, 0.0) == 1.0
                # monotonically decreasing with distance
                @test CLM.ProbabilityDensity(ipft, 1.0) > CLM.ProbabilityDensity(ipft, 2.0)
                # bounded in (0, 1]
                @test 0.0 < CLM.ProbabilityDensity(ipft, 3.0) <= 1.0
            end
        end

        @testset "exppower kernel" begin
            CLM.fates_dispersal_kernel_mode[] = CLM.fates_dispersal_kernel_exppower
            for ipft in 1:numpft
                a = pft.seed_dispersal_pdf_scale[ipft]
                b = pft.seed_dispersal_pdf_shape[ipft]
                dist = 2.0
                expected = (b / (2 * CLM.pi_const * CLM.fates_gamma(2 / b))) *
                           exp(-(dist^b) / (a^b))
                @test CLM.ProbabilityDensity(ipft, dist) ≈ expected
                @test CLM.ProbabilityDensity(ipft, dist) > 0.0
            end
            # local gamma matches known values
            @test CLM.fates_gamma(1.0) ≈ 1.0
            @test CLM.fates_gamma(2.0) ≈ 1.0
            @test CLM.fates_gamma(5.0) ≈ 24.0      # 4!
            @test CLM.fates_gamma(0.5) ≈ sqrt(CLM.pi_const)
        end

        @testset "logsech kernel" begin
            CLM.fates_dispersal_kernel_mode[] = CLM.fates_dispersal_kernel_logsech
            for ipft in 1:numpft
                a = pft.seed_dispersal_pdf_scale[ipft]
                b = pft.seed_dispersal_pdf_shape[ipft]
                dist = 2.0
                expected = (1 / (CLM.pi_const^2 * b * dist^2)) /
                           ((dist / a)^(1 / b) + (dist / a)^(-1 / b))
                @test CLM.ProbabilityDensity(ipft, dist) ≈ expected
                @test CLM.ProbabilityDensity(ipft, dist) > 0.0
            end
        end

        @testset "undefined kernel errors" begin
            CLM.fates_dispersal_kernel_mode[] = 999
            @test_throws Exception CLM.ProbabilityDensity(1, 1.0)
            # restore a valid kernel
            CLM.fates_dispersal_kernel_mode[] = CLM.fates_dispersal_kernel_exponential
        end

        @testset "dispersal_type init buffers" begin
            numprocs     = 2
            numgc_global = 5
            numgc_local  = 3

            # Set up the cadence + current date so init's GetCadenceDate works.
            CLM.hlm_seeddisp_cadence[] = CLM.fates_dispersal_cadence_yearly
            CLM.hlm_current_year[]     = 2026

            d = CLM.dispersal_type()
            CLM.init!(d, numprocs, numgc_global, numgc_local, numpft)

            @test size(d.outgoing_local)  == (numpft, numgc_local)
            @test size(d.outgoing_global) == (numpft, numgc_global)
            @test size(d.incoming_global) == (numpft, numgc_global)
            @test length(d.begg_array)    == numprocs
            @test length(d.ncells_array)  == numprocs

            @test all(d.outgoing_local  .== 0.0)
            @test all(d.outgoing_global .== 0.0)
            @test all(d.incoming_global .== 0.0)
            @test all(d.begg_array   .== CLM.fates_unset_int)
            @test all(d.ncells_array .== CLM.fates_unset_int)

            # init set the dispersal date to the current cadence date.
            @test CLM.dispersal_date[] == 2026
        end

        @testset "GetCadenceDate honors the configured cadence" begin
            CLM.hlm_current_day[]   = 12
            CLM.hlm_current_month[] = 6
            CLM.hlm_current_year[]  = 2026

            CLM.hlm_seeddisp_cadence[] = CLM.fates_dispersal_cadence_daily
            @test CLM.GetCadenceDate() == 12
            CLM.hlm_seeddisp_cadence[] = CLM.fates_dispersal_cadence_monthly
            @test CLM.GetCadenceDate() == 6
            CLM.hlm_seeddisp_cadence[] = CLM.fates_dispersal_cadence_yearly
            @test CLM.GetCadenceDate() == 2026
        end

        @testset "IsItDispersalTime calendar logic" begin
            CLM.hlm_seeddisp_cadence[] = CLM.fates_dispersal_cadence_yearly
            CLM.hlm_current_year[]     = 2026

            # Same date as last dispersal, flag clear -> not time.
            CLM.dispersal_date[] = 2026
            CLM.dispersal_flag[] = false
            @test CLM.IsItDispersalTime() == false

            # New date, flag clear, setdispersedflag not set -> time, no flag set.
            CLM.dispersal_date[] = 2025
            CLM.dispersal_flag[] = false
            @test CLM.IsItDispersalTime() == true
            @test CLM.dispersal_flag[] == false
            # date NOT advanced because setflag defaulted false
            @test CLM.dispersal_date[] == 2025

            # New date, setdispersedflag = true -> time, flag set, date advanced.
            CLM.dispersal_date[] = 2025
            CLM.dispersal_flag[] = false
            @test CLM.IsItDispersalTime(setdispersedflag = true) == true
            @test CLM.dispersal_flag[] == true
            @test CLM.dispersal_date[] == 2026

            # Now flag is set -> next call returns true and clears the flag,
            # regardless of date.
            CLM.dispersal_date[] = 2026   # same date
            @test CLM.IsItDispersalTime() == true
            @test CLM.dispersal_flag[] == false
        end

        @testset "neighbor / neighborhood topology + conservation" begin
            # Build a small synthetic neighbor topology: source cell with 3
            # neighbors, exponential kernel.
            CLM.fates_dispersal_kernel_mode[] = CLM.fates_dispersal_kernel_exponential

            nbrs = CLM.neighborhood_type()
            dists = [0.0, 1.0, 4.0]
            gidx  = [10, 11, 12]

            prev = nothing
            for (k, (gi, dd)) in enumerate(zip(gidx, dists))
                dp = [CLM.ProbabilityDensity(ip, dd) for ip in 1:numpft]
                n = CLM.neighbor_type(gindex = gi, gc_dist = dd, density_prob = dp)
                if prev === nothing
                    nbrs.first_neighbor = n
                else
                    prev.next_neighbor = n
                end
                nbrs.last_neighbor = n
                push!(nbrs.neighbor_indices, gi)
                prev = n
            end
            nbrs.neighbor_count = length(gidx)

            @test nbrs.neighbor_count == 3
            @test nbrs.neighbor_indices == gidx
            @test nbrs.first_neighbor.gindex == 10
            @test nbrs.last_neighbor.gindex  == 12

            # Walk the linked list, accumulate per-pft density.
            walked = 0
            totals = zeros(numpft)
            node = nbrs.first_neighbor
            while node !== nothing
                walked += 1
                totals .+= node.density_prob
                # density bounded in (0, 1] for exponential kernel
                @test all(0.0 .< node.density_prob .<= 1.0)
                node = node.next_neighbor
            end
            @test walked == nbrs.neighbor_count
            # nearest neighbor (dist 0) contributes density 1 per pft, so totals > 1
            @test all(totals .> 1.0)
        end
    end
end
