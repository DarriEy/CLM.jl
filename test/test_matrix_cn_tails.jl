using LinearAlgebra

# =============================================================================
# Matrix-CN deferred tails:
#   1. SASU spinup steady-state accelerator   (sasu_steady_state)
#   2. C13/C14 isotope matrix advance         (cn_veg_matrix_solve_iso!)
#   3. N-pool veg matrix solve                (cn_veg_matrix_solve_n!)
#   4. Matrix-state restart round-trip        (write/read_matrixcn_restart!)
#
# Each is gated behind use_matrixcn / spinup. The asserts are:
#   - SASU: the iterated pool advance X(t+1)=X(t)+A·X(t)+B converges to the
#     analytic steady state −A^{-1}B (to 1e-8), which is exactly what the SASU
#     accelerator computes in closed form.
#   - isotope: the isotope advance rides the SAME A operator as bulk C → for a
#     uniform isotope ratio r, Xiso == r·Xc (the operator is mass-proportional).
#   - N: the matrix N advance == a hand-coded sequential N pool update.
#   - restart: the accumulators/counters round-trip bit-for-bit.
# =============================================================================

@testset "Matrix-CN deferred tails" begin

    # ------------------------------------------------------------------ #
    # 1. SASU steady-state accelerator                                    #
    # ------------------------------------------------------------------ #
    @testset "SASU steady state == iterated matrix fixpoint" begin
        # A small stable linear pool system dX/dt = A·X + B.
        # A has negative diagonal (turnover) and off-diagonal transfers; B is
        # the constant input. Steady state X* satisfies A·X* + B = 0.
        A = [ -0.30   0.00   0.00;
               0.20  -0.50   0.10;
               0.05   0.30  -0.40 ]
        B = [ 0.040, 0.010, 0.025 ]

        # Analytic steady state via the SASU solver.
        Xss = CLM.sasu_steady_state(A, B)

        # Reference: closed-form −A^{-1}B (no near-zero diagonal here).
        Xref = -(A \ B)
        @test all(Xref .>= 0)
        @test isapprox(Xss, Xref; atol = 1e-12, rtol = 1e-12)

        # Iterate the forward pool advance X ← X + (A·X + B) and confirm it
        # converges to the same steady state. (Forward Euler with dt folded into
        # A/B; A is stable so the iteration contracts.)
        X = zeros(3)
        for _ in 1:200000
            X = X .+ (A * X .+ B)
        end
        @test isapprox(X, Xss; atol = 1e-8, rtol = 1e-6)

        # Residual A·X* + B ≈ 0 (definition of steady state).
        @test maximum(abs.(A * Xss .+ B)) < 1e-10
    end

    @testset "SASU near-zero-diagonal guard" begin
        # A pool with no turnover (zero diagonal) → guarded to 1e36 → capacity ~0.
        A = [ -0.30  0.00   0.00;
               0.20  0.00   0.00;    # pool 2 has zero turnover
               0.05  0.00  -0.40 ]
        B = [ 0.040, 0.010, 0.025 ]
        Xss = CLM.sasu_steady_state(A, B)
        @test all(Xss .>= 0)
        @test Xss[2] < 1e-30     # guarded pool driven to ~0
        @test isfinite(Xss[1]) && isfinite(Xss[3])
    end

    @testset "SASU clamps negative capacities" begin
        # A configuration whose −A^{-1}B has a negative component → clamped to 0.
        A = [ -0.10   0.00;
               0.50  -0.10 ]
        B = [ -0.05, 0.02 ]
        Xss = CLM.sasu_steady_state(A, B)
        @test all(Xss .>= 0)
    end

    @testset "veg SASU capacity accumulator → steady state" begin
        # Build a small veg-like matrix operator, accumulate it through the SASU
        # accumulator, and confirm the period-end capacity matches the analytic
        # steady state of the (per-X0-normalized) operator.
        nveg = 4
        bounds = 1:1
        mask = trues(1)
        s = CLM.CNVegMatrixSASU()

        # begin-of-period pools X0 and the per-step operator A (the I+AKall − I
        # part) + input B. We accumulate ONE step then take capacity.
        X0 = reshape([10.0, 5.0, 8.0, 3.0], 1, nveg)
        CLM.cn_veg_matrix_sasu_save_x0!(s, X0; mask_soilp = mask,
            bounds_patch = bounds, nveg = nveg)

        Arate = [ -0.20   0.00   0.00   0.00;
                   0.15  -0.30   0.00   0.00;
                   0.00   0.20  -0.25   0.00;
                   0.05   0.05   0.10  -0.40 ]
        # transfer flux matrix accumulated = Arate * diag(X0); we pass A_dense=Arate.
        A_dense = reshape(permutedims(Arate, (1, 2)), 1, nveg, nveg)
        # build the (1 × nveg × nveg) form
        A3 = Array{Float64}(undef, 1, nveg, nveg)
        A3[1, :, :] = Arate
        B = reshape([0.30, 0.10, 0.05, 0.02], 1, nveg)

        CLM.cn_veg_matrix_sasu_accumulate!(s, A3, B;
            mask_soilp = mask, bounds_patch = bounds, nveg = nveg)
        cap = CLM.cn_veg_matrix_sasu_capacity!(s;
            mask_soilp = mask, bounds_patch = bounds, nveg = nveg)

        # After one step the per-X0-normalized accumulated matrix == Arate, so
        # capacity == −Arate^{-1}·B.
        Xref = -(Arate \ B[1, :])
        Xref = max.(Xref, 0.0)
        @test isapprox(cap[1, :], Xref; atol = 1e-9, rtol = 1e-9)
        # accumulators were reset
        @test all(s.alloc_acc .== 0) && all(s.transfer_acc .== 0)
    end

    # ------------------------------------------------------------------ #
    # 2 + 3. veg N solve and isotope solve share _cn_veg_matrix_advance!  #
    # ------------------------------------------------------------------ #
    # Common small veg setup (no crop, no fire) reused by N and isotope tests.
    np = 2
    dt = 1800.0
    counts = CLM.veg_matrix_transfer_counts(false)
    nveg = CLM.NVEGPOOL_NATVEG
    npcropmin = 15
    ivt = [1, 2]
    bounds = 1:np
    mask = trues(np)

    # in-veg transfers (doner -> receiver), small rates (avoid the cap):
    #   1: leaf_st -> leaf_xf  2: leaf_xf -> leaf  3: froot_st -> froot_xf
    #   4: froot_xf -> froot   5: livestem -> deadstem
    doner    = [CLM.ILEAF_ST, CLM.ILEAF_XF, CLM.IFROOT_ST, CLM.IFROOT_XF, CLM.ILIVESTEM]
    receiver = [CLM.ILEAF_XF, CLM.ILEAF,    CLM.IFROOT_XF, CLM.IFROOT,    CLM.IDEADSTEM]
    nnon = length(doner)
    leaf_litter_rate  = 3.0e-7
    froot_litter_rate = 2.0e-7
    ph_rate = [ 5.0e-6  4.0e-6  3.0e-6  2.0e-6  1.0e-6;
                2.0e-6  6.0e-6  1.0e-6  5.0e-6  7.0e-7 ]

    @testset "veg N matrix advance == sequential" begin
        ns_input = [4.0e-6, 2.0e-6]    # gN/m2/s
        nalloc = zeros(np, nveg)
        for p in 1:np
            nalloc[p, CLM.ILEAF]     = 0.30
            nalloc[p, CLM.ILEAF_ST]  = 0.10
            nalloc[p, CLM.IFROOT]    = 0.25
            nalloc[p, CLM.IFROOT_ST] = 0.10
            nalloc[p, CLM.ILIVESTEM] = 0.25
        end

        # initial N pools
        function fresh_nstate()
            ns = CLM.CNVegNitrogenStateData()
            CLM.cnveg_nitrogen_state_init!(ns, np, 1, 1; use_matrixcn = false, nrepr = 1)
            ns.leafn_patch          .= [5.0, 3.0]
            ns.leafn_storage_patch  .= [2.0, 1.5]
            ns.leafn_xfer_patch     .= [1.0, 0.8]
            ns.frootn_patch         .= [4.0, 2.5]
            ns.frootn_storage_patch .= [1.8, 1.2]
            ns.frootn_xfer_patch    .= [0.9, 0.7]
            ns.livestemn_patch          .= [6.0, 4.5]
            ns.livestemn_storage_patch  .= [2.2, 1.6]
            ns.livestemn_xfer_patch     .= [1.1, 0.9]
            ns.deadstemn_patch          .= [10.0, 8.0]
            ns.deadstemn_storage_patch  .= [0.5, 0.4]
            ns.deadstemn_xfer_patch     .= [0.3, 0.25]
            ns.livecrootn_patch          .= [3.5, 2.8]
            ns.livecrootn_storage_patch  .= [1.2, 0.9]
            ns.livecrootn_xfer_patch     .= [0.6, 0.5]
            ns.deadcrootn_patch          .= [7.0, 5.5]
            ns.deadcrootn_storage_patch  .= [0.4, 0.3]
            ns.deadcrootn_xfer_patch     .= [0.2, 0.18]
            ns.retransn_patch            .= [2.0, 1.5]
            return ns
        end

        # ---- index <-> N pool accessor ----
        npget(ns, i, p) =
            i == CLM.ILEAF        ? ns.leafn_patch[p] :
            i == CLM.ILEAF_ST     ? ns.leafn_storage_patch[p] :
            i == CLM.ILEAF_XF     ? ns.leafn_xfer_patch[p] :
            i == CLM.IFROOT       ? ns.frootn_patch[p] :
            i == CLM.IFROOT_ST    ? ns.frootn_storage_patch[p] :
            i == CLM.IFROOT_XF    ? ns.frootn_xfer_patch[p] :
            i == CLM.ILIVESTEM    ? ns.livestemn_patch[p] :
            i == CLM.IDEADSTEM    ? ns.deadstemn_patch[p] :
            error("npget: $i")
        function npadd!(ns, i, p, v)
            if     i == CLM.ILEAF;        ns.leafn_patch[p] += v
            elseif i == CLM.ILEAF_ST;     ns.leafn_storage_patch[p] += v
            elseif i == CLM.ILEAF_XF;     ns.leafn_xfer_patch[p] += v
            elseif i == CLM.IFROOT;       ns.frootn_patch[p] += v
            elseif i == CLM.IFROOT_ST;    ns.frootn_storage_patch[p] += v
            elseif i == CLM.IFROOT_XF;    ns.frootn_xfer_patch[p] += v
            elseif i == CLM.ILIVESTEM;    ns.livestemn_patch[p] += v
            elseif i == CLM.IDEADSTEM;    ns.deadstemn_patch[p] += v
            else; error("npadd!: $i"); end
        end

        # ===================== (a) SEQUENTIAL reference =====================
        ns_seq = fresh_nstate()
        for p in 1:np
            t = Dict{Int,Float64}()
            for k in 1:nnon
                t[k] = ph_rate[p, k] * npget(ns_seq, doner[k], p)
            end
            litt_leaf  = leaf_litter_rate  * ns_seq.leafn_patch[p]
            litt_froot = froot_litter_rate * ns_seq.frootn_patch[p]
            ns_seq.leafn_patch[p]          += nalloc[p, CLM.ILEAF]     * ns_input[p] * dt
            ns_seq.leafn_storage_patch[p]  += nalloc[p, CLM.ILEAF_ST]  * ns_input[p] * dt
            ns_seq.frootn_patch[p]         += nalloc[p, CLM.IFROOT]    * ns_input[p] * dt
            ns_seq.frootn_storage_patch[p] += nalloc[p, CLM.IFROOT_ST] * ns_input[p] * dt
            ns_seq.livestemn_patch[p]      += nalloc[p, CLM.ILIVESTEM] * ns_input[p] * dt
            for k in 1:nnon
                npadd!(ns_seq, receiver[k], p,  t[k] * dt)
                npadd!(ns_seq, doner[k],    p, -t[k] * dt)
            end
            ns_seq.leafn_patch[p]  -= litt_leaf  * dt
            ns_seq.frootn_patch[p] -= litt_froot * dt
        end

        # ===================== (b) MATRIX advance =====================
        ns_mat = fresh_nstate()
        nf = CLM.CNVegNitrogenFluxData()
        CLM.cnveg_nitrogen_flux_init!(nf, np, 1, 1; use_matrixcn = true, nvegnpool = nveg)

        nph = counts.nnphtrans; nphout = counts.nnphouttrans
        @assert nnon <= nph - nphout
        full_doner    = fill(CLM.ILEAF, nph)
        full_receiver = fill(CLM.ILEAF, nph)
        for k in 1:nnon
            full_doner[k] = doner[k]; full_receiver[k] = receiver[k]
        end
        ioutn = nveg + 1
        # last nphout entries: leaf, froot out-transfers, pad rest with zero-rate livestem-out
        full_doner[nph-nphout+1] = CLM.ILEAF;     full_receiver[nph-nphout+1] = ioutn
        full_doner[nph-nphout+2] = CLM.IFROOT;    full_receiver[nph-nphout+2] = ioutn
        for k in (nph-nphout+3):nph
            full_doner[k] = CLM.ILIVESTEM; full_receiver[k] = ioutn
        end

        nf.matrix_nphtransfer_doner_patch    = full_doner
        nf.matrix_nphtransfer_receiver_patch = full_receiver
        nf.matrix_nphtransfer_patch = zeros(np, nph)
        nf.matrix_nphturnover_patch = zeros(np, nveg)
        # gap-mortality + fire: zero (inactive), but need valid doner/receiver of
        # full length so _build_ak_process! has a consistent in-veg block.
        ngm = counts.nngmtrans; ngmout = counts.nngmouttrans
        nf.matrix_ngmtransfer_doner_patch    = fill(CLM.ILEAF, ngm)
        nf.matrix_ngmtransfer_receiver_patch = fill(CLM.ILEAF, ngm)
        nf.matrix_ngmtransfer_patch = zeros(np, ngm)
        nf.matrix_ngmturnover_patch = zeros(np, nveg)
        nfi = counts.nnfitrans; nfiout = counts.nnfiouttrans
        nf.matrix_nfitransfer_doner_patch    = fill(CLM.ILEAF, nfi)
        nf.matrix_nfitransfer_receiver_patch = fill(CLM.ILEAF, nfi)
        nf.matrix_nfitransfer_patch = zeros(np, nfi)
        nf.matrix_nfiturnover_patch = zeros(np, nveg)
        nf.matrix_nalloc_patch  = nalloc
        nf.matrix_Ninput_patch  = ns_input

        # populate phenology transfer/turnover via the accumulator helper analog:
        # transfer[p,k] = rate; turnover[p,doner] += rate*dt.
        for p in 1:np
            for k in 1:nnon
                nf.matrix_nphtransfer_patch[p, k] = ph_rate[p, k]
                nf.matrix_nphturnover_patch[p, doner[k]] += ph_rate[p, k] * dt
            end
            # litterfall out-transfers contribute only to turnover
            nf.matrix_nphtransfer_patch[p, nph-nphout+1] = leaf_litter_rate
            nf.matrix_nphturnover_patch[p, CLM.ILEAF]  += leaf_litter_rate  * dt
            nf.matrix_nphtransfer_patch[p, nph-nphout+2] = froot_litter_rate
            nf.matrix_nphturnover_patch[p, CLM.IFROOT] += froot_litter_rate * dt
        end

        CLM.cn_veg_matrix_solve_n!(ns_mat, nf; mask_soilp = mask,
            bounds_patch = bounds, ivt = ivt, npcropmin = npcropmin,
            nvegnpool = nveg, counts = counts, dt = dt)

        # compare every loaded pool
        for p in 1:np
            @test isapprox(ns_mat.leafn_patch[p],          ns_seq.leafn_patch[p];          atol = 1e-10)
            @test isapprox(ns_mat.leafn_storage_patch[p],  ns_seq.leafn_storage_patch[p];  atol = 1e-10)
            @test isapprox(ns_mat.leafn_xfer_patch[p],     ns_seq.leafn_xfer_patch[p];     atol = 1e-10)
            @test isapprox(ns_mat.frootn_patch[p],         ns_seq.frootn_patch[p];         atol = 1e-10)
            @test isapprox(ns_mat.frootn_storage_patch[p], ns_seq.frootn_storage_patch[p]; atol = 1e-10)
            @test isapprox(ns_mat.frootn_xfer_patch[p],    ns_seq.frootn_xfer_patch[p];    atol = 1e-10)
            @test isapprox(ns_mat.livestemn_patch[p],      ns_seq.livestemn_patch[p];      atol = 1e-10)
            @test isapprox(ns_mat.deadstemn_patch[p],      ns_seq.deadstemn_patch[p];      atol = 1e-10)
            # untouched pools unchanged
            @test isapprox(ns_mat.retransn_patch[p],       ns_seq.retransn_patch[p];       atol = 1e-10)
        end
    end

    @testset "isotope matrix advance rides the bulk-C operator" begin
        # Build a bulk-C flux struct with the same phenology operator, advance the
        # bulk C AND an isotope vector seeded as r·Xc with matching r·B input. The
        # isotope advance must return exactly r·(bulk advance) since A is linear.
        cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cf, np, 1, 1; use_matrixcn = true)

        nph = counts.ncphtrans; nphout = counts.ncphouttrans
        full_doner    = fill(CLM.ILEAF, nph)
        full_receiver = fill(CLM.ILEAF, nph)
        for k in 1:nnon
            full_doner[k] = doner[k]; full_receiver[k] = receiver[k]
        end
        ioutc = nveg + 1
        full_doner[nph-nphout+1] = CLM.ILEAF;     full_receiver[nph-nphout+1] = ioutc
        full_doner[nph-nphout+2] = CLM.IFROOT;    full_receiver[nph-nphout+2] = ioutc
        for k in (nph-nphout+3):nph
            full_doner[k] = CLM.ILIVESTEM; full_receiver[k] = ioutc
        end
        cf.matrix_phtransfer_doner_patch    = full_doner
        cf.matrix_phtransfer_receiver_patch = full_receiver
        cf.matrix_phtransfer_patch = zeros(np, nph)
        cf.matrix_phturnover_patch = zeros(np, nveg)
        ngm = counts.ncgmtrans
        cf.matrix_gmtransfer_doner_patch    = fill(CLM.ILEAF, ngm)
        cf.matrix_gmtransfer_receiver_patch = fill(CLM.ILEAF, ngm)
        cf.matrix_gmtransfer_patch = zeros(np, ngm)
        cf.matrix_gmturnover_patch = zeros(np, nveg)
        nfi = counts.ncfitrans
        cf.matrix_fitransfer_doner_patch    = fill(CLM.ILEAF, nfi)
        cf.matrix_fitransfer_receiver_patch = fill(CLM.ILEAF, nfi)
        cf.matrix_fitransfer_patch = zeros(np, nfi)
        cf.matrix_fiturnover_patch = zeros(np, nveg)
        for p in 1:np
            for k in 1:nnon
                cf.matrix_phtransfer_patch[p, k] = ph_rate[p, k]
                cf.matrix_phturnover_patch[p, doner[k]] += ph_rate[p, k] * dt
            end
        end

        # Build a bulk-C pool vector Xc and input Bc, advance via the shared core.
        begp, endp = 1, np
        Xc = CLM.VectorType(); CLM.init_v!(Xc, nveg, begp, endp)
        Bc = CLM.VectorType(); CLM.init_v!(Bc, nveg, begp, endp)
        filter_soilp = collect(1:np)
        Xc0 = zeros(np, nveg)
        for p in 1:np, i in 1:nveg
            Xc0[p, i] = 10.0 + i + 0.5p
            Xc.V[p, i] = Xc0[p, i]
            Bc.V[p, i] = (i == CLM.ILEAF || i == CLM.IFROOT) ? 0.02 * p : 0.0
        end
        CLM._cn_veg_matrix_advance!(Xc, Bc, begp, endp, np, filter_soilp, nveg,
            cf.matrix_phtransfer_patch, cf.matrix_phturnover_patch,
            cf.matrix_phtransfer_doner_patch, cf.matrix_phtransfer_receiver_patch,
            counts.ncphtrans, counts.ncphouttrans,
            cf.matrix_gmtransfer_patch, cf.matrix_gmturnover_patch,
            cf.matrix_gmtransfer_doner_patch, cf.matrix_gmtransfer_receiver_patch,
            counts.ncgmtrans, counts.ncgmouttrans,
            cf.matrix_fitransfer_patch, cf.matrix_fiturnover_patch,
            cf.matrix_fitransfer_doner_patch, cf.matrix_fitransfer_receiver_patch,
            counts.ncfitrans, counts.ncfiouttrans,
            dt, 0)

        # Isotope: seed Xiso = r·Xc0, Biso = r·Bc, advance via the iso solver.
        r = 0.0112    # e.g. a 13C-like ratio
        Xiso = zeros(np, nveg)
        Biso = zeros(np, nveg)
        for p in 1:np, i in 1:nveg
            Xiso[p, i] = r * Xc0[p, i]
            Biso[p, i] = r * Bc.V[p, i]   # Bc was consumed? no — _advance! added it; Bc.V unchanged
        end
        # Bc.V was used read-only in the add; reconstruct from the same seed
        for p in 1:np, i in 1:nveg
            Biso[p, i] = r * ((i == CLM.ILEAF || i == CLM.IFROOT) ? 0.02 * p : 0.0)
        end

        CLM.cn_veg_matrix_solve_iso!(Xiso, Biso, cf; mask_soilp = mask,
            bounds_patch = bounds, nvegcpool = nveg, counts = counts, dt = dt)

        # Xiso must equal r · (advanced bulk Xc).
        for p in 1:np, i in 1:nveg
            @test isapprox(Xiso[p, i], r * Xc.V[p, i]; atol = 1e-10, rtol = 1e-9)
        end
    end

    # ------------------------------------------------------------------ #
    # 4. Matrix-state restart round-trip                                  #
    # ------------------------------------------------------------------ #
    @testset "matrix-CN restart round-trip" begin
        # Build a soil matrix state with populated AKX accumulators + counters.
        ndecomp_pools = 3; nlevdecomp = 2
        ndecomp_pools_vr = ndecomp_pools * nlevdecomp
        begc, endc = 1, 2
        nu = endc - begc + 1

        ms = CLM.CNSoilMatrixState()
        CLM.init_sm!(ms.AKXcacc, ndecomp_pools_vr, begc, endc)
        CLM.init_sm!(ms.AKXnacc, ndecomp_pools_vr, begc, endc)
        # populate a small known structure (3 entries)
        NE = 3
        ms.AKXcacc.RI[1:NE] = [1, 2, 3]
        ms.AKXcacc.CI[1:NE] = [1, 1, 2]
        ms.AKXcacc.NE = NE
        ms.AKXnacc.RI[1:NE] = [1, 2, 3]
        ms.AKXnacc.CI[1:NE] = [1, 1, 2]
        ms.AKXnacc.NE = NE
        for u in 1:nu, k in 1:NE
            ms.AKXcacc.M[u, k] = 0.1 * u + 0.01 * k
            ms.AKXnacc.M[u, k] = 0.2 * u + 0.02 * k
        end

        in_acc  = reshape(collect(1.0:(nu*ndecomp_pools_vr)), nu, ndecomp_pools_vr)
        in_nacc = in_acc .* 0.5
        decomp0_c = rand(endc, nlevdecomp, ndecomp_pools) .+ 1.0
        decomp0_n = decomp0_c .* 0.3
        iyr = 7; iloop = 2

        tmp = tempname() * ".nc"
        CLM.write_matrixcn_restart(tmp; ms = ms, in_acc = in_acc, in_nacc = in_nacc,
            decomp0_cpools_vr = decomp0_c, decomp0_npools_vr = decomp0_n,
            iyr = iyr, iloop = iloop)

        # read into fresh containers
        ms2 = CLM.CNSoilMatrixState()
        in_acc2  = zeros(nu, ndecomp_pools_vr)
        in_nacc2 = zeros(nu, ndecomp_pools_vr)
        decomp0_c2 = zeros(endc, nlevdecomp, ndecomp_pools)
        decomp0_n2 = zeros(endc, nlevdecomp, ndecomp_pools)
        (iyr2, iloop2) = CLM.read_matrixcn_restart!(tmp; ms = ms2,
            in_acc = in_acc2, in_nacc = in_nacc2,
            decomp0_cpools_vr = decomp0_c2, decomp0_npools_vr = decomp0_n2,
            begc = begc, endc = endc)

        @test iyr2 == iyr && iloop2 == iloop
        @test in_acc2  ≈ in_acc
        @test in_nacc2 ≈ in_nacc
        @test decomp0_c2 ≈ decomp0_c
        @test decomp0_n2 ≈ decomp0_n
        @test ms2.AKXcacc.NE == NE
        @test ms2.AKXcacc.RI[1:NE] == ms.AKXcacc.RI[1:NE]
        @test ms2.AKXcacc.CI[1:NE] == ms.AKXcacc.CI[1:NE]
        @test ms2.AKXcacc.M[:, 1:NE] ≈ ms.AKXcacc.M[:, 1:NE]
        @test ms2.AKXnacc.M[:, 1:NE] ≈ ms.AKXnacc.M[:, 1:NE]

        rm(tmp; force = true)
    end

end
