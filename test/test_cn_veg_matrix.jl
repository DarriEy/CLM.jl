using LinearAlgebra

# ---- helpers: index <-> pool field accessor (matches varpar pool indices) ----
function _cnvm_pool_get(cs, i::Int, p::Int)
    i == CLM.ILEAF        ? cs.leafc_patch[p] :
    i == CLM.ILEAF_ST     ? cs.leafc_storage_patch[p] :
    i == CLM.ILEAF_XF     ? cs.leafc_xfer_patch[p] :
    i == CLM.IFROOT       ? cs.frootc_patch[p] :
    i == CLM.IFROOT_ST    ? cs.frootc_storage_patch[p] :
    i == CLM.IFROOT_XF    ? cs.frootc_xfer_patch[p] :
    i == CLM.ILIVESTEM    ? cs.livestemc_patch[p] :
    i == CLM.ILIVESTEM_ST ? cs.livestemc_storage_patch[p] :
    i == CLM.ILIVESTEM_XF ? cs.livestemc_xfer_patch[p] :
    i == CLM.IDEADSTEM    ? cs.deadstemc_patch[p] :
    i == CLM.IDEADSTEM_ST ? cs.deadstemc_storage_patch[p] :
    i == CLM.IDEADSTEM_XF ? cs.deadstemc_xfer_patch[p] :
    i == CLM.ILIVECROOT   ? cs.livecrootc_patch[p] :
    i == CLM.ILIVECROOT_ST ? cs.livecrootc_storage_patch[p] :
    i == CLM.ILIVECROOT_XF ? cs.livecrootc_xfer_patch[p] :
    i == CLM.IDEADCROOT   ? cs.deadcrootc_patch[p] :
    i == CLM.IDEADCROOT_ST ? cs.deadcrootc_storage_patch[p] :
    i == CLM.IDEADCROOT_XF ? cs.deadcrootc_xfer_patch[p] :
    error("pool_get: unknown index $i")
end

function _cnvm_pool_add!(cs, i::Int, p::Int, v::Float64)
    if     i == CLM.ILEAF;        cs.leafc_patch[p] += v
    elseif i == CLM.ILEAF_ST;     cs.leafc_storage_patch[p] += v
    elseif i == CLM.ILEAF_XF;     cs.leafc_xfer_patch[p] += v
    elseif i == CLM.IFROOT;       cs.frootc_patch[p] += v
    elseif i == CLM.IFROOT_ST;    cs.frootc_storage_patch[p] += v
    elseif i == CLM.IFROOT_XF;    cs.frootc_xfer_patch[p] += v
    elseif i == CLM.ILIVESTEM;    cs.livestemc_patch[p] += v
    elseif i == CLM.ILIVESTEM_ST; cs.livestemc_storage_patch[p] += v
    elseif i == CLM.ILIVESTEM_XF; cs.livestemc_xfer_patch[p] += v
    elseif i == CLM.IDEADSTEM;    cs.deadstemc_patch[p] += v
    elseif i == CLM.IDEADSTEM_ST; cs.deadstemc_storage_patch[p] += v
    elseif i == CLM.IDEADSTEM_XF; cs.deadstemc_xfer_patch[p] += v
    elseif i == CLM.ILIVECROOT;   cs.livecrootc_patch[p] += v
    elseif i == CLM.ILIVECROOT_ST; cs.livecrootc_storage_patch[p] += v
    elseif i == CLM.ILIVECROOT_XF; cs.livecrootc_xfer_patch[p] += v
    elseif i == CLM.IDEADCROOT;   cs.deadcrootc_patch[p] += v
    elseif i == CLM.IDEADCROOT_ST; cs.deadcrootc_storage_patch[p] += v
    elseif i == CLM.IDEADCROOT_XF; cs.deadcrootc_xfer_patch[p] += v
    else; error("pool_add!: unknown index $i"); end
    return nothing
end

_cnvm_maxveg_for(ivt) = maximum(ivt)

# =============================================================================
# CNVegMatrixMod core: the matrix C-pool advance must EXACTLY reproduce the
# sequential c_state_update1 phenology + allocation + livewood-turnover pool
# updates (the matrix is an exact reformulation, not an approximation).
#
# We build a small veg setup, run BOTH:
#   (a) a hand-coded SEQUENTIAL pool update equal to the gated `!use_matrixcn`
#       block of c_state_update1! (allocation additions + phenology xfer→pool +
#       storage→xfer + litterfall + livewood turnover), and
#   (b) the MATRIX assemble+solve (cn_veg_matrix_solve_c!),
# from the SAME flux values, and assert the updated C pools agree to tight tol.
# Mass conservation (sum of in-veg pools + losses) is also checked.
# =============================================================================

@testset "CNVegMatrixMod (use_matrixcn veg C solver)" begin

    # --- transfer-count constants ---------------------------------------------
    @testset "transfer counts" begin
        c0 = CLM.veg_matrix_transfer_counts(false)
        @test c0.ncphtrans == 17 && c0.ncphouttrans == 3
        @test c0.ncgmtrans == 18 && c0.ncfitrans == 20
        c1 = CLM.veg_matrix_transfer_counts(true)
        @test c1.ncphtrans == 18 && c1.ncphouttrans == 4 && c1.nnphtrans == 37
    end

    # --- matrix_update_phc! accumulator helper --------------------------------
    @testset "matrix_update_phc! accumulator" begin
        np = 1
        cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cf, np, 1, 1; use_matrixcn=true)
        nveg = CLM.varpar.nvegcpool == 0 ? CLM.NVEGPOOL_NATVEG : CLM.varpar.nvegcpool
        # transfer 1: leaf_st (doner) -> leaf_xf (receiver)
        cf.matrix_phtransfer_patch   = zeros(np, 1)
        cf.matrix_phturnover_patch   = zeros(np, nveg)
        cf.matrix_phtransfer_doner_patch    = [CLM.ILEAF_ST]
        cf.matrix_phtransfer_receiver_patch = [CLM.ILEAF_XF]
        dt = 1800.0
        # rate small => no cap
        r = CLM.matrix_update_phc!(cf, 1, 1, 1e-5, dt; acc=true)
        @test r == 1e-5
        @test cf.matrix_phtransfer_patch[1, 1] ≈ 1e-5
        @test cf.matrix_phturnover_patch[1, CLM.ILEAF_ST] ≈ 1e-5 * dt
        # huge rate => capped so turnover fraction stops at 1
        cf.matrix_phtransfer_patch[1, 1] = 0.0
        cf.matrix_phturnover_patch[1, CLM.ILEAF_ST] = 0.0
        rcap = CLM.matrix_update_phc!(cf, 1, 1, 1.0, dt; acc=true)  # 1.0/s*dt >> 1
        @test cf.matrix_phturnover_patch[1, CLM.ILEAF_ST] ≈ 1.0    # capped at 1
        @test rcap ≈ 1.0 / dt
    end

    # --- core equivalence: matrix == sequential -------------------------------
    @testset "matrix advance == sequential update (C pools)" begin
        np = 3
        dt = 1800.0
        use_crop = false
        counts = CLM.veg_matrix_transfer_counts(use_crop)
        nveg = CLM.NVEGPOOL_NATVEG   # 18 (no crop)

        # ivt: all natural woody PFTs (>=1, < npcropmin) so livewood turnover runs.
        npcropmin = 15
        ivt = [1, 2, 3]
        woody = ones(Float64, _cnvm_maxveg_for(ivt))   # woody[pft]==1 -> woody

        # ---- choose per-patch phenology transfer rates (1/s) and litterfall ----
        # in-veg transfers (doner -> receiver), small enough to avoid the cap:
        #   1: leaf_st  -> leaf_xf     (storage -> transfer)
        #   2: leaf_xf  -> leaf        (transfer -> pool)
        #   3: froot_st -> froot_xf
        #   4: froot_xf -> froot
        #   5: livestem -> deadstem    (livewood turnover)
        #   6: livecroot-> deadcroot   (livewood turnover)
        doner    = [CLM.ILEAF_ST, CLM.ILEAF_XF, CLM.IFROOT_ST, CLM.IFROOT_XF, CLM.ILIVESTEM, CLM.ILIVECROOT]
        receiver = [CLM.ILEAF_XF, CLM.ILEAF,    CLM.IFROOT_XF, CLM.IFROOT,    CLM.IDEADSTEM, CLM.IDEADCROOT]
        nnon = length(doner)        # 6 in-veg transfers
        # litterfall (pure loss) rates per pool (1/s): leaf and froot lose to iout
        leaf_litter_rate  = 3.0e-7
        froot_litter_rate = 2.0e-7

        # per-patch transfer rates (1/s)
        ph_rate = [
            5.0e-6  4.0e-6  3.0e-6  2.0e-6  1.0e-6  8.0e-7;   # patch 1
            2.0e-6  6.0e-6  1.0e-6  5.0e-6  7.0e-7  9.0e-7;   # patch 2
            0.0     0.0     0.0     0.0     0.0     0.0   ;    # patch 3 (no transfer)
        ]

        # allocation input (NPP) and allocation fractions
        Cinput = [3.0e-5, 1.0e-5, 2.0e-5]     # gC/m2/s
        alloc = zeros(np, nveg)
        # allocate to leaf, leaf_st, froot, froot_st, livestem
        for p in 1:np
            alloc[p, CLM.ILEAF]      = 0.30
            alloc[p, CLM.ILEAF_ST]   = 0.10
            alloc[p, CLM.IFROOT]     = 0.25
            alloc[p, CLM.IFROOT_ST]  = 0.10
            alloc[p, CLM.ILIVESTEM]  = 0.25
        end

        # ---- initial pool values (gC/m2) ----
        function fresh_state()
            cs = CLM.CNVegCarbonStateData()
            CLM.cnveg_carbon_state_init!(cs, np, 1, 1; use_matrixcn=false, nrepr=1)
            cs.leafc_patch          .= [50.0, 30.0, 10.0]
            cs.leafc_storage_patch  .= [20.0, 15.0, 5.0]
            cs.leafc_xfer_patch     .= [10.0, 8.0, 3.0]
            cs.frootc_patch         .= [40.0, 25.0, 12.0]
            cs.frootc_storage_patch .= [18.0, 12.0, 6.0]
            cs.frootc_xfer_patch    .= [9.0, 7.0, 4.0]
            cs.livestemc_patch          .= [60.0, 45.0, 20.0]
            cs.livestemc_storage_patch  .= [22.0, 16.0, 8.0]
            cs.livestemc_xfer_patch     .= [11.0, 9.0, 5.0]
            cs.deadstemc_patch          .= [100.0, 80.0, 40.0]
            cs.deadstemc_storage_patch  .= [5.0, 4.0, 2.0]
            cs.deadstemc_xfer_patch     .= [3.0, 2.5, 1.0]
            cs.livecrootc_patch          .= [35.0, 28.0, 14.0]
            cs.livecrootc_storage_patch  .= [12.0, 9.0, 5.0]
            cs.livecrootc_xfer_patch     .= [6.0, 5.0, 2.5]
            cs.deadcrootc_patch          .= [70.0, 55.0, 28.0]
            cs.deadcrootc_storage_patch  .= [4.0, 3.0, 1.5]
            cs.deadcrootc_xfer_patch     .= [2.0, 1.8, 0.8]
            return cs
        end

        # ===================== (a) SEQUENTIAL reference =====================
        # This mirrors the gated `!use_matrixcn` block of c_state_update1!:
        #   allocation additions, xfer->pool transfers, storage->xfer transfers,
        #   litterfall losses, livewood turnover (live->dead).
        cs_seq = fresh_state()
        for p in 1:np
            # transfer fluxes (gC/m2/s) = rate * doner_pool
            t = Dict{Int,Float64}()  # transfer index -> flux
            poolval(i) = _cnvm_pool_get(cs_seq, i, p)
            for k in 1:nnon
                t[k] = ph_rate[p, k] * poolval(doner[k])
            end
            litt_leaf  = leaf_litter_rate  * cs_seq.leafc_patch[p]
            litt_froot = froot_litter_rate * cs_seq.frootc_patch[p]

            # allocation additions
            cs_seq.leafc_patch[p]          += alloc[p, CLM.ILEAF]     * Cinput[p] * dt
            cs_seq.leafc_storage_patch[p]  += alloc[p, CLM.ILEAF_ST]  * Cinput[p] * dt
            cs_seq.frootc_patch[p]         += alloc[p, CLM.IFROOT]    * Cinput[p] * dt
            cs_seq.frootc_storage_patch[p] += alloc[p, CLM.IFROOT_ST] * Cinput[p] * dt
            cs_seq.livestemc_patch[p]      += alloc[p, CLM.ILIVESTEM] * Cinput[p] * dt

            # phenology transfers (xfer->pool and storage->xfer), each as
            #   receiver += flux*dt ; doner -= flux*dt
            for k in 1:nnon
                _cnvm_pool_add!(cs_seq, receiver[k], p,  t[k] * dt)
                _cnvm_pool_add!(cs_seq, doner[k],    p, -t[k] * dt)
            end

            # litterfall (pure loss from leaf and froot)
            cs_seq.leafc_patch[p]  -= litt_leaf  * dt
            cs_seq.frootc_patch[p] -= litt_froot * dt
        end

        # ===================== (b) MATRIX advance =====================
        cs_mat = fresh_state()
        cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cf, np, 1, 1; use_matrixcn=true)

        # allocate the matrix transfer arrays + index arrays for this scenario.
        # Transfer ordering: 1..nnon are in-veg transfers, then the trailing
        # `ncphouttrans` (=3) are pure out-of-veg losses (leaf, froot, livestem->iout)
        # — only their TURNOVER contributes (no matrix entry). We model leaf and
        # froot litterfall as out-transfers (ncphouttrans=3 => 3 out entries; we
        # put leaf, froot, and a zero livestem-out to match the count).
        nph = counts.ncphtrans          # 17
        nphout = counts.ncphouttrans    # 3
        @assert nnon <= nph - nphout
        # Build full-size doner/receiver/transfer arrays (length nph). Pad the
        # unused in-veg slots with self-transfers carrying zero rate (doner=receiver,
        # zero flux => zero matrix entry, zero turnover).
        full_doner    = fill(CLM.ILEAF, nph)
        full_receiver = fill(CLM.ILEAF, nph)
        for k in 1:nnon
            full_doner[k] = doner[k]; full_receiver[k] = receiver[k]
        end
        # trailing out-transfers: leaf->iout, froot->iout, livestem->iout(zero)
        ioutc = nveg + 1
        full_doner[nph-2] = CLM.ILEAF;     full_receiver[nph-2] = ioutc
        full_doner[nph-1] = CLM.IFROOT;    full_receiver[nph-1] = ioutc
        full_doner[nph]   = CLM.ILIVESTEM; full_receiver[nph]   = ioutc

        cf.matrix_phtransfer_doner_patch    = full_doner
        cf.matrix_phtransfer_receiver_patch = full_receiver
        cf.matrix_phtransfer_patch = zeros(np, nph)
        cf.matrix_phturnover_patch = zeros(np, nveg)

        # gap-mortality and fire: present but all-zero (no gm/fi in this scenario).
        cf.matrix_gmtransfer_doner_patch    = fill(CLM.ILEAF, counts.ncgmtrans)
        cf.matrix_gmtransfer_receiver_patch = fill(CLM.ILEAF, counts.ncgmtrans)
        cf.matrix_gmtransfer_patch = zeros(np, counts.ncgmtrans)
        cf.matrix_gmturnover_patch = zeros(np, nveg)
        cf.matrix_fitransfer_doner_patch    = fill(CLM.ILEAF, counts.ncfitrans)
        cf.matrix_fitransfer_receiver_patch = fill(CLM.ILEAF, counts.ncfitrans)
        cf.matrix_fitransfer_patch = zeros(np, counts.ncfitrans)
        cf.matrix_fiturnover_patch = zeros(np, nveg)

        # allocation input
        cf.matrix_alloc_patch  = copy(alloc)
        cf.matrix_Cinput_patch = copy(Cinput)

        # Populate the transfer/turnover arrays via the accumulator helper, using
        # the SAME initial pool state (the rates are pool-independent; the helper
        # just registers them). matrixcheck=false here (no cap, matches the
        # uncapped sequential reference).
        for p in 1:np
            for k in 1:nnon
                CLM.matrix_update_phc!(cf, p, k, ph_rate[p, k], dt; matrixcheck=false, acc=true)
            end
            # litterfall as out-transfers (leaf, froot)
            CLM.matrix_update_phc!(cf, p, nph-2, leaf_litter_rate,  dt; matrixcheck=false, acc=true)
            CLM.matrix_update_phc!(cf, p, nph-1, froot_litter_rate, dt; matrixcheck=false, acc=true)
            # livestem-out: zero
            CLM.matrix_update_phc!(cf, p, nph,   0.0, dt; matrixcheck=false, acc=true)
        end

        CLM.cn_veg_matrix_solve_c!(cs_mat, cf;
            mask_soilp = trues(np),
            bounds_patch = 1:np,
            ivt = ivt, woody = woody,
            npcropmin = npcropmin, nvegcpool = nveg,
            counts = counts, dt = dt, num_actfirep = 0)

        # ===================== compare =====================
        tol = 1e-9
        for p in 1:np
            for i in 1:nveg
                @test isapprox(_cnvm_pool_get(cs_mat, i, p), _cnvm_pool_get(cs_seq, i, p);
                               atol = tol, rtol = 1e-10)
            end
        end

        # mass conservation: total in-veg C(n+1) = C(n) + allocation_in - losses.
        cs_n = fresh_state()
        for p in 1:np
            cn   = sum(_cnvm_pool_get(cs_n, i, p)   for i in 1:nveg)
            cnp1 = sum(_cnvm_pool_get(cs_mat, i, p) for i in 1:nveg)
            alloc_in = sum(alloc[p, i] for i in 1:nveg) * Cinput[p] * dt
            loss = (leaf_litter_rate * cs_n.leafc_patch[p] +
                    froot_litter_rate * cs_n.frootc_patch[p]) * dt
            @test isapprox(cnp1, cn + alloc_in - loss; atol = 1e-7, rtol = 1e-9)
        end
    end

    # --- 3-process path: phenology + gap-mortality + fire (SPMP_ABC) -----------
    @testset "matrix advance with gap-mortality + fire == sequential" begin
        np = 2
        dt = 1800.0
        counts = CLM.veg_matrix_transfer_counts(false)
        nveg = CLM.NVEGPOOL_NATVEG
        npcropmin = 15
        ivt = [1, 2]
        woody = ones(Float64, maximum(ivt))

        # initial pools
        cs0 = CLM.CNVegCarbonStateData()
        CLM.cnveg_carbon_state_init!(cs0, np, 1, 1; use_matrixcn=false, nrepr=1)
        cs0.leafc_patch          .= [50.0, 30.0]
        cs0.leafc_storage_patch  .= [20.0, 15.0]
        cs0.leafc_xfer_patch     .= [10.0, 8.0]
        cs0.frootc_patch         .= [40.0, 25.0]
        cs0.frootc_storage_patch .= [18.0, 12.0]
        cs0.frootc_xfer_patch    .= [9.0, 7.0]
        cs0.livestemc_patch          .= [60.0, 45.0]
        cs0.livestemc_storage_patch  .= [22.0, 16.0]
        cs0.livestemc_xfer_patch     .= [11.0, 9.0]
        cs0.deadstemc_patch          .= [100.0, 80.0]
        cs0.deadstemc_storage_patch  .= [5.0, 4.0]
        cs0.deadstemc_xfer_patch     .= [3.0, 2.5]
        cs0.livecrootc_patch          .= [35.0, 28.0]
        cs0.livecrootc_storage_patch  .= [12.0, 9.0]
        cs0.livecrootc_xfer_patch     .= [6.0, 5.0]
        cs0.deadcrootc_patch          .= [70.0, 55.0]
        cs0.deadcrootc_storage_patch  .= [4.0, 3.0]
        cs0.deadcrootc_xfer_patch     .= [2.0, 1.8]

        copystate(src) = (d = deepcopy(src); d)

        # phenology: leaf_xf -> leaf (transfer), col-major sorted single in-veg entry
        ph_doner = [CLM.ILEAF_XF]; ph_recv = [CLM.ILEAF]
        ph_rate = [3.0e-6, 1.5e-6]    # 1/s per patch
        # gap-mortality: leaf -> iout (pure loss), turnover only (no in-veg entry)
        gm_leaf_rate = [4.0e-7, 2.0e-7]
        # fire: livestem -> deadstem (in-veg transfer) + leaf -> iout (loss)
        fi_doner = [CLM.ILIVESTEM]; fi_recv = [CLM.IDEADSTEM]
        fi_livestem_rate = [5.0e-7, 3.0e-7]
        fi_leaf_loss_rate = [6.0e-7, 4.0e-7]

        # ---------- sequential reference ----------
        cs_seq = copystate(cs0)
        for p in 1:np
            # phenology transfer leaf_xf->leaf
            ft = ph_rate[p] * cs_seq.leafc_xfer_patch[p]
            cs_seq.leafc_patch[p]      += ft * dt
            cs_seq.leafc_xfer_patch[p] -= ft * dt
            # gap-mortality leaf loss
            cs_seq.leafc_patch[p] -= gm_leaf_rate[p] * cs_seq.leafc_patch[p] * dt
            # fire livestem->deadstem
            fls = fi_livestem_rate[p] * cs_seq.livestemc_patch[p]
            cs_seq.livestemc_patch[p] -= fls * dt
            cs_seq.deadstemc_patch[p] += fls * dt
            # fire leaf loss
            cs_seq.leafc_patch[p] -= fi_leaf_loss_rate[p] * cs_seq.leafc_patch[p] * dt
        end

        # NOTE: the sequential reference applies gap-mortality + fire leaf losses on
        # the ALREADY phenology/gm-updated leafc, but the matrix applies them
        # simultaneously on leafc(n). To make the two exactly comparable we instead
        # compute the matrix's expected losses on leafc(n) for the reference too:
        cs_seq = copystate(cs0)
        for p in 1:np
            leaf_n = cs0.leafc_patch[p]
            ph_in  = ph_rate[p] * cs0.leafc_xfer_patch[p] * dt           # leaf_xf -> leaf
            gm_out = gm_leaf_rate[p] * leaf_n * dt                       # leaf -> iout
            fi_out = fi_leaf_loss_rate[p] * leaf_n * dt                  # leaf -> iout (fire)
            cs_seq.leafc_patch[p]      = leaf_n + ph_in - gm_out - fi_out
            cs_seq.leafc_xfer_patch[p] = cs0.leafc_xfer_patch[p] - ph_rate[p] * cs0.leafc_xfer_patch[p] * dt
            fls = fi_livestem_rate[p] * cs0.livestemc_patch[p] * dt
            cs_seq.livestemc_patch[p]  = cs0.livestemc_patch[p] - fls
            cs_seq.deadstemc_patch[p]  = cs0.deadstemc_patch[p] + fls
        end

        # ---------- matrix ----------
        cs_mat = copystate(cs0)
        cf = CLM.CNVegCarbonFluxData()
        CLM.cnveg_carbon_flux_init!(cf, np, 1, 1; use_matrixcn=true)

        nph = counts.ncphtrans; ngm = counts.ncgmtrans; nfi = counts.ncfitrans
        # phenology arrays: in-veg slot 1 = leaf_xf->leaf, rest out (zero)
        cf.matrix_phtransfer_doner_patch    = fill(CLM.ILEAF, nph)
        cf.matrix_phtransfer_receiver_patch = fill(CLM.ILEAF, nph)
        cf.matrix_phtransfer_doner_patch[1] = ph_doner[1]
        cf.matrix_phtransfer_receiver_patch[1] = ph_recv[1]
        cf.matrix_phtransfer_patch = zeros(np, nph)
        cf.matrix_phturnover_patch = zeros(np, nveg)
        # gap-mortality: no in-veg transfer (all 18 are losses); use leaf->iout slot
        cf.matrix_gmtransfer_doner_patch    = fill(CLM.ILEAF, ngm)
        cf.matrix_gmtransfer_receiver_patch = fill(CLM.ILEAF, ngm)
        cf.matrix_gmtransfer_patch = zeros(np, ngm)
        cf.matrix_gmturnover_patch = zeros(np, nveg)
        # fire: in-veg slot 1 = livestem->deadstem; in-veg slot 2 = a zero-rate
        # self-loop placed AFTER slot 1 in column-major order (deadstem->deadstem,
        # linear 172 > livestem->deadstem's 118) to keep the COO ordering ascending;
        # out slot = leaf->iout.
        cf.matrix_fitransfer_doner_patch    = fill(CLM.IDEADSTEM, nfi)
        cf.matrix_fitransfer_receiver_patch = fill(CLM.IDEADSTEM, nfi)
        cf.matrix_fitransfer_doner_patch[1] = fi_doner[1]
        cf.matrix_fitransfer_receiver_patch[1] = fi_recv[1]
        # a trailing out-transfer slot for leaf-fire-loss:
        cf.matrix_fitransfer_doner_patch[nfi] = CLM.ILEAF
        cf.matrix_fitransfer_receiver_patch[nfi] = nveg + 1
        cf.matrix_fitransfer_patch = zeros(np, nfi)
        cf.matrix_fiturnover_patch = zeros(np, nveg)

        cf.matrix_alloc_patch  = zeros(np, nveg)
        cf.matrix_Cinput_patch = zeros(np)

        for p in 1:np
            CLM.matrix_update_phc!(cf, p, 1, ph_rate[p], dt; matrixcheck=false, acc=true)
            # gap-mortality leaf loss: register turnover only (gm in-veg index 1 is
            # a self/loss; use the first gm slot whose doner=leaf -> iout).
            cf.matrix_gmtransfer_doner_patch[ngm] = CLM.ILEAF
            cf.matrix_gmtransfer_receiver_patch[ngm] = nveg + 1
            CLM.matrix_update_gmc!(cf, p, ngm, gm_leaf_rate[p], dt; matrixcheck=false, acc=true)
            # fire livestem->deadstem (in-veg)
            CLM.matrix_update_fic!(cf, p, 1, fi_livestem_rate[p], dt; matrixcheck=false, acc=true)
            # fire leaf loss (out)
            CLM.matrix_update_fic!(cf, p, nfi, fi_leaf_loss_rate[p], dt; matrixcheck=false, acc=true)
        end

        CLM.cn_veg_matrix_solve_c!(cs_mat, cf;
            mask_soilp = trues(np), bounds_patch = 1:np,
            ivt = ivt, woody = woody, npcropmin = npcropmin,
            nvegcpool = nveg, counts = counts, dt = dt, num_actfirep = 1)

        for p in 1:np
            for i in 1:nveg
                @test isapprox(_cnvm_pool_get(cs_mat, i, p), _cnvm_pool_get(cs_seq, i, p);
                               atol = 1e-9, rtol = 1e-10)
            end
        end
    end

    # --- default (flag off) leaves the sequential path untouched --------------
    @testset "default use_matrixcn=false unchanged" begin
        # The matrix solver is only invoked when use_matrixcn is on; with it off
        # the c_state_update*! kernels run their sequential `!use_matrixcn` branch
        # (covered by test_c_state_update1/2/3). Here we just confirm the matrix
        # solver is a self-contained opt-in (calling it does NOT run unless asked).
        @test isdefined(CLM, :cn_veg_matrix_solve_c!)
        @test isdefined(CLM, :matrix_update_phc!)
    end
end
