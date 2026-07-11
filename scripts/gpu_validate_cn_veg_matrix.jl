# ==========================================================================
# gpu_validate_cn_veg_matrix.jl — focused REAL-Metal parity for the matrix-CN
# VEGETATION solves (cn_veg_matrix_solve_c! / cn_veg_matrix_solve_n!) driven the
# way the CN driver drives them: accumulate the transfer/turnover matrices from
# the LIVE per-flux fields (cn_veg_matrix_accumulate_{ph,gm,fi}_{c,n}!), then
# advance the pools in one matrix solve.
#
# This is a stronger / different proof than gpu_validate_matrixcn_e2e.jl: that one
# compares a matrix-CPU solve to a matrix-Metal solve on SYNTHETIC uniform flux
# fills. Here we
#   (C1) build a realistic 3-patch veg state + realistic phenology/gap/fire/alloc
#        fluxes, and prove the HOST Float64 matrix solve reproduces the SEQUENTIAL
#        (ODE) c_state_update result bit-close — i.e. matrix == sequential, which is
#        the whole point of the matrix method (byte-identity of the KA.CPU path);
#   (C2) move the same state + accumulated flux matrices + warm-up solve-structure
#        to the device (Metal, Float32) and prove the DEVICE matrix solve matches
#        the HOST matrix solve within ~1e-6 — for BOTH veg-C and veg-N.
# Transitively the device matrix solve reproduces the sequential result.
#
#   julia --project=scripts scripts/gpu_validate_cn_veg_matrix.jl
# ==========================================================================
using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))

# Float32-converting device adaptor for whole state/flux/solve-state structs: floats
# -> Float32 on device, ints/bools -> device (type preserved). Metal has no Float64.
struct DevF32{F}; f::F; end
CLM.Adapt.adapt_storage(d::DevF32, x::AbstractArray{<:AbstractFloat}) = d.f(Float32.(x))
CLM.Adapt.adapt_storage(d::DevF32, x::AbstractArray{<:Integer}) = d.f(collect(x))
CLM.Adapt.adapt_storage(d::DevF32, x::AbstractArray{Bool}) = d.f(collect(x))

reldiff(a, b) = begin
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    m
end

const NP = 3
const DT = 1800.0
const NVEG = CLM.NVEGPOOL_NATVEG
const COUNTS = CLM.veg_matrix_transfer_counts(false)
const IVT = [1, 2, 3]
const WOODY = ones(Float64, 5)
const NPCROPMIN = 15

# ---- realistic veg-C fixture (mirrors test_cn_veg_matrix_wiring) ---------------
const _DONER = CLM._PHC_DONER_NATVEG
const _RECV  = CLM._PHC_RECEIVER_NATVEG
const FLDLIST = [:leafc_storage_to_xfer_patch, :leafc_xfer_to_leafc_patch,
    :frootc_storage_to_xfer_patch, :frootc_xfer_to_frootc_patch,
    :livestemc_storage_to_xfer_patch, :livestemc_xfer_to_livestemc_patch,
    :deadstemc_storage_to_xfer_patch, :deadstemc_xfer_to_deadstemc_patch,
    :livecrootc_storage_to_xfer_patch, :livecrootc_xfer_to_livecrootc_patch,
    :deadcrootc_storage_to_xfer_patch, :deadcrootc_xfer_to_deadcrootc_patch,
    :livestemc_to_deadstemc_patch, :livecrootc_to_deadcrootc_patch,
    :leafc_to_litter_patch, :frootc_to_litter_patch, :livestemc_to_litter_patch]
const RATES = [
    5e-6 4e-6 3e-6 2e-6 1e-6 8e-7 6e-7 4e-7 3e-7 2e-7 1e-7 9e-8 5e-7 4e-7 3e-7 2e-7 0.0;
    2e-6 6e-6 1e-6 5e-6 7e-7 9e-7 5e-7 3e-7 2e-7 4e-7 8e-8 6e-8 6e-7 3e-7 2e-7 3e-7 0.0;
    1e-6 1e-6 1e-6 1e-6 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 0.0]
const GM_RATE = [1e-6,2e-6,3e-6,1.5e-6,2.5e-6,1e-6,8e-7,6e-7,4e-7,3e-7,2e-7,1e-7,
                 7e-7,5e-7,3e-7,2e-7,1e-7,9e-8]
const FI_OUT_RATE = [1e-6,1.5e-6,2e-6,8e-7,1e-6,1.2e-6,7e-7,5e-7,3e-7,2e-7,1e-7,9e-8,
                     6e-7,4e-7,2e-7,1e-7,8e-8,6e-8]
const FI_LS_RATE = 3e-7; const FI_LC_RATE = 2e-7
const ALLOC_FLUX = Dict(CLM.ILEAF=>3e-5, CLM.ILEAF_ST=>1e-5, CLM.IFROOT=>2.5e-5, CLM.IFROOT_ST=>8e-6,
    CLM.ILIVESTEM=>2e-5, CLM.ILIVESTEM_ST=>6e-6, CLM.IDEADSTEM=>1e-5, CLM.IDEADSTEM_ST=>4e-6,
    CLM.ILIVECROOT=>1.5e-5, CLM.ILIVECROOT_ST=>5e-6, CLM.IDEADCROOT=>8e-6, CLM.IDEADCROOT_ST=>3e-6)

const POOLFIELD = Dict(
    CLM.ILEAF=>:leafc_patch, CLM.ILEAF_ST=>:leafc_storage_patch, CLM.ILEAF_XF=>:leafc_xfer_patch,
    CLM.IFROOT=>:frootc_patch, CLM.IFROOT_ST=>:frootc_storage_patch, CLM.IFROOT_XF=>:frootc_xfer_patch,
    CLM.ILIVESTEM=>:livestemc_patch, CLM.ILIVESTEM_ST=>:livestemc_storage_patch, CLM.ILIVESTEM_XF=>:livestemc_xfer_patch,
    CLM.IDEADSTEM=>:deadstemc_patch, CLM.IDEADSTEM_ST=>:deadstemc_storage_patch, CLM.IDEADSTEM_XF=>:deadstemc_xfer_patch,
    CLM.ILIVECROOT=>:livecrootc_patch, CLM.ILIVECROOT_ST=>:livecrootc_storage_patch, CLM.ILIVECROOT_XF=>:livecrootc_xfer_patch,
    CLM.IDEADCROOT=>:deadcrootc_patch, CLM.IDEADCROOT_ST=>:deadcrootc_storage_patch, CLM.IDEADCROOT_XF=>:deadcrootc_xfer_patch)
poolget(cs, i, p) = CLM._vegc_pool_val(cs, i, p)
pooladd!(cs, i, p, v) = (haskey(POOLFIELD, i) && (getfield(cs, POOLFIELD[i])[p] += v))

function fresh_cs()
    cs = CLM.CNVegCarbonStateData(); CLM.cnveg_carbon_state_init!(cs, NP, 1, 1; use_matrixcn=false, nrepr=1)
    cs.leafc_patch .= [50.,30.,10.];  cs.leafc_storage_patch .= [20.,15.,5.];  cs.leafc_xfer_patch .= [10.,8.,3.]
    cs.frootc_patch .= [40.,25.,12.]; cs.frootc_storage_patch .= [18.,12.,6.]; cs.frootc_xfer_patch .= [9.,7.,4.]
    cs.livestemc_patch .= [60.,45.,20.]; cs.livestemc_storage_patch .= [22.,16.,8.]; cs.livestemc_xfer_patch .= [11.,9.,5.]
    cs.deadstemc_patch .= [100.,80.,40.]; cs.deadstemc_storage_patch .= [5.,4.,2.]; cs.deadstemc_xfer_patch .= [3.,2.5,1.]
    cs.livecrootc_patch .= [35.,28.,14.]; cs.livecrootc_storage_patch .= [12.,9.,5.]; cs.livecrootc_xfer_patch .= [6.,5.,2.5]
    cs.deadcrootc_patch .= [70.,55.,28.]; cs.deadcrootc_storage_patch .= [4.,3.,1.5]; cs.deadcrootc_xfer_patch .= [2.,1.8,.8]
    return cs
end
function set_c_fluxes!(cf, cs)
    for (k, f) in enumerate(FLDLIST)
        length(getfield(cf, f)) == 0 && setfield!(cf, f, zeros(NP))
        arr = getfield(cf, f); for p in 1:NP; arr[p] = poolget(cs, _DONER[k], p) * RATES[p, k]; end
    end
    for i in 1:NVEG
        f = CLM._GMC_FLUX[i]; length(getfield(cf, f)) == 0 && setfield!(cf, f, zeros(NP))
        arr = getfield(cf, f); for p in 1:NP; arr[p] = poolget(cs, i, p) * GM_RATE[i]; end
    end
    for i in 1:NVEG
        f = CLM._FIC_TO_FIRE[i]; length(getfield(cf, f)) == 0 && setfield!(cf, f, zeros(NP))
        fl = CLM._FIC_TO_LITTER[i]; length(getfield(cf, fl)) == 0 && setfield!(cf, fl, zeros(NP))
        fill!(getfield(cf, fl), 0.0)
        arr = getfield(cf, f); for p in 1:NP; arr[p] = poolget(cs, i, p)*FI_OUT_RATE[i]; end
    end
    length(cf.m_livestemc_to_deadstemc_fire_patch)==0 && (cf.m_livestemc_to_deadstemc_fire_patch=zeros(NP))
    length(cf.m_livecrootc_to_deadcrootc_fire_patch)==0 && (cf.m_livecrootc_to_deadcrootc_fire_patch=zeros(NP))
    for p in 1:NP
        cf.m_livestemc_to_deadstemc_fire_patch[p]   = poolget(cs, CLM.ILIVESTEM, p)*FI_LS_RATE
        cf.m_livecrootc_to_deadcrootc_fire_patch[p] = poolget(cs, CLM.ILIVECROOT, p)*FI_LC_RATE
    end
    for (pool, fld) in CLM._ALLOC_TARGET
        length(getfield(cf, fld)) == 0 && setfield!(cf, fld, zeros(NP))
        arr = getfield(cf, fld); for p in 1:NP; arr[p] = ALLOC_FLUX[pool]; end
    end
    return nothing
end

# Sequential (ODE) reference: apply allocation + phenology transfers + gap-mortality
# + fire the direct way (== c_state_update at matrixcheck=false).
function sequential_c_reference()
    cs = fresh_cs()
    cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf, NP, 1, 1; use_matrixcn=true)
    set_c_fluxes!(cf, cs)
    for p in 1:NP, (pool, fld) in CLM._ALLOC_TARGET
        pooladd!(cs, pool, p, getfield(cf, fld)[p] * DT)
    end
    for p in 1:NP, k in 1:17
        flux = getfield(cf, FLDLIST[k])[p]
        pooladd!(cs, _RECV[k], p,  flux*DT); pooladd!(cs, _DONER[k], p, -flux*DT)
    end
    for p in 1:NP, i in 1:NVEG
        pooladd!(cs, i, p, -getfield(cf, CLM._GMC_FLUX[i])[p] * DT)
    end
    for p in 1:NP
        fls = cf.m_livestemc_to_deadstemc_fire_patch[p]; flc = cf.m_livecrootc_to_deadcrootc_fire_patch[p]
        pooladd!(cs, CLM.ILIVESTEM, p, -fls*DT);  pooladd!(cs, CLM.IDEADSTEM, p, fls*DT)
        pooladd!(cs, CLM.ILIVECROOT, p, -flc*DT); pooladd!(cs, CLM.IDEADCROOT, p, flc*DT)
        for i in 1:NVEG
            pooladd!(cs, i, p, -(getfield(cf, CLM._FIC_TO_FIRE[i])[p]+getfield(cf, CLM._FIC_TO_LITTER[i])[p])*DT)
        end
    end
    return cs
end

# Build cf with the transfer/turnover matrices accumulated from the live flux fields,
# exactly as the driver does (cn_veg_matrix_accumulate_*!).
function accumulate_c_flux()
    cs0 = fresh_cs()
    cf = CLM.CNVegCarbonFluxData(); CLM.cnveg_carbon_flux_init!(cf, NP, 1, 1; use_matrixcn=true)
    set_c_fluxes!(cf, cs0)
    CLM.cn_veg_matrix_c_topology!(cf; use_crop=false, nvegcpool=NVEG)
    CLM.cn_veg_matrix_alloc_c!(cf, trues(NP), 1:NP)
    CLM.cn_veg_matrix_accumulate_ph_c!(cf, cs0, trues(NP), 1:NP; dt=DT, matrixcheck=false)
    CLM.cn_veg_matrix_accumulate_gm_c!(cf, cs0, trues(NP), 1:NP; dt=DT, matrixcheck=false)
    CLM.cn_veg_matrix_accumulate_fi_c!(cf, cs0, trues(NP), 1:NP; dt=DT, matrixcheck=false)
    return cf
end

Cargs(ivt_) = (; mask_soilp=trues(NP), bounds_patch=1:NP, ivt=ivt_, woody=WOODY,
               npcropmin=NPCROPMIN, nvegcpool=NVEG, counts=COUNTS, dt=DT, num_actfirep=NP)
cpools(cs) = Float64[Array(cs.leafc_patch); Array(cs.leafc_storage_patch); Array(cs.leafc_xfer_patch);
    Array(cs.frootc_patch); Array(cs.livestemc_patch); Array(cs.deadstemc_patch);
    Array(cs.livecrootc_patch); Array(cs.deadcrootc_patch)]

# ---- realistic veg-N fixture (mirrors test_cn_veg_matrix_wiring_n's flux setup) --
const NVN = CLM.IRETRANSN_NATVEG
const NFLDS = (:leafn_patch, :leafn_storage_patch, :leafn_xfer_patch, :frootn_patch, :frootn_storage_patch,
    :frootn_xfer_patch, :livestemn_patch, :livestemn_storage_patch, :livestemn_xfer_patch,
    :deadstemn_patch, :deadstemn_storage_patch, :deadstemn_xfer_patch, :livecrootn_patch,
    :livecrootn_storage_patch, :livecrootn_xfer_patch, :deadcrootn_patch, :deadcrootn_storage_patch,
    :deadcrootn_xfer_patch, :retransn_patch)
function fresh_ns()
    ns = CLM.CNVegNitrogenStateData(); CLM.cnveg_nitrogen_state_init!(ns, NP, 1, 1; nrepr=1)
    for (k, f) in enumerate(NFLDS); getfield(ns, f) .= Float64[5.0k + p for p in 1:NP]; end
    return ns
end
function accumulate_n_flux()
    nf = CLM.CNVegNitrogenFluxData(); CLM.cnveg_nitrogen_flux_init!(nf, NP, 1, 1)
    CLM.cn_veg_matrix_n_topology!(nf; use_crop=false, nvegnpool=NVN)
    # Per-pool-varying transfer/turnover fills so the solve is nontrivial (each pool has a
    # distinct rate; the device solve must reproduce the host solve entry-for-entry).
    for k in 1:size(nf.matrix_nphtransfer_patch, 2), p in 1:NP
        nf.matrix_nphtransfer_patch[p, k] = 1.0e-6 * (1 + 0.1k) * (1 + 0.05p)
    end
    for k in 1:size(nf.matrix_ngmtransfer_patch, 2), p in 1:NP
        nf.matrix_ngmtransfer_patch[p, k] = 5.0e-7 * (1 + 0.1k) * (1 + 0.05p)
    end
    fill!(nf.matrix_nfitransfer_patch, 2.0e-7)
    for i in 1:NVN, p in 1:NP
        nf.matrix_nphturnover_patch[p, i] = 0.05 * (1 + 0.02i)
        nf.matrix_ngmturnover_patch[p, i] = 0.02 * (1 + 0.02i)
        nf.matrix_nfiturnover_patch[p, i] = 0.01
    end
    fill!(nf.matrix_nalloc_patch, 0.0); fill!(nf.matrix_Ninput_patch, 0.0)
    return nf
end
Nargs(ivt_) = (; mask_soilp=trues(NP), bounds_patch=1:NP, ivt=ivt_, npcropmin=NPCROPMIN,
               nvegnpool=NVN, counts=COUNTS, dt=DT, num_actfirep=NP)
npools(ns) = Float64[Array(ns.leafn_patch); Array(ns.frootn_patch); Array(ns.livestemn_patch);
    Array(ns.deadstemn_patch); Array(ns.livecrootn_patch); Array(ns.deadcrootn_patch); Array(ns.retransn_patch)]

# ==========================================================================
function main(backend)
    println("=" ^ 68)
    println("Matrix-CN VEG solve — realistic live-flux parity (host + device)")
    println("=" ^ 68)
    if backend === nothing
        println("  No GPU backend available — running host matrix==sequential only.")
    end
    name, dev, FT = backend === nothing ? ("none", identity, Float32) : backend
    @printf("  Backend: %-8s (precision: %s)\n", name, FT)
    mvf(x) = CLM.Adapt.adapt(DevF32(dev), x)
    checks = Tuple{String,Float64,Bool}[]   # (name, metric, pass)

    # ---------- veg-C ----------
    cf_c = accumulate_c_flux()
    # (C1) HOST matrix solve == sequential ODE reference (Float64, matrix == sequential).
    cs_seq = sequential_c_reference()
    cs_host = fresh_cs(); CLM.cn_veg_matrix_solve_c!(cs_host, cf_c; Cargs(IVT)...)
    d_seq = 0.0
    for (a, b) in zip(cpools(cs_host), cpools(cs_seq)); d_seq = max(d_seq, abs(a - b)); end
    push!(checks, ("veg-C host matrix == sequential (abs)", d_seq, d_seq < 1e-8))

    # ---------- veg-N host solve (needed for the device compare below) ----------
    nf = accumulate_n_flux()
    ns_host = fresh_ns(); CLM.cn_veg_matrix_solve_n!(ns_host, nf; Nargs(IVT)...)

    if backend !== nothing
        # (C2) DEVICE veg-C matrix solve == host matrix solve (Float32, within ~1e-6).
        st_c = CLM.CNVegMatrixSolveState()
        CLM.cn_veg_matrix_solve_c!(fresh_cs(), cf_c; Cargs(IVT)..., state=st_c)   # host warm-up: build structure
        cs_d = mvf(fresh_cs()); cf_cd = mvf(cf_c); st_cd = mvf(st_c)
        CLM.cn_veg_matrix_solve_c!(cs_d, cf_cd; Cargs(dev(IVT))..., ref=dev(zeros(FT, 1, 1)), FT=FT, state=st_cd)
        r_c = reldiff(cpools(cs_host), cpools(cs_d))
        push!(checks, ("veg-C device == host matrix (rel)", r_c, r_c < 1e-6))

        # (N2) DEVICE veg-N matrix solve == host matrix solve.
        st_n = CLM.CNVegMatrixSolveState()
        CLM.cn_veg_matrix_solve_n!(fresh_ns(), nf; Nargs(IVT)..., state=st_n)     # host warm-up
        ns_d = mvf(fresh_ns()); nf_d = mvf(nf); st_nd = mvf(st_n)
        CLM.cn_veg_matrix_solve_n!(ns_d, nf_d; Nargs(dev(IVT))..., ref=dev(zeros(FT, 1, 1)), FT=FT, state=st_nd)
        r_n = reldiff(npools(ns_host), npools(ns_d))
        push!(checks, ("veg-N device == host matrix (rel)", r_n, r_n < 1e-6))
    end

    println()
    allpass = true
    for (name, metric, pass) in checks
        allpass &= pass
        @printf("  [%s] %-38s %.3e\n", pass ? "PASS" : "FAIL", name, metric)
    end
    println()
    if allpass
        println("  matrix-CN VEG solve MATCHES on $(name): host matrix==sequential + device==host (veg-C & veg-N)")
    else
        println("  FAILURES above")
        exit(1)
    end
    return 0
end

exit(main(detect_backend()))
