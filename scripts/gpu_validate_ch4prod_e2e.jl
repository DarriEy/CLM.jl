# ==========================================================================
# gpu_validate_ch4prod_e2e.jl — end-to-end GPU parity for the WHOLE
# ch4_prod! methane-production driver (the C/N decomposition -> CH4 production
# cascade). Two kernels run on the device:
#   1. _ch4prod_rrvr_kernel!  — per-PATCH scatter of root respiration into the
#      column-resolved rr_vr[c,j] (many patches -> one column, _scatter_add!).
#   2. _ch4prod_main_kernel!  — per-(column, level) production reading rr_vr,
#      with all CH4Data views bundled into _Ch4ProdState and scalar params/flags
#      bundled into _Ch4ProdScal / _Ch4ProdFlags.
#
# Builds a small instance, runs ch4_prod! on the CPU (Float64) and on Metal
# (Float32, fields adapted field-by-field — CH4Data is not @adapt_structure'd in
# this base), and compares the three mutated outputs:
#   ch4_prod_depth_unsat_col / o2_decomp_depth_unsat_col / sif_col
# across several config branches (sat 0/1, use_cn, anoxia, usephfact,
# ch4rmcnlim, anoxicmicrosites, use_nitrif_denitrif) with lake=false.
#
#   julia --project=scripts scripts/gpu_validate_ch4prod_e2e.jl
# ==========================================================================

using CLM
using Printf
include(joinpath(@__DIR__, "gpu_backends.jl"))

# reldiff: NaN-aware (both-NaN agrees; one-sided NaN flags divergence).
function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    return m
end
cpu_all_finite(a) = all(isfinite, Array(a))

const NC = 4
const NP = 6
const NLEV = 6

# Build the column/patch arrays + a CH4Data with the fields ch4_prod! touches, at
# precision FT.
function build(::Type{FT}) where {FT}
    nc = NC; np = NP; nlevsoi = NLEV

    ch4 = CLM.CH4Data{FT}(
        ch4_prod_depth_sat_col   = zeros(FT, nc, nlevsoi),
        ch4_prod_depth_unsat_col = zeros(FT, nc, nlevsoi),
        o2_decomp_depth_sat_col  = zeros(FT, nc, nlevsoi),
        o2_decomp_depth_unsat_col = zeros(FT, nc, nlevsoi),
        conc_o2_sat_col          = fill(FT(0.01), nc, nlevsoi),
        conc_o2_unsat_col        = fill(FT(0.05), nc, nlevsoi),
        lake_soilc_col           = fill(FT(100.0), nc, nlevsoi),
        layer_sat_lag_col        = fill(FT(0.5), nc, nlevsoi),
        sif_col                  = ones(FT, nc),
        annavg_finrw_col         = fill(FT(0.1), nc),
        finundated_col           = fill(FT(0.12), nc),
        finundated_lag_col       = fill(FT(0.08), nc),
        pH_col                   = FT[6.5, 5.0, 7.2, 8.4],
    )

    mask_soil  = trues(nc)
    mask_soilp = trues(np)
    patch_column = [1, 1, 2, 3, 4, 4]
    patch_itype  = [1, 2, 1, 0, 1, 1]   # itype 0 == noveg for patch 4 (skipped scatter)
    patch_wtcol  = FT[0.5, 0.5, 1.0, 1.0, 0.3, 0.7]
    is_fates     = falses(nc)

    crootfr   = fill(FT(0.2), np, nlevsoi)
    rootfr_col = fill(FT(0.2), nc, nlevsoi)
    watsat    = fill(FT(0.45), nc, nlevsoi)
    h2osoi_vol = fill(FT(0.3), nc, nlevsoi)
    t_soisno  = fill(FT(CLM.TFRZ + 15.0), nc, nlevsoi)
    somhr     = fill(FT(1.0e-6), nc)
    lithr     = fill(FT(1.0e-7), nc)
    hr_vr     = fill(FT(1.0e-8), nc, nlevsoi)
    o_scalar  = fill(FT(0.8), nc, nlevsoi)
    fphr      = fill(FT(0.9), nc, nlevsoi)
    pot_f_nit_vr = fill(FT(1.0e-9), nc, nlevsoi)
    rr        = fill(FT(1.0e-6), np)
    jwt       = [2, 3, 1, 4]

    dz = fill(FT(0.1), nc, nlevsoi)
    z  = zeros(FT, nc, nlevsoi)
    zi = zeros(FT, nc, nlevsoi)
    for j in 1:nlevsoi, c in 1:nc
        z[c, j]  = FT((j - 0.5) * 0.1)
        zi[c, j] = FT(j * 0.1)
    end

    return (; ch4, nc, np, nlevsoi, mask_soil, mask_soilp, patch_column, patch_itype,
            patch_wtcol, is_fates, crootfr, rootfr_col, watsat, h2osoi_vol, t_soisno,
            somhr, lithr, hr_vr, o_scalar, fphr, pot_f_nit_vr, rr, jwt, dz, z, zi)
end

# Run ch4_prod! for a given config tuple on data bundle B.
function run!(B, params, ch4vc, cfg)
    CLM.ch4_prod!(B.ch4, params, ch4vc, B.mask_soil, B.mask_soilp,
        B.patch_column, B.patch_itype, B.patch_wtcol, B.is_fates,
        B.crootfr, B.rootfr_col, B.watsat, B.h2osoi_vol, B.t_soisno,
        B.somhr, B.lithr, B.hr_vr, B.o_scalar, B.fphr, B.pot_f_nit_vr,
        B.rr, B.jwt, cfg.sat, false,
        B.dz, B.z, B.zi, B.nlevsoi, cfg.nlevdecomp, B.nlevsoi, 5, 0.01,
        1800.0, 0, cfg.use_cn, cfg.use_nitrif_denitrif, cfg.anoxia)
end

# Configurations exercising the Bool/Int branches (lake=false throughout).
const CONFIGS = [
    (; name="sat0 cn nd0 anox0",  sat=0, nlevdecomp=1, use_cn=true,  use_nitrif_denitrif=false, anoxia=false,
       usephfact=false, ch4rmcnlim=true,  anoxicmicrosites=true),
    (; name="sat1 cn nd1 anox1",  sat=1, nlevdecomp=NLEV, use_cn=true,  use_nitrif_denitrif=true,  anoxia=true,
       usephfact=true,  ch4rmcnlim=true,  anoxicmicrosites=true),
    (; name="sat0 cn ph rmcn0",   sat=0, nlevdecomp=NLEV, use_cn=true,  use_nitrif_denitrif=false, anoxia=true,
       usephfact=true,  ch4rmcnlim=false, anoxicmicrosites=false),
    (; name="sat1 nocn",          sat=1, nlevdecomp=NLEV, use_cn=false, use_nitrif_denitrif=true,  anoxia=false,
       usephfact=false, ch4rmcnlim=true,  anoxicmicrosites=true),
]

function main(backend)
    println("=" ^ 72)
    println("END-TO-END GPU parity for ch4_prod!  (CH4 production driver)")
    println("=" ^ 72)
    if backend === nothing
        println("  No GPU backend — nothing to validate (CPU path exercised by the suite).")
        return 0
    end
    name, to, FT = backend
    @printf("  Backend: %s   (working precision: %s)\n\n", name, FT)

    nfail = 0
    for cfg in CONFIGS
        params = CLM.CH4Params()
        ch4vc  = CLM.CH4VarCon(usephfact=cfg.usephfact, ch4rmcnlim=cfg.ch4rmcnlim,
                               anoxicmicrosites=cfg.anoxicmicrosites)

        # CPU reference (Float64)
        Bc = build(Float64)
        run!(Bc, params, ch4vc, cfg)

        # Device run (FT). Adapt the WHOLE CH4Data to the device (it is
        # @adapt_structure'd on this branch — the foundation), so every field the
        # device-view bundle reads is a device array. Building it via CH4Data{FT}(...)
        # would re-pin V=Vector{FT}/M=Matrix{FT} and pull the arrays back to the host.
        Bf = build(FT)
        ad(x) = _dev(to, x)
        dch4 = CLM.Adapt.adapt(device_array_type(), Bf.ch4)
        D = (; ch4=dch4, nc=Bf.nc, np=Bf.np, nlevsoi=Bf.nlevsoi,
             mask_soil=ad(Bf.mask_soil), mask_soilp=ad(Bf.mask_soilp),
             patch_column=ad(Bf.patch_column), patch_itype=ad(Bf.patch_itype),
             patch_wtcol=ad(Bf.patch_wtcol), is_fates=ad(Bf.is_fates),
             crootfr=ad(Bf.crootfr), rootfr_col=ad(Bf.rootfr_col),
             watsat=ad(Bf.watsat), h2osoi_vol=ad(Bf.h2osoi_vol), t_soisno=ad(Bf.t_soisno),
             somhr=ad(Bf.somhr), lithr=ad(Bf.lithr), hr_vr=ad(Bf.hr_vr),
             o_scalar=ad(Bf.o_scalar), fphr=ad(Bf.fphr), pot_f_nit_vr=ad(Bf.pot_f_nit_vr),
             rr=ad(Bf.rr), jwt=ad(Bf.jwt), dz=ad(Bf.dz), z=ad(Bf.z), zi=ad(Bf.zi))
        run!(D, params, ch4vc, cfg)

        # outputs depend on sat: pick the right arrays
        if cfg.sat == 0
            cp_cpu = Bc.ch4.ch4_prod_depth_unsat_col;  cp_dev = D.ch4.ch4_prod_depth_unsat_col
            o2_cpu = Bc.ch4.o2_decomp_depth_unsat_col; o2_dev = D.ch4.o2_decomp_depth_unsat_col
        else
            cp_cpu = Bc.ch4.ch4_prod_depth_sat_col;  cp_dev = D.ch4.ch4_prod_depth_sat_col
            o2_cpu = Bc.ch4.o2_decomp_depth_sat_col; o2_dev = D.ch4.o2_decomp_depth_sat_col
        end

        for (nm, a) in (("ch4_prod_depth", cp_cpu), ("o2_decomp_depth", o2_cpu),
                        ("sif_col", Bc.ch4.sif_col))
            cpu_all_finite(a) || @printf("  WARN[%s]: CPU ref %s not finite\n", cfg.name, nm)
        end

        checks = (
            ("ch4_prod_depth", cp_cpu, cp_dev),
            ("o2_decomp_depth", o2_cpu, o2_dev),
            ("sif_col",         Bc.ch4.sif_col, D.ch4.sif_col),
        )
        @printf("  --- config: %s ---\n", cfg.name)
        for (nm, a, b) in checks
            d = reldiff(a, b)
            ok = d < 5e-3
            @printf("    [%s] %-16s reldiff = %.3e\n", ok ? "PASS" : "FAIL", nm, d)
            ok || (nfail += 1)
        end
    end

    println()
    println(nfail == 0 ?
        "  WHOLE ch4_prod! MATCHES CPU ON $name ($FT)" :
        "  DIVERGENCE — investigate.")
    return nfail == 0 ? 0 : 1
end

const BACKEND = detect_backend()
exit(main(BACKEND))
