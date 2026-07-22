# ==========================================================================
# validate_mimics_live_init.jl — LIVE end-to-end activation of MIMICS through
# the production initializer clm_initialize!(decomp_method=2).
#
# Complements scripts/validate_mimics_fortran_parity.jl (the CI-safe oracle +
# config-derivation harness) by exercising the REAL init path on the Bow fixture:
# proves clm_initialize! sizes the soil-BGC state to 8 pools / 15 transitions,
# runs init_decompcascade_mimics! (populating the microbial cascade + texture
# params + the nue vector on the soil-BGC state), and reads the clm50_params
# mimics_* values — none of which the old dead-wired path did.
#
# Requires the local Bow fixture data (surfdata + clm5_params.nc); skips cleanly
# if absent (CI). Run:
#   julia +1.12 --project=. scripts/validate_mimics_live_init.jl
# ==========================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

if !isfile(FSURDAT) || !isfile(FPARAM)
    println("SKIP: Bow fixture data not present ($(FSURDAT)) — live-init check needs local data.")
    exit(0)
end

(inst, bounds, filt, tm) = build_bow_inst(; dtime=3600, use_cn=true, use_luna=true,
                                          decomp_method=2)
cs   = inst.soilbiogeochem_carbonstate
st   = inst.decomp_mimics_state
casc = inst.decomp_cascade
nue  = inst.soilbiogeochem_state.nue_decomp_cascade_col
prm  = inst.decomp_mimics_params

ok = true
check(name, got, want) = (p = got == want; global ok &= p;
    println("  ", rpad(name, 34), " = ", got, p ? "  PASS" : "  FAIL (want $want)"))

println("── LIVE clm_initialize!(decomp_method=2) MIMICS activation ──")
check("decomp_cpools_vr pools (dim3)", size(cs.decomp_cpools_vr_col, 3), 8)
check("i_cop_mic", st.i_cop_mic, 6)
check("i_oli_mic", st.i_oli_mic, 7)
check("cascade transitions", length(casc.cascade_donor_pool), 15)
check("pool short names", casc.decomp_pool_name_short,
      ["L1","L2","S1","S2","S3","M1","M2","CWD"])
check("nue length", length(nue), 15)
check("nue[1] (into mic)", nue[1], 0.85)
check("nue[7] (S2S1=1)", nue[7], 1.0)
check("initial_stock", casc.initial_stock, [1.0,1.0,200.0,200.0,200.0,1.0,1.0,1.0])
check("param densdep (clm50)", prm.mimics_densdep, 1.2)
check("param vint[1] (clm50)", prm.mimics_vint[1], 6.6)

# Seeded cold-start pools must be finite in the ACTIVE soil columns × decomp
# levels (init_cold_biogeochem seeds 8 pools there). Inactive columns and the
# padded levels nlevdecomp+1..nlevdecomp_full stay at the allocator NaN — this is
# IDENTICAL to the CENTURY cold start (verified: CENTURY 140/140, MIMICS 160/160
# active-finite), so checking the whole padded array would be a false failure.
let dk = cs.decomp_cpools_vr_col, nlev = CLM.varpar.nlevdecomp, msoil = filt.bgc_soilc
    nfin = 0; ntot = 0
    for c in 1:size(dk,1), j in 1:nlev, k in 1:size(dk,3)
        msoil[c] || continue
        ntot += 1; isfinite(dk[c,j,k]) && (nfin += 1)
    end
    finite = (ntot > 0 && nfin == ntot)
    global ok &= finite
    println("  ", rpad("active decomp_cpools_vr finite", 34), " = ", "$nfin/$ntot",
            finite ? "  PASS" : "  FAIL")
end

println()
if ok
    println("LIVE MIMICS INIT ACTIVATION OK")
    exit(0)
else
    println("LIVE MIMICS INIT ACTIVATION FAILED")
    exit(1)
end
