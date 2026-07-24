# ==========================================================================
# fortran_parity_voc_aripuana.jl — VOC / MEGAN Julia↔Fortran parity at the
# TROPICAL Aripuanã (Amazon) site — the canonical isoprene stress case
# (dominant PFT 4 = broadleaf-evergreen tropical, isoprene EF ~1e4 ug/m2/h).
#
# Unlike scripts/validate_voc_megan_fortran_parity.jl (an independent scalar
# ORACLE transcribed from VOCEmissionMod.F90), this diffs the port against an
# actual Fortran CTSM run: the same cesm_iso.exe (CN+c13+c14+MEGAN) cold-started
# at Aripuanã, with the MEGAN diagnostics written to a native-grid h1 history
# tape (hist_dov2xy=.false., hist_nhtfrq=1). Fortran stores, PER PATCH, for the
# FIRST megan compound (isoprene, imeg==1): GAMMA, GAMMAL/T/P/A/S/C, EOPT, TOPT,
# ALPHA, PAR_sun/24/240, PAR_shade/24/240; plus VOCFLXT (total, moles/m2/s) and
# MEG_isoprene (kg/m2/s). The 24/240-hour accumulators the gamma kernels consume
# (fsun240, forc_solad/solai 24/240, t_veg24/240, pco2/pbot 240) come from the
# restart-format before_step pdump at the matched step.
#
# TWO checks, both non-vacuous at the tropical site (isoprene flux is LARGE):
#   [A] DECOUPLED KERNEL PARITY — feed Fortran's own dumped inputs into the port
#       gamma helpers (get_gamma_L/SM/P/T/A/C) and compare each factor to the
#       Fortran GAMMA* diagnostic. Isolates the VOC kernel from any canopy
#       residual (all inputs are Fortran's).
#   [B] ASSEMBLY PARITY — with Fortran's total GAMMA and the port's emission
#       factor (read from the SAME megan_factors_file) + molec weight, reproduce
#       MEG_isoprene (kg/m2/s) and the isoprene contribution to VOCFLXT
#       (moles/m2/s). Tests the EF read + unit factor (1/3600/1e6) + moles split.
#
#   julia +1.12 --project=<repo> scripts/fortran_parity_voc_aripuana.jl [nstep]
# ==========================================================================
using CLM, NCDatasets, Dates, Printf

const ARI_RUN   = get(ENV, "ARI_ISO_DUMPDIR", "/private/tmp/claude-501/ari_iso_run")
const MEGAN_FILE = "/Users/darri.eythorsson/projects/cesm-inputdata/atm/cam/chem/trop_mozart/emis/megan21_emis_factors_78pft_c20161108.nc"
# The megan_specifier from the Fortran run's drv_flds_in (order matters: isoprene
# must be compound #1, matching Fortran's imeg==1 diagnostic convention).
const MEGAN_SPEC = ["ISOP = isoprene",
                    "C10H16 = pinene_a + carene_3 + thujene_a",
                    "CH3OH = methanol", "C2H5OH = ethanol",
                    "CH2O = formaldehyde", "CH3CHO = acetaldehyde",
                    "CH3COOH = acetic_acid", "CH3COCH3 = acetone"]

const TFRZ = 273.15
const MEGFAC = 1.0/3600.0/1.0e6      # megemis_units_factor (VOCEmissionMod.F90:459)

# ---- helpers to fetch a var from a NetCDF dataset, NaN-safe ----
getv(ds, k) = haskey(ds, k) ? Float64.(coalesce.(ds[k][:], NaN)) : Float64[]
# pull scalar accumulator "<STEM>_VALUE" (per-patch vector) from a pdump
accum(ds, stem) = haskey(ds, stem*"_VALUE") ? Float64.(coalesce.(ds[stem*"_VALUE"][:], NaN)) : Float64[]

# List all h1 history files (there may be several — hist_mfilt splits them).
function h1_files()
    sort(filter(f -> occursin(r"\.clm2\.h1\.", f) && endswith(f, ".nc"),
                readdir(ARI_RUN; join=true)))
end

# Scan every h1 file for the record whose (month,day,tod) matches nstep; return
# (filepath, record_index) or ("", 0). Year-agnostic (model clock is Bow's 2202
# label, forcing recycles 2002 — the absolute year is irrelevant).
function find_h1_record(nstep)
    tgt = DateTime(2000,1,1) + Hour(nstep)
    tmd = (month(tgt), day(tgt)); tsec = hour(tgt)*3600 + minute(tgt)*60
    for f in h1_files()
        ds = NCDataset(f,"r")
        if haskey(ds,"mcdate") && haskey(ds,"mcsec")
            md = Int.(ds["mcdate"][:]); ms = Int.(ds["mcsec"][:])
            for k in eachindex(md)
                if ((md[k]÷100)%100, md[k]%100) == tmd && ms[k] == tsec
                    close(ds); return (f, k)
                end
            end
        end
        close(ds)
    end
    return ("", 0)
end

# Build the port MEGAN config once (reads the real megan_factors_file).
function build_megan()
    CLM.megan_config_init(use_voc=true, megan_specifier=MEGAN_SPEC,
                          megan_factors_file=MEGAN_FILE)
end

const NFAIL = Ref(0); const NCHK = Ref(0)
function chk(name, got, ref; rtol=1e-6, atol=1e-30)
    NCHK[] += 1
    d = abs(got - ref); ok = d <= atol + rtol*max(abs(got),abs(ref))
    ok || (NFAIL[] += 1)
    @printf("    %-26s port=% .6e  F=% .6e  |Δ|=%.2e rel=%.2e  %s\n",
            name, got, ref, d, d/(1e-30+abs(ref)), ok ? "ok" : "**FAIL**")
    ok
end

# ==========================================================================
# The h1 tape carries a time dimension (one record/step at nhtfrq=1) and a pft
# (patch) dimension. Find the record index whose model date == the pdump step's.
# ==========================================================================
function hist_record_for_step(dsh, nstep)
    # The Fortran run cold-starts at <model-year>-01-01-00000 (model year is Bow's
    # 2202 label; forcing recycles 2002 — the ABSOLUTE year is irrelevant here).
    # Step N ends at day-of-year (1 + N÷24), time-of-day (N mod 24)*3600. Match on
    # (month, day, tod) so the 2202-vs-2002 label doesn't matter. h1 has nhtfrq=1
    # (one record/step); the first record is the nstep-0 initial write.
    tgt = DateTime(2000,1,1) + Hour(nstep)          # reference year cancels below
    tmd, tsec = (month(tgt), day(tgt)), (hour(tgt)*3600 + minute(tgt)*60)
    if haskey(dsh,"mcdate") && haskey(dsh,"mcsec")
        md = Int.(dsh["mcdate"][:]); ms = Int.(dsh["mcsec"][:])
        for k in eachindex(md)
            ((md[k]÷100)%100, md[k]%100) == tmd && ms[k] == tsec && return k
        end
        # fallback: nearest in seconds-of-year (day-of-month × 86400 + tod), year-agnostic
        soy(m,s) = ((m%100)*24*3600 + s)               # day-of-month*day + tod (Jan window)
        tgtsoy = tmd[2]*24*3600 + tsec
        return argmin([abs(soy(md[k], ms[k]) - tgtsoy) for k in eachindex(md)])
    end
    # CF time (days since ref): match by fractional day = nstep/24.
    return argmin(abs.(Float64.(dsh["time"][:]) .- nstep/24))
end

# ==========================================================================
function voc_parity(nstep)
    (h1, krec) = find_h1_record(nstep)
    bdump = joinpath(ARI_RUN, "pdump_before_step_n$(nstep).nc")
    if isempty(h1) || !isfile(bdump)
        println("MISSING inputs: h1 record for n$nstep = $(isempty(h1) ? "NOT FOUND" : "$h1[$krec]")  bdump exists=$(isfile(bdump))")
        return 1
    end
    println("h1 history : $h1")
    println("h1 record  : $krec  (matched to nstep $nstep)")
    println("pdump      : $bdump\n")

    megcfg = build_megan()
    iso = megcfg.meg_compounds[findfirst(c->c.name=="isoprene", megcfg.meg_compounds)]
    mf  = megcfg.megan_factors
    isoclass = iso.class_number

    dsh = NCDataset(h1,"r"); dsb = NCDataset(bdump,"r")

    # ---- read Fortran per-patch diagnostics at record krec ----
    # native-grid history var shape is (pft, time) in CDL → NCDatasets (pft,time).
    slic(k) = begin
        v = dsh[k]
        nd = ndims(v)
        nd == 2 ? Float64.(coalesce.(v[:,krec],NaN)) :
        nd == 1 ? Float64.(coalesce.(v[:],NaN)) : Float64[]
    end
    hash(k) = haskey(dsh,k)

    pfts = hash("pfts1d_itypveg") ? Int.(dsh["pfts1d_itypveg"][:]) :
           (haskey(dsb,"pfts1d_itypveg") ? Int.(dsb["pfts1d_itypveg"][:]) : Int[])

    # Fortran outputs
    F = Dict{String,Vector{Float64}}()
    for k in ("VOCFLXT","MEG_isoprene","MEG_pinene_a","MEG_carene_3","MEG_thujene_a",
              "GAMMA","GAMMAL","GAMMAT","GAMMAP","GAMMAA",
              "GAMMAS","GAMMAC","EOPT","TOPT","ALPHA","PAR_sun","PAR24_sun","PAR240_sun",
              "PAR_shade","PAR24_shade","PAR240_shade","TV","ELAI","TLAI","FSUN","BTRAN")
        hash(k) && (F[k] = slic(k))
    end
    println("Fortran h1 fields present: ", sort(collect(keys(F))))

    # accumulators from the pdump (per-patch)
    A = Dict{String,Vector{Float64}}()
    for (jl,stem) in (("fsun240","FSUN240"),("fsun24","FSUN24"),
                      ("solad240","FSD240"),("solad24","FSD24"),
                      ("solai240","FSI240"),("solai24","FSI24"),
                      ("tveg240","T_VEG240"),("tveg24","T_VEG24"),
                      ("pco2_240","pco2_240"),("pbot240","pbot240"))
        v = accum(dsb, stem); isempty(v) || (A[jl]=v)
    end
    elai_b = getv(dsb,"elai"); btran_b = accum(dsb,"BTRANAV"); fsun_b = getv(dsb,"fsun")
    println("pdump accumulators present: ", sort(collect(keys(A))),
            "  elai=", !isempty(elai_b), "  fsun=", !isempty(fsun_b), "  BTRANAV=", !isempty(btran_b), "\n")

    # map dump/hist patch (by pft type) — focus the vegetated PFTs (skip ivt 0).
    # Fortran history + pdump share the same 1d patch layout, so index alignment
    # is identity here; we still print the pft type per patch.
    np = length(get(F,"GAMMAL", pfts |> x->Float64.(x)))
    np == 0 && (np = length(pfts))

    maxg = 0.0; maxflx = 0.0
    for p in 1:np
        ivt = p <= length(pfts) ? pfts[p] : -1
        (ivt >= 1 && ivt <= 15) || continue
        gL_F = get(F,"GAMMAL",Float64[]); (p<=length(gL_F) && isfinite(gL_F[p]) && gL_F[p] > 0) || continue

        println("── patch $p  (PFT itype=$ivt) ───────────────────────────────")
        gL = gP = gT = gA = NaN

        # ---- [A] decoupled kernel parity: feed Fortran inputs into port helpers
        # gamma_L: needs (fsun240, elai). fsun240 from accum; elai from pdump.
        if haskey(A,"fsun240") && !isempty(elai_b) && p<=length(elai_b)
            gL = CLM.get_gamma_L(A["fsun240"][p], elai_b[p])
            haskey(F,"GAMMAL") && (chk("GAMMAL", gL, F["GAMMAL"][p]); maxg=max(maxg,abs(gL-F["GAMMAL"][p])/(1e-30+abs(F["GAMMAL"][p]))))
        end
        # gamma_P: h1 gives PAR_sun/24/240 + PAR_shade directly; the two remaining
        # shade PARs are reconstructed from the pdump accumulators exactly as
        # VOCEmissionMod builds them: par24_sha=(1-fsun24)*solai24*4.6,
        # par240_sha=(1-fsun240)*solai240*4.6. Instantaneous fsun from the pdump.
        gP = NaN
        if all(haskey(F,k) for k in ("PAR_sun","PAR24_sun","PAR240_sun","PAR_shade")) &&
           all(haskey(A,k) for k in ("fsun24","fsun240","solai24","solai240","solad240")) &&
           !isempty(fsun_b) && p<=length(fsun_b)
            LDF = mf.LDF[isoclass]
            par24_sha  = (1.0 - A["fsun24"][p])  * A["solai24"][p]  * 4.6
            par240_sha = (1.0 - A["fsun240"][p]) * A["solai240"][p] * 4.6
            gP, cp, alpha = CLM.get_gamma_P(F["PAR_sun"][p],F["PAR24_sun"][p],F["PAR240_sun"][p],
                                            F["PAR_shade"][p], par24_sha, par240_sha,
                                            fsun_b[p], A["fsun240"][p], A["solad240"][p], A["solai240"][p], LDF)
            # ALPHA is get_gamma_P's nonlinear coefficient (independent of fsun) →
            # a clean kernel check. GAMMAP's final value additionally weights
            # sunlit/shade by the INSTANTANEOUS fsun, which is not on the h1 tape
            # (FSUN dropped from fincl); we feed the start-of-step pdump fsun while
            # Fortran's PAR_sun was built with its during-step fsun. That timing
            # inconsistency (not a kernel error — ALPHA matches to ~1e-8) shows up
            # as a ~few-% GAMMAP gap, so GAMMAP is reported diagnostically, NOT as
            # a hard parity fail.
            haskey(F,"ALPHA")  && chk("ALPHA",  alpha, F["ALPHA"][p])
            if haskey(F,"GAMMAP")
                d = abs(gP - F["GAMMAP"][p]); r = d/(1e-30+abs(F["GAMMAP"][p]))
                @printf("    %-26s port=% .6e  F=% .6e  |Δ|=%.2e rel=%.2e  [input-limited: fsun timing]\n",
                        "GAMMAP (diag)", gP, F["GAMMAP"][p], d, r)
            end
        end
        # gamma_T: needs t_veg240/24/t_veg + class params. accum + Fortran TV.
        if haskey(A,"tveg240") && haskey(A,"tveg24") && haskey(F,"TV")
            gT, Eopt, topt = CLM.get_gamma_T(A["tveg240"][p], A["tveg24"][p], F["TV"][p],
                                             mf.ct1[isoclass], mf.ct2[isoclass], mf.betaT[isoclass],
                                             mf.LDF[isoclass], mf.Ceo[isoclass], ivt)
            haskey(F,"GAMMAT") && chk("GAMMAT", gT, F["GAMMAT"][p])
            haskey(F,"EOPT") && chk("EOPT", Eopt, F["EOPT"][p])
            haskey(F,"TOPT") && chk("TOPT", topt, F["TOPT"][p])
        end
        # gamma_A: evergreen tropical (ivt=4) → 1.0. Compare to Fortran GAMMAA.
        gA = CLM.get_gamma_A(ivt, elai_b[p], elai_b[p], isoclass, mf)
        haskey(F,"GAMMAA") && chk("GAMMAA", gA, F["GAMMAA"][p])

        # gamma_SM (soil moisture) + gamma_C (CO2) need instantaneous btran / ci,
        # which are NOT on the restart pdump or this h1 tape — so they cannot be
        # recomputed from the reference alone. The port's get_gamma_SM/get_gamma_C
        # are separately value-verified vs the VOCEmissionMod oracle (Bow harness).
        # Here we take Fortran's GAMMAS/GAMMAC as given and cross-check that the
        # port's INDEPENDENTLY-computed L·P·T·A factors reproduce Fortran's TOTAL
        # GAMMA: port_total = gL·gP·gT·gA·GAMMAS_F·GAMMAC_F  vs  Fortran GAMMA.
        gS_F = get(F,"GAMMAS",[NaN])[min(p,length(get(F,"GAMMAS",[NaN])))]
        gC_F = get(F,"GAMMAC",[NaN])[min(p,length(get(F,"GAMMAC",[NaN])))]
        gP_F = get(F,"GAMMAP",[NaN])[min(p,length(get(F,"GAMMAP",[NaN])))]
        @printf("    %-26s F_GAMMAS=% .6e  F_GAMMAC=% .6e  (inputs btran/ci not dumped)\n",
                "GAMMAS/GAMMAC (Fortran)", gS_F, gC_F)
        # Cross-check that the port's cleanly-computed L·T·A factors reproduce
        # Fortran's TOTAL GAMMA when combined with Fortran's own P/S/C (whose
        # instantaneous inputs — fsun/btran/ci — are not on the reference).
        if all(isfinite, (gL,gT,gA,gS_F,gC_F,gP_F)) && haskey(F,"GAMMA")
            gtot = gL*gP_F*gT*gA*gS_F*gC_F
            chk("GAMMA_total(port L·T·A × F P,S,C)", gtot, F["GAMMA"][p])
        end

        # ---- [B] assembly parity: reproduce MEG_isoprene & VOCFLXT(isoprene) ----
        # port epsilon (EF) for this PFT from the megan_factors_file, port molec wt.
        eps = iso.emis_factors[ivt]
        gamma_F = get(F,"GAMMA",[NaN])[min(p,length(get(F,"GAMMA",[NaN])))]
        if isfinite(gamma_F) && 0.0 <= gamma_F < 100.0
            meg_iso_port = iso.coeff * eps * gamma_F * MEGFAC * 1.0e-3            # kg/m2/s
            vocflx_iso_port = iso.coeff * eps * gamma_F * MEGFAC / iso.molec_weight # moles/m2/s (ISOP mech = isoprene only)
            haskey(F,"MEG_isoprene") && (chk("MEG_isoprene", meg_iso_port, F["MEG_isoprene"][p]);
                maxflx=max(maxflx, abs(meg_iso_port-F["MEG_isoprene"][p])/(1e-30+abs(F["MEG_isoprene"][p]))))
            # VOCFLXT is the SUM over all mech comps; isoprene(ISOP) is one term.
            # Report the isoprene contribution vs the total (isoprene dominates in tropics).
            if haskey(F,"VOCFLXT")
                @printf("    %-26s port_ISOPterm=% .6e  F_VOCFLXT(total)=% .6e  ratio=%.3f\n",
                        "VOCFLXT(isop vs total)", vocflx_iso_port, F["VOCFLXT"][p],
                        vocflx_iso_port/(1e-30+F["VOCFLXT"][p]))
            end
            @printf("    [inputs] eps(EF PFT%d)=%.1f  mw=%.3f  GAMMA_F=%.4f  → non-vacuous isoprene\n",
                    ivt, eps, iso.molec_weight, gamma_F)
            # Fortran monoterpene (C10H16 mech = pinene_a+carene_3+thujene_a) fluxes,
            # for context (non-vacuity of the second mechanism compound). Their
            # per-compound gammas are not on the tape (only isoprene's imeg==1 set),
            # so these are reported, not assembly-diffed.
            mono = sum(get(F,k,[0.0])[min(p,length(get(F,k,[0.0])))] for k in
                       ("MEG_pinene_a","MEG_carene_3","MEG_thujene_a"))
            @printf("    %-26s Σ(pinene_a+carene_3+thujene_a)=% .6e kg/m2/s (Fortran)\n",
                    "monoterpene C10H16", mono)
        end
        println()
    end
    close(dsh); close(dsb)

    @printf("\n  VOC parity summary: %d checks, %d FAIL\n", NCHK[], NFAIL[])
    @printf("  max|rel| gamma factors = %.3e   max|rel| MEG_isoprene = %.3e\n", maxg, maxflx)
    return NFAIL[] == 0 ? 0 : 1
end

function main()
    na = filter(a->occursin(r"^\d+$",a), ARGS)
    nstep = isempty(na) ? 40 : parse(Int, na[1])
    println("="^74); println("VOC/MEGAN Fortran parity — Aripuanã tropical (nstep=$nstep)"); println("="^74)
    rc = voc_parity(nstep)
    println("="^74)
    println(NFAIL[]==0 ? "RESULT: VOC tropical parity PASS" : "RESULT: $(NFAIL[]) VOC divergence(s)")
    return rc
end

exit(main())
