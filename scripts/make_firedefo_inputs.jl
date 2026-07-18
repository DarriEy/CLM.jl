# ==========================================================================
# Synthetic inputs for the TROPICAL DEFORESTATION-FIRE (`dtrotr`) parity case.
#
# WHY THIS EXISTS
#   The Li2016 deforestation-fire branch (`dtrotr_col` / `lfc` / `lfc2` /
#   `fbac1`) is `F≡0` vacuous at EVERY migrated SYMFLUENCE domain, for two
#   independent reasons:
#     (1) `dtrotr_col` accumulates `-dwt_smoothed(p)` ONLY when
#         `run_has_transient_landcover()` is true, which requires a
#         `flanduse_timeseries` file + `do_transient_pfts`. No migrated domain
#         has one, so `dwt ≡ 0` ⇒ `dtrotr ≡ 0`.
#     (2) The branch is additionally gated on `trotr1+trotr2 > 0.6`, i.e. the
#         column must be >60% tropical broadleaf tree (natpft index 4 = BET
#         tropical, index 6 = BDT tropical). Bow is boreal: natpft = [bare 5,
#         NET temperate 60, C3 arctic grass 35] ⇒ trotr1 = trotr2 = 0.
#
# THE UNBLOCK (same technique as the synthetic lnfm/hdm streams and the
# synthetic `peatf=0.5` surfdata of docs/CH4_FIRE_PARITY.md §2/§12)
#   Inject a SYNTHETIC SHARED input that BOTH codes read. This script writes
#   two files into the reference case directory:
#
#     surfdata_defo.nc            — a COPY of Bow's surfdata whose PCT_NAT_PFT
#                                   is re-composed to 40% BET-tropical + 30%
#                                   BDT-tropical (⇒ trotr1+trotr2 = 0.70 > 0.6)
#                                   with a small nonzero weight on EVERY other
#                                   natpft. Nothing else in the surfdata is
#                                   touched.
#     landuse_timeseries_defo.nc  — a minimal 2-slice `flanduse_timeseries`
#                                   (YEAR = 2202, 2203). Slice 1 is BYTE-EQUAL
#                                   to the surfdata composition (so CTSM's
#                                   `dynpft_check_consistency` passes without
#                                   having to disable it); slice 2 DECREASES
#                                   BET 40 → 38.5 and BDT 30 → 28.9 (total
#                                   −2.6 pp, redistributed to C3 arctic grass)
#                                   ⇒ `dwt_smoothed < 0` for both tropical
#                                   patches ⇒ `dtrotr_col > 0`.
#
# WHY EVERY natpft GETS A NONZERO WEIGHT
#   `subgridMod::natveg_patch_exists` creates patches for ALL 15 natural PFTs
#   when `do_transient_pfts` is on, but only for nonzero-weight PFTs otherwise.
#   Giving all 15 a nonzero weight makes the Fortran (transient) and CLM.jl
#   (non-transient) patch vectors the SAME length and the SAME order, so the
#   `pdump`/`bgcdump` patch arrays index-align with CLM.jl's — no remapping.
#
# WHY THE −2.6 pp MAGNITUDE
#   The deforestation climate term contains
#       max(0.0005, min(1, 19*dtrotr*dayspyr*secspday/dt - 0.001))
#   which SATURATES at 1 for any annual tropical loss above ~0.53 pp. A 2.6 pp
#   annual loss lands the term at 19*0.026 - 0.001 = 0.493 — in the INTERIOR of
#   both clamps, so the diff actually probes the formula rather than a clamp.
#   (Same reasoning keeps `cri` interior: with BET/BDT thresholds 4.0/1.8 the
#   area-weighted cri = 3.06 mm/day, well above Bow's January precipitation.)
#
# Usage: julia +1.12 --project=. scripts/make_firedefo_inputs.jl [outdir]
# ==========================================================================

using NCDatasets

const SRC_SURFDAT = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/" *
                    "optimization/CLM/dds_run_1/final_evaluation/settings/CLM/parameters/surfdata_clm.nc"
const DEFAULT_OUT = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bow_ref_firedefo"

# natpft index (0-based on file) → CLM PFT number. 4 = broadleaf evergreen
# tropical tree (nbrdlf_evr_trp_tree), 6 = broadleaf deciduous tropical tree
# (nbrdlf_dcd_trp_tree), 12 = C3 arctic grass.
const I_BET  = 4
const I_BDT  = 6
const I_GRAS = 12

"""PCT_NAT_PFT composition at the two flanduse time slices (1-based Julia idx = natpft+1)."""
function compositions()
    y1 = fill(2.0, 15)          # every natpft nonzero -> patch vectors align
    y1[1]        = 6.0          # bare soil
    y1[I_BET+1]  = 40.0
    y1[I_BDT+1]  = 30.0
    @assert sum(y1) ≈ 100.0 "year-1 PCT_NAT_PFT must sum to 100, got $(sum(y1))"

    y2 = copy(y1)
    y2[I_BET+1]  = 38.5         # −1.5 pp
    y2[I_BDT+1]  = 28.9         # −1.1 pp
    y2[I_GRAS+1] = y1[I_GRAS+1] + 2.6
    @assert sum(y2) ≈ 100.0 "year-2 PCT_NAT_PFT must sum to 100, got $(sum(y2))"
    @assert y2[I_BET+1] < y1[I_BET+1] && y2[I_BDT+1] < y1[I_BDT+1]
    return y1, y2
end

function write_surfdata(outdir, y1)
    dst = joinpath(outdir, "surfdata_defo.nc")
    cp(SRC_SURFDAT, dst; force = true)
    chmod(dst, 0o644)
    NCDataset(dst, "a") do ds
        v = ds["PCT_NAT_PFT"]
        # file layout is PCT_NAT_PFT(lsmlat, lsmlon, natpft); NCDatasets gives
        # the reversed (Fortran) order (natpft, lsmlon, lsmlat).
        old = Array(v)
        @assert size(old, 1) == 15 "expected natpft fastest-varying, got $(size(old))"
        for k in 1:15
            old[k, :, :] .= y1[k]
        end
        v[:, :, :] = old
    end
    return dst
end

function write_flanduse(outdir, y1, y2; years = (2202, 2203))
    dst = joinpath(outdir, "landuse_timeseries_defo.nc")
    isfile(dst) && rm(dst)
    NCDataset(dst, "c") do ds
        defDim(ds, "lsmlon", 1)
        defDim(ds, "lsmlat", 1)
        defDim(ds, "natpft", 15)
        defDim(ds, "time", 2)

        yv = defVar(ds, "YEAR", Int32, ("time",))
        yv.attrib["long_name"] = "year"
        yv[:] = Int32[years[1], years[2]]

        tv = defVar(ds, "time", Int32, ("time",))
        tv.attrib["long_name"] = "year"
        tv.attrib["units"] = "common_year"
        tv[:] = Int32[years[1], years[2]]

        pv = defVar(ds, "PCT_NAT_PFT", Float64, ("natpft", "lsmlon", "lsmlat", "time"))
        pv.attrib["long_name"] = "percent plant functional type on the natural veg landunit"
        pv.attrib["units"] = "unitless"
        pv[:, 1, 1, 1] = y1
        pv[:, 1, 1, 2] = y2

        ds.attrib["title"]  = "SYNTHETIC flanduse_timeseries for the CLM.jl dtrotr (deforestation-fire) parity case"
        ds.attrib["source"] = "CLM.jl scripts/make_firedefo_inputs.jl — NOT a physical land-use dataset"
    end
    return dst
end

function main()
    outdir = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_OUT
    mkpath(outdir)
    y1, y2 = compositions()

    trotr1, trotr2 = y1[I_BET+1] / 100, y1[I_BDT+1] / 100
    dtrotr_annual  = (y1[I_BET+1] - y2[I_BET+1] + y1[I_BDT+1] - y2[I_BDT+1]) / 100

    s = write_surfdata(outdir, y1)
    f = write_flanduse(outdir, y1, y2)

    println("wrote $s")
    println("wrote $f")
    println()
    println("  trotr1 + trotr2 = $(trotr1 + trotr2)  (gate: > 0.6)  -> ",
            trotr1 + trotr2 > 0.6 ? "OPEN" : "CLOSED")
    println("  annual tropical loss = $dtrotr_annual  (dtrotr accumulates to this over the year)")
    println("  19*dtrotr_annual - 0.001 = $(19 * dtrotr_annual - 0.001)  (must be in (0.0005, 1) to be interior)")
    cri = (4.0 * trotr1 + 1.8 * trotr2) / (trotr1 + trotr2)
    println("  cri = $cri mm/day  (area-weighted defo_fire_precip_thresh_bet/bdt)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
