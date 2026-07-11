#!/usr/bin/env julia
# ==========================================================================
# fates_singlepft.jl — single-PFT (tropical broadleaf-evergreen) FATES stand.
#
# Tests the "boom-bust is a STAND-COMPOSITION artifact, not a bug" hypothesis.
# The FATES natural-bare-ground cold start (EDInitMod.init_cohorts!) seeds ONE
# cohort per PFT for ALL ~14 FATES default PFTs — including COLD-DECIDUOUS ones
# (prt_params.season_decid[pft]==1) that don't belong at a tropical Amazon site.
# Over demographic time those cold-deciduous cohorts grow, dominate the leaf
# area, then FATES correctly forces a stand-clearing leaf-off (phenology
# nevercold reset) and they starve — the boom->bust die-back.
#
# This harness reuses fates_longhorizon.jl's build()/main()/census helpers, then
# RESTRICTS the stand to a single tropical broadleaf-evergreen PFT (default
# FATES PFT 1) two ways:
#   1. site.use_this_pft[pft] = (pft==target) for every site  — gates ALL
#      run-time recruitment (EDPhysiologyMod.recruitment, line ~1356) AND the
#      seed-rain distribution (line ~1174) to the target PFT only, so no
#      cold-deciduous cohort is ever re-recruited.
#   2. Prunes the already-seeded non-target cold-start cohorts via the supported
#      terminate_cohort path (biomass -> litter, list relinked), then
#      ed_update_site to rebuild canopy structure and TotalBalanceCheck(final)
#      to reseat the site mass-balance old_stock on the pruned stand.
# It does NOT touch core src/ — only the live in-memory FATES state after
# clm_initialize!.
#
#   FATES_NDAYS=3650 julia +1.12 --project=. scripts/fates_singlepft.jl
#   FATES_PFT=1      (target PFT index; default 1 = tropical broadleaf evergreen)
# ==========================================================================
include(joinpath(@__DIR__, "fates_longhorizon.jl"))

# ---- restrict the live FATES state to a single evergreen PFT --------------
function restrict_to_single_pft!(inst, fates, bounds)
    npft   = _C.numpft[]
    target = parse(Int, get(ENV, "FATES_PFT", "1"))

    # Show the deciduous taxonomy of the default PFT set so the choice is auditable.
    println("  PFT taxonomy (season_decid / stress_decid / evergreen):")
    for p in 1:npft
        @printf("    pft %2d : season_decid=%d  stress_decid=%d  evergreen=%d\n",
                p, _C.prt_params.season_decid[p], _C.prt_params.stress_decid[p],
                _C.prt_params.evergreen[p])
    end
    if _C.prt_params.season_decid[target] != _C.ifalse
        error("FATES_PFT=$target is season-deciduous (season_decid=1); pick an evergreen PFT")
    end
    @printf("  restricting stand to PFT %d (season_decid=%d, evergreen=%d)\n\n",
            target, _C.prt_params.season_decid[target], _C.prt_params.evergreen[target])

    for s in 1:fates.nsites
        site  = fates.sites[s]
        bc_in = fates.bc_in[s]

        # (1) Gate run-time recruitment + seed rain to the target PFT only.
        for p in 1:npft
            site.use_this_pft[p] = (p == target ? _C.itrue : _C.ifalse)
        end

        # (2) Prune the already-seeded non-target cold-start cohorts.
        pruned = 0
        cp = site.oldest_patch
        while cp !== nothing
            cc = cp.tallest
            while cc !== nothing
                nxt = cc.shorter                    # save link before unlink
                if cc.pft != target
                    _C.terminate_cohort(site, cp, cc, bc_in, _C.i_term_mort_type_numdens)
                    pruned += 1
                end
                cc = nxt
            end
            cp = cp.younger
        end
        # Rebuild canopy layering / patch structure on the pruned stand and
        # reseat the mass-balance baseline so the daily check starts clean.
        _C.ed_update_site(site, bc_in, fates.bc_out[s], false)
        _C.TotalBalanceCheck(site, _C._edmain_final_check_id)
        @printf("  site %d: pruned %d non-target cohort(s)\n", s, pruned)
    end
    return target
end

# ---- override build() so main() (from fates_longhorizon.jl) uses it -------
# Redefining `build` here REPLACES the all-PFT builder included from
# fates_longhorizon.jl (same generic function, identical zero-arg signature);
# the included main() resolves `build` by name at call time and so picks this
# up. We reproduce the all-PFT cold-start body verbatim (a wrapper cannot call
# the shadowed original — it is the same function object) and append the
# single-PFT restriction before returning.
function build()
    fsurdat = get(ENV, "FATES_FSURDAT",
        "$DATA/domain_Aripuana_Amazon/settings/CLM/parameters/surfdata_clm.nc")
    paramfile = get(ENV, "FATES_PARAMFILE",
        "$DATA/domain_Aripuana_Amazon/settings/CLM/parameters/clm5_params.nc")
    fyr = parse(Int, get(ENV, "FATES_YEAR", "2004"))
    inst, bounds, filt, _tm = _C.clm_initialize!(; fsurdat=fsurdat, paramfile=paramfile,
        use_fates=true, start_date=DateTime(fyr,1,1), dtime=1800)
    for g in 1:bounds.endg
        isfinite(inst.gridcell.latdeg[g]) || (inst.gridcell.latdeg[g] = inst.gridcell.lat[g]*180/π)
        isfinite(inst.gridcell.londeg[g]) || (inst.gridcell.londeg[g] = inst.gridcell.lon[g]*180/π)
    end
    let col=inst.column, lun=inst.landunit, nflag=0
        for c in 1:bounds.endc
            if lun.itype[col.landunit[c]] == _C.ISTSOIL
                col.is_fates[c] = true; nflag += 1
            end
        end
        @printf("  flagged %d FATES column(s); fates.nsites=%d\n", nflag, inst.fates.nsites)
    end
    if get(ENV,"FATES_WARMSOIL","1")=="1"
        col=inst.column; temp=inst.temperature; joff=_C.varpar.nlevsno; ngr=_C.varpar.nlevgrnd
        for c in 1:bounds.endc
            col.is_fates[c] || continue
            temp.t_grnd_col[c]=299.0
            for j in 1:ngr; temp.t_soisno_col[c, joff+j]=299.0; end
        end
    end
    config = _C.CLMDriverConfig(use_fates=true)
    filt_ia = _C.clump_filter_inactive_and_active

    fates = inst.fates
    target = restrict_to_single_pft!(inst, fates, bounds)

    # Cold-start probe: confirm ONLY the evergreen PFT now has cohorts.
    site = fates.sites[1]; bad = 0; pfts = Int[]
    cp = site.oldest_patch
    while cp !== nothing
        cc = cp.tallest
        while cc !== nothing
            push!(pfts, cc.pft)
            (_C.prt_params.season_decid[cc.pft] == _C.ifalse) || (bad += 1)
            cc = cc.shorter
        end
        cp = cp.younger
    end
    @printf("  cold-start cohorts after restriction: pfts=%s  (season_decid!=0 count=%d)\n",
            string(sort(unique(pfts))), bad)
    bad == 0 || error("restriction failed: $bad season-deciduous cohort(s) remain")
    all(==(target), pfts) || error("restriction failed: non-target pfts present: $(unique(pfts))")
    println("  ✔ single-evergreen-PFT stand confirmed\n")
    return inst, fates, config, bounds, filt, filt_ia
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
