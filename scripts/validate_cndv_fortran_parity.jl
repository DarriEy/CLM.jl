#!/usr/bin/env julia
# =============================================================================
# validate_cndv_fortran_parity.jl
#
# NON-VACUOUS Fortran-parity harness for CNDV (dynamic global vegetation).
#
# WHY THIS EXISTS
# ---------------
# The existing test_cndv.jl / test_cndv_wiring.jl assert that the CNDV call
# sequence "moves" patch weights and stays finite/in-[0,1]. That is structurally
# BLIND to two failure modes this project has actually shipped before:
#   1. A subsystem that runs but computes the WRONG VALUE (a mis-transcribed
#      Fortran formula / a wrong index offset) — a "moved and finite" test passes.
#   2. The CNDV establishment path being effectively DEAD because prec365 was a
#      hardcoded ZERO buffer, so `prec365 >= prec_min_estab` was false forever and
#      no PFT could EVER establish (MEMORY: fortran-ground-truth-unblocked). The
#      wiring test feeds prec365 = zeros, so it can NEVER exercise establishment.
#
# WHAT THIS HARNESS DOES (each is a VALUE assertion, not a finiteness assertion)
#   A. LIVE prec365 feed: drive the real atm2lnd_update_acc_vars! with use_cndv=true
#      and NONZERO precip forcing; assert prec365_col becomes NONZERO and equals the
#      analytic running-mean value. Directly guards the historical zero-buffer bug.
#   B. Light competition value-parity: run the port's cndv_light! and an INDEPENDENT
#      plain-scalar re-implementation of CNDVLightMod.F90::Light; assert bit-close.
#   C. Establishment value-parity + a REAL establishment event: run cndv_establishment!
#      and an independent re-implementation of CNDVEstablishmentMod.F90::Establishment
#      on a scenario where a bare (present=false) climate-adapted PFT MUST become
#      present with fpcgrid>0 given NONZERO prec365. Assert bit-close AND that the
#      establishment actually fired (the case the zero-prec365 test can't reach).
#
# The oracles are transcribed independently from the Fortran source at
#   $CLM_FORTRAN/src/biogeochem/CNDVLightMod.F90
#   $CLM_FORTRAN/src/biogeochem/CNDVEstablishmentMod.F90
# so an index-offset / threshold / formula transcription bug in the KERNELS is
# caught by disagreement with the oracle.
#
# Run:  julia +1.12 --project=. scripts/validate_cndv_fortran_parity.jl
# Exit: 0 = all parity checks pass, 1 = a divergence was found.
# =============================================================================

using CLM
using Test
using Printf

const C = CLM
const TFRZ = C.TFRZ
const PI_ = C.RPI
const SECSPDAY = C.SECSPDAY

# -----------------------------------------------------------------------------
# Scenario builder: one gridcell, one column, one istsoil landunit.
# PFT layout (Julia 1-based patch index -> Fortran ivt (0-based) via itype):
#   p1: bare ground   itype=0  (noveg)
#   p2: tree A        itype=1  woody, is_tree     -- present, dense
#   p3: tree B        itype=2  woody, is_tree     -- present, dense (forces >95% comp)
#   p4: grass         itype=13 non-woody          -- present
#   p5: tree C        itype=3  woody, is_tree     -- ABSENT, climate-adapted -> establishes
# -----------------------------------------------------------------------------
function build_scenario()
    np, nc, ng, nl = 5, 1, 1, 1

    pft = C.PftconType()
    C.pftcon_allocate!(pft)
    trees   = [2, 3, 4]      # Julia PFT indices for itype 1,2,3
    grasses = [14]           # Julia PFT index for itype 13
    for i in trees
        pft.woody[i] = 1.0; pft.is_tree[i] = true; pft.is_shrub[i] = false
        pft.dwood[i] = 2.5e5; pft.slatop[i] = 0.012; pft.dsladlai[i] = 0.0
    end
    for i in grasses
        pft.woody[i] = 0.0; pft.is_tree[i] = false; pft.is_shrub[i] = false
        pft.dwood[i] = 0.0; pft.slatop[i] = 0.030; pft.dsladlai[i] = 0.0
    end
    # eco constants (dgv_ecophyscon reads pftpar20/28/29/30/31). Pick values that
    # ALLOW establishment for the adapted tree: tcmin<=tmomin20-TFRZ<=tcmax, agdd20>=gddmin,
    # twmax>999 so the boreal agddtw gate is bypassed.
    pft.pftpar20 .= 15.0     # crownarea_max
    pft.pftpar28 .= -35.0    # tcmin
    pft.pftpar29 .= 28.0     # tcmax
    pft.pftpar30 .= 350.0    # gddmin
    pft.pftpar31 .= 1000.0   # twmax (>999 -> establishment agddtw gate bypassed)

    inst = C.CLMInstances()
    C.patch_init!(inst.patch, np)
    inst.patch.itype    .= [C.noveg, 1, 2, 13, 3]
    inst.patch.gridcell .= 1
    inst.patch.column   .= 1
    inst.patch.landunit .= 1
    inst.patch.wtcol    .= [0.05, 0.35, 0.35, 0.20, 0.05]

    C.landunit_init!(inst.landunit, nl)
    inst.landunit.itype[1]   = C.ISTSOIL
    inst.landunit.wtgcell[1] = 1.0

    C.dgvs_init_cold!(inst.dgvs, np)
    C.dgv_ecophyscon_init!(inst.dgv_ecophyscon, pft)

    return (; inst, pft, np, nc, ng, nl)
end

# -----------------------------------------------------------------------------
# ORACLE 1 — CNDVLightMod.F90::Light, plain-scalar, single gridcell.
# Mutates fpcgrid/nind/crownarea in place. `mask` selects natvegp patches.
# -----------------------------------------------------------------------------
function oracle_light!(fpcgrid, nind, crownarea, mask, itype, pft, eco,
                       deadstemc, leafcmax)
    taper = 200.0
    fpc_tree_max = 0.95
    np = length(fpcgrid)
    fpc_tree_total = 0.0; fpc_inc_tree = 0.0
    fpc_grass_total = 0.0; fpc_shrub_total = 0.0
    numtrees = 0
    fpc_inc = zeros(np)

    for p in 1:np
        mask[p] || continue
        iv = itype[p] + 1               # Fortran ivt(p) -> Julia array index
        if pft.woody[iv] == 1.0
            if fpcgrid[p] > 0.0 && nind[p] > 0.0
                stocking = nind[p] / fpcgrid[p]
                stemdiam = (24.0 * deadstemc[p] / (PI_ * stocking * pft.dwood[iv] * taper))^(1.0/3.0)
            else
                stemdiam = 0.0
            end
            crownarea[p] = min(eco.crownarea_max[iv], eco.allom1[iv] * stemdiam^eco.reinickerp[iv])
        end
        if crownarea[p] > 0.0 && nind[p] > 0.0
            lm_ind = leafcmax[p] * fpcgrid[p] / nind[p]
            if pft.dsladlai[iv] > 0.0
                lai_ind = max(0.001, ((exp(lm_ind*pft.dsladlai[iv] + log(pft.slatop[iv])) -
                                       pft.slatop[iv]) / pft.dsladlai[iv]) / crownarea[p])
            else
                lai_ind = lm_ind * pft.slatop[iv] / crownarea[p]
            end
        else
            lai_ind = 0.0
        end
        fpc_ind = 1.0 - exp(-0.5*lai_ind)
        fpcgrid_old = fpcgrid[p]
        fpcgrid[p] = crownarea[p] * nind[p] * fpc_ind
        fpc_inc[p] = max(0.0, fpcgrid[p] - fpcgrid_old)
        if pft.woody[iv] == 1.0
            if pft.is_tree[iv]
                numtrees += 1
                fpc_tree_total += fpcgrid[p]
                fpc_inc_tree += fpc_inc[p]
            elseif pft.is_shrub[iv]
                fpc_shrub_total += fpcgrid[p]
            end
        else
            fpc_grass_total += fpcgrid[p]
        end
    end

    fpc_grass_max = 1.0 - min(fpc_tree_total, fpc_tree_max)
    fpc_shrub_max = max(0.0, fpc_grass_max - fpc_grass_total)

    for p in 1:np
        mask[p] || continue
        iv = itype[p] + 1
        if pft.woody[iv] == 1.0 && pft.is_tree[iv]
            if fpc_tree_total > fpc_tree_max
                if fpc_inc_tree > 0.0
                    excess = (fpc_tree_total - fpc_tree_max) * fpc_inc[p] / fpc_inc_tree
                else
                    excess = (fpc_tree_total - fpc_tree_max) / numtrees
                end
                if fpcgrid[p] > 0.0
                    nind_kill = nind[p] * excess / fpcgrid[p]
                    nind[p] = max(0.0, nind[p] - nind_kill)
                    fpcgrid[p] = max(0.0, fpcgrid[p] - excess)
                else
                    nind[p] = 0.0; fpcgrid[p] = 0.0
                end
            end
        elseif pft.woody[iv] == 0.0
            if fpc_grass_total > fpc_grass_max
                excess = (fpc_grass_total - fpc_grass_max) * fpcgrid[p] / fpc_grass_total
                fpcgrid[p] = max(0.0, fpcgrid[p] - excess)
            end
        elseif pft.woody[iv] == 1.0 && pft.is_shrub[iv]
            if fpc_shrub_total > fpc_shrub_max
                excess = 1.0 - fpc_shrub_max / fpc_shrub_total
                if fpcgrid[p] > 0.0
                    nind_kill = nind[p] * excess / fpcgrid[p]
                    nind[p] = max(0.0, nind[p] - nind_kill)
                    fpcgrid[p] = max(0.0, fpcgrid[p] - excess)
                else
                    nind[p] = 0.0; fpcgrid[p] = 0.0
                end
            end
        end
    end
    return nothing
end

# -----------------------------------------------------------------------------
# ORACLE 2 — CNDVEstablishmentMod.F90::Establishment, plain-scalar, single gridcell.
# Mutates present/nind/fpcgrid/crownarea/greffic/heatstress. Returns fpc_total.
# -----------------------------------------------------------------------------
function oracle_establishment!(present, pftmayexist, nind, fpcgrid, crownarea,
                               greffic, heatstress, itype, column, landunit,
                               lun_itype, pft, eco, prec365, annsum_npp,
                               annsum_litfall, deadstemc, leafcmax, agddtw,
                               agdd20, tmomin20)
    ramp_agddtw = 300.0
    nind_min = 1.0e-10
    prec_min_estab = 100.0 / (365.0 * SECSPDAY)
    estab_max = 0.24
    taper = 200.0
    noveg = C.noveg
    nc3_arctic_grass = C.nc3_arctic_grass
    istsoil = C.ISTSOIL
    np = length(fpcgrid)

    survive = falses(np); estab = falses(np); dstemc = zeros(np)
    fpc_tree_total = 0.0; npft_estab = 0; ngrass = 0
    fpc_total_new = 0.0; fpc_total = 0.0

    for p in 1:np
        if nind[p] == 0.0; present[p] = false; end
        if !present[p]; nind[p] = 0.0; fpcgrid[p] = 0.0; end
        survive[p] = false; estab[p] = false; dstemc[p] = deadstemc[p]
    end
    for p in 1:np
        iv = itype[p] + 1
        if tmomin20[p] >= eco.tcmin[iv] + TFRZ
            if tmomin20[p] <= eco.tcmax[iv] + TFRZ && agdd20[p] >= eco.gddmin[iv]
                estab[p] = true
            end
            survive[p] = true
            if !pftmayexist[p]
                survive[p] = false; estab[p] = false; pftmayexist[p] = true
            end
        end
    end
    for p in 1:np
        c = column[p]; l = landunit[p]; iv = itype[p] + 1
        if present[p] && (!survive[p] || nind[p] < nind_min)
            present[p] = false; fpcgrid[p] = 0.0; nind[p] = 0.0
        end
        if lun_itype[l] == istsoil
            if !present[p] && prec365[c] >= prec_min_estab && estab[p]
                if eco.twmax[iv] > 999.0 || agddtw[p] == 0.0
                    present[p] = true; nind[p] = 0.0; fpcgrid[p] = 0.000844
                    if pft.woody[iv] < 1.0; fpcgrid[p] = 0.05; end
                    leafcmax[p] = 1.0
                    if dstemc[p] <= 0.0; dstemc[p] = 0.1; end
                end
            end
        end
    end
    for p in 1:np
        iv = itype[p] + 1
        if present[p]
            if pft.woody[iv] == 1.0
                fpc_tree_total += fpcgrid[p]
                if estab[p]; npft_estab += 1; end
            elseif pft.woody[iv] < 1.0 && itype[p] > noveg
                ngrass += 1
            end
        end
    end
    for p in 1:np
        iv = itype[p] + 1
        if present[p] && pft.woody[iv] == 1.0 && estab[p]
            estab_rate = estab_max * (1.0 - exp(5.0*(fpc_tree_total-1.0))) / npft_estab
            estab_grid = estab_rate * (1.0 - fpc_tree_total)
            nind[p] += estab_grid
            lm_ind = leafcmax[p] * fpcgrid[p] / nind[p]
            if fpcgrid[p] > 0.0 && nind[p] > 0.0
                stocking = nind[p] / fpcgrid[p]
                stemdiam = (24.0 * dstemc[p] / (PI_ * stocking * pft.dwood[iv] * taper))^(1.0/3.0)
            else
                stemdiam = 0.0
            end
            crownarea[p] = min(eco.crownarea_max[iv], eco.allom1[iv]*stemdiam^eco.reinickerp[iv])
            if crownarea[p] > 0.0
                if pft.dsladlai[iv] > 0.0
                    lai_ind = max(0.001, ((exp(lm_ind*pft.dsladlai[iv]+log(pft.slatop[iv])) -
                                           pft.slatop[iv])/pft.dsladlai[iv])/crownarea[p])
                else
                    lai_ind = lm_ind * pft.slatop[iv] / crownarea[p]
                end
            else
                lai_ind = 0.0
            end
            fpc_ind = 1.0 - exp(-0.5*lai_ind)
            fpcgrid[p] = crownarea[p] * nind[p] * fpc_ind
        end
        if present[p] && pft.woody[iv] == 1.0
            fpc_total_new += fpcgrid[p]
        end
    end
    for p in 1:np
        iv = itype[p] + 1
        if fpc_total_new > 0.95
            if pft.woody[iv] == 1.0 && present[p]
                nind[p] = nind[p] * 0.95 / fpc_total_new
                fpcgrid[p] = fpcgrid[p] * 0.95 / fpc_total_new
            end
            fpc_total = 0.95
        else
            fpc_total = fpc_total_new
        end
    end
    for p in 1:np
        iv = itype[p] + 1
        if present[p] && pft.woody[iv] < 1.0
            if leafcmax[p] <= 0.0 || fpcgrid[p] <= 0.0
                present[p] = false; nind[p] = 0.0
            else
                nind[p] = 1.0; crownarea[p] = 1.0
                lm_ind = leafcmax[p] * fpcgrid[p] / nind[p]
                if pft.dsladlai[iv] > 0.0
                    lai_ind = max(0.001, ((exp(lm_ind*pft.dsladlai[iv]+log(pft.slatop[iv])) -
                                           pft.slatop[iv])/pft.dsladlai[iv])/crownarea[p])
                else
                    lai_ind = lm_ind * pft.slatop[iv] / crownarea[p]
                end
                fpc_ind = 1.0 - exp(-0.5*lai_ind)
                fpcgrid[p] = crownarea[p] * nind[p] * fpc_ind
                fpc_total += fpcgrid[p]
            end
        end
    end
    for p in 1:np
        if fpc_total > 1.0
            if itype[p] >= nc3_arctic_grass && fpcgrid[p] > 0.0
                fpcgridtemp = fpcgrid[p]
                fpcgrid[p] = max(0.0, fpcgrid[p] - (fpc_total-1.0))
                fpc_total = fpc_total - fpcgridtemp + fpcgrid[p]
            end
        end
        if fpcgrid[p] < 1.0e-15
            fpc_total -= fpcgrid[p]; fpcgrid[p] = 0.0; present[p] = false; nind[p] = 0.0
        end
        if fpc_total < 1.0 && itype[p] == noveg
            fpcgrid[p] = 1.0 - fpc_total
            fpc_total += fpcgrid[p]
        end
    end
    for p in 1:np
        iv = itype[p] + 1
        if pft.woody[iv] == 1.0 && nind[p] > 0.0 && leafcmax[p] > 0.0 && fpcgrid[p] > 0.0
            if eco.twmax[iv] < 999.0
                heatstress[p] = max(0.0, min(1.0, agddtw[p]/ramp_agddtw))
            else
                heatstress[p] = 0.0
            end
            bm_delta = max(0.0, annsum_npp[p] - annsum_litfall[p])
            lm_ind = leafcmax[p] * fpcgrid[p] / nind[p]
            if pft.dsladlai[iv] > 0.0
                greffic[p] = bm_delta / max(0.001, (exp(lm_ind*pft.dsladlai[iv]+log(pft.slatop[iv])) -
                                                    pft.slatop[iv])/pft.dsladlai[iv])
            else
                greffic[p] = bm_delta / (lm_ind * pft.slatop[iv])
            end
        else
            greffic[p] = 0.0; heatstress[p] = 0.0
        end
    end
    return fpc_total
end

# =============================================================================
# Checks
# =============================================================================
const RTOL = 1.0e-10

@testset "CNDV Fortran-parity harness" begin

# -------------------------------------------------------------------------
# A. LIVE prec365 feed through atm2lnd_update_acc_vars! (use_cndv=true).
#    Regression guard for the historical ZERO-BUFFER bug that made CNDV
#    establishment impossible.
# -------------------------------------------------------------------------
@testset "A. prec365 is fed (not a zero buffer)" begin
    C.varpar_init!(C.varpar, 1, 14, 2, 5)
    prev_cndv = C.varctl.use_cndv
    prev_cn   = C.varctl.use_cn
    C.varctl.use_cndv = true
    C.varctl.use_cn   = true
    try
        ng, nc, np = 1, 1, 1
        a2l = C.Atm2LndData()
        C.atm2lnd_init!(a2l, ng, nc, np)          # allocates prec365_col etc.
        @test length(a2l.prec365_col) == nc
        # constant precip forcing: rain+snow = P (mm H2O/s)
        P = 2.0e-5
        fill!(a2l.forc_rain_downscaled_col, P)
        fill!(a2l.forc_snow_downscaled_col, 0.0)
        fill!(a2l.forc_t_downscaled_col, TFRZ + 10.0)
        pg = [1]; pc = [1]
        dtime = 1800
        # step the accumulator many times; a constant-forcing running mean -> P.
        nsteps = 400
        for n in 1:nsteps
            C.atm2lnd_update_acc_vars!(a2l, 1:np, pg, pc;
                bounds_c = 1:nc, nstep = n, dtime = dtime)
        end
        pv = a2l.prec365_col[1]
        @test isfinite(pv)
        @test pv > 0.0                    # THE guard: not stuck at zero
        @test isapprox(pv, P; rtol = 1e-6)  # constant forcing -> running mean == P
        prec_min_estab = 100.0 / (365.0 * SECSPDAY)
        @test pv >= prec_min_estab        # exceeds establishment threshold -> CNDV CAN establish
    finally
        C.varctl.use_cndv = prev_cndv
        C.varctl.use_cn   = prev_cn
    end
end

# -------------------------------------------------------------------------
# B. cndv_light! value-parity vs independent CNDVLightMod oracle.
# -------------------------------------------------------------------------
@testset "B. Light value-parity vs Fortran oracle" begin
    s = build_scenario()
    inst = s.inst; d = inst.dgvs; eco = inst.dgv_ecophyscon; pft = s.pft; np = s.np
    # seed a present, over-dense canopy so light competition (>95%) fires
    d.present_patch    .= [false, true, true, true, false]
    d.nind_patch       .= [0.0, 0.08, 0.06, 1.0, 0.0]
    d.fpcgrid_patch    .= [0.0, 0.55, 0.55, 0.20, 0.0]
    d.fpcgridold_patch .= copy(d.fpcgrid_patch)
    d.crownarea_patch  .= [0.0, 12.0, 12.0, 1.0, 0.0]
    # large deadstemc -> crownarea saturates at crownarea_max (15) so that the two
    # recomputed tree FPCs sum > 0.95 and the light-competition (excess) branch fires.
    deadstemc = [0.0, 5.0e6, 5.0e6, 0.0, 0.0]
    leafcmax  = [0.0, 250.0, 220.0, 180.0, 0.0]
    mask = collect(Bool, d.present_patch)   # natvegp = present (as cndv_driver! rebuilds)

    # oracle on copies
    o_fpc  = copy(d.fpcgrid_patch); o_nind = copy(d.nind_patch); o_ca = copy(d.crownarea_patch)
    oracle_light!(o_fpc, o_nind, o_ca, mask, inst.patch.itype, pft, eco, deadstemc, leafcmax)

    # port
    C.cndv_light!(d, eco, deadstemc, leafcmax, pft, inst.patch, mask, 1:np;
                  bounds_gridcell = 1:s.ng)

    @test isapprox(d.fpcgrid_patch,   o_fpc; rtol = RTOL, atol = 1e-14)
    @test isapprox(d.nind_patch,      o_nind; rtol = RTOL, atol = 1e-14)
    @test isapprox(d.crownarea_patch, o_ca; rtol = RTOL, atol = 1e-14)
    # non-vacuous: competition must actually have reduced total tree fpc to 0.95
    tree_total = d.fpcgrid_patch[2] + d.fpcgrid_patch[3]
    @test isapprox(tree_total, 0.95; rtol = 1e-6)
end

# -------------------------------------------------------------------------
# C. cndv_establishment! value-parity vs oracle + a REAL establishment event.
# -------------------------------------------------------------------------
@testset "C. Establishment value-parity + PFT establishes" begin
    s = build_scenario()
    inst = s.inst; d = inst.dgvs; eco = inst.dgv_ecophyscon; pft = s.pft
    np = s.np; nc = s.nc
    # climate that permits establishment for all adapted PFTs
    d.tmomin20_patch .= TFRZ + 0.0        # 0C: within [tcmin,tcmax]=[-35,28]
    d.agdd20_patch   .= 600.0             # >= gddmin 350
    d.agdd_patch     .= 600.0
    d.agddtw_patch   .= 0.0
    d.pftmayexist_patch .= true
    # present canopy: trees p2,p3 present; grass p4 present; p5 ABSENT (to establish)
    d.present_patch  .= [false, true, true, true, false]
    d.nind_patch     .= [0.0, 0.05, 0.04, 1.0, 0.0]
    d.fpcgrid_patch  .= [0.0, 0.30, 0.25, 0.15, 0.0]
    d.crownarea_patch.= [0.0, 10.0, 10.0, 1.0, 0.0]
    prec365 = [2.0e-5]                    # NONZERO, >> prec_min_estab
    annsum_npp      = [0.0, 300.0, 250.0, 120.0, 0.0]
    annsum_litfall  = [0.0, 80.0, 70.0, 40.0, 0.0]
    deadstemc       = [0.0, 5000.0, 4000.0, 0.0, 0.0]
    leafcmax        = [0.0, 250.0, 220.0, 180.0, 0.0]

    # oracle on independent copies of ALL mutated state
    o_present  = copy(d.present_patch); o_may = copy(d.pftmayexist_patch)
    o_nind     = copy(d.nind_patch);    o_fpc = copy(d.fpcgrid_patch)
    o_ca       = copy(d.crownarea_patch)
    o_greffic  = copy(d.greffic_patch); o_heat = copy(d.heatstress_patch)
    o_leafcmax = copy(leafcmax)
    fpc_total_oracle = oracle_establishment!(o_present, o_may, o_nind, o_fpc, o_ca,
        o_greffic, o_heat, inst.patch.itype, inst.patch.column, inst.patch.landunit,
        inst.landunit.itype, pft, eco, prec365, annsum_npp, annsum_litfall,
        deadstemc, o_leafcmax, d.agddtw_patch, d.agdd20_patch, d.tmomin20_patch)

    # port (mutates dgvs in place; leafcmax passed as its own array)
    port_leafcmax = copy(leafcmax)
    C.cndv_establishment!(d, eco, prec365, annsum_npp, annsum_litfall,
        deadstemc, port_leafcmax, pft, inst.patch, inst.landunit, 1:np;
        bounds_gridcell = 1:s.ng)

    @test d.present_patch    == o_present
    @test isapprox(d.nind_patch,      o_nind;    rtol = RTOL, atol = 1e-14)
    @test isapprox(d.fpcgrid_patch,   o_fpc;     rtol = RTOL, atol = 1e-14)
    @test isapprox(d.crownarea_patch, o_ca;      rtol = RTOL, atol = 1e-14)
    @test isapprox(d.greffic_patch,   o_greffic; rtol = RTOL, atol = 1e-12)
    @test isapprox(d.heatstress_patch,o_heat;    rtol = RTOL, atol = 1e-14)
    @test isapprox(port_leafcmax,     o_leafcmax;rtol = RTOL, atol = 1e-14)

    # NON-VACUOUS establishment event: the absent adapted tree p5 MUST now be present.
    # This is unreachable by any harness that feeds prec365 = zeros.
    @test o_present[5] == true            # oracle agrees establishment fired
    @test d.present_patch[5] == true      # port established the new PFT
    @test d.fpcgrid_patch[5] > 0.0        # with nonzero cover
    # fpc_total must close to ~1 (bare ground fills remainder) — Fortran's error check
    @test isapprox(fpc_total_oracle, 1.0; atol = 1e-6)
end

end # top testset
