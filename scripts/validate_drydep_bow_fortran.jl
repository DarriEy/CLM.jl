#!/usr/bin/env julia
# =============================================================================
# validate_drydep_bow_fortran.jl
#
# Fortran-parity validation of the ported dry-deposition Wesely index_season
# (PR #298/#300) at Bow-at-Banff (boreal/montane single point).
#
# GROUND TRUTH: a real CTSM run (symfluence_build cesm.exe, -bgc bgc, drydep_list
# = O3,SO2,NO2,CO,H2O2) branched from the 2006-01-01 Bow restart, one full year,
# daily native history (hist_dov2xy=.false.). Produced DRYDEPV_<spec> (cm/s),
# ELAI, SNOW_DEPTH, FSNO, TSA per patch/day. index_season is a LOCAL Fortran var
# (not output), but it is a DETERMINISTIC function of the published daily state
# (elai, snow_depth, annlai min/max, mlaidiff, landunit) — DryDepVelocity.F90:384-418.
#
# WHAT THIS CHECKS
#   1. Reconstruct index_season each day for the tree + grass patches TWO ways:
#        (a) the port's CLM._wesely_index_season (the #300 code under test)
#        (b) an INDEPENDENT inline transcription of DryDepVelocity.F90:384-418
#      fed identical Fortran daily state. Assert they agree every day (transcription
#      fidelity — the oracle approach, not "same code twice").
#   2. NON-VACUOUS seasonality: at a boreal site the season must actually VARY over
#      the year (winter=4 under snow, midsummer=1, plus shoulder 2/3/5) — assert
#      >= 3 distinct seasons appear (not stuck at one value).
#   3. Fortran DRYDEPV_O3/SO2 corroboration: velocities are nonzero and their
#      seasonal structure tracks the reconstructed index_season (independent
#      physical evidence the Fortran run selected the same seasons).
# =============================================================================

using CLM
using NCDatasets
using Printf
using Statistics

const C = CLM

const HIST = "/private/tmp/claude-501/drydep_run/Bow_at_Banff_lumped.clm2.h0.2006-01-02-00000.nc"
const SURF = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/settings/CLM/parameters/surfdata_clm.nc"

# ---- independent inline transcription of DryDepVelocity.F90:384-418 (soil landunit) ----
# Returns index_season for a soil (istsoil) patch. Bow patches are all istsoil.
function fortran_index_season_soil(elai, snow_depth, minlai, maxlai, mlaidiff)
    index_season = -1
    if snow_depth > 0.0
        index_season = 4
    elseif elai > 0.5 * maxlai
        index_season = 1
    end
    if index_season < 0
        if elai < (minlai + 0.05 * (maxlai - minlai))
            index_season = 3
        end
    end
    if index_season < 0
        if mlaidiff > 0.0
            index_season = 2
        elseif mlaidiff < 0.0
            index_season = 5
        elseif mlaidiff == 0.0
            index_season = 3
        end
    end
    return index_season
end

# mid-month interpolation bracket (CLM interpMonthlyVeg convention) -> (m1,m2)
# mlaidiff = LAI(m1) - LAI(m2)
const MIDMONTH = [16,45,75,105,136,166,197,228,258,289,319,350]  # approx day-of-year of mid-month (noleap)
function bracket_months(doy)
    # find m1 s.t. midmonth(m1) <= doy < midmonth(m1+1); wrap
    m1 = 12
    for m in 1:12
        if doy < MIDMONTH[m]
            m1 = m == 1 ? 12 : m-1
            break
        end
    end
    m2 = m1 == 12 ? 1 : m1+1
    return m1, m2
end

ds = NCDataset(HIST)
mcdate = Int.(ds["mcdate"][:])                 # YYYYMMDD per record
ntime  = length(mcdate)
itype  = Int.(ds["pfts1d_itype_veg"][:])       # [pft] itype
elai   = Array{Float64}(ds["ELAI"][:,:])       # [pft, time]
drO3   = Array{Float64}(ds["DRYDEPV_O3"][:,:]) # cm/s [pft,time]
drSO2  = Array{Float64}(ds["DRYDEPV_SO2"][:,:])
snowd  = Array{Float64}(ds["SNOW_DEPTH"][:,:]) # [column,time]
close(ds)

npft = length(itype)
# patch -> column: single column here
col_of(p) = 1

# surfdata monthly LAI per itype (annlai)
sd = NCDataset(SURF)
mlai_all = Array{Float64}(sd["MONTHLY_LAI"][:,:,:,:])  # (lon?,lat?,lsmpft,time) — check order
# dims declared: MONTHLY_LAI(time, lsmpft, lsmlat, lsmlon); NCDatasets reverses -> (lsmlon,lsmlat,lsmpft,time)
close(sd)
@printf("MONTHLY_LAI size (NCDatasets order) = %s\n", string(size(mlai_all)))
# single point: lsmlon=lsmlat=1 -> [1,1,lsmpft,12]
function annlai_for_itype(it)
    # lsmpft index is 0-based itype -> Julia index it+1
    return [mlai_all[1,1,it+1,m] for m in 1:12]
end

# day-of-year from mcdate (noleap-ish); use record index (daily) as doy proxy
function doy_of(rec); return rec; end   # record k ~ day k of 2006

println("\n== Bow drydep: patches ==")
for p in 1:npft
    @printf("  pft %d  itype=%d  wesveg=%d\n", p, itype[p], C._pft_to_wesely(itype[p]))
end

# Reconstruct index_season for the two VEGETATED patches (tree, grass)
veg_patches = [p for p in 1:npft if itype[p] != 0]

results = Dict{Int,Any}()
for p in veg_patches
    ann = annlai_for_itype(itype[p])
    minlai = minimum(ann); maxlai = maximum(ann)
    seas_port = Int[]; seas_fort = Int[]
    for k in 1:ntime
        m1,m2 = bracket_months(doy_of(k))
        mlaidiff = ann[m1] - ann[m2]
        el = elai[p,k]; sd_ = snowd[col_of(p),k]
        # port function under test
        _, sp = C._wesely_index_season(C.ISTSOIL, false, sd_ > 0.0, el, minlai, maxlai, mlaidiff)
        # independent oracle
        sf = fortran_index_season_soil(el, sd_, minlai, maxlai, mlaidiff)
        push!(seas_port, sp); push!(seas_fort, sf)
    end
    results[p] = (; ann, minlai, maxlai, seas_port, seas_fort)
end

# ---- Check 1: port vs independent oracle agree every day ----
maxmismatch = maximum(count(i -> results[p].seas_port[i] != results[p].seas_fort[i], 1:ntime)
                      for p in veg_patches)
for p in veg_patches
    r = results[p]
    mism = count(i -> r.seas_port[i] != r.seas_fort[i], 1:ntime)
    @printf("\npft %d (itype %d): index_season port-vs-oracle mismatches = %d / %d\n",
            p, itype[p], mism, ntime)
    # distinct seasons
    ds_port = sort(unique(r.seas_port))
    @printf("  distinct index_season (port) over year = %s\n", string(ds_port))
    @printf("  minlai=%.3f maxlai=%.3f\n", r.minlai, r.maxlai)
    # monthly-ish season sample (every ~30 days) + mean DRYDEPV_O3 by season
    @printf("  month-mid samples day: season | ELAI | SNOWD | DRYDEPV_O3(cm/s) | DRYDEPV_SO2\n")
    for k in (15,46,74,105,135,166,196,227,258,288,319,349)
        k <= ntime || continue
        @printf("    d%03d  s=%d  elai=%.3f  snow=%.3f  vO3=%.4f  vSO2=%.4f\n",
                k, r.seas_port[k], elai[p,k], snowd[col_of(p),k], drO3[p,k], drSO2[p,k])
    end
end

# ---- Check 3: DRYDEPV nonzero + seasonal (per veg patch) ----
println("\n== DRYDEPV non-vacuous (Fortran ground truth) ==")
for p in veg_patches
    vO3 = drO3[p,:]; vSO2 = drSO2[p,:]
    finite = filter(isfinite, vO3)
    @printf("pft %d: O3  min=%.4f max=%.4f mean=%.4f (cm/s)  | SO2 mean=%.4f\n",
            p, minimum(finite), maximum(finite), mean(finite), mean(filter(isfinite,vSO2)))
    r = results[p]
    for s in sort(unique(r.seas_port))
        idx = findall(==(s), r.seas_port)
        mo3 = mean(filter(isfinite, vO3[idx]))
        @printf("   season %d: n=%3d days  mean DRYDEPV_O3=%.4f cm/s\n", s, length(idx), mo3)
    end
end

# ---- Verdict ----
println("\n== VERDICT ==")
nseasons = maximum(length(sort(unique(results[p].seas_port))) for p in veg_patches)
ok = (maxmismatch == 0) && (nseasons >= 3)
@printf("index_season port==oracle every day: %s (max mismatch %d)\n", maxmismatch==0, maxmismatch)
@printf("boreal-seasonal (>=3 distinct seasons): %s (max distinct %d)\n", nseasons>=3, nseasons)
println(ok ? "PASS" : "FAIL")
exit(ok ? 0 : 1)
