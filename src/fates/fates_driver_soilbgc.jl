# fates_driver_soilbgc.jl
# FATES <-> CLM soil-biogeochemistry host coupling (use_fates_bgc).
#
# Host-side ports of the clmfates_interfaceMod.F90 routines that hand the FATES
# litter fluxes to CLM's decomposition cascade and hand the CLM soil carbon
# diagnostics back to FATES:
#
#   UpdateCLitterFluxes    (clmfates_interfaceMod.F90:1367-1441) -> fates_update_c_litter_fluxes!
#   UpdateNLitterFluxes    (clmfates_interfaceMod.F90:1294-1363) -> fates_update_n_litter_fluxes!
#   wrap_update_hifrq_hist (clmfates_interfaceMod.F90:2980-3034) -> fates_wrap_update_hifrq_hist!
#   wrap_update_hlmfates_dyn, litr_lig_c_to_n block (:1525-1545)  -> fates_update_litr_lig_c_to_n!
#
# The Fortran calls UpdateC/NLitterFluxes from CStateUpdate1 / NStateUpdate1 for every
# `col%is_fates` column of the bgc soil filter (CNCStateUpdate1Mod.F90:193-200,
# CNNStateUpdate1Mod.F90:149-156); the port calls them from cn_driver_no_leaching!
# immediately before c_state_update1! / n_state_update1!, which is the same point in
# the step (after the CNZero SetValues, before the cascade HR/transfer terms are added).
#
# Site <-> column map: the s-th `is_fates` column hosts FATES site s (the single-clump
# `f2hmap(nc)%fcolumn` / `hsites`), the convention every FATES driver hook uses.

"""
    fates_column_sites(col_is_fates, nsites) -> Vector{Int}

Per-column FATES site index (Fortran `f2hmap%hsites`): the s-th `is_fates` column maps
to site s; non-FATES columns (and FATES columns beyond `nsites`) map to 0.
"""
function fates_column_sites(col_is_fates::AbstractVector{Bool}, nsites::Integer)
    hs = zeros(Int, length(col_is_fates))
    s = 0
    for c in eachindex(col_is_fates)
        col_is_fates[c] || continue
        s += 1
        s <= nsites || break
        hs[c] = s
    end
    return hs
end

"""
    fates_update_c_litter_fluxes!(fates, cf_soil, col_is_fates, mask_soilc;
                                  dtime, nlevdecomp, decomp_method, i_met_lit=1)

Port of `UpdateCLitterFluxes`. For each FATES column of the bgc soil filter:
run `FluxIntoLitterPools` (fragmentation of the FATES litter pools into
labile/cellulose/lignin fluxes per decomposition layer, g/m3/s), add
`litt_flux_{lab,cel,lig}_c_si * dtime` to `decomp_cpools_sourcesink_col` of the
metabolic / cellulose / lignin litter pools (lignin == cellulose pool under MIMICS),
and diagnose the column litter input `fates_litter_flux` (gC/m2/s) for the CN
balance check.

`decomp_cpools_sourcesink_col` must already be zeroed for the step (CNZero
SetValues — `_zero_soilbgc_cflux!` in the port), exactly as the Fortran assumes.
"""
function fates_update_c_litter_fluxes!(fates, cf_soil, col_is_fates::AbstractVector{Bool},
                                       mask_soilc::AbstractVector{Bool};
                                       dtime::Real, nlevdecomp::Int, decomp_method::Int,
                                       i_met_lit::Int = 1)
    fates === nothing && return nothing
    hsites = fates_column_sites(col_is_fates, fates.nsites)
    ss = cf_soil.decomp_cpools_sourcesink_col
    i_cel_lit = i_met_lit + 1
    # CENTURY carries a separate lignin litter pool; MIMICS's structural pool holds
    # both cellulose and lignin.
    i_lig_lit = decomp_method == MIMICS_DECOMP ? i_cel_lit : i_cel_lit + 1
    for c in eachindex(mask_soilc)
        mask_soilc[c] || continue
        s = hsites[c]
        s == 0 && continue
        site = fates.sites[s]
        bc_in = fates.bc_in[s]
        bc_out = fates.bc_out[s]

        FluxIntoLitterPools(site, bc_in, bc_out)

        flit = 0.0
        for j in 1:nlevdecomp
            # (gC/m3/timestep)
            ss[c, j, i_met_lit] += bc_out.litt_flux_lab_c_si[j] * dtime
            # Used for mass balance checking (gC/m2/s)
            flit += bc_out.litt_flux_lab_c_si[j] * bc_in.dz_decomp_sisl[j]
        end
        for j in 1:nlevdecomp
            ss[c, j, i_cel_lit] += bc_out.litt_flux_cel_c_si[j] * dtime
        end
        flit += sum(bc_out.litt_flux_cel_c_si[j] * bc_in.dz_decomp_sisl[j] for j in 1:nlevdecomp)
        for j in 1:nlevdecomp
            ss[c, j, i_lig_lit] += bc_out.litt_flux_lig_c_si[j] * dtime
        end
        flit += sum(bc_out.litt_flux_lig_c_si[j] * bc_in.dz_decomp_sisl[j] for j in 1:nlevdecomp)
        if length(cf_soil.fates_litter_flux) >= c
            cf_soil.fates_litter_flux[c] = flit
        end
    end
    return nothing
end

"""
    fates_update_n_litter_fluxes!(fates, nf_soil, col_is_fates, mask_soilc)

Port of `UpdateNLitterFluxes`. The Fortran (ctsm5.3) zeroes the FATES column's
`decomp_npools_sourcesink_col` and sets `fates_litter_flux = 0`: the N litter
transfer is commented out upstream (FATES N coupling not yet active), so FATES
litter enters the CLM litter pools as carbon only.
"""
function fates_update_n_litter_fluxes!(fates, nf_soil, col_is_fates::AbstractVector{Bool},
                                       mask_soilc::AbstractVector{Bool})
    fates === nothing && return nothing
    hsites = fates_column_sites(col_is_fates, fates.nsites)
    ss = nf_soil.decomp_npools_sourcesink_col
    for c in eachindex(mask_soilc)
        mask_soilc[c] || continue
        hsites[c] == 0 && continue
        ss[c, :, :] .= 0.0
        # Fortran assigns the whole array (`nf_soil%fates_litter_flux = 0._r8`).
        fill!(nf_soil.fates_litter_flux, 0.0)
    end
    return nothing
end

"""
    fates_wrap_update_hifrq_hist!(inst; dtime, decomp_active)

Port of `wrap_update_hifrq_hist`: pass the host's heterotrophic respiration and
soil/litter carbon totals to FATES (`bc_in.tot_het_resp/tot_somc/tot_litc`, zero
without soil decomposition) and fill the FATES per-timestep history buffers
(FATES_HET_RESP = HR, FATES_NEP = -HR + FATES NEP terms). The Fortran calls it
after EcosystemDynamicsPostDrainage, so HR is the current step's.
"""
function fates_wrap_update_hifrq_hist!(inst::CLMInstances; dtime::Real,
                                       decomp_active::Bool)
    fates = inst.fates
    fates === nothing && return inst
    cf = inst.soilbiogeochem_carbonflux
    cs = inst.soilbiogeochem_carbonstate
    fcol = fates_site_columns_from_flags(inst.column.is_fates, fates.nsites)
    for s in 1:fates.nsites
        bc = fates.bc_in[s]
        c = fcol[s]
        if decomp_active && c != 0
            bc.tot_het_resp = cf.hr_col[c]
            bc.tot_somc     = cs.totsomc_col[c]
            bc.tot_litc     = cs.totlitc_col[c]
        else
            bc.tot_het_resp = 0.0
            bc.tot_somc     = 0.0
            bc.tot_litc     = 0.0
        end
    end
    fates_hifrq_history_step!(inst; dt_tstep = dtime)
    return inst
end

"""
    fates_site_columns_from_flags(col_is_fates, nsites) -> Vector{Int}

Host column of each FATES site (Fortran `f2hmap%fcolumn`): the s-th `is_fates`
column; 0 where a site has no flagged column.
"""
function fates_site_columns_from_flags(col_is_fates::AbstractVector{Bool}, nsites::Integer)
    fc = zeros(Int, nsites)
    s = 0
    for c in eachindex(col_is_fates)
        col_is_fates[c] || continue
        s += 1
        s <= nsites || break
        fc[s] = c
    end
    return fc
end

"""
    fates_update_litr_lig_c_to_n!(inst)

The `wrap_update_hlmfates_dyn` litter-quality hand-off (clmfates_interfaceMod.F90,
"FATES calculation of ligninNratio"): with soil decomposition active, copy
`bc_out.litt_flux_ligc_per_n` into the host column's `litr_lig_c_to_n_col` (read
by MIMICS; CENTURY carries but never reads it). The first-restart-step skip is
moot here (no restart path).
"""
function fates_update_litr_lig_c_to_n!(inst::CLMInstances)
    fates = inst.fates
    fates === nothing && return inst
    cf = inst.soilbiogeochem_carbonflux
    isempty(cf.litr_lig_c_to_n_col) && return inst
    fcol = fates_site_columns_from_flags(inst.column.is_fates, fates.nsites)
    for s in 1:fates.nsites
        c = fcol[s]
        c == 0 && continue
        cf.litr_lig_c_to_n_col[c] = fates.bc_out[s].litt_flux_ligc_per_n
    end
    return inst
end
