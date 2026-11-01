# ==========================================================================================
# fates_driver_restart.jl — the host side of a FATES restart (clmfates_interfaceMod.F90
# `restart`, :1745-2045): build the restart IO interface for the attached sites, write the
# packed restart vectors to a NetCDF file, and rebuild every site from such a file.
#
# The host (CLM) restart is separate (src/driver/restart_io.jl); these routines restart
# the FATES site -> patch -> cohort state only.
# ==========================================================================================

"""
    fates_restart_interface(fates) -> fates_restart_interface_type

The restart IO object for `fates`'s sites (one clump): the cohort dimension spans
`nsites * fates_maxElementsPerSite`, the column dimension the sites. Mirrors the Fortran
`restart` PART I-III (Init / SetThreadBoundsEach / restart_map / assemble types /
initialize_restart_vars), with site `s` at column `s` and its first cohort slot at
`(s-1)*fates_maxElementsPerSite + 1` (the Fortran `(g-1)*fates_maxElementsPerSite + 1`
for one site per gridcell).
"""
function fates_restart_interface(fates)
    nsites = fates.nsites
    nper = fates_maxElementsPerSite[]
    nper > 0 || fates_endrun("fates_restart_interface: fates_maxElementsPerSite is unset")
    fb = fates_bounds_type(; cohort_begin = 1, cohort_end = nsites * nper,
                             column_begin = 1, column_end = nsites)
    rio = fates_restart_interface_type()
    Init!(rio, 1, fb)
    SetThreadBoundsEach!(rio, 1, fb)
    rio.restart_map[1].site_index    = collect(1:nsites)
    rio.restart_map[1].cohort1_index = [(s - 1) * nper + 1 for s in 1:nsites]
    assemble_restart_output_types!(rio)
    initialize_restart_vars!(rio)
    return rio
end

_fates_rst_dim(v) = (v.vtype == cohort_r8 || v.vtype == cohort_int) ? "cohort" : "column"

"""
    write_fates_restart(path, fates)

Pack the FATES state (`set_restart_vectors!`) and write every restart variable to
`path` on the `cohort` / `column` dimensions (the `restartvar(flag='write')` loop).
"""
function write_fates_restart(path::AbstractString, fates)
    rio = fates_restart_interface(fates)
    set_restart_vectors!(rio, 1, fates.nsites, fates.sites)
    NCDataset(path, "c") do ds
        defDim(ds, "cohort", fates.nsites * fates_maxElementsPerSite[])
        defDim(ds, "column", fates.nsites)
        for v in rio.rvars
            data = v.is_int ? Int32.(v.int1d) : v.r81d
            defVar(ds, v.vname, data, (_fates_rst_dim(v),);
                   attrib = ["units" => v.units, "long_name" => v.long])
        end
    end
    return path
end

"""
    read_fates_restart!(path, inst)

Rebuild every FATES site of `inst.fates` from the restart file `path`, as the Fortran
`restart(flag='read')` does: read the vectors, rebuild the site skeleton
(`create_patchcohort_structure!`, which re-runs `init_site_vars!` + `zero_site!`),
unpack (`get_restart_vectors!`), then `ed_update_site(is_restarting=true)`, the
`wrap_update_hlmfates_dyn` hand-off to the HLM, `update_3dpatch_radiation!`, the
dynamics history, and the restart call of the global seed dispersal.
"""
function read_fates_restart!(path::AbstractString, inst)
    fates = inst.fates
    nsites = fates.nsites
    hlm_use_planthydro[] == itrue &&
        fates_endrun("read_fates_restart!: plant-hydraulics restarts (RestartHydrStates) are not wired")
    rio = fates_restart_interface(fates)
    NCDataset(path, "r") do ds
        for v in rio.rvars
            a = Array(ds[v.vname])
            if v.is_int
                v.int1d .= Int.(a)
            else
                v.r81d .= Float64.(coalesce.(a, NaN))
            end
        end
    end

    # A restart run starts from freshly allocated sites.
    for s in 1:nsites
        fates.sites[s] = ed_site_type()
        fates.sites[s].h_gid = s
    end
    create_patchcohort_structure!(rio, 1, nsites, fates.sites;
                                  current_tod = hlm_current_tod[],
                                  bc_in = fates.bc_in, bc_out = fates.bc_out)
    get_restart_vectors!(rio, 1, nsites, fates.sites)

    col = inst.column
    site_col = zeros(Int, nsites)
    let s = 0
        for c in 1:length(col.is_fates)
            col.is_fates[c] || continue
            s += 1
            s <= nsites || break
            site_col[s] = c
        end
    end
    for s in 1:nsites
        c = site_col[s]
        bc = fates.bc_in[s]
        if c > 0
            bc.max_rooting_depth_index_col =
                min(bc.nlevsoil, inst.active_layer.altmax_lastyear_indx_col[c])
        end
        ed_update_site(fates.sites[s], bc, fates.bc_out[s], true)
    end

    # wrap_update_hlmfates_dyn(..., is_initing_from_restart = .true.): snow bc packed,
    # the restarted site%snow_depth kept.
    fates_update_snow_burial!(inst, site_col; is_initing_from_restart = true)
    canopy_summarization!(nsites, fates.sites, fates.bc_in)
    fates_wrap_seed_outgoing!(fates)
    update_hlm_dynamics!(nsites, fates.sites, site_col, fates.bc_out)
    varctl.use_fates_bgc && fates_update_litr_lig_c_to_n!(inst)
    for s in 1:nsites
        c = site_col[s]
        c == 0 && continue
        fates_flush_hlm_canopy!(inst; c = c)
        for (ifp, p) in fates_veg_patches(fates.sites[s], c, col)
            fates_unpack_bcout_canopy_structure!(inst; s = s, c = c, p = p, ifp = ifp)
        end
        fates_set_filters!(inst; c = c, s = s, set_active = false)
    end

    update_3dpatch_radiation!(rio, nsites, fates.sites, fates.bc_out)
    if fates.hist !== nothing
        update_history_dyn!(fates.hist, 1, nsites, fates.sites, fates.bc_in)
    end
    fates_wrap_global_seed_dispersal!(fates; is_restart = true)
    return inst
end
