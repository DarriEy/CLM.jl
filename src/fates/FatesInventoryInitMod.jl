# FatesInventoryInitMod.jl
# FATES (Tier F) Batch 16 — inventory-based site initialization.
#
# Faithful port of CTSM FATES `main/FatesInventoryInitMod.F90` (1226 lines,
# Ryan Knox June 2017; Jessica Needham Oct 2023 column trim — FATES issue #1062).
#
# This module initializes the FATES patch/cohort demographic state from external
# forest-inventory files (the PSS patch file + CSS cohort file, "type1" format):
# it reads patch records (year, patch-id, land-use track, age, area) and cohort
# records (year, patch-id, dbh, height, pft, nplant), creates the corresponding
# patch and cohort objects via the already-ported demographic constructors, and
# assembles them into each site's age-ordered patch linked list / per-patch
# height-ordered cohort linked list.
#
# Notes on this port (vs. Fortran):
#   * The Fortran reads HLM file paths through fortran file units (shr_file_*).
#     We abstract the PSS/CSS records as *iterables of lines* (Vector{AbstractString})
#     so the parsing is fully testable from in-memory strings.  A thin wrapper
#     `initialize_sites_by_inventory!` keeps the public Fortran entry-point shape
#     but accepts an `inventory` provider that maps a site -> (format, lat, lon,
#     pss_lines, css_lines); a convenience file-reading provider is also given.
#   * `set_inventory_patch_type1!` / `set_inventory_cohort_type1!` read ONE record
#     from a stateful line cursor (a `LineCursor`) — mirroring the Fortran read of
#     one line per call and the `ios /= 0 -> exit` loop control.
#   * Reuses Batch 12–15 constructors: `Create!` (patch), `create_cohort`,
#     `sort_cohorts`, `count_cohorts`, `fuse_cohorts`, `fuse_patches`,
#     `InitPRTObject!`, `SetState!`, `StorageNutrientTarget`, allometry, litter
#     `InitConditions!` — none of which are reimplemented here.
#
# Upstream-Fortran quirks preserved (see inline comments):
#   * pft==0 is a SPECIAL CASE: one cohort is created per PFT (1..numpft), each
#     getting n = c_nplant*cpatch.area/numpft.
#   * cohort %n is scaled by cpatch.area (the absolute patch area in m2), NOT by a
#     fraction — `temp_cohort%n = c_nplant * cpatch%area / ncohorts_to_create`.
#   * patch area = p_area * AREA (fraction -> m2 of the notional 10000 m2 forest).
#   * the insertion-sort middle-of-list branch in the Fortran has a known quirk:
#     it does NOT `exit` after inserting, it keeps walking — preserved here.
#   * nutrient pools (N/P) initialized at the *p1* (max) stoichiometry, storage via
#     StorageNutrientTarget — comment in Fortran says "half way between max and
#     min" but the code uses nitr_stoich_p1/phos_stoich_p1 directly.
#   * `write_inventory_type1` lat/lon filename mangling kept faithfully.

# ----------------------------------------------------------------------------------
# Module-level parameters (Fortran module constants)
# ----------------------------------------------------------------------------------

const inv_debug_inv = false   # Fortran `debug_inv` — dev debug flag

# The maximum distance in degrees allowed between a site's coordinate defined in
# model memory and a physical site listed in the file.
const max_site_adjacency_deg = 0.05

const inv_do_inventory_out = false  # Fortran `do_inventory_out`

# Catch-bad-value thresholds from set_inventory_cohort_type1
const inv_abnormal_large_nplant = 1000.0   # Used to catch bad values
const inv_abnormal_large_dbh    = 500.0    # I've never heard of a tree > 3m (sic)
const inv_abnormal_large_height = 500.0    # I've never heard of a tree > 500m tall

# ----------------------------------------------------------------------------------
# Line cursor — a small stateful reader over an iterable of lines, used to mirror
# the Fortran "read one record per call, ios/=0 -> stop" file-unit semantics
# without any real file I/O.  The first data line of each file is the header.
# ----------------------------------------------------------------------------------

"""
    LineCursor

Stateful cursor over a vector of text lines.  `pos` is the index of the NEXT line
to be read (1-based).  `read_line!` returns `(line, ok)` where `ok=false` signals
end-of-data (the Fortran `ios /= 0`).
"""
mutable struct LineCursor
    lines::Vector{String}
    pos::Int
end
LineCursor(lines) = LineCursor(String.(collect(lines)), 1)

"""Read the next line; returns (line::String, ok::Bool). ok=false at EOF."""
function read_line!(cur::LineCursor)
    if cur.pos > length(cur.lines)
        return ("", false)
    end
    ln = cur.lines[cur.pos]
    cur.pos += 1
    return (ln, true)
end

"""Skip the header line (Fortran `read(...,fmt=*) header_str`)."""
function skip_header!(cur::LineCursor)
    read_line!(cur)
    return nothing
end

# Split a whitespace-delimited record into fields (Fortran list-directed read).
_inv_fields(line::AbstractString) = split(strip(line))

# ==================================================================================
# count_inventory_sites  (Fortran function)
# ==================================================================================

"""
    count_inventory_sites(sitelist_lines) -> Int

Count the number of inventory sites listed in the descriptor file (every line
after the single header line).  Mirrors the Fortran `count_inventory_sites`.
"""
function count_inventory_sites(sitelist_lines)
    cur = LineCursor(sitelist_lines)
    skip_header!(cur)   # header line
    nsites = 0
    while true
        _, ok = read_line!(cur)
        ok || break
        nsites += 1
    end
    return nsites
end

# ==================================================================================
# assess_inventory_sites  (Fortran subroutine)
# ==================================================================================

"""
    assess_inventory_sites(sitelist_lines, nsites)
        -> (inv_format_list, inv_pss_list, inv_css_list, inv_lat_list, inv_lon_list)

Parse the inventory descriptor file, line by line, extracting per-site:
format-id, latitude, longitude, pss-file-path, css-file-path.  Performs the same
sanity checks as the Fortran (lat in [-90,90]; lon wrapped to [0,360]).  The file
existence checks (`inquire`) are intentionally NOT performed here — the records
are abstracted away from the filesystem so this is testable in memory; the caller
(file provider) is responsible for resolving paths.

File format (one header line + one line per site):
    type  latitude  longitude  pss-name  css-name
"""
function assess_inventory_sites(sitelist_lines, nsites::Integer)
    inv_format_list = Vector{Int}(undef, nsites)
    inv_pss_list    = Vector{String}(undef, nsites)
    inv_css_list    = Vector{String}(undef, nsites)
    inv_lat_list    = Vector{Float64}(undef, nsites)
    inv_lon_list    = Vector{Float64}(undef, nsites)

    cur = LineCursor(sitelist_lines)
    skip_header!(cur)   # header line

    for isite in 1:nsites
        site_str, ok = read_line!(cur)
        if !ok
            fates_endrun("inventory descriptor file ended before nsites lines were read")
        end

        flds = _inv_fields(site_str)
        if length(flds) < 5
            fates_endrun("inventory site line has fewer than 5 fields: $(site_str)")
        end

        file_format = parse(Int, flds[1])
        site_lat    = parse(Float64, flds[2])
        site_lon    = parse(Float64, flds[3])
        pss_file    = String(flds[4])
        css_file    = String(flds[5])

        if site_lat < -90.0 || site_lat > 90.0
            fates_endrun("read invalid latitude: $(site_lat) from inventory site list")
        end

        # Longitude should be converted to positive coordinate
        if site_lon < 0.0
            site_lon = 360.0 + site_lon
        end
        if site_lon < 0.0 || site_lon > 360.0
            fates_endrun("read invalid longitude: $(site_lon) from inventory site list")
        end

        inv_format_list[isite] = file_format
        inv_pss_list[isite]    = pss_file
        inv_css_list[isite]    = css_file
        inv_lat_list[isite]    = site_lat
        inv_lon_list[isite]    = site_lon
    end

    return inv_format_list, inv_pss_list, inv_css_list, inv_lat_list, inv_lon_list
end

# ==================================================================================
# set_inventory_patch_type1!  (Fortran subroutine)
# ==================================================================================

"""
    set_inventory_patch_type1!(newpatch, pss_cur, ipa) -> (patch_name, ios)

Read one "type 1" PSS record from the line cursor and populate `newpatch`.
Returns the patch name string and an `ios` flag (0 = ok, nonzero = EOF/error,
mirroring the Fortran out-arg).  Record fields:
    time(year)  patch(string)  trk(LU index)  age(years)  area(fraction)
"""
function set_inventory_patch_type1!(newpatch::fates_patch_type, pss_cur::LineCursor,
                                    ipa::Integer)
    line, ok = read_line!(pss_cur)
    if !ok
        return ("", 1)   # ios /= 0
    end

    flds = _inv_fields(line)
    if length(flds) < 5
        return ("", 1)
    end

    p_time = parse(Float64, flds[1])   # read but unused downstream (kept for fidelity)
    p_name = String(flds[2])
    p_trk  = parse(Int, flds[3])       # LU track index (read but not used in type1)
    p_age  = parse(Float64, flds[4])
    p_area = parse(Float64, flds[5])

    patch_name = strip(p_name)

    # Fill in the patch's memory structures.
    newpatch.age       = p_age
    newpatch.age_class = get_age_class_index(newpatch.age)
    # QUIRK: patch fraction -> absolute m2 of the notional `area` (=10000 m2) forest.
    newpatch.area      = p_area * area

    # The litter and CWD pools are initialized to zero here.  The Fortran has a
    # long comment (RGK 06-2017) noting that estimating them from the cohort data
    # under a steady-state assumption is significant science modeling with no
    # simple first-hack solution — so they are left at zero.
    for el in 1:num_elements[]
        InitConditions!(newpatch.litter[el], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end

    return (String(patch_name), 0)
end

# ==================================================================================
# set_inventory_cohort_type1!  (Fortran subroutine)
# ==================================================================================

"""
    set_inventory_cohort_type1!(csite, bc_in, css_cur, patch_pointer_vec,
                                patch_name_vec) -> ios

Read one "type 1" CSS record from the line cursor, build the corresponding cohort
(or, if pft==0, one cohort per PFT) using FATES allometry + PARTEH, and insert it
into the height-ordered cohort list of the matching patch.  Returns `ios`
(0 = ok, nonzero = EOF) mirroring the Fortran out-arg.  Record fields:
    time(year)  patch(string)  dbh(cm)  height(m)  pft(int)  n(/m2)

`patch_pointer_vec` / `patch_name_vec` are parallel arrays (the Fortran `pp_array`
+ `patch_name_vec`) used to match a cohort to its patch by name.
"""
function set_inventory_cohort_type1!(csite::ed_site_type, bc_in,
                                     css_cur::LineCursor,
                                     patch_pointer_vec::Vector{fates_patch_type},
                                     patch_name_vec::Vector{String})

    line, ok = read_line!(css_cur)
    if !ok
        return 1   # ios /= 0
    end
    flds = _inv_fields(line)
    if length(flds) < 6
        return 1
    end

    c_time   = parse(Float64, flds[1])   # year, read but unused
    p_name   = String(flds[2])
    c_dbh    = parse(Float64, flds[3])
    c_height = parse(Float64, flds[4])
    c_pft    = parse(Int, flds[5])
    c_nplant = parse(Float64, flds[6])

    # Identify the patch based on the patch_name.
    matched_patch = false
    cpatch = patch_pointer_vec[1]   # placeholder; reassigned on match
    for ipa in 1:length(patch_name_vec)
        if strip(p_name) == strip(patch_name_vec[ipa])
            cpatch = patch_pointer_vec[ipa]
            matched_patch = true
        end
    end
    if !matched_patch
        fates_endrun("could not match a cohort with a patch (patch name: $(p_name))")
    end

    # Sanity checks on the input data (pft, nplant and dbh are the critical ones).
    if c_pft > numpft[]
        fates_endrun("inventory pft $(c_pft) greater than numpft $(numpft[])")
    end
    if c_pft < 0
        fates_endrun("inventory produced a cohort with <0 pft index: $(c_pft)")
    end
    if c_dbh < nearzero && c_height < nearzero
        fates_endrun("inventory dbh $(c_dbh) and height $(c_height) both <= 0; one must be positive")
    end
    if c_dbh > nearzero && c_height > nearzero
        fates_endrun("inventory dbh $(c_dbh) and height $(c_height) both positive; one must be <= 0")
    end
    if c_dbh > inv_abnormal_large_dbh
        fates_endrun("inventory produced a cohort with very large diameter [cm]: $(c_dbh)")
    end
    if c_height > inv_abnormal_large_height
        fates_endrun("inventory produced a cohort with very large height [m]: $(c_height)")
    end
    if c_nplant <= 0
        fates_endrun("inventory produced a cohort with <= 0 density /m2: $(c_nplant)")
    end
    if c_nplant > inv_abnormal_large_nplant
        fates_endrun("inventory produced a cohort with very large density /m2: $(c_nplant)")
    end

    # QUIRK (pft==0 SPECIAL CASE): create one identical cohort per PFT, each with
    # n = n_orig/numpft.
    if c_pft == 0
        ncohorts_to_create = numpft[]
    else
        ncohorts_to_create = 1
    end

    rstatus = 0   # recruit status (Fortran `rstatus`/`recruitstatus`)

    for i_pft in 1:ncohorts_to_create
        # temp_cohort is a scratch cohort needed by the allometry funcs.
        temp_cohort = fates_cohort_type()

        if c_pft != 0
            temp_cohort.pft = c_pft               # normal case
        else
            temp_cohort.pft = i_pft               # special case: one per PFT
        end

        # QUIRK: density scaled by absolute patch area (m2), split across the
        # cohorts created.
        temp_cohort.n = c_nplant * cpatch.area / float(ncohorts_to_create)

        temp_cohort.crowndamage = 1   # assume undamaged
        ft = temp_cohort.pft

        if c_dbh > 0.0
            temp_cohort.dbh = c_dbh
            h, _ = h_allom(c_dbh, ft)
            temp_cohort.height = h
        else
            temp_cohort.height = c_height
            d, _ = h2d_allom(c_height, ft)
            temp_cohort.dbh = d
        end

        temp_cohort.canopy_trim = 1.0

        # Determine the phenology status and the elongation factors.
        fnrt_drop_fraction = prt_params.phen_fnrt_drop_fraction[ft]
        stem_drop_fraction = prt_params.phen_stem_drop_fraction[ft]

        if prt_params.season_decid[ft] == itrue &&
           (csite.cstatus == phen_cstat_nevercold || csite.cstatus == phen_cstat_iscold)
            # Cold deciduous and the season is for leaves-off.
            temp_cohort.efleaf_coh = 0.0
            temp_cohort.effnrt_coh = 1.0 - fnrt_drop_fraction
            temp_cohort.efstem_coh = 1.0 - stem_drop_fraction
            temp_cohort.status_coh = leaves_off

        elseif prt_params.stress_decid[ft] == ihard_stress_decid ||
               prt_params.stress_decid[ft] == isemi_stress_decid
            # Drought deciduous: elongation factor from the site; tissues other than
            # leaves use a combination of elongation factor (e) and drop fraction (x)
            # so the remaining biomass = e when x=1, and original when x=0.
            temp_cohort.efleaf_coh = csite.elong_factor[ft]
            temp_cohort.effnrt_coh = 1.0 - (1.0 - temp_cohort.efleaf_coh) * fnrt_drop_fraction
            temp_cohort.efstem_coh = 1.0 - (1.0 - temp_cohort.efleaf_coh) * stem_drop_fraction
            if temp_cohort.efleaf_coh > 0.0
                temp_cohort.status_coh = leaves_on   # growing even if not fully flushed
            else
                temp_cohort.status_coh = leaves_off  # abscissing
            end
        else
            # Evergreen, or deciduous PFT during the growing season: fully flushed.
            temp_cohort.efleaf_coh = 1.0
            temp_cohort.effnrt_coh = 1.0
            temp_cohort.efstem_coh = 1.0
            temp_cohort.status_coh = leaves_on
        end

        # --- Allometric biomass pools (carbon) ---
        c_agw, _    = bagw_allom(temp_cohort.dbh, ft, temp_cohort.crowndamage, temp_cohort.efstem_coh)
        c_bgw, _    = bbgw_allom(temp_cohort.dbh, ft, temp_cohort.efstem_coh)
        c_leaf, _   = bleaf(temp_cohort.dbh, ft, temp_cohort.crowndamage,
                            temp_cohort.canopy_trim, temp_cohort.efleaf_coh)

        temp_cohort.l2fr = prt_params.allom_l2fr[ft]
        c_fnrt, _   = bfineroot(temp_cohort.dbh, ft, temp_cohort.canopy_trim,
                                temp_cohort.l2fr, temp_cohort.effnrt_coh)

        _, c_sapw, _ = bsap_allom(temp_cohort.dbh, ft, temp_cohort.crowndamage,
                                  temp_cohort.canopy_trim, temp_cohort.efstem_coh)
        c_struct, _  = bdead_allom(c_agw, c_bgw, c_sapw, ft)
        c_store, _   = bstore_allom(temp_cohort.dbh, ft, temp_cohort.crowndamage,
                                    temp_cohort.canopy_trim)

        # --- Build the PARTEH allocation object and set the per-element states ---
        prt_obj = InitPRTObject!()

        for el in 1:num_elements[]
            element_id = element_list[el]

            if element_id == carbon12_element
                m_struct = c_struct
                m_leaf   = c_leaf
                m_fnrt   = c_fnrt
                m_sapw   = c_sapw
                m_store  = c_store
                m_repro  = 0.0
            elseif element_id == nitrogen_element
                # QUIRK: comment says "half way between max and min" but code uses p1.
                m_struct = c_struct * prt_params.nitr_stoich_p1[ft, prt_params.organ_param_id[struct_organ]]
                m_leaf   = c_leaf   * prt_params.nitr_stoich_p1[ft, prt_params.organ_param_id[leaf_organ]]
                m_fnrt   = c_fnrt   * prt_params.nitr_stoich_p1[ft, prt_params.organ_param_id[fnrt_organ]]
                m_sapw   = c_sapw   * prt_params.nitr_stoich_p1[ft, prt_params.organ_param_id[sapw_organ]]
                m_repro  = 0.0
                m_store  = StorageNutrientTarget(ft, element_id, m_leaf, m_fnrt, m_sapw, m_struct)
            elseif element_id == phosphorus_element
                m_struct = c_struct * prt_params.phos_stoich_p1[ft, prt_params.organ_param_id[struct_organ]]
                m_leaf   = c_leaf   * prt_params.phos_stoich_p1[ft, prt_params.organ_param_id[leaf_organ]]
                m_fnrt   = c_fnrt   * prt_params.phos_stoich_p1[ft, prt_params.organ_param_id[fnrt_organ]]
                m_sapw   = c_sapw   * prt_params.phos_stoich_p1[ft, prt_params.organ_param_id[sapw_organ]]
                m_repro  = 0.0
                m_store  = StorageNutrientTarget(ft, element_id, m_leaf, m_fnrt, m_sapw, m_struct)
            else
                m_struct = 0.0; m_leaf = 0.0; m_fnrt = 0.0
                m_sapw = 0.0; m_store = 0.0; m_repro = 0.0
            end

            if hlm_parteh_mode[] == prt_carbon_allom_hyp ||
               hlm_parteh_mode[] == prt_cnp_flex_allom_hyp
                # Equally distribute leaf mass into available leaf age-bins.
                for iage in 1:nleafage[]
                    SetState!(prt_obj, leaf_organ, element_id, m_leaf / float(nleafage[]), iage)
                end
                SetState!(prt_obj, fnrt_organ,   element_id, m_fnrt)
                SetState!(prt_obj, sapw_organ,   element_id, m_sapw)
                SetState!(prt_obj, store_organ,  element_id, m_store)
                SetState!(prt_obj, struct_organ, element_id, m_struct)
                SetState!(prt_obj, repro_organ,  element_id, m_repro)
            else
                fates_endrun("Unspecified PARTEH module during inventory initialization")
            end
        end

        CheckInitialConditions(prt_obj)

        create_cohort(csite, cpatch, temp_cohort.pft, temp_cohort.n, temp_cohort.height,
                      temp_cohort.coage, temp_cohort.dbh, prt_obj,
                      temp_cohort.efleaf_coh, temp_cohort.effnrt_coh, temp_cohort.efstem_coh,
                      temp_cohort.status_coh, rstatus, temp_cohort.canopy_trim,
                      temp_cohort.c_area, 1, temp_cohort.crowndamage, csite.spread, bc_in)
        # temp_cohort goes out of scope (Fortran `deallocate(temp_cohort)`)
    end

    return 0
end

# ==================================================================================
# write_inventory_type1  (Fortran subroutine)
# ==================================================================================

"""
    write_inventory_type1(currentSite) -> (pss_name_out, css_name_out, pss_lines, css_lines)

Build the "type 1" PSS/CSS text for a site from its in-memory patch/cohort lists.
The Fortran writes these to files in the run directory with lat/long-tagged names;
here we return the generated file names and the line vectors so the caller can
write or test them.  The lat/lon filename mangling is preserved faithfully.
"""
function write_inventory_type1(currentSite::ed_site_type)
    # Generate pss/css file name based on the location of the site.
    ilat_int = abs(trunc(Int, currentSite.lat))
    ilat_dec = trunc(Int, 100000 * (abs(currentSite.lat) - float(ilat_int)))
    ilon_int = abs(trunc(Int, currentSite.lon))
    ilon_dec = trunc(Int, 100000 * (abs(currentSite.lon) - float(ilon_int)))

    ilat_sign = currentSite.lat >= 0.0 ? "N" : "S"
    ilon_sign = currentSite.lon >= 0.0 ? "E" : "W"

    latstr = string(lpad(ilat_int, 2, '0'), ".", lpad(ilat_dec, 5, '0'), ilat_sign)
    lonstr = string(lpad(ilon_int, 3, '0'), ".", lpad(ilon_dec, 5, '0'), ilon_sign)
    pss_name_out = string("pss_out_", latstr, "_", lonstr, ".txt")
    css_name_out = string("css_out_", latstr, "_", lonstr, ".txt")

    pss_lines = String[]
    css_lines = String[]
    push!(pss_lines, "time patch trk age area")
    push!(css_lines, "time patch dbh height pft nplant")

    ipatch = 0
    currentpatch = currentSite.youngest_patch
    while currentpatch !== nothing
        ipatch += 1
        patch_str = string("<patch_", lpad(ipatch, 4, '0'), ">")
        push!(pss_lines, string("0000 ", patch_str, " 2 ",
                                currentpatch.age, " ", currentpatch.area / area))

        currentcohort = currentpatch.tallest
        while currentcohort !== nothing
            push!(css_lines, string("0000 ", patch_str, " ",
                                    currentcohort.dbh, " ", -3.0, " ",
                                    currentcohort.pft, " ",
                                    currentcohort.n / currentpatch.area))
            currentcohort = currentcohort.shorter
        end
        currentpatch = currentpatch.older
    end

    return (pss_name_out, css_name_out, pss_lines, css_lines)
end

# ==================================================================================
# initialize_site_by_inventory!  (per-site core; the body of the Fortran site loop)
# ==================================================================================

"""
    initialize_site_by_inventory!(site, bc_in, file_format, pss_lines, css_lines)

Initialize ONE site's FATES demographic state from its PSS/CSS line records.
This is the per-site body of the Fortran `initialize_sites_by_inventory` loop:
build the age-ordered patch list from the PSS, insert cohorts from the CSS,
re-number patches, fuse cohorts/patches.  Returns the total cohort count.
"""
function initialize_site_by_inventory!(site::ed_site_type, bc_in,
                                       file_format::Integer,
                                       pss_lines, css_lines)

    # ---- PSS: build the age-ordered patch linked list ----
    pss_cur = LineCursor(pss_lines)
    skip_header!(pss_cur)   # header line

    # One quick pass to count patch records.
    npatches = 0
    while true
        _, ok = read_line!(pss_cur)
        ok || break
        npatches += 1
    end
    # Rewind and re-skip header.
    pss_cur = LineCursor(pss_lines)
    skip_header!(pss_cur)

    patch_name_vec    = Vector{String}(undef, npatches)
    patch_pointer_vec = Vector{fates_patch_type}(undef, npatches)

    for ipa in 1:npatches
        # Create the patch with nominal values; the PSS record fills in the rest.
        newpatch = fates_patch_type()
        Create!(newpatch, 0.0, 0.0, primaryland, fates_unset_int, num_swb,
                numpft[], site.nlevsoil, hlm_current_tod[], EDParams[].regeneration_model)
        newpatch.patchno = ipa
        newpatch.younger = nothing
        newpatch.older   = nothing

        patch_name = ""
        if file_format == 1
            patch_name, _ = set_inventory_patch_type1!(newpatch, pss_cur, ipa)
        end

        patch_name_vec[ipa]    = String(strip(patch_name))
        patch_pointer_vec[ipa] = newpatch

        if ipa == 1
            # First patch: starts as both oldest and youngest.
            site.youngest_patch = newpatch
            site.oldest_patch   = newpatch
        else
            # Insert this patch into the age-ordered LL (youngest..oldest).
            if newpatch.age <= site.youngest_patch.age
                # Youngest patch.
                newpatch.older = site.youngest_patch
                newpatch.younger = nothing
                site.youngest_patch.younger = newpatch
                site.youngest_patch = newpatch
            elseif newpatch.age > site.oldest_patch.age
                # Oldest patch.
                newpatch.older = nothing
                newpatch.younger = site.oldest_patch
                site.oldest_patch.older = newpatch
                site.oldest_patch = newpatch
            else
                # Somewhere in the middle.
                # QUIRK: the Fortran does NOT break after inserting — it keeps
                # walking the whole list.  Preserved here.
                currentpatch = site.youngest_patch
                while currentpatch !== nothing
                    olderpatch = currentpatch.older
                    if currentpatch.older !== nothing
                        if newpatch.age >= currentpatch.age &&
                           newpatch.age < olderpatch.age
                            newpatch.older   = currentpatch.older
                            newpatch.younger = currentpatch
                            currentpatch.older = newpatch
                            olderpatch.younger = newpatch
                        end
                    end
                    currentpatch = olderpatch
                end
            end
        end
    end

    # ---- CSS: insert cohorts, each matched to a patch by name ----
    css_cur = LineCursor(css_lines)
    skip_header!(css_cur)   # header line

    while true
        ios = 0
        if file_format == 1
            ios = set_inventory_cohort_type1!(site, bc_in, css_cur,
                                              patch_pointer_vec, patch_name_vec)
        end
        ios != 0 && break
    end

    # ---- Re-number patches, fuse cohorts, count cohorts ----
    ipa = 1
    total_cohorts = 0
    currentpatch = site.youngest_patch
    while currentpatch !== nothing
        currentpatch.patchno = ipa
        ipa += 1

        fuse_cohorts(site, currentpatch, bc_in)
        sort_cohorts(currentpatch)
        count_cohorts(currentpatch)
        total_cohorts += currentpatch.countcohorts

        currentpatch = currentpatch.older
    end

    if total_cohorts == 0
        fates_endrun("The inventory initialization produced no cohorts.")
    end

    # ---- Fuse patches ----
    fuse_patches(site, bc_in)

    if inv_do_inventory_out
        write_inventory_type1(site)
    end

    return total_cohorts
end

# ==================================================================================
# initialize_sites_by_inventory!  (public entry point)
# ==================================================================================

"""
    initialize_sites_by_inventory!(nsites, sites, bc_in;
                                   sitelist_lines, pss_provider, css_provider)

Public entry point — a faithful, testable analogue of the Fortran
`initialize_sites_by_inventory(nsites, sites, bc_in)`.

Because the Fortran reads HLM file paths through fortran file units, the file I/O
is abstracted via *providers*:
  * `sitelist_lines`  — iterable of lines of the inventory descriptor file.
  * `pss_provider(pss_path)` — returns the PSS file's lines for a given path.
  * `css_provider(css_path)` — returns the CSS file's lines for a given path.

The default providers (`read_lines_file`) read real files; tests pass closures over
in-memory dictionaries.  For each site this:
  1. Parses the descriptor, finds the most-proximal inventory site (min lat/lon
     distance), checking it is within `max_site_adjacency_deg`.
  2. Loads that site's PSS/CSS lines and calls `initialize_site_by_inventory!`.
"""
function initialize_sites_by_inventory!(nsites::Integer, sites::AbstractVector,
                                        bc_in::AbstractVector;
                                        sitelist_lines,
                                        pss_provider = read_lines_file,
                                        css_provider = read_lines_file)

    nfilesites = count_inventory_sites(sitelist_lines)
    if nfilesites < 1
        fates_endrun("The inventory file does not contain at least one site.")
    end

    inv_format_list, inv_pss_list, inv_css_list, inv_lat_list, inv_lon_list =
        assess_inventory_sites(sitelist_lines, nfilesites)

    for s in 1:nsites
        site = sites[s]

        # Identify the most proximal PSS/CSS couplet (Fortran minloc).
        invsite = 1
        best = Inf
        for i in 1:nfilesites
            d2 = (site.lat - inv_lat_list[i])^2 + (site.lon - inv_lon_list[i])^2
            if d2 < best
                best = d2
                invsite = i
            end
        end

        if sqrt(best) > max_site_adjacency_deg
            fates_endrun(string("Model site at lat:", site.lat, " lon:", site.lon,
                " has no reasonably proximal site in the inventory site list. Closest is at lat:",
                inv_lat_list[invsite], " lon:", inv_lon_list[invsite],
                ". Separation must be less than ", max_site_adjacency_deg, " degrees."))
        end

        pss_lines = pss_provider(inv_pss_list[invsite])
        css_lines = css_provider(inv_css_list[invsite])

        initialize_site_by_inventory!(site, bc_in[s], inv_format_list[invsite],
                                      pss_lines, css_lines)
    end

    return nothing
end

"""Default file provider: read all lines of `path` into a Vector{String}."""
read_lines_file(path::AbstractString) = readlines(path)
