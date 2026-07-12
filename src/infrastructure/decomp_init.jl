# ==========================================================================
# Ported from: src/main/decompInitMod.F90
#
# Domain-decomposition initialization: given a global land-gridcell count
# (plus an optional clumps-per-proc and process count), assign gridcells to
# clumps round-robin, clumps to processes round-robin, and build the per-clump
# / per-proc subgrid bounds (begg/endg/...) and the gindex global-index maps.
#
# This is the index/bookkeeping layer ONLY — no MPI, no state distribution.
# It makes nclumps>1 and multi-process decompositions representable + correct.
#
# Single-process single-clump (npes=1, clump_pproc=1) reproduces the bounds
# implied by the existing get_proc_bounds / get_clump_bounds accessors exactly,
# so the existing single-clump driver path is unaffected.
#
# Fortran provenance:
#   - decompInit_lnd    : nclumps setup, clump->proc round robin (pid=mod(n-1,npes)),
#                         gridcell->clump assignment, per-clump/per-proc begg/endg,
#                         gindex_global (gdc2glo) build.
#   - decompInit_clumps : accumulate per-gridcell landunit/column/patch/cohort
#                         counts into per-clump/per-proc begl/endl/begc/endc/...
#                         and numg/numl/numc/nump/numCohort.
#   - decompInit_glcp   : gindex_grc / gindex_lun / gindex_col / gindex_patch
#                         (compressed, land-only, in task order).
#
# subgridMod is partly ported (its add_landunit!/add_column!/add_patch! half lives
# in infrastructure/init_subgrid.jl); its `subgrid_get_gcellinfo` counter is NOT,
# so the per-gridcell
# subgrid counts are supplied as arguments (defaulting to 1 landunit / 1 column /
# 1 patch / 0 cohorts per gridcell). The bookkeeping is otherwise faithful to
# the Fortran.
# ==========================================================================

"""
    decompInit!(numg::Integer;
                clump_pproc::Integer = 1,
                npes::Integer = 1,
                iam::Integer = 0,
                nsegspc::Integer = 20,
                nlunits_per_g::Union{Nothing,AbstractVector{<:Integer}} = nothing,
                ncols_per_g::Union{Nothing,AbstractVector{<:Integer}} = nothing,
                npatches_per_g::Union{Nothing,AbstractVector{<:Integer}} = nothing,
                ncohorts_per_g::Union{Nothing,AbstractVector{<:Integer}} = nothing,
                decomp_data::DecompData = decomp) -> DecompData

Initialize the global decomposition `decomp_data` from a land-gridcell count.

Assigns the `numg` land gridcells to `nclumps = clump_pproc * npes` clumps, the
clumps to `npes` processes round-robin (`owner = mod(cid-1, npes)`, matching the
Fortran `pid = mod(n-1, npes)`), and builds the per-clump and per-process subgrid
bounds plus the gindex global-index maps. `iam` selects which process's
`procinfo` is populated (the local process); all clumps' owners/bounds/counts are
populated globally so any process is representable.

Gridcells are distributed to clumps in global-gridcell order. When
`numg/nclumps < nsegspc` (the common case, and the default `nsegspc=20`,
matching CLM's `nsegspc` default) this is a plain round-robin
(`cid = mod(ng-1, nclumps) + 1`); otherwise it follows the Fortran segment
formula (contiguous segments). Gridcell global indices here are simply
`1:numg` (i.e. the land grid is already compressed — no ocean mask), so
`gindex_global == gindex_grc` and both list the local process's gridcells in
global order.

`nlunits_per_g` / `ncols_per_g` / `npatches_per_g` / `ncohorts_per_g`, if given,
are length-`numg` per-gridcell subgrid counts (indexed by global gridcell);
otherwise each gridcell has 1 landunit / 1 column / 1 patch / 0 cohorts.

Returns `decomp_data` (also mutated in place).
"""
function decompInit!(numg::Integer;
                     clump_pproc::Integer = 1,
                     npes::Integer = 1,
                     iam::Integer = 0,
                     nsegspc::Integer = 20,
                     nlunits_per_g::Union{Nothing,AbstractVector{<:Integer}} = nothing,
                     ncols_per_g::Union{Nothing,AbstractVector{<:Integer}} = nothing,
                     npatches_per_g::Union{Nothing,AbstractVector{<:Integer}} = nothing,
                     ncohorts_per_g::Union{Nothing,AbstractVector{<:Integer}} = nothing,
                     decomp_data::DecompData = decomp)

    # --- validate inputs (decompInit_lnd guards) ---
    clump_pproc > 0 || error("decompInit!: clump_pproc=$clump_pproc must be greater than 0")
    npes > 0        || error("decompInit!: npes=$npes must be greater than 0")
    numg > 0        || error("decompInit!: numg=$numg must be greater than 0")
    nsegspc > 0     || error("decompInit!: nsegspc=$nsegspc must be greater than 0")
    (0 <= iam <= npes - 1) || error("decompInit!: iam=$iam out of range 0:$(npes-1)")

    nclumps = clump_pproc * npes
    if nclumps < npes
        error("decompInit!: number of clumps=$nclumps is less than the number of processes=$npes")
    end
    if npes > numg
        error("decompInit!: number of processes ($npes) exceeds number of land gridcells ($numg)")
    end
    if nclumps > numg
        error("decompInit!: number of clumps ($nclumps) exceeds number of land gridcells ($numg)")
    end

    # per-gridcell subgrid counts (default 1 landunit/column/patch, 0 cohorts)
    lun_g = nlunits_per_g  === nothing ? fill(1, numg) : nlunits_per_g
    col_g = ncols_per_g    === nothing ? fill(1, numg) : ncols_per_g
    pch_g = npatches_per_g === nothing ? fill(1, numg) : npatches_per_g
    coh_g = ncohorts_per_g === nothing ? fill(0, numg) : ncohorts_per_g
    length(lun_g) == numg || error("decompInit!: nlunits_per_g length $(length(lun_g)) != numg $numg")
    length(col_g) == numg || error("decompInit!: ncols_per_g length $(length(col_g)) != numg $numg")
    length(pch_g) == numg || error("decompInit!: npatches_per_g length $(length(pch_g)) != numg $numg")
    length(coh_g) == numg || error("decompInit!: ncohorts_per_g length $(length(coh_g)) != numg $numg")

    # --- allocate / reset procinfo and clumps ---
    procinfo = ProcessorType()
    procinfo.nclumps   = clump_pproc
    procinfo.cid       = fill(-1, clump_pproc)
    procinfo.ncells    = 0
    procinfo.nlunits   = 0
    procinfo.ncols     = 0
    procinfo.npatches  = 0
    procinfo.nCohorts  = 0
    # beg initialized to 1, end to 0 for simple incremental accumulation
    procinfo.begg = 1; procinfo.begl = 1; procinfo.begc = 1; procinfo.begp = 1; procinfo.begCohort = 1
    procinfo.endg = 0; procinfo.endl = 0; procinfo.endc = 0; procinfo.endp = 0; procinfo.endCohort = 0

    clumps = [ClumpType() for _ in 1:nclumps]
    for cl in clumps
        cl.owner = -1
        cl.ncells = 0; cl.nlunits = 0; cl.ncols = 0; cl.npatches = 0; cl.nCohorts = 0
        cl.begg = 1; cl.begl = 1; cl.begc = 1; cl.begp = 1; cl.begCohort = 1
        cl.endg = 0; cl.endl = 0; cl.endc = 0; cl.endp = 0; cl.endCohort = 0
    end

    # --- assign clumps to procs round robin (Fortran: pid = mod(n-1,npes)) ---
    cid_local = 0
    for n in 1:nclumps
        pid = mod(n - 1, npes)
        clumps[n].owner = pid
        if iam == pid
            cid_local += 1
            (1 <= cid_local <= clump_pproc) ||
                error("decompInit!: round robin cid error n=$n pid=$pid npes=$npes")
            procinfo.cid[cid_local] = n
        end
    end

    # --- assign gridcells to clumps (in global gridcell order) ---
    # Fortran seglen1 branch: numg/nclumps < nsegspc  => plain round robin.
    seglen1 = (Float64(numg) / Float64(nclumps)) < Float64(nsegspc)

    lcid = Vector{Int}(undef, numg)  # clump id owning each global gridcell
    for ng in 1:numg
        if seglen1
            cid = mod(ng - 1, nclumps) + 1
        else
            rcid = (Float64(ng - 1) / Float64(numg)) * Float64(nsegspc) * Float64(nclumps)
            cid = mod(trunc(Int, rcid), nclumps) + 1
        end
        lcid[ng] = cid

        owner_cid = clumps[cid].owner

        # per-proc gridcell bookkeeping (relative to local process iam)
        if iam == owner_cid
            procinfo.ncells += 1
        end
        if iam > owner_cid
            procinfo.begg += 1
        end
        if iam >= owner_cid
            procinfo.endg += 1
        end

        # per-clump gridcell bookkeeping
        clumps[cid].ncells += 1
        for m in 1:nclumps
            owner_m = clumps[m].owner
            if (owner_m > owner_cid) || (owner_m == owner_cid && m > cid)
                clumps[m].begg += 1
            end
            if (owner_m > owner_cid) || (owner_m == owner_cid && m >= cid)
                clumps[m].endg += 1
            end
        end
    end

    # ----------------------------------------------------------------------
    # gindex_global (decompInit_lnd: gdc2glo). The land grid here is already
    # compressed (gridcell global indices 1:numg), so gdc2glo maps each task-
    # order land slot back to its global gridcell index. clumpcnt(cid) = start
    # gdc index of clump cid, walked in (pid, cid) order.
    # ----------------------------------------------------------------------
    clumpcnt = zeros(Int, nclumps)
    ag = 1
    for pid in 0:(npes - 1)
        for cid in 1:nclumps
            if clumps[cid].owner == pid
                clumpcnt[cid] = ag
                ag += clumps[cid].ncells
            end
        end
    end

    gdc2glo = zeros(Int, numg)  # task-order slot -> global gridcell index
    cnt = copy(clumpcnt)
    for an in 1:numg               # an = global gridcell index (compressed land grid)
        cid = lcid[an]
        slot = cnt[cid]
        gdc2glo[slot] = an
        cnt[cid] += 1
    end

    # gindex_global covers this process's gridcells (begg:endg), in task order.
    # The local array is sized by the PROC-RELATIVE endg (= ncells = endg-begg+1,
    # matching Fortran's bounds%endg from get_proc_bounds), not the absolute
    # global endg. procinfo.begg/endg are absolute task-order slot positions used
    # to index gdc2glo.
    endg_proc = procinfo.endg - procinfo.begg + 1
    gindex_global = Vector{Int}(undef, endg_proc)
    for n in procinfo.begg:procinfo.endg
        gindex_global[n - procinfo.begg + 1] = gdc2glo[n]
    end

    # ----------------------------------------------------------------------
    # decompInit_clumps: accumulate per-gridcell subgrid counts into per-clump
    # and per-proc begl/endl/begc/endc/begp/endp/begCohort/endCohort and the
    # global numg/numl/numc/nump/numCohort totals.
    #
    # allvecg(cid,k) = total subgrid count of kind k over all gridcells of clump
    # cid (summed over all processes). Build it directly from lcid + per-g counts.
    # ----------------------------------------------------------------------
    allvec_l = zeros(Int, nclumps)
    allvec_c = zeros(Int, nclumps)
    allvec_p = zeros(Int, nclumps)
    allvec_co = zeros(Int, nclumps)
    for g in 1:numg
        cid = lcid[g]
        allvec_l[cid]  += lun_g[g]
        allvec_c[cid]  += col_g[g]
        allvec_p[cid]  += pch_g[g]
        allvec_co[cid] += coh_g[g]
    end

    numl = 0; numc = 0; nump = 0; numCohort = 0
    for cid in 1:nclumps
        ilunits  = allvec_l[cid]
        icols    = allvec_c[cid]
        ipatches = allvec_p[cid]
        icohorts = allvec_co[cid]

        numl      += ilunits
        numc      += icols
        nump      += ipatches
        numCohort += icohorts

        clumps[cid].nlunits  += ilunits
        clumps[cid].ncols    += icols
        clumps[cid].npatches += ipatches
        clumps[cid].nCohorts += icohorts

        owner_cid = clumps[cid].owner
        for m in 1:nclumps
            owner_m = clumps[m].owner
            if (owner_m > owner_cid) || (owner_m == owner_cid && m > cid)
                clumps[m].begl      += ilunits
                clumps[m].begc      += icols
                clumps[m].begp      += ipatches
                clumps[m].begCohort += icohorts
            end
            if (owner_m > owner_cid) || (owner_m == owner_cid && m >= cid)
                clumps[m].endl      += ilunits
                clumps[m].endc      += icols
                clumps[m].endp      += ipatches
                clumps[m].endCohort += icohorts
            end
        end

        if iam == owner_cid
            procinfo.nlunits  += ilunits
            procinfo.ncols    += icols
            procinfo.npatches += ipatches
            procinfo.nCohorts += icohorts
        end
        if iam > owner_cid
            procinfo.begl      += ilunits
            procinfo.begc      += icols
            procinfo.begp      += ipatches
            procinfo.begCohort += icohorts
        end
        if iam >= owner_cid
            procinfo.endl      += ilunits
            procinfo.endc      += icols
            procinfo.endp      += ipatches
            procinfo.endCohort += icohorts
        end
    end

    # ----------------------------------------------------------------------
    # decompInit_glcp: build gindex_grc / gindex_lun / gindex_col / gindex_patch
    # (compressed, land-only, GLOBAL subgrid indices for this process's points).
    # No MPI: the Fortran gather-to-master / running-count / scatter-back collapses
    # to a single in-process walk over all numg land gridcells in compressed global
    # order, accumulating the running global subgrid index and recording it at the
    # local slots.
    #
    # gindex_grc: the compressed (land-only) global gridcell index of each local
    # gridcell. With a fully-compressed land grid (gdc2glo[n] in 1:numg, ordered),
    # the land-only global index of task-order slot n is gdc2glo[n] itself.
    # ----------------------------------------------------------------------
    lsize_g = endg_proc
    gindex_grc = Vector{Int}(undef, lsize_g)
    for n in 1:lsize_g
        gindex_grc[n] = gindex_global[n]   # land-only global gridcell index
    end

    # Per-gridcell subgrid counts on this process, in local (task) order.
    lcount = Vector{Int}(undef, lsize_g)
    ccount = Vector{Int}(undef, lsize_g)
    pcount = Vector{Int}(undef, lsize_g)
    cocount = Vector{Int}(undef, lsize_g)
    for n in 1:lsize_g
        g = gindex_global[n]   # global gridcell index of local slot n
        lcount[n]  = lun_g[g]
        ccount[n]  = col_g[g]
        pcount[n]  = pch_g[g]
        cocount[n] = coh_g[g]
    end

    # start[gi] = first GLOBAL subgrid index for local gridcell gi. Walk ALL land
    # gridcells in compressed global order, accumulating the global per-gridcell
    # counts (so the subgrid index is globally unique across processes, matching
    # the Fortran `count = count + temp` running index); record `start` only at
    # the local slots. local_of_global maps a land-only global gridcell index back
    # to the local (task-order) slot (0 if not on this process).
    local_of_global = zeros(Int, numg)
    for n in 1:lsize_g
        local_of_global[gindex_grc[n]] = n
    end

    # count_global[gg] = per-gridcell subgrid count, indexed by land-only global
    # gridcell index gg (gindex_grc here equals the global gridcell id, so global
    # gridcell index == global gridcell id on the compressed land grid).
    function _build_start(count_global::AbstractVector{<:Integer})
        start = zeros(Int, lsize_g)
        running = 1
        for gg in 1:numg            # gg = land-only global gridcell index
            li = local_of_global[gg]
            if li != 0
                start[li] = running
            end
            running += count_global[gg]
        end
        return start
    end

    start_l  = _build_start(lun_g)
    start_c  = _build_start(col_g)
    start_p  = _build_start(pch_g)
    start_co = _build_start(coh_g)

    # Expand per-gridcell starts into per-subgrid gindex (offset within gridcell).
    function _expand(start::Vector{Int}, count_local::Vector{Int}, total::Int)
        gindex = Vector{Int}(undef, total)
        k = 0
        for gi in 1:lsize_g
            for off in 0:(count_local[gi] - 1)
                k += 1
                gindex[k] = start[gi] + off
            end
        end
        return gindex
    end

    # proc-relative sizes (Fortran lsize_* = bounds%end* from get_proc_bounds);
    # these equal the per-process subgrid counts (== sum of local count arrays).
    lsize_l  = procinfo.endl      - procinfo.begl      + 1
    lsize_c  = procinfo.endc      - procinfo.begc      + 1
    lsize_p  = procinfo.endp      - procinfo.begp      + 1
    lsize_co = procinfo.endCohort - procinfo.begCohort + 1
    gindex_lun    = _expand(start_l,  lcount,  lsize_l)
    gindex_col    = _expand(start_c,  ccount,  lsize_c)
    gindex_patch  = _expand(start_p,  pcount,  lsize_p)
    gindex_cohort = _expand(start_co, cocount, lsize_co)

    # --- commit into decomp_data ---
    decomp_data.procinfo      = procinfo
    decomp_data.clumps        = clumps
    decomp_data.nclumps       = nclumps
    decomp_data.numg          = numg
    decomp_data.numl          = numl
    decomp_data.numc          = numc
    decomp_data.nump          = nump
    decomp_data.numCohort     = numCohort
    decomp_data.gindex_global = gindex_global
    decomp_data.gindex_grc    = gindex_grc
    decomp_data.gindex_lun    = gindex_lun
    decomp_data.gindex_col    = gindex_col
    decomp_data.gindex_patch  = gindex_patch
    decomp_data.gindex_cohort = gindex_cohort

    return decomp_data
end
