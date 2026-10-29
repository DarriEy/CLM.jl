# Shared FATES live-run helpers, included by test_fates_live_modes.jl and
# test_fates_spinup.jl (each file is self-contained; neither relies on the other
# having run first).
#
# _fates_save_globals / _fates_restore_globals snapshot and restore every FATES +
# varpar global clm_fates_init! mutates, plus the host `pftcon` (the live-run
# setup writes single pftcon fields by hand, which used to leave a partially
# written pftcon behind for later files).

const _C = CLM
# ---------------------------------------------------------------------------
# Shared helpers (build a driver-ready FATES-tagged single-column inst, run the
# real clm_drv! loop, walk the cohort/patch census).
# ---------------------------------------------------------------------------

# Snapshot of every FATES + varpar global clm_fates_init! mutates.
function _fates_save_globals()
    return (
        prt_global = _C.prt_global[], numel = _C.num_elements[],
        ellist = copy(_C.element_list), parteh = _C.hlm_parteh_mode[],
        use_sp = _C.hlm_use_sp[], hydro = _C.hlm_use_planthydro[],
        agetrk = _C.hlm_use_cohort_age_tracking[], nleafage = _C.nleafage[],
        numpft = _C.numpft[], tod = _C.hlm_current_tod[],
        nlevsc = _C.nlevsclass[], nlevca = _C.nlevcoage[], nlevage = _C.nlevage[],
        nlevdam = _C.nlevdamage[], nlevhgt = _C.nlevheight[],
        sfp = _C.SFParams[], edparams = _C.EDParams[], edpft = _C.EDPftvarcon_inst[],
        paramd = _C.ParamDerived[], restart = _C.hlm_is_restart[],
        nocomp = _C.hlm_use_nocomp[], biogeog = _C.hlm_use_fixed_biogeog[],
        luh = _C.hlm_use_luh[], inv = _C.hlm_use_inventory_init[],
        damage = _C.hlm_use_tree_damage[], doy = _C.hlm_day_of_year[],
        dpy = _C.hlm_days_per_year[], freqday = _C.hlm_freq_day[],
        nharv = _C.hlm_num_lu_harvest_cats[], ch4 = _C.hlm_use_ch4[],
        vert = _C.hlm_use_vertsoilc[], spit = _C.hlm_spitfire_mode[],
        paramfp = _C.fates_maxPatchesPerSite[],
        maxcomp = _C.max_comp_per_site[], npcomp = _C.fates_np_comp_scaling[],
        nspec = _C.hlm_nitrogen_spec[], pspec = _C.hlm_phosphorus_spec[],
        nlevsno = _C.varpar.nlevsno, nlevsoi = _C.varpar.nlevsoi,
        nlevgrnd = _C.varpar.nlevgrnd, nlevmaxu = _C.varpar.nlevmaxurbgrnd,
        pftcon = deepcopy(_C.pftcon))
end

function _fates_restore_globals(g)
    _C.prt_global[] = g.prt_global
    _C.num_elements[] = g.numel
    empty!(_C.element_list); append!(_C.element_list, g.ellist)
    _C.hlm_parteh_mode[] = g.parteh
    _C.hlm_use_sp[] = g.use_sp
    _C.hlm_use_planthydro[] = g.hydro
    _C.hlm_use_cohort_age_tracking[] = g.agetrk
    _C.nleafage[] = g.nleafage
    _C.numpft[] = g.numpft
    _C.hlm_current_tod[] = g.tod
    _C.nlevsclass[] = g.nlevsc
    _C.nlevcoage[] = g.nlevca
    _C.nlevage[] = g.nlevage
    _C.nlevdamage[] = g.nlevdam
    _C.nlevheight[] = g.nlevhgt
    _C.SFParams[] = g.sfp
    _C.EDParams[] = g.edparams
    _C.EDPftvarcon_inst[] = g.edpft
    _C.ParamDerived[] = g.paramd
    _C.hlm_is_restart[] = g.restart
    _C.hlm_use_nocomp[] = g.nocomp
    _C.hlm_use_fixed_biogeog[] = g.biogeog
    _C.hlm_use_luh[] = g.luh
    _C.hlm_use_inventory_init[] = g.inv
    _C.hlm_use_tree_damage[] = g.damage
    _C.hlm_day_of_year[] = g.doy
    _C.hlm_days_per_year[] = g.dpy
    _C.hlm_freq_day[] = g.freqday
    _C.hlm_num_lu_harvest_cats[] = g.nharv
    _C.hlm_use_ch4[] = g.ch4
    _C.hlm_use_vertsoilc[] = g.vert
    _C.hlm_spitfire_mode[] = g.spit
    _C.fates_maxPatchesPerSite[] = g.paramfp
    _C.max_comp_per_site[] = g.maxcomp
    _C.fates_np_comp_scaling[] = g.npcomp
    _C.hlm_nitrogen_spec[] = g.nspec
    _C.hlm_phosphorus_spec[] = g.pspec
    _C.varpar.nlevsno = g.nlevsno
    _C.varpar.nlevsoi = g.nlevsoi
    _C.varpar.nlevgrnd = g.nlevgrnd
    _C.varpar.nlevmaxurbgrnd = g.nlevmaxu
    for f in fieldnames(typeof(_C.pftcon))
        setfield!(_C.pftcon, f, getfield(g.pftcon, f))
    end
end

# Walk the cohort/patch linked lists: census + finiteness + total patch area.
function _fates_census(inst)
    site = inst.fates.sites[1]
    npatch = 0; ncoh = 0; totarea = 0.0; bad = String[]
    cp = site.oldest_patch
    while cp !== nothing
        npatch += 1
        isfinite(cp.area) || push!(bad, "patch.area")
        totarea += cp.area
        cc = cp.tallest
        while cc !== nothing
            ncoh += 1
            (isfinite(cc.n) && cc.n >= 0.0) || push!(bad, "cohort.n")
            (isfinite(cc.dbh) && cc.dbh > 0.0) || push!(bad, "cohort.dbh")
            cc = cc.shorter
        end
        cp = cp.younger
    end
    return (npatch=npatch, ncoh=ncoh, totarea=totarea, bad=bad)
end
