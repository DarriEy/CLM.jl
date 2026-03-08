# ==========================================================================
# Ported from: src/main/readParamsMod.F90
# Read parameters from clm5_params.nc — populates pftcon and scalar params
# ==========================================================================

"""
    readParameters!(paramfile::String)

Master parameter reader. Opens clm5_params.nc, reads PFT parameters into the
global `pftcon`, and reads scalar parameters used by various modules.
"""
function readParameters!(paramfile::String)
    ds = NCDataset(paramfile, "r")
    try
        # Read PFT parameters
        pftcon_read!(pftcon, ds)

        # Read initVertical scalar parameters
        readParams_initVertical!(ds)

        # Read saturated excess runoff parameters (SaturatedExcessRunoffMod)
        if haskey(ds, "fff")
            sat_excess_runoff_params.fff = ds["fff"][1]
        end

        # Read photosynthesis parameters (PhotosynthesisMod)
        readParams_photosynthesis!(ds)
    finally
        close(ds)
    end
    nothing
end

# --- Helper: read a 1D PFT variable from NetCDF ---

"""
    _read_pft_var!(ds, varname, dest; offset=0)

Read a 1D variable dimensioned on `pft` from the dataset and store it into
`dest`. NetCDF arrays are 0-based PFT (size npft_file), Julia arrays are
1-based (size MXPFT+1). `offset` shifts the start position in `dest`.
"""
function _read_pft_var!(ds::NCDataset, varname::String, dest::Vector{Float64};
                         offset::Int=0)
    if haskey(ds, varname)
        data = Array(ds[varname])
        npft_file = length(data)
        n = min(npft_file, length(dest) - offset)
        for i in 1:n
            val = data[i]
            dest[i + offset] = ismissing(val) ? 0.0 : Float64(val)
        end
    end
    nothing
end

function _read_pft_var_int!(ds::NCDataset, varname::String, dest::Vector{Int};
                             offset::Int=0)
    if haskey(ds, varname)
        data = Array(ds[varname])
        npft_file = length(data)
        n = min(npft_file, length(dest) - offset)
        for i in 1:n
            val = data[i]
            dest[i + offset] = ismissing(val) ? 0 : Int(round(val))
        end
    end
    nothing
end

"""
    _read_scalar(ds, varname, default) -> Float64

Read a scalar variable from the dataset with a fallback default.
"""
function _read_scalar(ds::NCDataset, varname::String, default::Float64)
    if haskey(ds, varname)
        v = Array(ds[varname])
        val = v isa AbstractArray ? v[1] : v
        return ismissing(val) ? default : Float64(val)
    end
    return default
end

"""
    pftcon_read!(p::PftconType, ds::NCDataset)

Populate PftconType arrays from an open clm5_params.nc dataset.
Corresponds to Fortran pftconMod::InitRead.
"""
function pftcon_read!(p::PftconType, ds::NCDataset)
    # Ensure pftcon is allocated
    if isempty(p.slatop)
        pftcon_allocate!(p)
    end

    # --- 1D float PFT variables ---
    # Mapping: NetCDF varname => pftcon field
    float_vars = [
        ("z0mr",           p.z0mr),
        ("displar",        p.displar),
        ("dleaf",          p.dleaf),
        ("c3psn",          p.c3psn),
        ("xl",             p.xl),
        ("roota_par",      p.roota_par),
        ("rootb_par",      p.rootb_par),
        ("slatop",         p.slatop),
        ("dsladlai",       p.dsladlai),
        ("leafcn",         p.leafcn),
        ("biofuel_harvfrac", p.biofuel_harvfrac),
        ("flnr",           p.flnr),
        ("smpso",          p.smpso),
        ("smpsc",          p.smpsc),
        ("fnitr",          p.fnitr),
        ("woody",          p.woody),
        ("lflitcn",        p.lflitcn),
        ("frootcn",        p.frootcn),
        ("livewdcn",       p.livewdcn),
        ("deadwdcn",       p.deadwdcn),
        ("grperc",         p.grperc),
        ("grpnow",         p.grpnow),
        ("froot_leaf",     p.froot_leaf),
        ("stem_leaf",      p.stem_leaf),
        ("croot_stem",     p.croot_stem),
        ("flivewd",        p.flivewd),
        ("fcur",           p.fcur),
        ("fcurdv",         p.fcurdv),
        ("lf_flab",        p.lf_flab),
        ("lf_fcel",        p.lf_fcel),
        ("lf_flig",        p.lf_flig),
        ("fr_flab",        p.fr_flab),
        ("fr_fcel",        p.fr_fcel),
        ("fr_flig",        p.fr_flig),
        ("leaf_long",      p.leaf_long),
        ("evergreen",      p.evergreen),
        ("stress_decid",   p.stress_decid),
        ("season_decid",   p.season_decid),
        ("season_decid_temperate", p.season_decid_temperate),
        ("crit_onset_gdd_sf", p.crit_onset_gdd_sf),
        ("ndays_on",       p.ndays_on),
        ("pftpar20",       p.pftpar20),
        ("pftpar28",       p.pftpar28),
        ("pftpar29",       p.pftpar29),
        ("pftpar30",       p.pftpar30),
        ("pftpar31",       p.pftpar31),
        ("a_fix",          p.a_fix),
        ("b_fix",          p.b_fix),
        ("c_fix",          p.c_fix),
        ("s_fix",          p.s_fix),
        ("akc_active",     p.akc_active),
        ("akn_active",     p.akn_active),
        ("ekc_active",     p.ekc_active),
        ("ekn_active",     p.ekn_active),
        ("kc_nonmyc",      p.kc_nonmyc),
        ("kn_nonmyc",      p.kn_nonmyc),
        ("kr_resorb",      p.kr_resorb),
        ("perecm",         p.perecm),
        ("fun_cn_flex_a",  p.fun_cn_flex_a),
        ("fun_cn_flex_b",  p.fun_cn_flex_b),
        ("fun_cn_flex_c",  p.fun_cn_flex_c),
        ("FUN_fracfixers", p.FUN_fracfixers),
        ("manunitro",      p.manunitro),
        ("fleafcn",        p.fleafcn),
        ("ffrootcn",       p.ffrootcn),
        ("fstemcn",        p.fstemcn),
        ("pconv",          p.pconv),
        ("pprod10",        p.pprod10),
        ("pprod100",       p.pprod100),
        ("pprodharv10",    p.pprodharv10),
        ("graincn",        p.graincn),
        ("mxtmp",          p.mxtmp),
        ("baset",          p.baset),
        ("declfact",       p.declfact),
        ("bfact",          p.bfact),
        ("aleaff",         p.aleaff),
        ("arootf",         p.arootf),
        ("astemf",         p.astemf),
        ("arooti",         p.arooti),
        ("fleafi",         p.fleafi),
        ("allconsl",       p.allconsl),
        ("allconss",       p.allconss),
        ("crop",           p.crop),
        ("irrigated",      p.irrigated),
        ("ztopmx",         p.ztopmx),
        ("laimx",          p.laimx),
        ("gddmin",         p.gddmin),
        ("hybgdd",         p.hybgdd),
        ("lfemerg",        p.lfemerg),
        ("grnfill",        p.grnfill),
        ("mbbopt",         p.mbbopt),
        ("medlynslope",    p.medlynslope),
        ("medlynintercept", p.medlynintercept),
        ("cc_leaf",        p.cc_leaf),
        ("cc_lstem",       p.cc_lstem),
        ("cc_dstem",       p.cc_dstem),
        ("cc_other",       p.cc_other),
        ("fm_leaf",        p.fm_leaf),
        ("fm_lstem",       p.fm_lstem),
        ("fm_dstem",       p.fm_dstem),
        ("fm_other",       p.fm_other),
        ("fm_root",        p.fm_root),
        ("fm_lroot",       p.fm_lroot),
        ("fm_droot",       p.fm_droot),
        ("fsr_pft",        p.fsr_pft),
        ("fd_pft",         p.fd_pft),
        ("rswf_min",       p.rswf_min),
        ("rswf_max",       p.rswf_max),
        ("nstem",          p.nstem),
        ("taper",          p.taper),
        # Raupach92 roughness parameters
        ("z0v_Cr",         p.z0v_Cr),
        ("z0v_Cs",         p.z0v_Cs),
        ("z0v_c",          p.z0v_c),
        ("z0v_cw",         p.z0v_cw),
        ("z0v_LAImax",     p.z0v_LAImax),
        ("z0v_LAIoff",     p.z0v_LAIoff),
        # Planting temperature
        ("planting_temp",       p.planttemp),
        ("min_planting_temp",   p.minplanttemp),
        # Flexible CN
        ("i_vcad",         p.i_vcad),
        ("s_vcad",         p.s_vcad),
        ("i_flnr",         p.i_flnr),
        ("s_flnr",         p.s_flnr),
    ]

    for (varname, dest) in float_vars
        _read_pft_var!(ds, varname, dest)
    end

    # --- 1D integer PFT variables ---
    _read_pft_var_int!(ds, "mergetoclmpft", p.mergetoclmpft)
    _read_pft_var_int!(ds, "mxmat", p.mxmat)
    _read_pft_var_int!(ds, "min_NH_planting_date", p.mnNHplantdate)
    _read_pft_var_int!(ds, "max_NH_planting_date", p.mxNHplantdate)
    _read_pft_var_int!(ds, "min_SH_planting_date", p.mnSHplantdate)
    _read_pft_var_int!(ds, "max_SH_planting_date", p.mxSHplantdate)

    # --- 2D optical properties (VIS=1, NIR=2) ---
    IVIS = 1; INIR = 2
    _read_pft_col!(ds, "rholvis", p.rhol, IVIS)
    _read_pft_col!(ds, "rholnir", p.rhol, INIR)
    _read_pft_col!(ds, "rhosvis", p.rhos, IVIS)
    _read_pft_col!(ds, "rhosnir", p.rhos, INIR)
    _read_pft_col!(ds, "taulvis", p.taul, IVIS)
    _read_pft_col!(ds, "taulnir", p.taul, INIR)
    _read_pft_col!(ds, "tausvis", p.taus, IVIS)
    _read_pft_col!(ds, "tausnir", p.taus, INIR)

    # --- 2D rootprof_beta [pft, variant] ---
    if haskey(ds, "rootprof_beta")
        data = Array(ds["rootprof_beta"])
        npft_file = size(data, 1)
        nvar = size(data, 2)
        for v in 1:min(nvar, size(p.rootprof_beta, 2))
            for i in 1:min(npft_file, size(p.rootprof_beta, 1))
                p.rootprof_beta[i, v] = Float64(data[i, v])
            end
        end
    end

    # --- Constants that don't come from file ---
    n = length(p.dwood)
    for i in 1:n
        p.dwood[i] = DWOOD_PARAM
        p.root_radius[i] = ROOT_RADIUS_PARAM
        p.root_density[i] = ROOT_DENSITY_PARAM
    end

    # --- Set vegetation type flags ---
    for i in 1:n
        p.is_tree[i]  = (p.woody[i] > 0.5) && (p.slatop[i] > 0.0)
        p.is_shrub[i] = false  # derived from type index later if needed
        p.is_grass[i] = (p.woody[i] < 0.5) && (p.crop[i] < 0.5) && (i > 1)
    end

    # Set mergetoclmpft-based flags
    set_is_pft_known_to_model!(p)

    nothing
end

"""
    _read_pft_col!(ds, varname, dest_matrix, col_idx)

Read a 1D PFT variable and store it into column `col_idx` of a 2D matrix.
"""
function _read_pft_col!(ds::NCDataset, varname::String, dest::Matrix{Float64}, col_idx::Int)
    if haskey(ds, varname)
        data = Array(ds[varname])
        npft_file = length(data)
        for i in 1:min(npft_file, size(dest, 1))
            val = data[i]
            dest[i, col_idx] = ismissing(val) ? 0.0 : Float64(val)
        end
    end
    nothing
end

# --- initVertical scalar parameters ---
# Module-level storage (replaces Fortran params_inst)
const _initvert_slopebeta  = Ref(0.0)
const _initvert_slopemax   = Ref(0.0)
const _initvert_zbedrock   = Ref(-1.0)   # negative means don't substitute
const _initvert_zbedrock_sf = Ref(1.0)

"""
    readParams_initVertical!(ds::NCDataset)

Read initVertical scalar parameters from the parameter file.
"""
function readParams_initVertical!(ds::NCDataset)
    _initvert_slopebeta[]  = _read_scalar(ds, "slopebeta", 0.0)
    _initvert_slopemax[]   = _read_scalar(ds, "slopemax", 0.0)
    _initvert_zbedrock[]   = _read_scalar(ds, "zbedrock", -1.0)
    _initvert_zbedrock_sf[] = _read_scalar(ds, "zbedrock_sf", 1.0)
    nothing
end

"""
    readParams_photosynthesis!(ds::NCDataset)

Read photosynthesis module parameters from the parameter file.
Initializes the global `params_inst` (PhotoParamsData).
"""
function readParams_photosynthesis!(ds::NCDataset)
    photo_params_init!(params_inst)

    # 1D PFT parameters
    _read_pft_var!(ds, "theta_cj", params_inst.theta_cj)
    _read_pft_var!(ds, "krmax", params_inst.krmax)
    _read_pft_var!(ds, "lmr_intercept_atkin", params_inst.lmr_intercept_atkin)

    # Scalar parameters (kinetics, activation/deactivation energies)
    params_inst.theta_ip     = _read_scalar(ds, "theta_ip", 0.95)
    params_inst.act25        = _read_scalar(ds, "act25", 72.0)
    params_inst.fnr          = _read_scalar(ds, "fnr", 7.16)
    params_inst.cp25_yr2000  = _read_scalar(ds, "cp25_yr2000", 42.75e-6)
    params_inst.kc25_coef    = _read_scalar(ds, "kc25_coef", 404.9e-6)
    params_inst.ko25_coef    = _read_scalar(ds, "ko25_coef", 278.4e-3)
    params_inst.fnps         = _read_scalar(ds, "fnps", 0.15)
    params_inst.theta_psii   = _read_scalar(ds, "theta_psii", 0.7)
    params_inst.vcmaxha      = _read_scalar(ds, "vcmaxha", 65330.0)
    params_inst.jmaxha       = _read_scalar(ds, "jmaxha", 43540.0)
    params_inst.tpuha        = _read_scalar(ds, "tpuha", 53100.0)
    params_inst.lmrha        = _read_scalar(ds, "lmrha", 46390.0)
    params_inst.kcha         = _read_scalar(ds, "kcha", 79430.0)
    params_inst.koha         = _read_scalar(ds, "koha", 36380.0)
    params_inst.cpha         = _read_scalar(ds, "cpha", 37830.0)
    params_inst.vcmaxhd      = _read_scalar(ds, "vcmaxhd", 149250.0)
    params_inst.jmaxhd       = _read_scalar(ds, "jmaxhd", 152040.0)
    params_inst.tpuhd        = _read_scalar(ds, "tpuhd", 150650.0)
    params_inst.lmrhd        = _read_scalar(ds, "lmrhd", 150650.0)
    params_inst.lmrse        = _read_scalar(ds, "lmrse", 490.0)
    params_inst.tpu25ratio   = _read_scalar(ds, "tpu25ratio", 0.167)
    params_inst.kp25ratio    = _read_scalar(ds, "kp25ratio", 20160.0)
    params_inst.vcmaxse_sf   = _read_scalar(ds, "vcmaxse_sf", 1.0)
    params_inst.jmaxse_sf    = _read_scalar(ds, "jmaxse_sf", 1.0)
    params_inst.tpuse_sf     = _read_scalar(ds, "tpuse_sf", 1.0)
    params_inst.jmax25top_sf = _read_scalar(ds, "jmax25top_sf", 1.0)

    # 2D parameters (pft × nvegwcs): kmax, psi50, ck
    if haskey(ds, "kmax")
        data = Array(ds["kmax"])
        data = replace(data, missing => NaN)
        n1 = min(size(data, 1), size(params_inst.kmax, 1))
        n2 = min(size(data, 2), size(params_inst.kmax, 2))
        params_inst.kmax[1:n1, 1:n2] .= Float64.(data[1:n1, 1:n2])
    end
    if haskey(ds, "psi50")
        data = Array(ds["psi50"])
        data = replace(data, missing => NaN)
        n1 = min(size(data, 1), size(params_inst.psi50, 1))
        n2 = min(size(data, 2), size(params_inst.psi50, 2))
        params_inst.psi50[1:n1, 1:n2] .= Float64.(data[1:n1, 1:n2])
    end
    if haskey(ds, "ck")
        data = Array(ds["ck"])
        data = replace(data, missing => NaN)
        n1 = min(size(data, 1), size(params_inst.ck, 1))
        n2 = min(size(data, 2), size(params_inst.ck, 2))
        params_inst.ck[1:n1, 1:n2] .= Float64.(data[1:n1, 1:n2])
    end

    nothing
end
