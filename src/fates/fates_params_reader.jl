# fates_params_reader.jl
# Real FATES parameter-file reader — replaces the synthetic in-memory param
# tables (`_fates_spike_setup_pft!`) with a reader of the official FATES default
# parameter file `data/fates/fates_params_default.cdl`.
#
# WHY a CDL text parser (vs ncgen+NCDatasets):
#   The source of truth is the CDL (text) form of the FATES parameter NetCDF
#   shipped with CTSM (src/fates/parameter_files/fates_params_default.cdl). We
#   copy it verbatim into the repo at `data/fates/fates_params_default.cdl` and
#   parse it directly. A self-contained CDL parser keeps this CI-portable: no
#   `ncgen` binary and no NetCDF file to track. NCDatasets is a dependency, but
#   reading the text avoids a build-time `ncgen` step and a committed binary.
#
# HOW it maps onto the existing port:
#   The Julia port ALREADY mirrors the Fortran name->field mapping exactly via the
#   register/receive infrastructure (FatesParametersInterface.jl):
#     PRTRegisterParams!/PRTReceiveParams!   (PRTParamsFATESMod.jl  -> prt_params)
#     FatesRegisterParams!/FatesReceiveParams! (EDParamsMod.jl      -> EDParams[])
#     Register!/Receive!  (EDPftvarcon.jl     -> EDPftvarcon_inst[])
#     SpitFireRegisterParams!/SpitFireReceiveParams! (SFParamsMod.jl -> SFParams[])
#   This reader therefore does what the Fortran `FatesReadParameters` /
#   `FatesUnitTestParamReaderMod` host glue does: it builds a
#   `fates_parameters_type`, calls every module's Register routine, fills the
#   registry's `data` arrays FROM THE FILE (matching by `name=` and dims), then
#   calls every module's Receive routine. The name->field mapping lives (and stays
#   tested) in the existing modules; this file only bridges the CDL file to the
#   registry. Derived params (FatesParameterDerivedMod) and dim Refs are then set.
#
# DIMENSION ORDER:
#   NetCDF/CDL declares multi-dim vars in C row-major order with the FASTEST-
#   varying dim LAST, e.g. `double fates_stoich_nitr(fates_plant_organs,
#   fates_pft)` lists all PFTs for organ 1, then all PFTs for organ 2, ...  Fortran
#   (column-major) reads that into a `(npft, norgan)` array — PFT first. The Julia
#   struct fields are likewise `[pft, organ]`. So we reshape the flat value stream
#   (PFT-fastest) COLUMN-MAJOR into `(dims reversed)` = `(npft, norgan)`, which is
#   exactly the registry/Julia `[pft, organ]` layout. The registry's
#   `dimension_sizes` are written in that same (Julia) order: dim1=pft, dim2=organ.

# Default repo path to the committed CDL parameter file.
const FATES_PARAMS_DEFAULT_CDL = joinpath(@__DIR__, "..", "..", "data", "fates",
                                          "fates_params_default.cdl")

# -------------------------------------------------------------------------------------
# CDL parser
# -------------------------------------------------------------------------------------

"""
    _CDLVar

One parsed CDL variable: its element type tag, its declared dimension names
(in CDL/file order — fastest-varying LAST), and either its flat numeric values
(`values`, file order) or its character-string rows (`strings`).
"""
struct _CDLVar
    name::String
    typ::Symbol                 # :double, :float, :int, :char
    dim_names::Vector{String}   # CDL declaration order (fastest-varying last)
    values::Vector{Float64}     # numeric vars: flat stream in file order
    strings::Vector{String}     # char vars: one string per leading-dim row
end

"""
    parse_cdl(path) -> (dims::Dict{String,Int}, vars::Dict{String,_CDLVar})

Parse a NetCDF CDL (text) file into a dimensions map and a name->`_CDLVar` map.
Self-contained: handles the `dimensions:`, `variables:` and `data:` sections,
multi-line comma/semicolon-delimited numeric data, and quoted char-array rows.
"""
function parse_cdl(path::AbstractString)
    isfile(path) || error("parse_cdl: file not found: $path")
    text = read(path, String)
    lines = split(text, '\n')

    dims = Dict{String,Int}()
    # variable declaration metadata: name -> (typ, dim_names)
    decl_typ = Dict{String,Symbol}()
    decl_dims = Dict{String,Vector{String}}()
    vars = Dict{String,_CDLVar}()

    section = :none   # :dimensions, :variables, :data
    i = 1
    n = length(lines)
    while i <= n
        raw = lines[i]
        line = strip(raw)
        i += 1

        if isempty(line)
            continue
        end
        # Section headers.
        if line == "dimensions:"
            section = :dimensions; continue
        elseif line == "variables:"
            section = :variables; continue
        elseif startswith(line, "data:")
            section = :data; continue
        elseif startswith(line, "//") || startswith(line, "netcdf ") || line == "{" || line == "}"
            continue
        end

        if section == :dimensions
            # e.g.  fates_pft = 14 ;
            m = match(r"^(\w+)\s*=\s*([0-9]+)\s*;", line)
            if m !== nothing
                dims[m.captures[1]] = parse(Int, m.captures[2])
            end

        elseif section == :variables
            # Variable declarations:  double fates_xxx(dimA, dimB) ;   or scalar:
            #   double fates_xxx ;
            # Attribute lines (fates_xxx:units = "...") are ignored.
            m = match(r"^(double|float|int|char)\s+(\w+)\s*(\(([^)]*)\))?\s*;", line)
            if m !== nothing
                typ = Symbol(m.captures[1])
                name = m.captures[2]
                dimstr = m.captures[4]
                dim_names = String[]
                if dimstr !== nothing
                    for d in split(dimstr, ',')
                        ds = strip(d)
                        isempty(ds) || push!(dim_names, String(ds))
                    end
                end
                decl_typ[name] = typ
                decl_dims[name] = dim_names
            end

        elseif section == :data
            # A data assignment may span several lines and ends at ';'.
            m = match(r"^(\w+)\s*=(.*)$", line)
            m === nothing && continue
            name = m.captures[1]
            buf = String(m.captures[2])
            while !occursin(';', buf) && i <= n
                buf *= "\n" * lines[i]
                i += 1
            end
            # Strip trailing ';'
            scolon = findlast(==(';'), buf)
            body = scolon === nothing ? buf : buf[1:prevind(buf, scolon)]

            typ = get(decl_typ, name, :double)
            dim_names = get(decl_dims, name, String[])

            if typ == :char
                strs = _parse_cdl_strings(body)
                vars[name] = _CDLVar(name, typ, dim_names, Float64[], strs)
            else
                vals = _parse_cdl_numbers(body)
                vars[name] = _CDLVar(name, typ, dim_names, vals, String[])
            end
        end
    end

    return dims, vars
end

# Parse the comma-separated numeric body of a data assignment into Float64s.
function _parse_cdl_numbers(body::AbstractString)
    out = Float64[]
    for tok in split(body, ',')
        t = strip(tok)
        isempty(t) && continue
        # CDL uses 'f'/'d'/'L' type suffixes and may use _ as a fill placeholder.
        t = replace(t, r"[fFdDlL]$" => "")
        if t == "_"
            push!(out, NaN)
        else
            push!(out, parse(Float64, t))
        end
    end
    return out
end

# Parse quoted string rows from a char-var data body (one entry per leading dim).
function _parse_cdl_strings(body::AbstractString)
    out = String[]
    for m in eachmatch(r"\"((?:[^\"\\]|\\.)*)\"", body)
        push!(out, String(m.captures[1]))
    end
    return out
end

# -------------------------------------------------------------------------------------
# Fill the parameter registry from the parsed CDL, honoring the registered dims.
# -------------------------------------------------------------------------------------

"""
    _set_registry_from_cdl!(fates_params, dims, vars) -> (n_set, missing_names)

For each parameter registered on `fates_params`, look up its values in the parsed
CDL `vars`, set its `dimension_sizes` (Julia order: as registered), and store the
data via `SetDataScalar!/SetData1D!/SetData2D!`. Returns the number of params set
and the list of registered names with no CDL variable.

Reshape rule: the CDL flat value stream lists the registry's dim1 fastest (the
file declares dims reversed), so a column-major reshape into the registered
`(size_dim1, size_dim2)` reproduces the Fortran/Julia `[dim1, dim2]` layout.
"""
function _set_registry_from_cdl!(fates_params::fates_parameters_type,
                                 dims::Dict{String,Int}, vars::Dict{String,_CDLVar})
    n_set = 0
    missing_names = String[]

    for idx in 1:fates_params.num_parameters
        p = fates_params.parameters[idx]
        name = strip(p.name)
        if !haskey(vars, name)
            push!(missing_names, String(name))
            continue
        end
        cv = vars[name]

        if p.dimension_shape == dimension_shape_scalar
            length(cv.values) >= 1 ||
                error("_set_registry_from_cdl!: scalar $name has no value")
            SetDataScalar!(fates_params, idx, cv.values[1])
            n_set += 1

        elseif p.dimension_shape == dimension_shape_1d
            # 1-D: size from the registered dimension name.
            dn = p.dimension_names[1]
            sz = get(dims, dn, length(cv.values))
            p.dimension_sizes[1] = sz
            length(cv.values) == sz ||
                error("_set_registry_from_cdl!: 1d $name expected $sz got $(length(cv.values))")
            SetData1D!(fates_params, idx, cv.values)
            n_set += 1

        elseif p.dimension_shape == dimension_shape_2d
            # Registered dim order is the Julia/Fortran order (dim1, dim2). The
            # CDL declares them reversed; flat stream has dim1 fastest.
            dn1 = p.dimension_names[1]
            dn2 = p.dimension_names[2]
            s1 = get(dims, dn1, 0)
            s2 = get(dims, dn2, 0)
            (s1 > 0 && s2 > 0) ||
                error("_set_registry_from_cdl!: 2d $name unknown dims $dn1/$dn2")
            p.dimension_sizes[1] = s1
            p.dimension_sizes[2] = s2
            length(cv.values) == s1 * s2 ||
                error("_set_registry_from_cdl!: 2d $name expected $(s1*s2) got $(length(cv.values))")
            # Column-major reshape: dim1 fastest => (s1, s2) == [dim1, dim2].
            mat = reshape(cv.values, s1, s2)
            SetData2D!(fates_params, idx, mat)
            n_set += 1
        end
    end

    return n_set, missing_names
end

# -------------------------------------------------------------------------------------
# Top-level reader
# -------------------------------------------------------------------------------------

"""
    read_fates_params!(; path = FATES_PARAMS_DEFAULT_CDL, verbose = false) -> summary

Read the FATES default parameter file and populate ALL FATES parameter globals
(`prt_params`, `EDParams[]`, `EDPftvarcon_inst[]`, `SFParams[]`, `ParamDerived[]`)
plus the dimension Refs (`numpft[]`, `nlevsclass[]`, `nlevcoage[]`, `nlevage[]`,
`nlevdamage[]`). Mirrors the Fortran host `FatesReadParameters` sequence:
register all parameters, read values from the file into the registry, then have
each module receive them.

This is the physically-meaningful replacement for `_fates_spike_setup_pft!`.

Returns a NamedTuple summary:
  * `npft`              — fates_pft dimension from the file
  * `n_registered`      — number of registered FATES parameters
  * `n_populated`       — registered params filled from the file
  * `n_cdl_vars`        — total `fates_*` variables in the file
  * `n_cdl_unmapped`    — CDL vars with no registered (Julia) consumer
  * `cdl_unmapped`      — their names (sorted)
  * `registered_missing`— registered params absent from the file (should be empty)

`hlm_parteh_mode[]` must already be set (carbon vs CNP) before calling, since the
organ-id derived map (`PRTDerivedParams!`) and the consistency checks branch on
it. `clm_fates_init!` sets it via `_fates_spike_set_ctrlparms!` before this call.
"""
function read_fates_params!(; path::AbstractString = FATES_PARAMS_DEFAULT_CDL,
                            verbose::Bool = false)
    dims, vars = parse_cdl(path)

    haskey(dims, "fates_pft") ||
        error("read_fates_params!: file is missing the fates_pft dimension")
    npft_file = dims["fates_pft"]

    # Two separate registry passes, mirroring the Fortran host exactly. The FATES
    # registry capacity (max_params=250) is per-instance, and the full FATES param
    # set (~300) does NOT fit in one registry — so the Fortran splits it:
    #   pass 1 (FatesReadParameters,  FatesInterfaceMod.F90:2508): ED + SPITFIRE
    #           + PRT + synchronized;
    #   pass 2 (FatesReadPFTs, clmfates_paraminterfaceMod.F90:70): EDPftvarcon.
    n_registered = 0
    n_populated = 0
    registered_missing = String[]
    registered_names = String[]

    # ---- Pass 1: ED + SPITFIRE + PRT + synchronized ----
    fp1 = fates_parameters_type()
    Init!(fp1)
    FatesRegisterParams!(fp1)                          # EDParamsMod
    SpitFireRegisterParams!(fp1)                       # SPITFIRE
    PRTRegisterParams!(fp1)                            # PARTEH / prt_params
    RegisterParams!(FatesSynchronizedParamsInst, fp1) # synchronized (no-op)
    for k in 1:fp1.num_parameters
        push!(registered_names, String(strip(fp1.parameters[k].name)))
    end
    np1, miss1 = _set_registry_from_cdl!(fp1, dims, vars)
    FatesReceiveParams!(fp1)
    SpitFireReceiveParams!(fp1)
    PRTReceiveParams!(fp1)
    ReceiveParams!(FatesSynchronizedParamsInst, fp1)
    n_registered += fp1.num_parameters
    n_populated  += np1
    append!(registered_missing, miss1)

    # ---- Pass 2: EDPftvarcon (its own registry, as in the host) ----
    EDpftconInit!(EDPftvarcon_inst[])
    fp2 = fates_parameters_type()
    Init!(fp2)
    Register!(EDPftvarcon_inst[], fp2)
    for k in 1:fp2.num_parameters
        push!(registered_names, String(strip(fp2.parameters[k].name)))
    end
    np2, miss2 = _set_registry_from_cdl!(fp2, dims, vars)
    Receive!(EDPftvarcon_inst[], fp2)
    n_registered += fp2.num_parameters
    n_populated  += np2
    append!(registered_missing, miss2)

    # ---- 4. Dimension Refs the cold-start chain reads ----
    numpft[]     = npft_file
    nleafage[]   = size(prt_params.leaf_long, 2)
    nlevsclass[] = length(EDParams[].ED_val_history_sizeclass_bin_edges)
    nlevage[]    = length(EDParams[].ED_val_history_ageclass_bin_edges)
    nlevcoage[]  = length(EDParams[].ED_val_history_coageclass_bin_edges)
    nlevdamage[] = length(EDParams[].ED_val_history_damage_bin_edges)
    nlevheight[] = length(EDParams[].ED_val_history_height_bin_edges)

    # VAI radiative-transfer bin edges. In Fortran these are computed inside
    # SetFatesGlobalElements2 (FatesInterfaceMod.F90:929) from the just-read
    # vai_top_bin_width / vai_width_increase_factor; the param file does not store
    # dinc_vai/dlower_vai directly. Compute them here so EDParams is fully usable
    # by the cold-start allometry (tree_sai reads sum(dinc_vai)).
    edp = EDParams[]
    edp.dinc_vai = [edp.ED_val_vai_top_bin_width *
                    edp.ED_val_vai_width_increase_factor^(i - 1) for i in 1:nlevleaf]
    edp.dlower_vai = [sum(@view edp.dinc_vai[1:i]) for i in 1:nlevleaf]

    # ---- 5. Derived parameters (jmax25/tpu25/kp25/branch_frac/damage matrix) ----
    # PRTDerivedParams! sets the global-organ -> param-file-organ reverse lookup
    # that the maintenance-respiration / mass-init paths read.
    PRTDerivedParams!()
    pd = param_derived_type()
    Init!(pd, npft_file)
    ParamDerived[] = pd

    # ---- 6. Coverage report: which CDL fates_* vars have no Julia consumer ----
    registered_set = Set(registered_names)
    cdl_unmapped = String[]
    n_cdl_vars = 0
    for (nm, _) in vars
        startswith(nm, "fates_") || continue
        n_cdl_vars += 1
        nm in registered_set || push!(cdl_unmapped, nm)
    end
    sort!(cdl_unmapped)
    sort!(registered_missing)

    summary = (npft = npft_file,
               n_registered = n_registered,
               n_populated = n_populated,
               n_cdl_vars = n_cdl_vars,
               n_cdl_unmapped = length(cdl_unmapped),
               cdl_unmapped = cdl_unmapped,
               registered_missing = registered_missing)

    if verbose
        @info "read_fates_params!" npft = summary.npft n_registered = summary.n_registered n_populated = summary.n_populated n_cdl_vars = summary.n_cdl_vars n_cdl_unmapped = summary.n_cdl_unmapped
        isempty(registered_missing) ||
            @warn "read_fates_params!: registered params absent from file" registered_missing
        isempty(cdl_unmapped) ||
            @info "read_fates_params!: CDL vars with no Julia consumer (informational)" cdl_unmapped
    end

    return summary
end
