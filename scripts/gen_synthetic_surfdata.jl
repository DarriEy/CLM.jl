# ==========================================================================
# gen_synthetic_surfdata.jl — FULLY SYNTHETIC single-point surfdata fixtures.
#
# Replaces the tracked byte-copy fixtures (redistribution audit D05-D07) with
# files containing NO upstream data: every value below is a constant, a
# standard pedotransfer formula evaluated on a synthetic soil texture, or an
# analytic seasonal cycle. The coordinates are a generic mid-latitude point
# chosen to be nobody's study site.
#
# Produces:
#   test_inputs/lake/surfdata_lake100.nc    — PCT_LAKE=100
#   test_inputs/lake/surfdata_mixed.nc      — PCT_NATVEG=50, PCT_LAKE=50
#   test_inputs/glacier/surfdata_glacier100.nc — PCT_GLACIER=100 (class-1 glc_mec)
#
# These are STRUCTURAL fixtures: they exercise landunit/subgrid construction,
# cold-start init, and device parity. They are not any real site. Workflows
# that need a real-site variant (Fortran reference parity) regenerate one
# locally with scripts/gen_lake_surfdata.jl / gen_glacier_surfdata.jl from
# SYMFLUENCE_DATA; those outputs are not tracked.
#
# Usage:
#   julia --project=. scripts/gen_synthetic_surfdata.jl
# Set CLM_FIXTURE_OUTDIR to write elsewhere (regeneration check without
# overwriting the tracked fixtures).
# ==========================================================================

using NCDatasets

const REPO = joinpath(@__DIR__, "..")
const OUTBASE = get(ENV, "CLM_FIXTURE_OUTDIR", "")

outpath(rel...) = isempty(OUTBASE) ? joinpath(REPO, rel...) : joinpath(OUTBASE, rel...)

# ---- synthetic soil: uniform loam, CLM5 pedotransfer (Clapp-Hornberger via
# the standard CLM sand/clay fits; see SoilStateInitTimeConstMod) -----------
const SAND = 40.0      # percent
const CLAY = 20.0      # percent
watsat_f(sand) = 0.489 - 0.00126 * sand
bsw_f(clay)    = 2.91 + 0.159 * clay
sucsat_f(sand) = 10.0 * 10.0^(1.88 - 0.0131 * sand)
hksat_f(sand)  = 0.0070556 * 10.0^(-0.884 + 0.0153 * sand)

# ---- synthetic phenology: analytic seasonal cycle ------------------------
# Evergreen-needleleaf-like PFT at natpft index 2 (0-based pft 1) and a C3
# grass at index 14 (0-based pft 13). Peak LAI mid-summer (month 7).
lai_tree(m)  = 2.0 + 1.0 * (1 + cospi(2 * (m - 7) / 12)) / 2      # 2.0..3.0
lai_grass(m) = 0.3 + 1.2 * (1 + cospi(2 * (m - 7) / 12)) / 2      # 0.3..1.5

function write_synthetic_surfdata(dst::AbstractString;
        pct_natveg::Float64, pct_lake::Float64, pct_glacier::Float64,
        lakedepth::Float64 = 10.0)
    mkpath(dirname(dst))
    isfile(dst) && rm(dst)
    NCDataset(dst, "c") do ds
        defDim(ds, "lsmlon", 1); defDim(ds, "lsmlat", 1)
        defDim(ds, "natpft", 15); defDim(ds, "cft", 2); defDim(ds, "lsmpft", 17)
        defDim(ds, "nlevsoi", 10); defDim(ds, "nlevurb", 5)
        defDim(ds, "nglcec", 10); defDim(ds, "nglcecp1", 11); defDim(ds, "time", 12)

        g2(name, val; T = Float64) = begin
            v = defVar(ds, name, T, ("lsmlon", "lsmlat"))
            v[:, :] .= T(val); v
        end
        gsoil(name, vals) = begin
            v = defVar(ds, name, Float64, ("nlevsoi", "lsmlon", "lsmlat"))
            v[:, 1, 1] = vals; v
        end

        # Generic synthetic mid-latitude point — not a study site.
        g2("LATIXY", 45.0); g2("LONGXY", 262.5)
        g2("LANDFRAC_PFT", 1.0); g2("PFTDATA_MASK", 1; T = Int32)

        g2("PCT_NATVEG", pct_natveg); g2("PCT_CROP", 0.0)
        g2("PCT_LAKE", pct_lake); g2("PCT_WETLAND", 0.0)
        g2("PCT_GLACIER", pct_glacier); g2("PCT_OCEAN", 0.0); g2("PCT_URBAN", 0.0)

        pnat = defVar(ds, "PCT_NAT_PFT", Float64, ("natpft", "lsmlon", "lsmlat"))
        pnat[:, 1, 1] = [20.0, 40.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 40.0, 0]
        pcft = defVar(ds, "PCT_CFT", Float64, ("cft", "lsmlon", "lsmlat"))
        pcft[:, 1, 1] = [100.0, 0.0]

        g2("FMAX", 0.38); g2("TOPO", 1000.0); g2("STD_ELEV", 100.0); g2("SLOPE", 0.05)

        gsoil("PCT_SAND", fill(SAND, 10)); gsoil("PCT_CLAY", fill(CLAY, 10))
        gsoil("watsat", fill(watsat_f(SAND), 10))
        gsoil("hksat", fill(hksat_f(SAND), 10))
        gsoil("sucsat", fill(sucsat_f(SAND), 10))
        gsoil("bsw", fill(bsw_f(CLAY), 10))
        gsoil("ORGANIC", [5.0, 4.0, 3.0, 2.0, 1.0, 0.5, 0, 0, 0, 0])
        g2("SOIL_COLOR", 15; T = Int32); g2("zbedrock", 8.0)

        glc = defVar(ds, "PCT_GLC_MEC", Float64, ("nglcecp1", "lsmlon", "lsmlat"))
        glc[:, 1, 1] = [100.0; zeros(10)]      # all glacier in class 0 (lowest)
        tglc = defVar(ds, "TOPO_GLC_MEC", Float64, ("nglcecp1", "lsmlon", "lsmlat"))
        tglc[:, 1, 1] = fill(1000.0, 11)
        g2("GLACIER_REGION", 0; T = Int32)

        mx = defVar(ds, "mxsoil_color", Int32, ())
        mx[] = Int32(20)
        g2("gdp", 10.0); g2("peatf", 0.0); g2("abm", 6; T = Int32)
        g2("LAKEDEPTH", lakedepth)
        for ef in ("EF1_BTR", "EF1_FET", "EF1_FDT", "EF1_SHR", "EF1_GRS", "EF1_CRP")
            g2(ef, 0.0)
        end

        for (name, tree_v, grass_v) in (
                ("MONTHLY_LAI", lai_tree, lai_grass),
                ("MONTHLY_SAI", m -> 0.2 + 0.1 * lai_tree(m), m -> 0.05 + 0.1 * lai_grass(m)),
                ("MONTHLY_HEIGHT_TOP", m -> 17.0, m -> 0.5),
                ("MONTHLY_HEIGHT_BOT", m -> 8.5, m -> 0.05))
            v = defVar(ds, name, Float64, ("lsmlon", "lsmlat", "lsmpft", "time"))
            v[:, :, :, :] .= 0.0
            for m in 1:12
                v[1, 1, 2, m] = tree_v(m)      # lsmpft slot 2 = natpft index 2
                v[1, 1, 14, m] = grass_v(m)    # lsmpft slot 14 = natpft index 14
            end
        end

        ds.attrib["Dataset_Version"] = "5.3"
        ds.attrib["title"] = "Fully synthetic single-point CLM surface data (structural test fixture)"
        ds.attrib["source"] = "scripts/gen_synthetic_surfdata.jl — constants, standard pedotransfer formulas, and analytic phenology only; no external data"
        ds.attrib["synthetic"] = "true"
    end
    return dst
end

function main()
    f1 = write_synthetic_surfdata(outpath("test_inputs", "lake", "surfdata_lake100.nc");
        pct_natveg = 0.0, pct_lake = 100.0, pct_glacier = 0.0)
    f2 = write_synthetic_surfdata(outpath("test_inputs", "lake", "surfdata_mixed.nc");
        pct_natveg = 50.0, pct_lake = 50.0, pct_glacier = 0.0)
    f3 = write_synthetic_surfdata(outpath("test_inputs", "glacier", "surfdata_glacier100.nc");
        pct_natveg = 0.0, pct_lake = 0.0, pct_glacier = 100.0)
    foreach(f -> println("wrote ", f), (f1, f2, f3))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
