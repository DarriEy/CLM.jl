# make_luh_1pt.jl — single-point subsets of the CESM inputdata LUH2 FATES files.
#
# Picks the 4x5 cell nearest the site (Bow at Banff surfdata LATIXY/LONGXY; the
# nearest cell is 50N 245E) and writes 1x1 copies that CTSM's dyn_file / ncd_io
# (grlnd) and the port's reader both read on a single-point grid:
#
#   LUH2 timeseries : YEAR, time, lon, lat, LONGXY, LATIXY + the 12 states, 108
#                     transitions and 10 harvest variables (time, lsmlat=1, lsmlon=1)
#   landuse x pft   : frac_primr/secnd/pastr/range/csurf (natpft, lat=1, lon=1),
#                     frac_brgnd, mask, natpft, lon, lat
#
# Source files (fetched from https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/
# lnd/clm2/surfdata_map/, 339 MB + 1.9 MB; the paths namelist_defaults_ctsm.xml gives):
#   fates-sci.1.68.3_api.31.0.0_tools.1.0.1/LUH2_states_transitions_management.timeseries_4x5_hist_simyr1850-2015_c231101.nc
#   fates-sci.1.77.0_api.36.0.0/fates_landuse_pft_map_4x5_240206.nc
#
# usage: julia --project=. make_luh_1pt.jl <luh_4x5.nc> <lupft_4x5.nc> <surfdata.nc> <outdir>
using NCDatasets, CLM
using CLM: landuse_state_varnames, landuse_transition_varnames, landuse_harvest_area_varnames, landuse_harvest_mass_varnames

function nearest(lon, lat, tlon, tlat)
    i = argmin(abs.(mod.(lon .- mod(tlon, 360) .+ 180, 360) .- 180))
    j = argmin(abs.(lat .- tlat))
    return i, j
end

function main(luh_src, lupft_src, surf, outdir)
    mkpath(outdir)
    tlon, tlat = NCDataset(surf) do s; (Float64(s["LONGXY"][1, 1]), Float64(s["LATIXY"][1, 1])); end
    luh_out = joinpath(outdir, "LUH2_1pt_Bow_hist_simyr1850-2015.nc")
    pft_out = joinpath(outdir, "fates_landuse_pft_map_1pt_Bow.nc")
    keep = vcat(landuse_state_varnames, landuse_transition_varnames, landuse_harvest_area_varnames,
                landuse_harvest_mass_varnames)
    NCDataset(luh_src) do d
        i, j = nearest(d["lon"][:], d["lat"][:], tlon, tlat)
        println("LUH2: site ($tlat, $tlon) -> cell i=$i j=$j ($(d["lat"][j]), $(d["lon"][i]))")
        nt = d.dim["time"]
        isfile(luh_out) && rm(luh_out)
        NCDataset(luh_out, "c") do o
            defDim(o, "time", nt); defDim(o, "lsmlon", 1); defDim(o, "lsmlat", 1)
            o.attrib["source"] = basename(luh_src)
            o.attrib["subset"] = "nearest 4x5 cell i=$i j=$j to the Bow at Banff site"
            defVar(o, "time", Float64.(d["time"].var[:]), ("time",); attrib = Dict(
                "units" => d["time"].attrib["units"], "calendar" => "noleap"))
            defVar(o, "YEAR", Float64.(d["YEAR"][:]), ("time",))
            defVar(o, "lon", [Float64(d["lon"][i])], ("lsmlon",))
            defVar(o, "lat", [Float64(d["lat"][j])], ("lsmlat",))
            defVar(o, "LONGXY", fill(Float64(d["lon"][i]), 1, 1), ("lsmlon", "lsmlat"))
            defVar(o, "LATIXY", fill(Float64(d["lat"][j]), 1, 1), ("lsmlon", "lsmlat"))
            for v in keep
                haskey(d, v) || error("$v not on $luh_src")
                x = Float32.(coalesce.(d[v][i, j, :], NaN32))
                defVar(o, v, reshape(x, 1, 1, nt), ("lsmlon", "lsmlat", "time");
                       fillvalue = NaN32)
            end
        end
    end
    NCDataset(lupft_src) do d
        i, j = nearest(d["lon"][:], d["lat"][:], tlon, tlat)
        println("landuse x pft: cell i=$i j=$j ($(d["lat"][j]), $(d["lon"][i]))")
        npft = d.dim["natpft"]
        isfile(pft_out) && rm(pft_out)
        NCDataset(pft_out, "c") do o
            defDim(o, "lon", 1); defDim(o, "lat", 1); defDim(o, "natpft", npft)
            o.attrib["source"] = basename(lupft_src)
            o.attrib["subset"] = "nearest 4x5 cell i=$i j=$j to the Bow at Banff site"
            defVar(o, "lon", [Float64(d["lon"][i])], ("lon",))
            defVar(o, "lat", [Float64(d["lat"][j])], ("lat",))
            defVar(o, "natpft", Int32.(d["natpft"][:]), ("natpft",))
            defVar(o, "mask", fill(Int64(d["mask"][i, j]), 1, 1), ("lon", "lat"))
            defVar(o, "frac_brgnd", fill(Float64(d["frac_brgnd"][i, j]), 1, 1), ("lon", "lat");
                   fillvalue = NaN)
            for v in ("frac_primr", "frac_secnd", "frac_pastr", "frac_range", "frac_csurf")
                x = Float64.(coalesce.(d[v][i, j, :], NaN))
                defVar(o, v, reshape(x, 1, 1, npft), ("lon", "lat", "natpft"); fillvalue = NaN)
            end
        end
    end
    println("wrote ", luh_out, "\nwrote ", pft_out)
end

main(ARGS...)
