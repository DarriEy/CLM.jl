# ==========================================================================
# Synthetic forcing file generator for CLM.jl integration tests
#
# Creates a minimal NetCDF forcing file for the Bow at Banff test site
# with realistic winter conditions.
# ==========================================================================

using NCDatasets
using Dates

"""
    generate_forcing(filepath; start_date, end_date, dtime, lat, lon)

Generate a synthetic atmospheric forcing file with realistic diurnal cycles.

# Arguments
- `filepath::String` — Output NetCDF file path
- `start_date::DateTime` — Start of forcing period
- `end_date::DateTime` — End of forcing period
- `dtime::Int` — Timestep in seconds (default 1800)
- `lat::Float64` — Latitude (default 51.17 for Bow at Banff)
- `lon::Float64` — Longitude (default -115.57 for Bow at Banff)
"""
function generate_forcing(filepath::String;
    start_date::DateTime = DateTime(2000, 1, 1),
    end_date::DateTime = DateTime(2000, 1, 2),
    dtime::Int = 1800,
    lat::Float64 = 51.17,
    lon::Float64 = -115.57)

    # Generate time steps
    times = collect(start_date:Second(dtime):end_date - Second(1))
    ntimes = length(times)

    # Allocate arrays
    tbot = zeros(ntimes)   # Temperature [K]
    psrf = zeros(ntimes)   # Surface pressure [Pa]
    wind = zeros(ntimes)   # Wind speed [m/s]
    flds = zeros(ntimes)   # Longwave radiation [W/m2]
    fsds = zeros(ntimes)   # Shortwave radiation [W/m2]
    prec = zeros(ntimes)   # Precipitation [mm/s]
    qbot = zeros(ntimes)   # Specific humidity [kg/kg]

    for (i, t) in enumerate(times)
        # Hour of day as fraction (0-1)
        hod = (hour(t) * 3600 + minute(t) * 60 + second(t)) / 86400.0

        # Day of year
        doy = dayofyear(t)

        # ---- Temperature: winter diurnal cycle ----
        # Base: 270 K (-3°C) with ±3K diurnal range
        tbot[i] = 270.0 + 3.0 * sin(2π * (hod - 0.25))  # min at 6am, max at 6pm

        # ---- Pressure: ~85 kPa (elevation ~1400m) ----
        psrf[i] = 85000.0 + 100.0 * sin(2π * doy / 365.0)  # slight seasonal variation

        # ---- Wind: 3 m/s base with gustiness ----
        wind[i] = 3.0 + 1.0 * sin(2π * (hod + 0.1))

        # ---- Longwave: ~250 W/m2 with diurnal variation ----
        flds[i] = 250.0 + 20.0 * sin(2π * (hod - 0.25))

        # ---- Shortwave: cosine zenith angle, peak ~400 W/m2 ----
        # Simple approximation: solar noon at ~12:00
        lat_rad = lat * π / 180.0
        declin = -23.44 * π / 180.0 * cos(2π * (doy + 10) / 365.0)
        hour_angle = 2π * (hod - 0.5)
        cosz = sin(lat_rad) * sin(declin) + cos(lat_rad) * cos(declin) * cos(hour_angle)
        cosz = max(cosz, 0.0)
        fsds[i] = 1360.0 * 0.75 * cosz  # top-of-atm * transmissivity

        # ---- Precipitation: light snowfall in early hours ----
        if hod > 0.0 && hod < 0.25  # 0:00 - 6:00
            prec[i] = 0.0001  # 0.1 mm/hr = 2.78e-5 mm/s
        end

        # ---- Specific humidity: low for cold winter ----
        # ~3 g/kg = 0.003 kg/kg
        qbot[i] = 0.003 + 0.001 * sin(2π * (hod - 0.25))
    end

    # ---- Write NetCDF ----
    ds = NCDataset(filepath, "c")

    defDim(ds, "time", Inf)  # unlimited

    # Time variable
    defVar(ds, "time", Float64, ("time",);
           attrib = Dict("units" => "days since 2000-01-01",
                         "calendar" => "noleap"))

    # Forcing variables (1D: time only, single gridcell)
    defVar(ds, "TBOT", Float64, ("time",); attrib = Dict("units" => "K", "long_name" => "air temperature"))
    defVar(ds, "PSRF", Float64, ("time",); attrib = Dict("units" => "Pa", "long_name" => "surface pressure"))
    defVar(ds, "WIND", Float64, ("time",); attrib = Dict("units" => "m/s", "long_name" => "wind speed"))
    defVar(ds, "FLDS", Float64, ("time",); attrib = Dict("units" => "W/m2", "long_name" => "longwave radiation"))
    defVar(ds, "FSDS", Float64, ("time",); attrib = Dict("units" => "W/m2", "long_name" => "shortwave radiation"))
    defVar(ds, "PRECTmms", Float64, ("time",); attrib = Dict("units" => "mm/s", "long_name" => "total precipitation"))
    defVar(ds, "QBOT", Float64, ("time",); attrib = Dict("units" => "kg/kg", "long_name" => "specific humidity"))

    # Write data
    ref_date = DateTime(2000, 1, 1)
    for i in 1:ntimes
        ds["time"][i] = Dates.value(times[i] - ref_date) / (1000 * 86400)  # ms to days
        ds["TBOT"][i] = tbot[i]
        ds["PSRF"][i] = psrf[i]
        ds["WIND"][i] = wind[i]
        ds["FLDS"][i] = flds[i]
        ds["FSDS"][i] = fsds[i]
        ds["PRECTmms"][i] = prec[i]
        ds["QBOT"][i] = qbot[i]
    end

    # Global attributes
    ds.attrib["title"] = "Synthetic forcing for CLM.jl integration test"
    ds.attrib["site"] = "Bow at Banff"
    ds.attrib["latitude"] = lat
    ds.attrib["longitude"] = lon

    close(ds)

    return filepath
end
