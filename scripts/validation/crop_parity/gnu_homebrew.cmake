string(APPEND LDFLAGS " -framework Accelerate -Wl,-rpath $(NETCDF)/lib")
# ParallelIO bundled with CTSM 5.3.012 expects the legacy netCDF-C
# attribute-name macro, which newer compatible netCDF releases no longer
# expose. Supply the equivalent literal at compile time without changing
# the frozen upstream source.
string(APPEND CPPDEFS " -DNETCDF_ENABLE_LEGACY_MACROS")
