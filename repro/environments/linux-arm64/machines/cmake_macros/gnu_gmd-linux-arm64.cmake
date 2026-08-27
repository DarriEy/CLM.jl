# ParallelIO bundled with CTSM 5.3.012 uses the netCDF legacy attribute-name macro.
string(APPEND CPPDEFS " -DNETCDF_ENABLE_LEGACY_MACROS")
set(LAPACK_LIBDIR "/opt/conda/lib")
