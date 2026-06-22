# FatesRadiationMemMod.jl
# Julia port of FATES src/fates/radiation/FatesRadiationMemMod.F90
#
# Radiation memory: data that defines how FATES uses its radiation schemes —
# solver flags, radiation-stream indices, shortwave-band (vis/nir) indices, and
# small albedo/transmittance constants. Deps: FatesConstantsMod (fates_r8 ->
# Float64).

# Radiation solver flags
const norman_solver = 1
const twostr_solver = 2

# Radiation streams (direct/diffuse)
const num_rad_stream_types = 2  # number of radiation streams used
const idirect = 1               # array index for direct radiation
const idiffuse = 2              # array index for diffuse radiation

# Shortwave bands. This needs to match what is used in the host model.
# This is visible (1) and near-infrared (2).
const num_swb = 2  # number of shortwave bands we use

const ivis = 1     # array index for visible-spectrum shortwave radiation
const inir = 2     # array index for near-infrared-spectrum shortwave radiation
const ipar = ivis  # photosynthetically active band ~= visible band

# albedo land ice by waveband (1=vis, 2=nir) — module variables (mutable Refs in
# Fortran); kept as plain vectors here.
const alb_ice = [0.80, 0.55]   # albedo, land ice
const rho_snow = [0.80, 0.55]  # snow reflectance
const tau_snow = [0.01, 0.01]  # snow transmittance
