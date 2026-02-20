# ==========================================================================
# Ported from: src/main/landunit_varcon.F90
# Landunit type constants and utility functions
# ==========================================================================

# Landunit type indices
const ISTSOIL     = 1   # vegetated or bare soil
const ISTCROP     = 2   # crop
const ISTOCN      = 3   # ocean (unused in land model)
const ISTICE      = 4   # land ice (glacier)
const ISTDLAK     = 5   # deep lake
const ISTWET      = 6   # wetland
const ISTURB_MIN  = 7   # minimum urban type index
const ISTURB_TBD  = 7   # urban tall building district
const ISTURB_HD   = 8   # urban high density
const ISTURB_MD   = 9   # urban medium density
const ISTURB_MAX  = 9   # maximum urban type index
const MAX_LUNIT   = 9   # maximum number of landunit types
const NUMURBL     = ISTURB_MAX - ISTURB_MIN + 1  # number of urban landunit types

const LANDUNIT_NAMES = [
    "vegetated_or_bare_soil",  # 1
    "crop",                     # 2
    "ocean",                    # 3
    "landice",                  # 4
    "deep_lake",                # 5
    "wetland",                  # 6
    "urban_tbd",                # 7
    "urban_hd",                 # 8
    "urban_md",                 # 9
]

"""
    landunit_is_special(ltype::Int) -> Bool

Returns true if the landunit type is "special" (not vegetated soil or crop).
Special landunits: ocean, ice, lake, wetland, urban.
"""
function landunit_is_special(ltype::Int)
    @assert 1 <= ltype <= MAX_LUNIT "Invalid landunit type: $ltype"
    return ltype != ISTSOIL && ltype != ISTCROP
end
