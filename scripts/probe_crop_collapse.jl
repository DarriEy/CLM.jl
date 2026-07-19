# Probe: does the STATIC surfdata path collapse the 64-CFT set the way CTSM does?
#
# CTSM calls collapse_crop_types TWICE — surfrdMod.F90:991 (static surfdata read)
# and dyncropFileMod.F90:176 (transient landuse). CLM.jl only ever called it from
# the transient path (dyn_pft_crop_file.jl:347), so a STATIC crop run kept all 22
# nonzero CFTs of the US-plains surfdata instead of collapsing to Fortran's 8.
#
# Run:
#   julia +1.12 --project=. scripts/probe_crop_collapse.jl
#
# Prints the built crop patch (itype, weight) vector next to the Fortran reference.

using CLM
using Printf

const SURFDATA = joinpath(
    homedir(),
    "Library/CloudStorage/GoogleDrive-dareyt@gmail.com/My Drive/data",
    "SYMFLUENCE_data/crop_cft_surfdata/surfdata_cropCFT_USplains_1pt.nc",
)

# The 8 crop CFT patches the Fortran reference builds, from the collapsed
# mergetoclmpft mapping with irrigate=.true. (so the irrigated CFTs survive).
const FORTRAN_ITYPES = [17, 18, 19, 20, 23, 24, 75, 76]

function main()
    isfile(SURFDATA) || error("crop surfdata not found: $SURFDATA")

    CLM.varctl.use_crop = true
    CLM.varctl.create_crop_landunit = true
    CLM.varctl.irrigate = true

    (numpft, numcft) = CLM.surfrd_get_num_patches(SURFDATA)
    nlevurb = CLM.surfrd_get_nlevurb(SURFDATA)
    @printf("surfdata dims: natpft=%d cft=%d\n", numpft, numcft)

    CLM.varpar_init!(CLM.varpar, 1, numpft, numcft, nlevurb)
    @printf("varpar: natpft_lb=%d natpft_ub=%d cft_lb=%d cft_ub=%d maxveg=%d\n",
            CLM.varpar.natpft_lb, CLM.varpar.natpft_ub,
            CLM.varpar.cft_lb, CLM.varpar.cft_ub, CLM.varpar.maxveg)
    @printf("CTSM expects:   natpft_lb=0 natpft_ub=%d cft_lb=%d cft_ub=%d maxveg=%d\n",
            numpft - 1, numpft, numpft + numcft - 1, numpft + numcft - 1)

    n_merge = length(CLM.pftcon.mergetoclmpft)
    @printf("mergetoclmpft length=%d; loop indexes up to maxveg+1=%d -> %s\n",
            n_merge, CLM.varpar.maxveg + 1,
            CLM.varpar.maxveg + 1 <= n_merge ? "IN BOUNDS" : "OUT OF BOUNDS (would throw)")

    CLM.varcon_init!()
    surf = CLM.SurfaceInputData()
    CLM.surfrd_get_data!(surf, 1, 1, SURFDATA)

    nz = [(CLM.varpar.cft_lb + m - 1, surf.wt_cft[1, m])
          for m in 1:size(surf.wt_cft, 2) if surf.wt_cft[1, m] > 0.0]
    @printf("\nCFTs with nonzero weight after the static read: %d\n", length(nz))
    for (t, w) in nz
        @printf("  itype %3d  wt %.6f\n", t, w)
    end

    julia_itypes = [t for (t, _) in nz]
    @printf("\nJulia itypes : %s\n", julia_itypes)
    @printf("Fortran itypes: %s\n", FORTRAN_ITYPES)
    @printf("MATCH: %s\n", julia_itypes == FORTRAN_ITYPES ? "yes" : "NO")
    @printf("weight sum: %.10f\n", sum(w for (_, w) in nz))
end

main()
