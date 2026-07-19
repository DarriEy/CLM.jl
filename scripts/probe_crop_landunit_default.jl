# =====================================================================
# create_crop_landunit default divergence probe
#
# CTSM keys create_crop_landunit on `.not. use_fates`
# (namelist_defaults_ctsm.xml:2377-2378), and CLMBuildNamelist.pm:2248-2250
# makes `.false.` a FATAL error for any non-FATES run. CLM.jl instead derives it
# from `use_crop` (clm_initialize.jl:138) and branches varpar_init! on `use_crop`
# (varpar.jl:178).
#
# This probe builds the REAL subgrid (initGridCells!) both ways for a DEFAULT
# non-crop, non-FATES run and reports what actually moves.
#
#   julia +1.12 --project=. scripts/probe_crop_landunit_default.jl
# =====================================================================

using CLM

const DATAROOT = joinpath(homedir(), "Library/CloudStorage",
                          "GoogleDrive-dareyt@gmail.com", "My Drive", "data",
                          "SYMFLUENCE_data")
const PARAMFILE = joinpath(DATAROOT, "domain_Aripuana_Amazon", "settings",
                           "CLM", "parameters", "clm5_params.nc")

# A standard scorecard domain (PCT_CROP = 0) and the crop-CFT point
# (PCT_CROP = 45.0%). The contrast between them IS the blast radius.
const CASES = [
    ("Bow_at_Banff (PCT_CROP=0)",
     joinpath(DATAROOT, "domain_Bow_at_Banff_lumped", "settings", "CLM",
              "parameters", "surfdata_clm.nc")),
    ("US_plains crop point (PCT_CROP=45%)",
     joinpath(DATAROOT, "crop_cft_surfdata", "surfdata_cropCFT_USplains_1pt.nc")),
]

"""Build the subgrid for a given create_crop_landunit setting, use_crop=false."""
function build(surfdata::String, create_crop_landunit::Bool)
    # Default non-crop, non-FATES run.
    CLM.varctl.use_crop = false
    CLM.varctl.use_fates = false
    CLM.varctl.irrigate = false
    CLM.varctl.create_crop_landunit = create_crop_landunit

    (numpft, numcft) = CLM.surfrd_get_num_patches(surfdata)
    nlevurb = CLM.surfrd_get_nlevurb(surfdata)
    CLM.varpar_init!(CLM.varpar, 1, numpft, numcft, nlevurb)
    CLM.varcon_init!()
    CLM.readParameters!(PARAMFILE)

    surf = CLM.SurfaceInputData()
    CLM.surfrd_get_data!(surf, 1, 1, surfdata)

    ng = 1
    (nl, nc, np) = CLM.count_subgrid_elements(surf, ng)

    # count_subgrid_elements and set_landunit_crop_noncompete! disagree when
    # create_crop_landunit is on but cft_size==0: the counter has an
    # `n_cft == 0` fallback that reserves a column+patch (surfdata.jl:740-747)
    # while the builder bails out entirely (init_gridcells.jl:147). That state is
    # exactly the half-applied fix (flag flipped, varpar_init! still on
    # use_crop), so surface it rather than dying inside initGridCells!.
    if create_crop_landunit && surf.wt_lunit[1, CLM.ISTCROP] > 0.0 &&
       CLM.varpar.cft_size == 0
        error("inconsistent: create_crop_landunit=true with cft_size==0 — " *
              "count says nl=$nl but the builder will skip the crop landunit")
    end
    bounds = CLM.BoundsType(begg=1, endg=ng, begl=1, endl=nl,
                            begc=1, endc=nc, begp=1, endp=np,
                            level=CLM.BOUNDS_LEVEL_CLUMP, clump_index=1)
    grc = CLM.GridcellData(); CLM.gridcell_init!(grc, ng)
    lun = CLM.LandunitData(); CLM.landunit_init!(lun, nl)
    col = CLM.ColumnData();   CLM.column_init!(col, nc)
    pch = CLM.PatchData();    CLM.patch_init!(pch, np)
    CLM.initGridCells!(bounds, surf, grc, lun, col, pch)

    crop_p = [p for p in 1:np
              if lun.itype[col.landunit[pch.column[p]]] == CLM.ISTCROP]

    return (; nl, nc, np,
            cft_size = CLM.varpar.cft_size,
            cft_lb = CLM.varpar.cft_lb,
            cft_ub = CLM.varpar.cft_ub,
            natpft_ub = CLM.varpar.natpft_ub,
            maxveg = CLM.varpar.maxveg,
            wt_soil = round(surf.wt_lunit[1, CLM.ISTSOIL], digits=6),
            wt_crop = round(surf.wt_lunit[1, CLM.ISTCROP], digits=6),
            lun_types = sort([lun.itype[l] for l in 1:nl]),
            crop_itypes = [pch.itype[p] for p in crop_p],
            crop_wts = [col.wtlunit[pch.column[p]] for p in crop_p])
end

function report(name, a, b)
    println("\n", "="^72)
    println(name)
    println("="^72)
    println(rpad("  field", 16), rpad("CLM.jl now (clu=false)", 26), "CTSM (clu=true)")
    println("-"^72)
    for k in (:nl, :nc, :np, :cft_size, :cft_lb, :cft_ub, :natpft_ub, :maxveg,
              :wt_soil, :wt_crop, :lun_types, :crop_itypes, :crop_wts)
        av, bv = getfield(a, k), getfield(b, k)
        mark = av == bv ? "  " : "* "
        println(mark, rpad(String(k), 14), rpad(repr(av), 26), repr(bv))
    end
    println("\n(* = differs)")
end

for (name, f) in CASES
    if !isfile(f)
        @info "skipping (surfdata absent)" name f
        continue
    end
    local a, b
    try
        a = build(f, false)   # what CLM.jl does today
        b = build(f, true)    # what CTSM does
    catch e
        println("\n", "="^72, "\n", name, "\n", "="^72)
        println("BUILD FAILED: ", sprint(showerror, e))
        continue
    end
    report(name, a, b)
end
