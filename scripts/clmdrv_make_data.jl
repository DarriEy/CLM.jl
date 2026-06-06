# clmdrv_make_data.jl — cloned from test/test_clm_driver.jl make_driver_data.
# Builds a minimal Float64 driver state and sets the shared module globals
# (varpar / varcon / pftcon / SNOW_DZ* / urban_ctrl). Kept in sync with the test.

function make_driver_data(; ng=2, nl=3, nc=4, np=6)
    CLM.varpar_init!(CLM.varpar, 1, 14, 2, 5)
    CLM.varcon_init!()

    nlevsno = CLM.varpar.nlevsno
    nlevgrnd = CLM.varpar.nlevgrnd
    nlevtot = nlevsno + nlevgrnd
    nlevcan = CLM.NLEVCAN

    inst = CLM.CLMInstances()
    CLM.clm_instInit!(inst; ng=ng, nl=nl, nc=nc, np=np,
                       nlevdecomp_full=CLM.varpar.nlevdecomp_full)

    for c in 1:nc
        inst.column.landunit[c] = 1
        inst.column.gridcell[c] = 1
        inst.column.snl[c] = 0
        inst.column.patchi[c] = 1
        inst.column.patchf[c] = 0
        inst.column.nbedrock[c] = CLM.varpar.nlevsoi
    end
    for l in 1:nl
        inst.landunit.itype[l] = CLM.ISTSOIL
        inst.landunit.urbpoi[l] = false
        inst.landunit.lakpoi[l] = false
        inst.landunit.active[l] = true
    end
    for p in 1:np
        inst.patch.active[p] = true
        inst.patch.landunit[p] = 1
        inst.patch.gridcell[p] = 1
        inst.patch.column[p] = min(p, nc)
        inst.patch.itype[p] = 1
    end

    ppercol = max(1, np ÷ nc)
    pidx = 1
    for c in 1:nc
        inst.column.patchi[c] = pidx
        inst.column.patchf[c] = min(pidx + ppercol - 1, np)
        npc = inst.column.patchf[c] - inst.column.patchi[c] + 1
        for p in inst.column.patchi[c]:inst.column.patchf[c]
            inst.patch.column[p] = c
            inst.patch.wtcol[p] = 1.0 / npc
        end
        pidx = inst.column.patchf[c] + 1
    end
    for p in pidx:np
        inst.patch.column[p] = nc
        inst.patch.wtcol[p] = 0.0
    end

    bounds = CLM.BoundsType(
        begg=1, endg=ng, begl=1, endl=nl, begc=1, endc=nc, begp=1, endp=np,
        begCohort=0, endCohort=0,
        level=CLM.BOUNDS_LEVEL_CLUMP, clump_index=1)

    filt = CLM.ClumpFilter()
    CLM.alloc_filters!(filt, nc, np, nl)
    filt.allc .= true
    filt.nolakec .= true
    filt.nolakep .= true
    filt.soilc .= true
    filt.soilp .= true
    filt.nolakeurbanp .= true
    filt.nourbanp .= true
    filt.nourbanc .= true
    filt.hydrologyc .= true
    filt.urbanc .= false
    filt.urbanl .= false
    filt.urbanp .= false
    filt.snowc .= false
    filt.nosnowc .= true
    filt.do_smb_c .= false
    filt.exposedvegp .= true
    filt.noexposedvegp .= false
    filt.lakec .= false
    filt.lakep .= false
    filt.lakesnowc .= false
    filt.lakenosnowc .= false
    filt.bgc_soilc .= false
    filt.bgc_vegp .= false
    filt.pcropp .= false
    filt.soilnopcropp .= true
    filt.actfirec .= false
    filt.actfirep .= false

    filt_ia = CLM.ClumpFilter()
    CLM.alloc_filters!(filt_ia, nc, np, nl)
    filt_ia.allc .= true
    filt_ia.nolakec .= true
    filt_ia.nolakep .= true
    filt_ia.nourbanc .= true
    filt_ia.nourbanp .= true

    CLM.urban_read_nml!(CLM.urban_ctrl)

    # Initialize the global photosynthesis params (the suite gets this as a
    # side effect of test_photosynthesis.jl running before the driver test).
    CLM.photo_params_init!(CLM.params_inst)

    dzmin = zeros(nlevsno)
    dzmax_u = zeros(nlevsno)
    dzmax_l = zeros(nlevsno)
    dzmin[1] = 0.010; dzmax_u[1] = 0.02; dzmax_l[1] = 0.03
    dzmin[2] = 0.015; dzmax_u[2] = 0.05; dzmax_l[2] = 0.07
    for j in 3:nlevsno
        dzmin[j] = dzmax_u[j-1] * 0.5
        dzmax_u[j] = 2.0 * dzmax_u[j-1] + 0.01
        dzmax_l[j] = dzmax_u[j] + dzmax_l[j-1]
        if j == nlevsno
            dzmax_u[j] = floatmax(Float64)
            dzmax_l[j] = floatmax(Float64)
        end
    end
    CLM.SNOW_DZMIN[] = dzmin
    CLM.SNOW_DZMAX_U[] = dzmax_u
    CLM.SNOW_DZMAX_L[] = dzmax_l

    config = CLM.CLMDriverConfig()

    npft = CLM.MXPFT + 1
    p = CLM.pftcon
    p.dleaf         = fill(0.04, npft)
    p.slatop        = fill(0.01, npft)
    p.leafcn        = fill(25.0, npft)
    p.flnr          = fill(0.1, npft)
    p.fnitr         = fill(0.1, npft)
    p.mbbopt        = fill(9.0, npft)
    p.c3psn         = fill(1.0, npft)
    p.woody         = fill(0.0, npft)
    p.smpso         = fill(-66000.0, npft)
    p.smpsc         = fill(-275000.0, npft)
    p.z0mr          = fill(0.055, npft)
    p.displar       = fill(0.67, npft)
    p.xl            = fill(0.1, npft)
    p.rhol          = fill(0.1, npft, 2)
    p.rhos          = fill(0.2, npft, 2)
    p.taul          = fill(0.05, npft, 2)
    p.taus          = fill(0.1, npft, 2)
    p.medlynintercept = fill(100.0, npft)
    p.medlynslope     = fill(6.0, npft)
    p.crop          = fill(0.0, npft)

    photosyns = CLM.PhotosynthesisData()
    CLM.photosynthesis_data_init!(photosyns, np)

    return inst, bounds, filt, filt_ia, config, photosyns
end
