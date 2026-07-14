# ==========================================================================
# probe_cbalance_terms.jl — is the CARBON balance check able to fail?
#
# Runs ONE real use_cn=true clm_drv! step from the BGC IC and prints every
# term the column/gridcell C balance reads. A term that is 0.0 or NaN across
# the board is a DEAD aggregation, not a physical zero.
#
#   julia +1.12 --project=. scripts/probe_cbalance_terms.jl
# ==========================================================================
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const BGC_DUMPDIR = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_dumps_2202"
const NSTEP = 1753153

function main()
    inst, bounds = run_one_parity_step!(NSTEP; use_cn=true, dumpdir=BGC_DUMPDIR,
                                        step_date=DateTime(2003, 1, 1))
    cvf = inst.bgc_vegetation.cnveg_carbonflux_inst
    scf = inst.soilbiogeochem_carbonflux
    snf = inst.soilbiogeochem_nitrogenflux
    scs = inst.soilbiogeochem_carbonstate

    function show1(label, arr)
        if isempty(arr)
            println(rpad(label, 34), "EMPTY (unallocated)")
        else
            v = Float64.(arr[1:min(end, 4)])
            println(rpad(label, 34), join([string(round(x, sigdigits=6)) for x in v], "  "))
        end
    end

    println("\n=== CARBON: column terms the C balance check reads (col 1..4) ===")
    show1("gpp_col  [INPUT]",            cvf.gpp_col)
    show1("er_col   [OUTPUT]",           cvf.er_col)
    show1("  ar_col",                    cvf.ar_col)
    show1("  rr_col",                    cvf.rr_col)
    show1("  hr_col (soilbgc)",          scf.hr_col)
    show1("npp_col",                     cvf.npp_col)
    show1("nep_col",                     cvf.nep_col)
    show1("fire_closs_col",              cvf.fire_closs_col)
    show1("hrv_xsmrpool_to_atm_col",     cvf.hrv_xsmrpool_to_atm_col)
    show1("xsmrpool_to_atm_col",         cvf.xsmrpool_to_atm_col)
    show1("gru_conv_cflux_col",          cvf.gru_conv_cflux_col)
    show1("wood_harvestc_col",           cvf.wood_harvestc_col)
    show1("som_c_leached_col",           scf.som_c_leached_col)
    show1("totc_col (state)",            scs.totc_col)

    println("\n=== CARBON: patch terms (the p2c SOURCES — are these alive?) ===")
    show1("gpp_patch",                   cvf.gpp_patch)
    show1("ar_patch",                    cvf.ar_patch)
    show1("npp_patch",                   cvf.npp_patch)

    println("\n=== CARBON: gridcell ===")
    show1("nbp_grc [INPUT]",             cvf.nbp_grc)
    show1("nee_grc",                     cvf.nee_grc)
    show1("landuseflux_grc",             cvf.landuseflux_grc)
    show1("totc_grc",                    scs.totc_grc)

    println("\n=== NITROGEN: the #221 1e-7 warning suspects ===")
    show1("f_n2o_nit_col",               snf.f_n2o_nit_col)
    show1("denit_col",                   snf.denit_col)
    show1("smin_no3_leached_col",        snf.smin_no3_leached_col)
    show1("smin_no3_runoff_col",         snf.smin_no3_runoff_col)
    show1("som_n_leached_col",           snf.som_n_leached_col)
    show1("f_n2o_nit_vr_col[:,1] (src)", isempty(snf.f_n2o_nit_vr_col) ? Float64[] : snf.f_n2o_nit_vr_col[:, 1])

    return 0
end

exit(main())
