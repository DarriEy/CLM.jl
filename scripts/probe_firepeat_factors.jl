# Attribution probe for the BAF_PEATF 7.7% divergence.
# baf_peatf (boreal) = C/secsphr * exp(-PI*max(wf2,0)/denom) * min(1,(tsoi17-tfrz)/10) * peatf * (1-fsat)
# Prints Julia's wf2/tsoi17/fsat + the three factors and Julia vs Fortran baf_peatf,
# and back-derives Fortran's implied product to attribute the residual per-factor.
include(joinpath(@__DIR__, "fortran_parity_firepeat.jl"))

using NCDatasets
const TFRZ = 273.15
const DENOM = 0.3
const C = 0.09e-4 / 3600.0   # boreal_peatfire_c /secsphr

for n in [4720, 4804, 4888]
    inst, _ = run_one_parity_step!(n;
        use_cn=true, use_lch4=true, use_hydrstress=true, dumpdir=DUMPDIR_PEAT,
        step_date=DATE0+Hour(n-N0), forcing_file=FFORCING2002, forcing_date=peat_forcing_date(n),
        fndep=FNDEP, cnfire_method=:li2016, flnfm=FLNFM, fhdm=FHDM,
        pre_step_hook=(i,b,df)->begin
            apply_bow_lifire_inparm!(i); inject_fire_accum!(i,b,df)
            fill!(i.cnfire_li2014.peatf_lf_col, PEATF_INJECT)
        end)
    wf2  = inst.water.waterdiagnosticbulk_inst.wf2_col[1]
    ts17 = inst.temperature.t_soi17cm_col[1]
    fsat = inst.sat_excess_runoff.fsat_col[1]
    bpj  = inst.bgc_vegetation.cnveg_state_inst.baf_peatf_col[1]
    f_exp  = exp(-pi*max(wf2,0.0)/DENOM)
    f_temp = max(0.0, min(1.0, (ts17-TFRZ)/10.0))
    f_wet  = 1.0 - fsat
    dump = joinpath(DUMPDIR_PEAT, "bgcdump_after_fire_n$(n).nc")
    bpf = NCDataset(dump,"r") do ds; Float64(ds["BAF_PEATF"][1]); end
    # implied Fortran product = bpf/(C*peatf) ; Julia product = f_exp*f_temp*f_wet
    prodJ = f_exp*f_temp*f_wet
    prodF = bpf/(C*PEATF_INJECT)
    println("n=$n")
    @printf("  wf2=%.6f tsoi17=%.4f (t-tfrz=%.4f) fsat=%.6g\n", wf2, ts17, ts17-TFRZ, fsat)
    @printf("  factors: exp(-PI*wf2/0.3)=%.6f  ftemp=%.6f  (1-fsat)=%.6f\n", f_exp, f_temp, f_wet)
    @printf("  baf_peatf Julia=%.6e  Fortran=%.6e  rel=%.4f\n", bpj, bpf, abs(bpj-bpf)/bpf)
    @printf("  product Julia=%.6e Fortran(implied)=%.6e  ratio=%.5f\n", prodJ, prodF, prodJ/prodF)
    # if only wf2 differs: implied Fortran wf2 = -0.3/PI*log(prodF/(f_temp*f_wet))
    wf2F = -DENOM/pi*log(prodF/(f_temp*f_wet))
    @printf("  if residual is ALL wf2: Fortran wf2=%.6f (Julia %.6f, dwf2=%.6f)\n", wf2F, wf2, wf2-wf2F)
end
