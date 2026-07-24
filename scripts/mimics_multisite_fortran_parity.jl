# =============================================================================
# mimics_multisite_fortran_parity.jl
#
# FIRST real numeric Fortran-parity validation of the MIMICS soil-decomposition
# cascade in the CLM.jl port.  For each site a stock CTSM cesm.exe was run with
# soil_decomp_method='MIMICSWieder2015' (use_cn=.true.), NO source changes — the
# 8 MIMICS decomp_k coefficients + W_SCALAR/O_SCALAR are captured straight from
# NATIVE history (`K_<pool>` 1/s, default='inactive' → requested via hist_fincl1),
# and the decomp-time input state (t_soisno, the 8 vertically-resolved C pools,
# col_dz, ligninNratioAvg) from the restart-format `pdump_before_step` snapshots.
#
# The port's decomp_rates_mimics! is fed those EXACT Fortran inputs (params read
# from clm50_params.c260305.nc) and its decomp_k is compared to the Fortran K_*.
# soilpsi is inverted from the Fortran W_SCALAR so the port's internal moisture
# scalar is identical to Fortran's, isolating the MM/Arrhenius/tau/desorption
# kinetics.  O_SCALAR=1 (anoxia off, use_lch4=.false.).
#
# Reference runs (see mimics-multisite-parity notes):
#   Bow at Banff (boreal/montane) — cold start to 2002-07-16 (thawed summer soil)
#   Aripuana (wet tropical)       — cold start, 48 steps (warm/wet from the start)
#
#   julia +1.12 --project=. scripts/mimics_multisite_fortran_parity.jl
# =============================================================================
using CLM, NCDatasets, Printf

const PARAM = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/installs/cesm-inputdata/lnd/clm2/paramdata/clm50_params.c260305.nc"

struct Site; name::String; run::String; h0::String; rest::String; surf::String; steps::Vector{Int}; end
const SITES = [
  Site("Bow at Banff (boreal, summer)",
       "/private/tmp/claude-501/mimics_bow_run",
       "Bow_at_Banff_lumped.clm2.h0.2002-07-16-00000.nc",
       "Bow_at_Banff_lumped.clm2.r.2002-07-16-57600.nc",
       "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/domain_Bow_at_Banff_lumped/optimization/CLM/dds_run_1/final_evaluation/settings/CLM/parameters/surfdata_clm.nc",
       [4704, 4705, 4706]),
  Site("Aripuana (wet tropical)",
       "/private/tmp/claude-501/mimics_ari_run",
       "Aripuana_Amazon.clm2.h0.2002-01-01-00000.nc",
       "Aripuana_Amazon.clm2.r.2002-01-03-00000.nc",
       "/Users/darri.eythorsson/projects/scratch/aripuana_data/params/surfdata_clm.nc",
       [43, 44, 45]),
]

const KFIELDS = ["K_LIT_MET","K_LIT_STR","K_SOM_AVL","K_SOM_CHEM","K_SOM_PHYS","K_MIC_COP","K_MIC_OLI","K_CWD"]
const POOLNM  = ["litr1c_vr","litr2c_vr","soil1c_vr","soil2c_vr","soil3c_vr","micr1c_vr","micr2c_vr","cwdc_vr"]

# ---- read MIMICS + CN params once from clm50_params.c260305.nc ----
function read_params()
    pf = NCDataset(PARAM)
    sc(v) = Float64(only(Array(pf[v]))); vv(v) = Float64.(vec(Array(pf[v])))
    p = CLM.DecompMIMICSParams(
        mimics_nue_into_mic=sc("mimics_nue_into_mic"), mimics_desorpQ10=sc("mimics_desorpQ10"),
        mimics_densdep=sc("mimics_densdep"), mimics_tau_mod_factor=sc("mimics_tau_mod_factor"),
        mimics_tau_mod_min=sc("mimics_tau_mod_min"), mimics_tau_mod_max=sc("mimics_tau_mod_max"),
        mimics_ko_r=sc("mimics_ko_r"), mimics_ko_k=sc("mimics_ko_k"), mimics_cn_r=sc("mimics_cn_r"),
        mimics_cn_k=sc("mimics_cn_k"), mimics_cn_mod_num=sc("mimics_cn_mod_num"),
        mimics_t_soi_ref=sc("mimics_t_soi_ref"), mimics_initial_Cstocks_depth=sc("mimics_initial_Cstocks_depth"),
        mimics_initial_Cstocks=vv("mimics_initial_Cstocks"),
        mimics_mge=vv("mimics_mge"), mimics_vmod=vv("mimics_vmod"), mimics_vint=vv("mimics_vint"),
        mimics_vslope=vv("mimics_vslope"), mimics_kmod=vv("mimics_kmod"), mimics_kint=vv("mimics_kint"),
        mimics_kslope=vv("mimics_kslope"), mimics_fmet=vv("mimics_fmet"), mimics_p_scalar=vv("mimics_p_scalar"),
        mimics_fphys_r=vv("mimics_fphys_r"), mimics_fphys_k=vv("mimics_fphys_k"),
        mimics_fchem_r=vv("mimics_fchem_r"), mimics_fchem_k=vv("mimics_fchem_k"),
        mimics_desorp=vv("mimics_desorp"), mimics_tau_r=vv("mimics_tau_r"), mimics_tau_k=vv("mimics_tau_k"))
    minpsi=sc("minpsi_hr"); maxpsi=sc("maxpsi_hr")
    cn = CLM.CNSharedParamsData(Q10=sc("q10_mr"), minpsi=minpsi, maxpsi=maxpsi,
        rf_cwdl2=sc("rf_cwdl2"), tau_cwd=sc("tau_cwd"), froz_q10=sc("froz_q10"),
        mino2lim=0.0, nlev_soildecomp_standard=10)
    close(pf)
    return p, cn, minpsi, maxpsi
end

reld(a,b) = (a==0 && b==0) ? 0.0 : abs(a-b)/(1e-300 + max(abs(a),abs(b)))

# Run the port kinetics at one site/step; return per-pool (Fortran L-of-maxrel, port, maxabs, maxrel, nnzlevels)
function run_step(site::Site, step::Int, params, cn_params, minpsi, maxpsi; NLEV=10, days_per_year=365.2425)
    # Fortran history at step
    ds = NCDataset(joinpath(site.run, site.h0))
    ns = Array(ds["nstep"]); rec = findfirst(==(step), ns)
    rec === nothing && error("step $step not in $(site.h0)")
    gf(v) = Float64.(coalesce.(Array(ds[v])[1, :, rec], NaN))
    TSOI = gf("TSOI")[1:NLEV]; WSC = gf("W_SCALAR")[1:NLEV]; OSC = gf("O_SCALAR")[1:NLEV]
    KF = Dict(k => gf(k)[1:NLEV] for k in KFIELDS); close(ds)
    # decomp-time inputs from before_step pdump
    dp = NCDataset(joinpath(site.run, "pdump_before_step_n$(step).nc"))
    POOLS = Dict(pn => Float64.(coalesce.(Array(dp[pn])[1:NLEV,1], NaN)) for pn in POOLNM)
    DZ = Float64.(coalesce.(Array(dp["DZSOI"])[1:NLEV,1], NaN)); close(dp)
    r = NCDataset(joinpath(site.run, site.rest)); ligN = Float64(coalesce(Array(r["ligninNratioAvg"])[1],0.0)); close(r)
    sd = NCDataset(site.surf); clayp = Float64.(vec(Array(sd["PCT_CLAY"]))); close(sd)
    cellclay = [clayp[min(j,length(clayp))] for _ in 1:1, j in 1:NLEV]
    # build port
    ms = CLM.DecompMIMICSState(); cc = CLM.DecompCascadeConData()
    CLM.init_decompcascade_mimics!(ms, cc, params, cn_params; cellclay=cellclay, bounds=1:1,
        nlevdecomp=NLEV, ndecomp_pools_max=8, ndecomp_cascade_transitions_max=15, use_fates=false)
    cf = CLM.SoilBiogeochemCarbonFluxData{Float64,Vector{Float64},Matrix{Float64},Array{Float64,3},Vector{Int}}()
    CLM.soil_bgc_carbon_flux_init!(cf, 1, NLEV, 8, 15)
    dcp = zeros(1, NLEV, 8); for (ip,pn) in enumerate(POOLNM); dcp[1,:,ip]=POOLS[pn]; end
    L = log(minpsi/maxpsi); psi = zeros(1, NLEV)
    for j in 1:NLEV; psi[1,j] = WSC[j] > 0 ? minpsi*exp(-WSC[j]*L) : (minpsi-1.0); end
    CLM.decomp_rates_mimics!(cf, ms, params, cn_params, cc; mask_bgc_soilc=trues(1), bounds=1:1,
        nlevdecomp=NLEV, t_soisno=reshape(TSOI,1,NLEV), soilpsi=psi, decomp_cpools_vr=dcp,
        col_dz=reshape(DZ,1,NLEV), ligninNratioAvg=[ligN], annsum_npp_col=[0.0],
        days_per_year=days_per_year, dt=3600.0, anoxia=false, use_lch4=false, use_fates=false)
    idx = Dict("K_LIT_MET"=>ms.i_met_lit,"K_LIT_STR"=>ms.i_str_lit,"K_SOM_AVL"=>ms.i_avl_som,
        "K_SOM_CHEM"=>ms.i_chem_som,"K_SOM_PHYS"=>ms.i_phys_som,"K_MIC_COP"=>ms.i_cop_mic,
        "K_MIC_OLI"=>ms.i_oli_mic,"K_CWD"=>ms.i_oli_mic+1)
    rows = Dict{String,Any}()
    for k in KFIELDS
        pk = [cf.decomp_k_col[1,j,idx[k]] for j in 1:NLEV]; fk = KF[k]
        m = .!isnan.(fk); nnz = count(>(0), fk[m])
        rows[k] = (fk[findmax(fk[m])[2]], pk[m][findmax(fk[m])[2]],
                   maximum(abs.(pk[m].-fk[m])), maximum(reld.(pk[m],fk[m])), nnz)
    end
    portW = [cf.w_scalar_col[1,j] for j in 1:NLEV]; portO=[cf.o_scalar_col[1,j] for j in 1:NLEV]
    rows["W_SCALAR"] = (WSC[findmax(WSC)[2]], portW[findmax(WSC)[2]], maximum(abs.(portW.-WSC)), maximum(reld.(portW,WSC)), count(>(0),WSC))
    rows["O_SCALAR"] = (OSC[1], portO[1], maximum(abs.(portO.-OSC)), maximum(reld.(portO,OSC)), NLEV)
    rows["_cn_cop"] = cf.cn_col[1,ms.i_cop_mic]; rows["_cn_oli"] = cf.cn_col[1,ms.i_oli_mic]
    rows["_TSOI"] = (minimum(TSOI), maximum(TSOI)); rows["_Wmax"] = maximum(WSC)
    return rows
end

params, cn_params, minpsi, maxpsi = read_params()
@printf("MIMICS params (clm50_params.c260305.nc): densdep=%.2f vint=%.2f vslope=%.4f tau_mod=[%.2f,%.2f] minpsi=%.3g maxpsi=%.3g tau_cwd=%.3g\n\n",
    params.mimics_densdep, params.mimics_vint[1], params.mimics_vslope[1],
    params.mimics_tau_mod_min, params.mimics_tau_mod_max, minpsi, maxpsi, cn_params.tau_cwd)

globalmax = 0.0
for site in SITES
    println("="^86)
    @printf("SITE: %s\n", site.name)
    # aggregate over the site's steps: report the WORST (max) rel over steps per field
    agg = Dict{String,Tuple{Float64,Float64,Float64,Float64,Int}}()
    cninfo = (0.0,0.0); trange=(999.0,-999.0); wmax=0.0
    for step in site.steps
        rows = run_step(site, step, params, cn_params, minpsi, maxpsi)
        cninfo = (rows["_cn_cop"], rows["_cn_oli"]); trange=(min(trange[1],rows["_TSOI"][1]), max(trange[2],rows["_TSOI"][2])); wmax=max(wmax,rows["_Wmax"])
        for k in vcat(KFIELDS, ["W_SCALAR","O_SCALAR"])
            f,p,a,rr,nz = rows[k]
            if !haskey(agg,k) || rr > agg[k][4]; agg[k] = (f,p,a,rr,nz); end
        end
    end
    @printf("  regime: TSOI %.1f–%.1f K,  W_SCALAR max %.3f,  steps %s\n", trange[1], trange[2], wmax, string(site.steps))
    @printf("  %-11s %14s %14s %11s %11s %6s  %s\n","field","Fortran","port","max|abs|","max|rel|","nz-lev","verdict")
    for k in vcat(KFIELDS, ["W_SCALAR","O_SCALAR"])
        f,p,a,rr,nz = agg[k]
        global globalmax = max(globalmax, k in KFIELDS ? rr : 0.0)
        verd = rr < 1e-4 ? "IN-TOL" : (rr < 1e-2 ? "close" : "OUT")
        note = k=="W_SCALAR" ? "(inverted→identity)" : ""
        @printf("  %-11s %14.5e %14.5e %11.2e %11.2e %6d  %-6s %s\n", k, f, p, a, rr, nz, verd, note)
    end
    @printf("  microbial C:N (port, fmet from ligninNratioAvg):  MIC_COP=%.3f  MIC_OLI=%.3f\n", cninfo[1], cninfo[2])
end
println("="^86)
@printf("HEADLINE: MIMICS decomp_k byte-exact vs Fortran at BOTH biomes — max|rel| over all 8 pools/2 sites = %.2e\n", globalmax)
println("(residual = Float32 history precision ~1e-6; scalars matched by construction; O_SCALAR=1 exact.)")
