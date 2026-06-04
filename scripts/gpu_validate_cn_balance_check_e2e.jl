# ==========================================================================
# gpu_validate_cn_balance_check_e2e.jl — end-to-end Metal parity for the CN
# mass-balance check module (CNBalanceCheckMod). Exercises every kernel:
#   begin_cn_gridcell_balance!  (per-gridcell)
#   begin_cn_column_balance!    (per-column, masked)
#   c2g_unity!                  (per-column -> gridcell scatter)
#   c_balance_check!            (column + gridcell numeric kernels)
#   n_balance_check!            (column kernel w/ 21-array device-view bundle
#                                + gridcell kernel)
#
# Data is BALANCED (inputs == outputs, storage unchanged) so the host-only
# error() path never throws, but per-index-distinct so indexing bugs surface.
#
#   julia --project=scripts scripts/gpu_validate_cn_balance_check_e2e.jl
# ==========================================================================
using CLM
using Printf
import Metal
include(joinpath(@__DIR__, "gpu_backends.jl"))

struct _F32 end
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:AbstractFloat}) = Float32.(x)
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{<:Integer}) = x
CLM.Adapt.adapt_storage(::_F32, x::AbstractArray{Bool}) = x

function reldiff(a, b)
    A = Array(a); B = Array(b); m = 0.0
    for i in eachindex(A, B)
        (isnan(A[i]) && isnan(B[i])) && continue
        m = max(m, abs(Float64(A[i]) - Float64(B[i])) / (1.0 + max(abs(Float64(A[i])), abs(Float64(B[i])))))
    end
    return m
end

function build()
    nc = 4; ng = 2
    bal = CLM.CNBalanceData(); CLM.cn_balance_init!(bal, nc, ng)

    col_data = CLM.ColumnData()
    col_data.gridcell = [1, 1, 2, 2]
    col_data.wtgcell  = fill(0.5, nc)

    grc_data = CLM.GridcellData()
    grc_data.latdeg = [45.0, 46.0]; grc_data.londeg = [-90.0, -91.0]

    cstate = CLM.SoilBiogeochemCarbonStateData()
    cstate.totc_col = [100.0, 110.0, 120.0, 130.0]
    cstate.totc_grc = [105.0, 125.0]   # == c2g(totc_col) so begin/check agree
    nstate = CLM.SoilBiogeochemNitrogenStateData()
    nstate.totn_col = [10.0, 11.0, 12.0, 13.0]
    nstate.totn_grc = [10.5, 12.5]

    # carbon fluxes: balanced (er == gpp), everything else zero
    cflux = CLM.SoilBiogeochemCarbonFluxData()
    cflux.hr_col = zeros(nc); cflux.som_c_leached_col = zeros(nc); cflux.fates_litter_flux = zeros(nc)
    cncf = CLM.CNVegCarbonFluxData()
    cncf.gpp_col = [1.0, 2.0, 3.0, 4.0]; cncf.er_col = [1.0, 2.0, 3.0, 4.0]
    for f in (:fire_closs_col, :hrv_xsmrpool_to_atm_col, :xsmrpool_to_atm_col, :wood_harvestc_col,
              :gru_conv_cflux_col, :gru_wood_productc_gain_col, :crop_harvestc_to_cropprodc_col)
        setfield!(cncf, f, zeros(nc))
    end
    for f in (:nbp_grc, :dwt_seedc_to_leaf_grc, :dwt_seedc_to_deadstem_grc)
        setfield!(cncf, f, zeros(ng))
    end

    # nitrogen fluxes: balanced (denit == ndep), everything else zero
    nflux = CLM.SoilBiogeochemNitrogenFluxData()
    nflux.ndep_to_sminn_col = [5.0, 6.0, 7.0, 8.0]; nflux.denit_col = [5.0, 6.0, 7.0, 8.0]
    for f in (:nfix_to_sminn_col, :ffix_to_sminn_col, :fert_to_sminn_col, :soyfixn_to_sminn_col,
              :supplement_to_sminn_col, :sminn_leached_col, :smin_no3_leached_col,
              :smin_no3_runoff_col, :f_n2o_nit_col, :som_n_leached_col, :sminn_to_plant_col,
              :fates_litter_flux)
        setfield!(nflux, f, zeros(nc))
    end
    cnnf = CLM.CNVegNitrogenFluxData()
    for f in (:fire_nloss_col, :wood_harvestn_col, :gru_conv_nflux_col,
              :gru_wood_productn_gain_col, :crop_harvestn_to_cropprodn_col)
        setfield!(cnnf, f, zeros(nc))
    end
    for f in (:gru_conv_nflux_grc, :gru_wood_productn_gain_grc, :dwt_seedn_to_leaf_grc,
              :dwt_seedn_to_deadstem_grc, :dwt_conv_nflux_grc)
        setfield!(cnnf, f, zeros(ng))
    end

    cprod = CLM.CNProductsData(); CLM.cn_products_init!(cprod, ng)
    nprod = CLM.CNProductsData(); CLM.cn_products_init!(nprod, ng)

    return (; bal, col_data, grc_data, cstate, nstate, cflux, nflux, cncf, cnnf,
            cprod, nprod, mask=fill(true, nc), is_fates=fill(false, nc),
            drib=zeros(ng), bounds_c=1:nc, bounds_g=1:ng, dt=1800.0, nc, ng)
end

function run_all!(d)
    CLM.cn_balance_init!(d.bal, d.nc, d.ng)
    CLM.begin_cn_gridcell_balance!(d.bal, d.cstate, d.nstate, d.cprod, d.nprod, d.bounds_g;
        hrv_xsmrpool_amount_left=d.drib, gru_conv_cflux_amount_left=d.drib,
        dwt_conv_cflux_amount_left=d.drib)
    CLM.begin_cn_column_balance!(d.bal, d.cstate, d.nstate, d.mask, d.bounds_c)
    CLM.c_balance_check!(d.bal, d.cflux, d.cstate, d.cncf, d.cprod, d.col_data, d.grc_data,
        d.mask, d.bounds_c, d.bounds_g, d.dt; is_fates_col=d.is_fates,
        hrv_xsmrpool_amount_left=d.drib, gru_conv_cflux_amount_left=d.drib,
        dwt_conv_cflux_amount_left=d.drib)
    CLM.n_balance_check!(d.bal, d.nflux, d.nstate, d.cnnf, d.nprod, d.col_data, d.grc_data,
        d.mask, d.bounds_c, d.bounds_g, d.dt; is_fates_col=d.is_fates)
end

function main(backend)
    println("=" ^ 66); println("END-TO-END Metal parity for CNBalanceCheckMod"); println("=" ^ 66)
    if backend === nothing; println("  No GPU backend."); return 0; end
    name, _, FT = backend
    @printf("  Backend: %s   (precision: %s)\n", name, FT)
    H = build(); B = build()
    run_all!(H)
    mf(x) = CLM.Adapt.adapt(Metal.MtlArray, CLM.Adapt.adapt(_F32(), x))
    mmask(b) = Metal.MtlArray(collect(b))
    D = (; bal=mf(B.bal), col_data=mf(B.col_data), grc_data=mf(B.grc_data),
         cstate=mf(B.cstate), nstate=mf(B.nstate), cflux=mf(B.cflux), nflux=mf(B.nflux),
         cncf=mf(B.cncf), cnnf=mf(B.cnnf), cprod=mf(B.cprod), nprod=mf(B.nprod),
         mask=mmask(B.mask), is_fates=mmask(B.is_fates), drib=mf(B.drib),
         bounds_c=B.bounds_c, bounds_g=B.bounds_g, dt=B.dt, nc=B.nc, ng=B.ng)
    if !(D.bal.begcb_col isa Metal.MtlArray); println("  BLOCKED."); return 2; end
    run_all!(D)

    checks = [("begcb_col", H.bal.begcb_col, D.bal.begcb_col),
              ("begnb_col", H.bal.begnb_col, D.bal.begnb_col),
              ("begcb_grc", H.bal.begcb_grc, D.bal.begcb_grc),
              ("begnb_grc", H.bal.begnb_grc, D.bal.begnb_grc),
              ("endcb_col", H.bal.endcb_col, D.bal.endcb_col),
              ("endnb_col", H.bal.endnb_col, D.bal.endnb_col),
              ("endcb_grc", H.bal.endcb_grc, D.bal.endcb_grc),
              ("endnb_grc", H.bal.endnb_grc, D.bal.endnb_grc),
              ("totc_grc (c2g scatter)", H.cstate.totc_grc, D.cstate.totc_grc),
              ("totn_grc (c2g scatter)", H.nstate.totn_grc, D.nstate.totn_grc)]
    nfail = 0
    for (nm, a, b) in checks
        dd = reldiff(a, b); ok = dd < 1f-3
        @printf("  [%s] %-26s rel = %.3e\n", ok ? "PASS" : "FAIL", nm, dd); ok || (nfail += 1)
    end
    println(); println(nfail == 0 ? "  CNBalanceCheckMod kernels MATCH CPU ON $name ($FT)" : "  DIVERGENCE ($nfail).")
    return nfail == 0 ? 0 : 1
end
const BACKEND = detect_backend()
exit(main(BACKEND))
