# fates_planthydro_replay.jl — teacher-forced replay of the FATES 1-D hydraulics solve
# against a Fortran `hydro_pdump.txt` (instrument_planthydro.py).
#
# At every hydraulics call the port's rhizosphere shells, cohort water contents and patch
# transpiration are overwritten with the Fortran pre-solve values. Hydraulics_BC then
# runs, and its post-solve change in cohort water content is compared with the Fortran
# one. This isolates the solver from the host: given the same inputs the changes should
# agree to round-off (1e-7..1e-10 relative on the fates_parity_1pt case). `old` swaps in
# the pre-fix copy semantics of GetImTaylorKAB's derivative arguments, for comparison.
#
# usage: julia --project=. scripts/validation/fates_fortran_parity/fates_planthydro_replay.jl \
#            <hydro_pdump.txt> <out.txt> <nsteps> [old]
const FDUMP = ARGS[1]
const OUT = ARGS[2]
const NSTEPS = parse(Int, ARGS[3])
const OLD = length(ARGS) >= 4 && ARGS[4] == "old"
include(joinpath(@__DIR__, "fates_planthydro_parity.jl"))   # run_julia
using Printf

const FREC = Dict{Any,Vector{Float64}}()
for line in eachline(FDUMP)
    t = split(line)
    if t[1] == "COH"
        FREC[("COH", parse(Int, t[2]), parse(Int, t[3]), parse(Int, t[4]))] = parse.(Float64, t[6:end])
    else
        FREC[(t[1], parse(Int, t[2]))] = parse.(Float64, t[4:end])
    end
end
const IO_ = open(OUT, "w")
const PD = Ref(0)

if OLD   # pre-fix semantics: GetImTaylorKAB's derivative zeroing does not reach the caller
    @eval CLM @inline function _imtaylor_kab_path!(k_eff, a_term, b_term, j::Integer, kmax_up, kmax_dn,
                                         i_up::Integer, i_dn::Integer, ftc_node, h_node,
                                         dftc_dtheta_node, dpsi_dtheta_node)
        k_eff[j], a_term[j], b_term[j] = GetImTaylorKAB(kmax_up, kmax_dn,
            ftc_node[i_up], ftc_node[i_dn], h_node[i_up], h_node[i_dn],
            dftc_dtheta_node[i_up], dftc_dtheta_node[i_dn],
            dpsi_dtheta_node[i_up], dpsi_dtheta_node[i_dn])
        return nothing
    end
end

@eval CLM function hydraulics_drive(nsites::Integer, sites::AbstractVector,
                                    bc_in::AbstractVector, bc_out::AbstractVector, dtime::Real)
    Main.PD[] += 1
    k = Main.PD[]
    FillDrainRhizShells(nsites, sites, bc_in, bc_out)
    s = 1
    ch = sites[s].si_hydr; nl = ch.nlevrhiz; nsl = bc_in[s].nlevsoil
    sp = Main.FREC[("SPRE", k)]
    off = 1 + nsl
    ch.h2osoi_liqvol_shell[1:nl, :] .= reshape(sp[off+1:off+nl*nshell], nl, nshell)
    ntr = length(sp) - (off + 2 * nl * nshell + nl)
    for i in 1:min(ntr, length(bc_in[s].qflx_transp_pa))
        bc_in[s].qflx_transp_pa[i] = sp[off + 2 * nl * nshell + nl + i]
    end
    pre = Dict{Any,Vector{Float64}}()
    ifp = 0
    cp = sites[s].oldest_patch
    while cp !== nothing
        if cp.nocomp_pft_label != nocomp_bareground
            ifp += 1
            cc = cp.tallest
            while cc !== nothing
                key = ("COH", k, ifp, cc.pft)
                if haskey(Main.FREC, key)
                    v = Main.FREC[key]; h = cc.co_hydr
                    h.th_ag[1] = v[7]; h.th_ag[2] = v[8]; h.th_troot = v[9]
                    h.th_aroot[1:nl] .= v[10:9+nl]
                end
                cc = cc.shorter
            end
        end
        cp = cp.younger
    end
    Hydraulics_BC(nsites, sites, bc_in, bc_out, dtime)
    ifp = 0
    worst = 0.0; wk = ""
    cp = sites[s].oldest_patch
    while cp !== nothing
        if cp.nocomp_pft_label != nocomp_bareground
            ifp += 1
            cc = cp.tallest
            while cc !== nothing
                key = ("COH", k, ifp, cc.pft)
                if haskey(Main.FREC, key)
                    v = Main.FREC[key]; h = cc.co_hydr
                    o = 9 + nl
                    post = [h.th_ag[1], h.th_ag[2], h.th_troot, h.th_aroot[1:nl]...]
                    fpost = v[o+1:o+3+nl]
                    # relative error of the CHANGE in water content (pre -> post)
                    dJ = post .- v[7:9+nl]; dF = fpost .- v[7:9+nl]
                    e = maximum(abs.(dJ .- dF) ./ max.(abs.(dF), 1e-14))
                    el = abs(post[1] - fpost[1]) / abs(fpost[1])
                    psiF = v[o+4+nl]; psiJ = h.psi_ag[1]
                    println(Main.IO_, Main.@sprintf("%4d %2d %2d qtopF %.4e qtopJ %.4e dthl F %.6e J %.6e  max_rel_dth %.3e  psil F %.6g J %.6g",
                                               k, ifp, cc.pft, v[3], h.qtop / dtime, dF[1], dJ[1], e, psiF, psiJ))
                    if e > worst; worst = e; wk = "pft $(cc.pft)"; end
                end
                cc = cc.shorter
            end
        end
        cp = cp.younger
    end
    println(Main.IO_, Main.@sprintf("CALL %4d worst rel(dth) %.3e (%s)", k, worst, wk))
    flush(Main.IO_)
    return nothing
end

run_julia(NSTEPS)
close(IO_)
