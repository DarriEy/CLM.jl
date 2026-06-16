# Probe: run summer CN steps with PHS+LUNA, compare converged grass vegwp
# (4 segments) + htop/tsai/qflx against the Fortran after_canopyfluxes dump.
include(joinpath(@__DIR__, "fortran_parity_common.jl"))

const DD = "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data/clm_bgc_spinup/bgc_ref_summer"

function fvec(ds, v)
    haskey(ds, v) || return Float64[]
    return Float64.(ds[v][:])
end

function run_step(nstep)
    inst, bounds = run_one_parity_step!(nstep; use_cn=true, dumpdir=DD,
        use_hydrstress=true, use_luna=true,
        step_date=DateTime(2002,1,1)+Hour(nstep-1753153),
        forcing_file=replace(FFORCING, "clmforc.2003.nc"=>"clmforc.2002.nc"))
    ds = NCDataset(joinpath(DD, "pdump_after_canopyfluxes_n$(nstep).nc"), "r")
    f_pfts = Int.(ds["pfts1d_itypveg"][:])
    jl_pfts = Int.(inst.patch.itype)
    vegwp_f = ds["vegwp"][:, :]   # dims (vegwcs, pft) in Julia col-major read
    htop_f = fvec(ds,"htop"); tsai_f = fvec(ds,"tsai")
    cs = inst.canopystate
    seg = ["SUN","SHA","XYL","ROOT"]
    println("\n===== nstep $nstep =====")
    for pj in 1:length(jl_pfts)
        it = jl_pfts[pj]; it == 0 && continue
        pf = findfirst(==(it), f_pfts); pf === nothing && continue
        tag = it == 1 ? "TREE" : "PFT$it"
        @printf("%-6s htop J=%.4g F=%.4g | tsai J=%.4g F=%.4g\n",
                tag, cs.htop_patch[pj], htop_f[pf], cs.tsai_patch[pj], tsai_f[pf])
        for s in 1:4
            jv = cs.vegwp_patch[pj, s]; fv = vegwp_f[s, pf]
            @printf("   vegwp[%-4s] J=%12.1f F=%12.1f  Δ=%10.1f\n", seg[s], jv, fv, jv-fv)
        end
        dj = cs.vegwp_patch[pj,1]-cs.vegwp_patch[pj,4]
        df = vegwp_f[1,pf]-vegwp_f[4,pf]
        @printf("   leaf-root drop: J=%.1f F=%.1f\n", dj, df)
    end
    close(ds)
end

for n in (1757852, 1757855)
    run_step(n)
end
