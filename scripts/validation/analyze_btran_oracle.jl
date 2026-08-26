#!/usr/bin/env julia

using NCDatasets
using Printf

const FIELDS = (
    "BTRAN", "VEGWP", "SMP", "HK", "ROOT_CONDUCTANCE", "PBOT240",
    "PCO2_240", "PO2_240", "ELAI240", "PAR240D_Z", "PAR240X_Z",
    "T_A10", "T_VEG10_DAY", "T_VEG10_NIGHT", "RH10_AF", "RB10",
    "VCMX25_Z", "JMX25_Z", "PNLC_Z", "ENZS_Z",
    "K_SOIL_ROOT", "SOIL_CONDUCTANCE", "ROOTFR", "BSW", "SUCSAT",
    "LAISUN", "LAISHA", "ELAI", "ESAI", "TSAI", "FDRY", "FORC_RHO",
    "FORC_PBOT", "THM", "QFLX_TRAN_VEG", "BSUN", "BSHA",
    "NRAD", "TLAI_Z", "FSUN_Z", "FABD_SUN_Z", "FABD_SHA_Z",
    "FABI_SUN_Z", "FABI_SHA_Z",
)

function read_oracle(path)
    NCDataset(path) do ds
        data = Dict{String,Any}(name => Array(ds[name]) for name in FIELDS if haskey(ds, name))
        data["nstep"] = Int(ds["nstep"][])
        data["ymd"] = Int(ds["timemgr_rst_curr_ymd"][])
        data["tod"] = Int(ds["timemgr_rst_curr_tod"][])
        data
    end
end

function main(run_dir)
    rx = r"^btran_oracle_(before_step|after_hydrologydrainage)_n(\d+)\.nc$"
    paths = filter(path -> match(rx, basename(path)) !== nothing,
                   readdir(run_dir; join=true))
    isempty(paths) && error("no BTRAN oracle files found in $run_dir")

    records = [(path=path, data=read_oracle(path)) for path in paths]
    sort!(records; by=r -> (r.data["nstep"], r.path))
    target = filter(r -> r.data["nstep"] >= 2832, records)
    isempty(target) && error("no target-window records found")

    @printf("records=%d target_records=%d steps=%d:%d\n", length(records),
            length(target), first(target).data["nstep"], last(target).data["nstep"])
    npft = length(first(target).data["BTRAN"])
    for p in 1:npft
        values = [(Float64(r.data["BTRAN"][p]), r) for r in target
                  if isfinite(r.data["BTRAN"][p])]
        value, record = values[argmin(first.(values))]
        d = record.data
        @printf("pft=%d min_BTRAN=%.17g nstep=%d ymd=%08d tod=%05d boundary=%s\n",
                p, value, d["nstep"], d["ymd"], d["tod"],
                occursin("before_step", record.path) ? "before_step" : "after_hydrologydrainage")
    end

    for name in FIELDS
        vals = reduce(vcat, vec(Float64.(r.data[name])) for r in target)
        finite = filter(isfinite, vals)
        @printf("field=%s finite=%d total=%d min=%.17g max=%.17g\n",
                name, length(finite), length(vals), minimum(finite), maximum(finite))
    end
end

main(length(ARGS) == 1 ? ARGS[1] : pwd())
