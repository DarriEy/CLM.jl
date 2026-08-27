#!/usr/bin/env julia

"""Extract the committed Bow 2009 variable subset from its frozen CTSM history file."""

using NCDatasets
using SHA

const SOURCE_SHA256 = "6793ab309c9da08a5e18e7acc79c488b9758f16cee4f21ccedfb20794cc9734b"
const VARIABLES = [
    "time", "TSA", "TG", "EFLX_LH_TOT", "FSH", "FSA", "QRUNOFF", "QOVER",
    "H2OSNO", "SNOWDP", "FCEV", "FCTR", "FGEV", "RAIN", "SNOW", "levsoi",
    "SOILLIQ", "SOILICE", "levgrnd", "TSOI", "FSNO", "FSAT", "ZWT", "BTRANMN",
    "QINFL", "QDRAI", "ELAI", "ESAI",
]

file_sha256(path) = bytes2hex(open(sha256, path))

function copy_variable!(dst, src, name)
    input = src[name]
    attributes = Dict{String,Any}(String(key) => value for (key, value) in input.attrib)
    fillvalue = pop!(attributes, "_FillValue", nothing)
    if name in ("time", "levsoi", "levgrnd")
        fillvalue = NaN32
    end
    if name == "time"
        attributes["units"] = "days since 2002-01-01"
    end
    datatype = eltype(input.var)
    output = if isnothing(fillvalue)
        defVar(dst, name, datatype, NCDatasets.dimnames(input); attrib = attributes)
    else
        defVar(dst, name, datatype, NCDatasets.dimnames(input);
               attrib = attributes, fillvalue = fillvalue)
    end
    indices = ntuple(axis -> 1:size(input.var, axis), ndims(input.var))
    output[indices...] = input.var[indices...]
end

function extract_reference(source, output)
    isfile(source) || error("source history file not found: $source")
    file_sha256(source) == SOURCE_SHA256 || error("source history checksum does not match")
    ispath(output) && error("refusing to overwrite output: $output")

    NCDataset(source, "r") do src
        size(src["time"], 1) == 365 || error("source does not contain 365 daily records")
        NCDataset(output, "c") do dst
            defDim(dst, "time", 365)
            defDim(dst, "lndgrid", size(src["TSA"], 1))
            defDim(dst, "levsoi", size(src["levsoi"], 1))
            defDim(dst, "levgrnd", size(src["levgrnd"], 1))
            for name in VARIABLES
                copy_variable!(dst, src, name)
            end
            dst.attrib["source"] = "Fortran CLM5 (CTSM release-clm5.0.30) via SYMFLUENCE DDS calibration"
            dst.attrib["site"] = "Bow at Banff (51.17N, 115.57W)"
            dst.attrib["period"] = "2009 (evaluation year from DDS run_1)"
            dst.attrib["params"] = "DDS best: 29 params, KGE=0.917"
        end
    end
end

length(ARGS) == 2 || error("usage: extract_bow_reference.jl SOURCE_HISTORY OUTPUT")
extract_reference(ARGS[1], ARGS[2])
println("wrote ", ARGS[2], " from source SHA-256 ", SOURCE_SHA256)
