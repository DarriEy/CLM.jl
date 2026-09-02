#!/usr/bin/env julia

"""Emit a machine-readable environment record for a GMD experiment run."""

using Dates
using Pkg
using SHA

function command_output(cmd::Cmd)
    strip(read(cmd, String))
end

function optional_command(cmd::Cmd)
    try
        command_output(cmd)
    catch
        nothing
    end
end

function file_sha256(path)
    isfile(path) || return nothing
    open(path, "r") do io
        bytes2hex(sha256(io))
    end
end

function git_record(root)
    head = command_output(`git -C $root rev-parse HEAD`)
    status = command_output(`git -C $root status --porcelain=v1 --untracked-files=all`)
    Dict(
        "commit" => head,
        "branch" => optional_command(`git -C $root branch --show-current`),
        "describe" => optional_command(`git -C $root describe --always --dirty --tags`),
        "clean" => isempty(status),
        "status" => isempty(status) ? String[] : split(status, '\n'),
        "remote" => optional_command(`git -C $root remote get-url origin`),
    )
end

function package_record()
    dependencies = Pkg.dependencies()
    rows = [Dict(
        "name" => pkg.name,
        "uuid" => string(uuid),
        "version" => isnothing(pkg.version) ? nothing : string(pkg.version),
        "tree_hash" => isnothing(pkg.tree_hash) ? nothing : string(pkg.tree_hash),
        "direct" => pkg.is_direct_dep,
    ) for (uuid, pkg) in dependencies]
    sort!(rows; by = row -> (row["name"], row["uuid"]))
end

function environment_record(root)
    Dict(
        "schema_version" => 1,
        "recorded_utc" => Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sssZ"),
        "git" => git_record(root),
        "julia" => Dict(
            "version" => string(VERSION),
            "commit" => Base.GIT_VERSION_INFO.commit,
            "threads" => Threads.nthreads(),
            "cpu_threads" => Sys.CPU_THREADS,
            "word_size" => Sys.WORD_SIZE,
            "machine" => Sys.MACHINE,
            "kernel" => Sys.KERNEL,
            "architecture" => string(Sys.ARCH),
            "bindir" => Sys.BINDIR,
        ),
        "host" => Dict(
            "hostname" => gethostname(),
            "cpu" => Sys.CPU_NAME,
            "os_version" => optional_command(`sw_vers -productVersion`),
        ),
        "environment" => Dict(
            "project_sha256" => file_sha256(joinpath(root, "Project.toml")),
            "manifest_sha256" => file_sha256(joinpath(root, "Manifest.toml")),
            "packages" => package_record(),
        ),
    )
end

root = normpath(joinpath(@__DIR__, "..", ".."))
# Stdlib-only JSON emission. The original `using JSON` resolved through the
# DEPOT DEFAULT environment, not the locked project Manifest — an
# undocumented dependency that broke the first frozen campaign preflight in
# the clean container (exactly the defect class this runner exists to
# catch). The environment record must not depend on anything outside the
# frozen project + stdlib.
function _json(io::IO, x; indent::Int = 0)
    pad = "  "^indent
    if x isa AbstractDict
        print(io, "{\n")
        for (i, k) in enumerate(sort!(collect(keys(x)); by = string))
            i > 1 && print(io, ",\n")
            print(io, pad, "  ", repr(string(k)), ": ")
            _json(io, x[k]; indent = indent + 1)
        end
        print(io, "\n", pad, "}")
    elseif x isa AbstractVector
        print(io, "[\n")
        for (i, v) in enumerate(x)
            i > 1 && print(io, ",\n")
            print(io, pad, "  ")
            _json(io, v; indent = indent + 1)
        end
        print(io, "\n", pad, "]")
    elseif x isa AbstractString || x isa Symbol || x isa VersionNumber
        print(io, repr(string(x)))
    elseif x === nothing
        print(io, "null")
    elseif x isa Bool || x isa Integer
        print(io, x)
    elseif x isa Real
        print(io, isfinite(x) ? string(x) : repr(string(x)))
    else
        print(io, repr(string(x)))
    end
end

_json(stdout, environment_record(root))
println()
println()
