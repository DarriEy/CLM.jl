#!/usr/bin/env julia

"""Emit a machine-readable environment record for a GMD experiment run."""

using Dates
using JSON
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
JSON.print(stdout, environment_record(root), 2)
println()
