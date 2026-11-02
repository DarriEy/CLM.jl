# Compile-time sweep of Enzyme's type-analysis lattice knobs on single CUDA reverse
# phases (scripts/cuda_reverse_phase_diff.jl), one fresh Julia process per setting, run
# one after another (one GPU job at a time). Usage, from the server run folder:
#   powershell -File scripts\cuda_reverse_typeknob_sweep.ps1 [-phases a,b] [-configs 6:512,4:512]
#   (wraps the server helper rr.ps1; logs rev\knob_<depth>_<offset>.log)
param([string]$phases = "cf_rev_friction", [string]$configs = "6:512,4:512,6:128,3:64")
foreach ($c in $configs.Split(",")) {
    $d, $o = $c.Split(":")
    & C:\Users\darey\clm_cuda_run\rr.ps1 scripts/cuda_reverse_phase_diff.jl "knob_${d}_${o}" `
        "PHASES=$phases" "ENZYME_MAXTYPEDEPTH=$d" "ENZYME_MAXTYPEOFFSET=$o"
}
