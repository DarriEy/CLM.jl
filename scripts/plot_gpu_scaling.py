#!/usr/bin/env python3
# ==========================================================================
# plot_gpu_scaling.py — publication log-log scaling figure for the CLM.jl
# clm_drv! biogeophys step, CPU (Float64) vs Apple Metal GPU (Float32), as the
# number of soil columns N grows. Data from scripts/gpu_scaling_bench.jl
# (min-of-5 trials per point). Writes paper/gpu_scaling.png.
# ==========================================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# --- benchmark data (one clm_drv! SP step, wall time in ms) ---
N   = np.array([1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144], float)
CPU = np.array([1.247, 2.657, 6.640, 23.815, 93.226, 446.301, 1987.041,
                8076.900, 32697.276, 144387.384])
GPU = np.array([395.807, 352.742, 387.017, 377.982, 390.383, 509.174, 496.658,
                520.149, 731.864, 1118.205])
speedup = CPU / GPU

# --- crossover (speedup == 1), interpolated in log-log space ---
ls = np.log(speedup); lN = np.log(N)
i = np.where(np.diff(np.sign(ls)))[0][0]           # bracket where speedup crosses 1
Ncross = np.exp(lN[i] + (0 - ls[i]) * (lN[i+1]-lN[i]) / (ls[i+1]-ls[i]))

C_CPU, C_GPU, C_X = "#2171b5", "#2ca02c", "0.35"
plt.rcParams.update({"font.size": 11, "axes.grid": True,
                     "grid.color": "0.9", "grid.linewidth": 0.6})
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.5))

# ---- Panel 1: wall time vs N (log-log) ----
ax.loglog(N, CPU, "o-", color=C_CPU, lw=1.8, ms=6, label="CPU (Float64)")
ax.loglog(N, GPU, "s-", color=C_GPU, lw=1.8, ms=6, label="GPU — Metal (Float32)")
ax.axvline(Ncross, color=C_X, ls="--", lw=1.1)
ax.axvspan(Ncross, N.max()*1.6, color=C_GPU, alpha=0.06)
ax.text(Ncross*1.15, 2.0, f"crossover\n≈ {Ncross:,.0f} cols", color=C_X, fontsize=9, va="bottom")
ax.text(N.max(), GPU[-1]*1.5, f"{speedup[-1]:.0f}× faster", color=C_GPU,
        fontsize=9.5, ha="right", va="bottom", fontweight="bold")
ax.text(1.2, CPU[0]*0.62, "single column\n(parity runs):\nGPU ~320× slower",
        color=C_CPU, fontsize=8.2, va="top")
ax.set_xlabel("Number of soil columns  $N$")
ax.set_ylabel("Wall time per clm_drv! step (ms)")
ax.set_title("clm_drv! (biogeophys) step — CPU vs GPU scaling", fontsize=11)
ax.legend(loc="upper left", fontsize=9.5, framealpha=0.95)
ax.set_xlim(0.7, N.max()*2.2)

# ---- Panel 2: speedup vs N (semilog-x) ----
ax2.semilogx(N, speedup, "D-", color="#6a3d9a", lw=1.8, ms=5.5)
ax2.axhline(1.0, color=C_X, ls="--", lw=1.1)
ax2.axvline(Ncross, color=C_X, ls="--", lw=1.1)
ax2.axhspan(1.0, speedup.max()*1.4, color=C_GPU, alpha=0.06)
ax2.text(Ncross*1.2, 1.15, "GPU wins  →", color=C_X, fontsize=9)
for xx, ss in ((16384, speedup[7]), (262144, speedup[-1])):
    ax2.annotate(f"{ss:.0f}×", (xx, ss), textcoords="offset points", xytext=(-2, 6),
                 fontsize=9, color="#6a3d9a", ha="right")
ax2.axhline(280, color="0.6", ls=":", lw=1.0)
ax2.text(2, 300, "asymptotic ceiling ~280×  (per-column cost ratio)", fontsize=8.2, color="0.4")
ax2.set_yscale("log")
ax2.set_xlabel("Number of soil columns  $N$")
ax2.set_ylabel("GPU speedup  (CPU time / GPU time)")
ax2.set_title("Speedup vs problem size", fontsize=11)
ax2.yaxis.set_major_formatter(mticker.FuncFormatter(lambda y, _: f"{y:g}×"))
ax2.set_ylim(0.002, 420)

fig.suptitle("CLM.jl on Apple Metal (Float32): single-point is overhead-bound, "
             "grid scale is throughput-bound", fontsize=12, y=1.00)
plt.tight_layout(rect=[0, 0, 1, 0.97])
out = __file__.rsplit("/", 1)[0] + "/gpu_scaling.png"
plt.savefig(out, dpi=200, bbox_inches="tight")
print("saved", out)
print(f"crossover N ≈ {Ncross:,.0f} columns; max measured speedup {speedup[-1]:.0f}× at N={int(N[-1]):,}")
