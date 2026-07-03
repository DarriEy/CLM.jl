#!/usr/bin/env python3
"""Multi-biome parity scorecard overview heatmap for the README.

Rows = biomes, columns = the 51 validated output variables (grouped by category).
Cell colour = signed relative error Julia-vs-Fortran (annual mean), diverging around
0. Near-zero quantities (surface water on dry/glacier sites) are judged by an
absolute noise floor and shown as agreement. Renders paper/parity_scorecard.png.

Paths point at the local SYMFLUENCE_data references; regenerate after a sweep.
"""
import os, numpy as np
from datetime import datetime
from netCDF4 import Dataset
import cftime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle

_D = os.environ.get("SYMFLUENCE_DATA", "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data")
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "paper", "data")

# (key, fortran_h0, julia_nc, short label, biome descriptor)
DOMAINS = [
    ("Bow",        f"{_D}/clm_parity_run/Bow_at_Banff_lumped.clm2.h0.2003-01-01-00000.nc", "julia_clm_bow_phs_2003.nc", "Bow at Banff", "Alpine (Canada)"),
    ("Stillwater", f"{_D}/domain_Stillwater_Oklahoma/optimization/CLM/dds_clm_dds_calibration/final_evaluation/Stillwater_Oklahoma.clm2.h0.2003-01-01-00000.nc", "julia_clm_stillwater_phs_2003.nc", "Stillwater", "Semi-arid prairie (USA)"),
    ("Aripuana",   f"{_D}/domain_Aripuana_Amazon/simulations/clm_validation/CLM/Aripuana_Amazon.clm2.h0.2004-01-01-00000.nc", "julia_clm_aripuana_phs_2004.nc", "Aripuanã", "Tropical rainforest (Brazil)"),
    ("HubbardBrook", f"{_D}/domain_Temperate_HubbardBrook_USA/simulations/clm_tempforest/CLM/Temperate_HubbardBrook_USA.clm2.h0.2016-12-31-00000.nc", "julia_clm_hubbardbrook_phs_2017.nc", "Hubbard Brook", "Temperate deciduous forest (USA)"),
    ("Krycklan",   f"{_D}/domain_Boreal_Krycklan_Sweden/simulations/clm_boreal/CLM/Boreal_Krycklan_Sweden.clm2.h0.2012-12-30-00000.nc", "julia_clm_krycklan_phs_2013.nc", "Krycklan", "Boreal forest (Sweden)"),
    ("Tagus",      f"{_D}/domain_Mediterranean_Tagus_Spain/simulations/clm_mediterranean/CLM/Mediterranean_Tagus_Spain.clm2.h0.2012-12-30-00000.nc", "julia_clm_tagus_phs_2013.nc", "Tagus", "Mediterranean (Spain)"),
    ("Abisko",     f"{_D}/domain_Arctic_Abisko_Sweden/simulations/clm_arctic/CLM/Arctic_Abisko_Sweden.clm2.h0.2012-12-30-00000.nc", "julia_clm_abisko_phs_2013.nc", "Abisko", "Arctic tundra (Sweden)"),
    ("Massa",      f"{_D}/domain_Alps_Massa_Aletsch_CH/simulations/clm_alpineglacier/CLM/Alps_Massa_Aletsch_CH.clm2.h0.2012-12-30-00000.nc", "julia_clm_massa_phs_2013.nc", "Massa Aletsch", "Alpine glacier (Switzerland)"),
    ("Baltimore",  f"{_D}/domain_Urban_DeadRun_Baltimore/simulations/clm_urban/CLM/Urban_DeadRun_Baltimore.clm2.h0.2012-12-30-00000.nc", "julia_clm_baltimore_phs_2013.nc", "Baltimore", "Urban (USA)"),
    ("Iceland",    f"{_D}/domain_Iceland_Jokulsa_Fjollum/simulations/clm_glacier/CLM/Iceland_Jokulsa_Fjollum.clm2.h0.2016-12-31-00000.nc", "julia_clm_iceland_phs_2017.nc", "Iceland Jökulsá", "Glacier outwash (Iceland)"),
]

# (julia, fortran, label, group, unit-floor). Mirrors plot_parity_full.py's 51 vars.
VARS = [
    ("FSA","FSA","Absorbed SW","Energy",0.1),("FSR","FSR","Reflected SW","Energy",0.1),
    ("FIRA","FIRA","Net LW","Energy",0.1),("EFLX_SH_TOT","FSH","Sensible H","Energy",0.1),
    ("EFLX_LH_TOT","EFLX_LH_TOT","Latent H","Energy",0.1),("FCTR","FCTR","Transpiration","Energy",0.1),
    ("FCEV","FCEV","Canopy evap","Energy",0.1),("FGEV","FGEV","Ground evap","Energy",0.1),
    ("FGR","FGR","Ground heat","Energy",0.1),("SABV","SABV","SW abs veg","Energy",0.1),
    ("SABG","SABG","SW abs grnd","Energy",0.1),("FSRVD","FSRVD","Refl vis-dir","Energy",0.1),
    ("FSRNI","FSRNI","Refl nir-dif","Energy",0.1),("FIRE","FIRE","Emitted LW","Energy",0.1),
    ("TAUX","TAUX","Wind stress x","Energy",0.002),("TAUY","TAUY","Wind stress y","Energy",0.002),
    ("QRUNOFF","QRUNOFF","Total runoff","Water",0.02),("QOVER","QOVER","Surf runoff","Water",0.02),
    ("QFLX_DRAIN","QDRAI","Drainage","Water",0.02),("QFLX_EVAP_TOT","QFLX_EVAP_TOT","Total ET","Water",0.02),
    ("QSOIL","QSOIL","Ground evap q","Water",0.02),("QVEGT","QVEGT","Transp q","Water",0.02),
    ("QVEGE","QVEGE","Canopy evap q","Water",0.02),("QFLX_SNOMELT","QSNOMELT","Snowmelt","Water",0.02),
    ("RAIN","RAIN","Rainfall","Water",2e-6*86400),("SNOW","SNOW","Snowfall","Water",2e-6*86400),
    ("FSAT","FSAT","Sat fraction","Water",0.002),("H2OSFC","H2OSFC","Surface water","Water",0.02),
    ("FH2OSFC","FH2OSFC","Frac sfc water","Water",0.002),("TOTSOILLIQ","TOTSOILLIQ","Soil liquid","Water",0.02),
    ("H2OSNO","H2OSNO","Snow water eq","Snow",0.02),("SNOW_DEPTH","SNOW_DEPTH","Snow depth","Snow",0.001),
    ("FRAC_SNO","FSNO","Snow cover","Snow",0.002),
    ("TSA","TSA","2-m air T","State",0.5),("T_GRND","TG","Ground T","State",0.5),
    ("TV","TV","Veg T","State",0.5),("TSOI_10CM","TSOI_10CM","Soil T 10cm","State",0.5),
    ("SOILWATER_10CM","SOILWATER_10CM","Soil water 10cm","State",0.02),("SOILICE","TOTSOILICE","Soil ice","State",0.02),
    ("ZWT","ZWT","Water table","State",0.001),("ELAI","ELAI","Exposed LAI","State",0.02),
    ("BTRAN","BTRANMN","BTRAN","State",0.002),("U10","U10","10-m wind","State",0.05),
    ("Q2M","Q2M","2-m humidity","State",5e-5),("RH2M","RH2M","2-m RH","State",0.5),
    ("TSAI","TSAI","Stem area idx","State",0.02),("TLAI","TLAI","Total LAI","State",0.02),
    ("LAISUN","LAISUN","Sunlit LAI","State",0.02),("LAISHA","LAISHA","Shaded LAI","State",0.02),
    ("SOILRESIS","SOILRESIS","Soil resist","State",5.0),
    ("FPSN","FPSN","Photosynthesis","Carbon",0.05),
]
SCALE = {"QRUNOFF":86400,"QOVER":86400,"QFLX_DRAIN":86400,"QFLX_EVAP_TOT":86400,"QSOIL":86400,
         "QVEGT":86400,"QVEGE":86400,"QFLX_SNOMELT":86400,"RAIN":86400,"SNOW":86400}

def date_ord(ds):
    t = ds["time"]; dts = cftime.num2date(t[:], t.units, getattr(t,"calendar","noleap"))
    return np.array([d.year*10000+d.month*100+d.day for d in dts])

def series(ds, v, idx):
    if v is None or v not in ds.variables: return None
    a = np.ma.filled(ds[v][:].astype("f8"), np.nan)
    if a.ndim > 1: a = a.reshape(a.shape[0], -1)[:, 0]
    a = a.flatten()[idx]
    return np.where(np.isfinite(a) & (a > -9990), a, np.nan)

grid = np.full((len(DOMAINS), len(VARS)), np.nan)
passed = np.zeros((len(DOMAINS), len(VARS)), bool)
for i, (key, fh0, jnc, lab, biome) in enumerate(DOMAINS):
    dj = Dataset(os.path.join(DATA_DIR, jnc)); df = Dataset(fh0)
    jo, fo = date_ord(dj), date_ord(df)
    common = np.intersect1d(jo, fo)
    ji = np.array([np.where(jo==c)[0][0] for c in common])
    fi = np.array([np.where(fo==c)[0][0] for c in common])
    for k, (jn, fn, vl, grp, floor) in enumerate(VARS):
        sc = SCALE.get(jn, 1)
        j = series(dj, jn, ji); f = series(df, fn, fi)
        if j is None or f is None: continue
        jm, fm = np.nanmean(j)*sc, np.nanmean(f)*sc
        d = jm - fm
        pct = d/abs(fm)*100 if abs(fm) > 1e-9 else 0.0
        grid[i, k] = pct
        ok = (abs(d) < 0.5) if grp == "State" and "T" == vl.split()[-1][:1] and floor == 0.5 and vl not in ("Sat fraction",) else False
        # temperature vars: |Δ|<0.5K; else |%|<10 OR |Δ|<unit floor
        if floor == 0.5 and ("T" in vl):
            ok = abs(d) < 0.5
        else:
            ok = (abs(pct) < 10) or (abs(d) < floor)
        passed[i, k] = ok

npass = passed.sum(1)
print("Scorecard (within 10% / 0.5K):")
for i,(key,_,_,lab,biome) in enumerate(DOMAINS):
    print(f"  {lab:16s} {npass[i]:2d}/{len(VARS)}  {biome}")
print(f"TOTAL: {passed.sum()}/{passed.size}")

# ---- Plot ----
# For cells that pass only via the absolute noise floor (a physically near-zero
# quantity, e.g. surface water on a dry prairie), the relative % is meaningless —
# neutralise the displayed value so it isn't a misleading saturated colour.
disp = grid.copy()
for i in range(len(DOMAINS)):
    for k in range(len(VARS)):
        if passed[i, k] and np.isfinite(grid[i, k]) and abs(grid[i, k]) >= 10:
            disp[i, k] = 0.0

n_b, n_v = len(DOMAINS), len(VARS)
fig, ax = plt.subplots(figsize=(19, 6.6))
clip = np.clip(disp, -10, 10)
norm = TwoSlopeNorm(vmin=-10, vcenter=0, vmax=10)
im = ax.imshow(clip, aspect="auto", cmap="RdBu_r", norm=norm)
# mark the (few) genuine out-of-tolerance cells with an ×
for i in range(n_b):
    for k in range(n_v):
        if not passed[i, k] and np.isfinite(grid[i, k]):
            ax.text(k, i, "×", ha="center", va="center", color="black", fontsize=11, fontweight="bold")

ax.set_yticks(range(n_b))
ax.set_yticklabels([f"{lab}\n{biome}" for _,_,_,lab,biome in DOMAINS], fontsize=8.5)
ax.set_xticks(range(n_v))
ax.set_xticklabels([vl for _,_,vl,_,_ in VARS], rotation=90, fontsize=7.2)
# group separators + headers
groups = [g for _,_,_,g,_ in VARS]
bounds_g = [k for k in range(1, n_v) if groups[k] != groups[k-1]]
for b in bounds_g:
    ax.axvline(b-0.5, color="white", lw=2.5)
start = 0
for gi, g in enumerate(["Energy","Water","Snow","State","Carbon"]):
    idxs = [k for k in range(n_v) if groups[k]==g]
    if not idxs: continue
    ax.text(np.mean(idxs), -0.85, g, ha="center", va="bottom", fontsize=10, fontweight="bold")
for i in range(n_b+1): ax.axhline(i-0.5, color="white", lw=0.6)

cb = fig.colorbar(im, ax=ax, pad=0.008, fraction=0.018, extend="both")
cb.set_label("Julia − Fortran, annual mean (% of Fortran; clipped ±10%)", fontsize=9)
total = passed.sum()
ax.set_title(f"CLM.jl vs Fortran CLM5 — {n_b} biomes × {n_v} variables:  "
             f"{total}/{passed.size} within tolerance (|Δ| < 10% or < 0.5 K)  ·  × = out of tolerance",
             fontsize=13, fontweight="bold", pad=26)
plt.tight_layout()
out = os.path.join(os.path.dirname(__file__), "..", "scripts", "parity_scorecard.png")
plt.savefig(out, dpi=140, bbox_inches="tight")
print("saved", out)
