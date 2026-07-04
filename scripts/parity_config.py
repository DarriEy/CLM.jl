"""Shared parity configuration: the domain registry, the validated-variable set,
and the two evaluation tiers.

TWO TIERS (deliberately separate):

1. COVERAGE scorecard (`plot_scorecard_overview.py`): |Δ| < 10 % (0.5 K for
   temperatures) or below a per-unit noise floor. This answers "does the port
   reproduce Fortran across many biomes and variables at all" — breadth, not
   exactness. It is NOT the parity bar.

2. SCIENTIFIC-PARITY gate (`parity_gate.py`, the real bar): every (domain,
   variable) cell must meet BOTH
     (a) an annual-mean gate — |Δ%| ≤ 1 %  (temperatures |ΔK| ≤ 0.05 K;
         near-zero quantities |Δ| ≤ the per-unit noise floor), AND
     (b) a daily-RMSE gate — normalized RMSE (RMSE / std of the Fortran daily
         series) ≤ 0.10  (near-constant series, std below floor: |RMSE| ≤ the
         per-unit noise floor; temperatures: daily RMSE ≤ 0.20 K).
   Scientific parity is declared ONLY when every cell that fails this gate is a
   DOCUMENTED EXCEPTION (below) — i.e. there are zero undocumented failures.
"""
import os
import numpy as np
import cftime

_D = os.environ.get("SYMFLUENCE_DATA", "/Users/darri.eythorsson/compHydro/SYMFLUENCE_data")
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "paper", "data")

# (key, fortran_h0, julia_nc, short label, biome descriptor)
DOMAINS = [
    ("Bow",        f"{_D}/clm_parity_run/Bow_at_Banff_lumped.clm2.h0.2003-01-01-00000.nc", "julia_clm_bow_phs_2003.nc", "Bow at Banff", "Alpine (Canada)"),
    ("Stillwater", f"{_D}/domain_Stillwater_Oklahoma/optimization/CLM/dds_clm_dds_calibration/final_evaluation/Stillwater_Oklahoma.clm2.h0.2003-01-01-00000.nc", "julia_clm_stillwater_phs_2003.nc", "Stillwater", "Semi-arid prairie (USA)"),
    ("Desert",     f"{_D}/domain_Desert_WalnutGulch_USA/simulations/clm_desert/CLM/Desert_WalnutGulch_USA.clm2.h0.2016-12-31-00000.nc", "julia_clm_desert_phs_2017.nc", "Walnut Gulch", "Hot desert (USA)"),
    ("Aripuana",   f"{_D}/domain_Aripuana_Amazon/simulations/clm_validation/CLM/Aripuana_Amazon.clm2.h0.2004-01-01-00000.nc", "julia_clm_aripuana_phs_2004.nc", "Aripuanã", "Tropical rainforest (Brazil)"),
    ("Savanna",    f"{_D}/domain_Savanna_Donga_Benin/simulations/clm_savanna/CLM/Savanna_Donga_Benin.clm2.h0.2016-12-31-00000.nc", "julia_clm_savanna_phs_2017.nc", "Donga", "Tropical savanna (Benin)"),
    ("Larch",      f"{_D}/domain_Yakutia_Larch_Russia/simulations/clm_larch/CLM/Yakutia_Larch_Russia.clm2.h0.2016-12-31-00000.nc", "julia_clm_larch_phs_2017.nc", "Yakutia", "Larch taiga (Russia)"),
    ("Peatland",   f"{_D}/domain_Peatland_MerBleue_Canada/simulations/clm_peatland/CLM/Peatland_MerBleue_Canada.clm2.h0.2016-12-31-00000.nc", "julia_clm_peatland_phs_2017.nc", "Mer Bleue", "Temperate bog (Canada)"),
    ("Steppe",     f"{_D}/domain_Steppe_Kherlen_Mongolia/simulations/clm_steppe/CLM/Steppe_Kherlen_Mongolia.clm2.h0.2016-12-31-00000.nc", "julia_clm_steppe_phs_2017.nc", "Kherlen", "Continental steppe (Mongolia)"),
    ("Eucalyptus", f"{_D}/domain_Eucalyptus_Tumbarumba_Australia/simulations/clm_eucalyptus/CLM/Eucalyptus_Tumbarumba_Australia.clm2.h0.2016-12-31-00000.nc", "julia_clm_eucalyptus_phs_2017.nc", "Tumbarumba", "Temperate evergreen broadleaf (Australia)"),
    ("Aspen",      f"{_D}/domain_Aspen_BOREAS_Canada/simulations/clm_aspen/CLM/Aspen_BOREAS_Canada.clm2.h0.2016-12-31-00000.nc", "julia_clm_aspen_phs_2017.nc", "BOREAS", "Boreal deciduous broadleaf (Canada)"),
    ("Paramo",     f"{_D}/domain_Paramo_Antisana_Ecuador/simulations/clm_paramo/CLM/Paramo_Antisana_Ecuador.clm2.h0.2016-12-31-00000.nc", "julia_clm_paramo_phs_2017.nc", "Antisana", "Tropical montane páramo (Ecuador)"),
    ("Cropland",   f"{_D}/domain_Cropland_Mead_USA/simulations/clm_crop/CLM/Cropland_Mead_USA.clm2.h0.2016-12-31-00000.nc", "julia_clm_cropland_phs_2017.nc", "Mead", "Cropland (USA)"),
    ("HubbardBrook", f"{_D}/domain_Temperate_HubbardBrook_USA/simulations/clm_tempforest/CLM/Temperate_HubbardBrook_USA.clm2.h0.2016-12-31-00000.nc", "julia_clm_hubbardbrook_phs_2017.nc", "Hubbard Brook", "Temperate deciduous forest (USA)"),
    ("Maritime",   f"{_D}/domain_Maritime_HJAndrews_USA/simulations/clm_maritime/CLM/Maritime_HJAndrews_USA.clm2.h0.2016-12-31-00000.nc", "julia_clm_maritime_phs_2017.nc", "HJ Andrews", "Maritime conifer (USA)"),
    ("Krycklan",   f"{_D}/domain_Boreal_Krycklan_Sweden/simulations/clm_boreal/CLM/Boreal_Krycklan_Sweden.clm2.h0.2012-12-30-00000.nc", "julia_clm_krycklan_phs_2013.nc", "Krycklan", "Boreal forest (Sweden)"),
    ("Tagus",      f"{_D}/domain_Mediterranean_Tagus_Spain/simulations/clm_mediterranean/CLM/Mediterranean_Tagus_Spain.clm2.h0.2012-12-30-00000.nc", "julia_clm_tagus_phs_2013.nc", "Tagus", "Mediterranean (Spain)"),
    ("Abisko",     f"{_D}/domain_Arctic_Abisko_Sweden/simulations/clm_arctic/CLM/Arctic_Abisko_Sweden.clm2.h0.2012-12-30-00000.nc", "julia_clm_abisko_phs_2013.nc", "Abisko", "Arctic tundra (Sweden)"),
    ("Massa",      f"{_D}/domain_Alps_Massa_Aletsch_CH/simulations/clm_alpineglacier/CLM/Alps_Massa_Aletsch_CH.clm2.h0.2012-12-30-00000.nc", "julia_clm_massa_phs_2013.nc", "Massa Aletsch", "Alpine glacier (Switzerland)"),
    ("Baltimore",  f"{_D}/domain_Urban_DeadRun_Baltimore/simulations/clm_urban/CLM/Urban_DeadRun_Baltimore.clm2.h0.2012-12-30-00000.nc", "julia_clm_baltimore_phs_2013.nc", "Baltimore", "Urban (USA)"),
    ("Iceland",    f"{_D}/domain_Iceland_Jokulsa_Fjollum/simulations/clm_glacier/CLM/Iceland_Jokulsa_Fjollum.clm2.h0.2016-12-31-00000.nc", "julia_clm_iceland_phs_2017.nc", "Iceland Jökulsá", "Glacier outwash (Iceland)"),
]

# (julia, fortran, label, group, unit-floor).
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
    ("FSRVI","FSRVI","Refl vis-dif","Energy",0.1),("FSRND","FSRND","Refl nir-dir","Energy",0.1),
    ("FSDSVD","FSDSVD","Inc vis-dir","Energy",0.1),("FSDSVI","FSDSVI","Inc vis-dif","Energy",0.1),
    ("FSDSND","FSDSND","Inc nir-dir","Energy",0.1),("FSDSNI","FSDSNI","Inc nir-dif","Energy",0.1),
    ("QFLX_INFL","QINFL","Infiltration","Water",0.02),("QH2OSFC","QH2OSFC","Sfc water flux","Water",0.02),
    ("ZWT_PERCH","ZWT_PERCH","Perched WT","State",0.001),("DSL","DSL","Dry surf layer","State",0.05),
    ("FSH_V","FSH_V","Sens H veg","Energy",0.1),("FSH_G","FSH_G","Sens H grnd","Energy",0.1),
    ("SABG_PEN","SABG_PEN","SW penetr","Energy",0.1),("QINTR","QINTR","Interception","Water",0.02),
    ("SNOWLIQ","SNOWLIQ","Snow liquid","Snow",0.02),("SNOWICE","SNOWICE","Snow ice","Snow",0.02),
    ("ESAI","ESAI","Exposed SAI","State",0.02),("H2OCAN","H2OCAN","Canopy water","State",0.02),
    ("FPSN","FPSN","Photosynthesis","Carbon",0.05),
]
# Keep columns grouped by category (stable-sort) for the heatmap.
_GORDER = {"Energy": 0, "Water": 1, "Snow": 2, "State": 3, "Carbon": 4}
VARS = sorted(VARS, key=lambda v: _GORDER[v[3]])
SCALE = {"QRUNOFF":86400,"QOVER":86400,"QFLX_DRAIN":86400,"QFLX_EVAP_TOT":86400,"QSOIL":86400,
         "QVEGT":86400,"QVEGE":86400,"QFLX_SNOMELT":86400,"RAIN":86400,"SNOW":86400,
         "QFLX_INFL":86400,"QH2OSFC":86400,"QINTR":86400}

# ── Coverage tier ──
COVERAGE_PCT = 10.0
COVERAGE_TEMP_K = 0.5

# ── Scientific-parity (strict) tier ──
STRICT_PCT = 1.0        # annual-mean relative gate (%)
STRICT_TEMP_K = 0.05    # annual-mean temperature gate (K)
STRICT_NRMSE = 0.10     # daily normalized-RMSE gate (RMSE / std_Fortran)
STRICT_TEMP_RMSE_K = 0.20   # daily temperature RMSE gate (K)
# Temperatures are judged in K (Δ and RMSE), not %.
TEMP_VARS = {"TSA", "T_GRND", "TV", "TSOI_10CM"}

# Documented exceptions: (domain_key, julia_var) -> reason. KNOWN, characterized
# residuals that do not meet the strict gate. They are exempt from the
# scientific-parity declaration but ALWAYS listed in the report. Add a cell here
# only with a written mechanism, never to hide an unexplained miss.
DOCUMENTED_EXCEPTIONS = {
    ("Savanna", "FSH_V"): (
        "Coupled Monin-Obukhov canopy-air-temperature (taf) partition on a sparse "
        "canopy; ~0.6 W/m2, and total FSH matches to 1.3%. Instrumented timestep-"
        "matched chase confirmed the resistances and temperatures reproduce Fortran; "
        "the veg-share sub-partition is a small residual of large conductance x "
        "temperature terms (coupled-solve floor)."),
    ("Baltimore", "SNOW_DEPTH"): (
        "Thin/intermittent trace-snow (snl=0) ghost-depth averaging present in BOTH "
        "models; over actual snow-present days the depth matches (Julia is if anything "
        "shallower). ~1.7 mm annual-mean artifact, not a density bug."),
}

def date_ord(ds):
    t = ds["time"]; dts = cftime.num2date(t[:], t.units, getattr(t, "calendar", "noleap"))
    return np.array([d.year*10000 + d.month*100 + d.day for d in dts])

def series(ds, v, idx):
    if v is None or v not in ds.variables:
        return None
    a = np.ma.filled(ds[v][:].astype("f8"), np.nan)
    if a.ndim > 1:
        a = a.reshape(a.shape[0], -1)[:, 0]
    a = a.flatten()[idx]
    return np.where(np.isfinite(a) & (a > -9990), a, np.nan)
