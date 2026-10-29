"""stage_niche_mode.py <src_rundir> <dst_rundir> <ndays> [key=value ...]

Stage a FATES "research mode" reference run (scripts/validation/fates_fortran_parity/README.md "RESEARCH MODES") from a
staged single-point FATES run dir (the fates_parity_1pt case, e.g. runs/fates_bgc_1pt_v2:
lnd_in, nuopc.runconfig, datm_*, fates_params_default.nc).

Every key=value pair is a lnd_in namelist override (the rewrite setup_case.sh does at
the lnd_in step), e.g. use_fates_tree_damage=.true. or fates_seeddisp_cadence=1.
Two keys are special:
  paramfile=<path.nc>   use this FATES parameter file instead of the source's copy
  fields=A,B,C          extra daily history fields (on top of the common set)
  restart_ndays=N       write a restart after N days (for restart round-trips)
  mesh=<path.nc>        use this ESMF mesh for atm/lnd/mask instead of the site mesh
                        (inventory init: a copy with longitudes +360, see the docs)

The history is DAILY (hist_nhtfrq=-24) with the demographic set used by
fates_niche_compare.jl. Run cesm.exe in <dst_rundir> afterwards.
"""
import os, re, shutil, sys

src, dst, ndays = sys.argv[1], sys.argv[2], int(sys.argv[3])
over, extra, paramfile, rst_ndays, mesh = {}, [], None, None, None
for kv in sys.argv[4:]:
    k, v = kv.split("=", 1)
    if k == "fields":
        extra = [f for f in v.split(",") if f]
    elif k == "paramfile":
        paramfile = v
    elif k == "restart_ndays":
        rst_ndays = int(v)
    elif k == "mesh":
        mesh = v
    else:
        over[k] = v

name = os.path.basename(dst.rstrip("/"))
os.makedirs(os.path.join(dst, "timing", "checkpoints"), exist_ok=True)
for f in ["nuopc.runseq", "fd.yaml", "datm_in", "datm.streams.xml", "drv_in", "drv_flds_in",
          "lnd_in", "nuopc.runconfig"]:
    shutil.copy(os.path.join(src, f), os.path.join(dst, f))
if mesh is not None:
    for f in ["nuopc.runconfig", "datm_in", "datm.streams.xml"]:
        q = os.path.join(dst, f)
        t = open(q).read()
        t, n = re.subn(r"/[^\s\"'<>]*/esmf_mesh\.nc", mesh, t)
        open(q, "w").write(t)
        print("  mesh ->", f, n)
fp = os.path.join(dst, "fates_params_default.nc")
shutil.copy(paramfile or os.path.join(src, "fates_params_default.nc"), fp)

fields = ['FATES_NCOHORTS', 'FATES_NPATCHES', 'FATES_RECRUITMENT_PF', 'FATES_NPLANT_PF',
          'FATES_NPLANT_SZPF', 'FATES_MORTALITY_PF', 'FATES_NPP', 'FATES_GPP', 'FATES_VEGC',
          'FATES_VEGC_PF', 'FATES_STRUCTC', 'FATES_LEAFC', 'FATES_STOREC', 'FATES_SEED_BANK',
          'FATES_LAI', 'FATES_SEEDS_IN_LOCAL', 'FATES_BASALAREA_SZPF', 'TSOI', 'H2OSOI', 'TV',
          'FSH', 'HR'] + extra

p = os.path.join(dst, "lnd_in")
s = open(p).read().replace(os.path.abspath(src).rstrip("/") + "/", os.path.abspath(dst).rstrip("/") + "/")
s = re.sub(r"fates_paramfile = '[^']*'", f"fates_paramfile = '{fp}'", s)
new = (" hist_nhtfrq = -24\n hist_mfilt = 3000\n hist_empty_htapes = .true.\n hist_fincl1 = "
       + ",".join(f"'{f}'" for f in fields) + "\n")
s, n = re.subn(r" hist_nhtfrq = [^\n]*\n hist_mfilt = [^\n]*\n(?: hist_empty_htapes = [^\n]*\n)? hist_fincl1 = [^\n]*\n",
               new, s)
assert n == 1, "history block not found"
for k, v in over.items():
    s, n = re.subn(rf"(?m)^(\s*){k}\s*=.*$", rf"\g<1>{k} = {v}", s)
    if n == 0:   # not in the generated namelist: append to &clm_inparm
        s, n = re.subn(r"(&clm_inparm\n)", rf"\g<1> {k} = {v}\n", s, count=1)
    assert n == 1, k
open(p, "w").write(s)

p = os.path.join(dst, "nuopc.runconfig")
s = open(p).read().replace(os.path.abspath(src).rstrip("/"), os.path.abspath(dst).rstrip("/"))
s = re.sub(r"case_name = \S+", f"case_name = {name}", s)
s = re.sub(r"stop_option = \S+", "stop_option = ndays", s)
s = re.sub(r"stop_n = \S+", f"stop_n = {ndays}", s)
if rst_ndays is not None:
    s = re.sub(r"restart_option = \S+", "restart_option = ndays", s)
    s = re.sub(r"restart_n = \S+", f"restart_n = {rst_ndays}", s)
open(p, "w").write(s)
print("ok", dst, over)
