"""stage_multiyear.py <rundir> <nyears>

Turn a staged single-point FATES run dir (a copy of the fates_parity_1pt run dir:
lnd_in, nuopc.runconfig, datm_*, fates_params_default.nc) into the MULTI-YEAR
daily-history reference run used by docs/FATES_MULTIYEAR_ENERGY.md "Defect 2":
hist_nhtfrq=-24 with the FATES demography fields, stop_option=nyears.
The source run dir must carry `hist_nhtfrq = 1 / hist_mfilt = 1000 / hist_fincl1`
lines (as the fates_bgc_1pt template does). Then run cesm.exe in the dir and
extract with fates_multiyear_extract.jl.
"""
import re, sys
R = sys.argv[1]
nyears = sys.argv[2]
p = f"{R}/lnd_in"; s = open(p).read()
s = s.replace("runs/fates_bgc_1pt_v2/", "runs/" + R.rstrip("/").split("/")[-1] + "/")
fields = ['FATES_NCOHORTS','FATES_NPATCHES','FATES_RECRUITMENT_PF','FATES_NPLANT_PF','FATES_NPLANT_SZPF',
          'FATES_DDBH_SZPF','FATES_MORTALITY_PF','FATES_NPP','FATES_GPP','FATES_NPP_PF','FATES_VEGC',
          'FATES_VEGC_PF','FATES_STRUCTC','FATES_LEAFC','FATES_STOREC','FATES_SEED_BANK','FATES_SEED_ALLOC',
          'FATES_STEM_ALLOC','FATES_LEAF_ALLOC','FATES_FROOT_ALLOC','FATES_STORE_ALLOC','FATES_CROOT_ALLOC',
          'FATES_LAI','FATES_SEED_GERM_EL','FATES_SEEDS_IN_LOCAL','TSOI','H2OSOI','FSH','TV']
new = (" hist_nhtfrq = -24\n hist_mfilt = 3000\n hist_empty_htapes = .true.\n hist_fincl1 = "
       + ",".join(f"'{f}'" for f in fields) + "\n")
s, n = re.subn(r" hist_nhtfrq = 1\n hist_mfilt = 1000\n hist_fincl1 = [^\n]*\n", new, s)
assert n == 1
open(p, "w").write(s)
p = f"{R}/nuopc.runconfig"; s = open(p).read()
for a, b in [("case_name = fates_bgc_1pt", "case_name = fates_recruit_my"),
             ("stop_n = 240", f"stop_n = {nyears}"), ("stop_option = nsteps", "stop_option = nyears")]:
    assert a in s, a
    s = s.replace(a, b)
open(p, "w").write(s)
print("ok")
