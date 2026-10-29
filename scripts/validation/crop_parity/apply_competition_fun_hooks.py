#!/usr/bin/env python3
"""Insert the CNFUN parity dump hooks into a stock ctsm5.3.012 SoilBiogeochemCompetitionMod.F90.

Usage: apply_competition_fun_hooks.py <stock SoilBiogeochemCompetitionMod.F90> <out SourceMods/src.clm/...>

Adds cropdump_fun_write(... 'before_cnfun' ...) immediately before, and
cropdump_fun_write(... 'after_cnfun' ...) immediately after (i.e. after CNFUN's
p2c calls), BOTH CNFUN call sites (non-nitrif and nitrif branches). Each anchor
must match exactly twice, so a CTSM version drift fails loudly.
"""
import sys

ARGS = ("(bounds, '{label}', num_bgc_soilc, filter_bgc_soilc, num_bgc_vegp, filter_bgc_vegp, &\n"
        "                 waterstatebulk_inst, waterfluxbulk_inst, temperature_inst, soilstate_inst, &\n"
        "                 cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst, &\n"
        "                 cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &\n"
        "                 soilbiogeochem_nitrogenflux_inst, soilbiogeochem_carbonflux_inst, &\n"
        "                 canopystate_inst, soilbiogeochem_nitrogenstate_inst)   ! CLM.jl soyfix-fun parity\n")


def hook(label):
    return "            call cropdump_fun_write" + ARGS.format(label=label)


src, dst = sys.argv[1], sys.argv[2]
text = open(src).read()

edits = [
    ("    use perf_mod         , only : t_startf, t_stopf\n",
     "    use perf_mod         , only : t_startf, t_stopf\n"
     "    use cropdumpMod      , only : cropdump_fun_write   ! CLM.jl soyfix-fun parity\n", 1),
    ("            call t_startf( 'CNFUN' )\n",
     hook('before_cnfun') + "            call t_startf( 'CNFUN' )\n", 2),
    ("            call t_stopf( 'CNFUN' )\n",
     "            call t_stopf( 'CNFUN' )\n" + hook('after_cnfun'), 2),
]

for old, new, want in edits:
    n = text.count(old)
    if n != want:
        sys.exit(f"anchor matched {n} times (need {want}):\n{old}")
    text = text.replace(old, new)

open(dst, "w").write(text)
print(f"wrote {dst}")
