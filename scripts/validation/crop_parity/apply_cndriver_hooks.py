#!/usr/bin/env python3
"""Insert the four cropdump_write hooks into a stock ctsm5.3.012 CNDriverMod.F90.

Usage: apply_cndriver_hooks.py <stock CNDriverMod.F90> <out SourceMods/src.clm/CNDriverMod.F90>

Each anchor must match exactly once, so a CTSM version drift fails loudly
instead of silently dropping a hook.
"""
import sys

ARGS = (", &\n         cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst, &\n"
        "         cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, crop_inst, &\n"
        "         soilbiogeochem_state_inst, soilbiogeochem_nitrogenstate_inst, &\n"
        "         soilbiogeochem_nitrogenflux_inst, waterdiagnosticbulk_inst, &\n"
        "         temperature_inst, canopystate_inst)")


# In-model CNSoyfix unit sweep (no-op unless CROPDUMP_SOYFIX_SWEEP=1); restores
# every array it touches, so the trajectory is unchanged. use_crop is in scope
# in CNDriverNoLeaching (it gates CNNFert/CNSoyfix just above the hook).
SWEEP = ("    if (use_crop) call cropdump_soyfix_sweep(bounds, num_bgc_soilc, filter_bgc_soilc, &\n"
         "         num_bgc_vegp, filter_bgc_vegp, waterdiagnosticbulk_inst, crop_inst, &\n"
         "         cnveg_state_inst, cnveg_nitrogenflux_inst, soilbiogeochem_state_inst, &\n"
         "         soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)   ! CLM.jl CNSoyfix sweep\n")


def hook(label, indent="    "):
    return f"{indent}call cropdump_write(bounds, '{label}'{ARGS}   ! CLM.jl crop-tail parity\n"

src, dst = sys.argv[1], sys.argv[2]
text = open(src).read()

hooks = [
    # use statement
    ("  use perf_mod                        , only : t_startf, t_stopf\n",
     "  use perf_mod                        , only : t_startf, t_stopf\n"
     "  use cropdumpMod                     , only : cropdump_write, cropdump_soyfix_sweep   ! CLM.jl crop-tail parity\n"),
    # after CNSoyfix (end of the use_crop block, before CNMResp)
    ("    end if\n\n    call t_startf('CNMResp')\n",
     "    end if\n\n" + hook('after_soyfix') + SWEEP + "\n"
     "    call t_startf('CNMResp')\n"),
    # after calc_allometry
    ("       call calc_allometry(num_bgc_vegp, filter_bgc_vegp, &\n"
     "            cnveg_carbonflux_inst, cnveg_state_inst)\n"
     "       call t_stopf('cnalloc')\n    end if\n",
     "       call calc_allometry(num_bgc_vegp, filter_bgc_vegp, &\n"
     "            cnveg_carbonflux_inst, cnveg_state_inst)\n"
     "       call t_stopf('cnalloc')\n    end if\n"
     + hook('after_cropalloc')),
    # after calc_plant_nutrient_competition
    ("    call t_stopf('calc_plant_nutrient_competition')\n",
     "    call t_stopf('calc_plant_nutrient_competition')\n"
     + hook('after_cnalloc')),
    # end of CNDriverNoLeaching
    ("    end associate\n\n  end subroutine CNDriverNoLeaching\n",
     "    end associate\n\n" + hook('end_noleach') + "\n"
     "  end subroutine CNDriverNoLeaching\n"),
]

for old, new in hooks:
    n = text.count(old)
    if n != 1:
        sys.exit(f"anchor matched {n} times (need 1):\n{old}")
    text = text.replace(old, new)

open(dst, "w").write(text)
print(f"wrote {dst} with {len(hooks)} insertions")
