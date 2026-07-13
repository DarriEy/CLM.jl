#!/usr/bin/env python3
"""
instrument.py — inject the schema-v2, phase-tagged FATES parity dump into a CTSM
source tree, writing the patched files into a case's SourceMods/src.clm.

Three files are patched (all exact-anchor replacements; every anchor is asserted,
so an upstream source change fails loudly instead of silently producing an
un-instrumented executable):

  clmfates_interfaceMod.F90
      * `use FatesParityDumpMod`
      * a new public type-bound procedure `parity_dump(nc, phase)` that dumps every
        site of a clump at (parity_nstep, phase)
      * init_coldstart -> phase 0 dump (the nstep=0 ground truth of PR #213)
      * dynamics_driv  -> set `parity_site = s` / `parity_daily = .true.` around the
        per-site ed_ecosystem_dynamics + ed_update_site calls, so the EDMainMod
        sub-phase hooks know which site they are in and do NOT fire during the
        cold-start / restart calls to ed_update_site.

  clm_driver.F90
      * `use FatesParityDumpMod, only : parity_nstep`
      * in the use_fates block, immediately before the `if (is_beg_curr_day())`
        FATES-dynamics branch: set parity_nstep and emit the phase-1 ("fast") dump.
        That point is AFTER the photosynthesis/flux solve, the FATES running-mean
        update and the hifrq history fill, and BEFORE any daily dynamics — i.e. it
        isolates the fast (sub-daily) FATES thread.

  EDMainMod.F90
      * `use FatesParityDumpMod`
      * dumps at every phase boundary of ed_ecosystem_dynamics + ed_update_site
        (phases 10..18; see FatesParityDumpMod.F90 for the code table).

usage: instrument.py <CLMROOT> <SOURCEMODS_DIR>
"""
import os
import shutil
import sys

CLMROOT, SM = sys.argv[1], sys.argv[2]
HERE = os.path.dirname(os.path.abspath(__file__))
os.makedirs(SM, exist_ok=True)


def patch(rel, edits, dst=None):
    """Apply (anchor, replacement) exact-text edits to CLMROOT/rel -> SM/<basename>."""
    src = os.path.join(CLMROOT, rel)
    text = open(src).read()
    for anchor, repl in edits:
        assert text.count(anchor) == 1, \
            f"{rel}: anchor found {text.count(anchor)}x (want 1):\n{anchor[:160]}"
        text = text.replace(anchor, repl)
    out = os.path.join(SM, dst or os.path.basename(src))
    open(out, "w").write(text)
    print(f"[ok] instrumented {os.path.basename(src)} -> {out}")


# --------------------------------------------------------------------------------
# 0. the dump module itself
# --------------------------------------------------------------------------------
shutil.copy(os.path.join(HERE, "FatesParityDumpMod.F90"),
            os.path.join(SM, "FatesParityDumpMod.F90"))
print("[ok] FatesParityDumpMod.F90 ->", SM)

# --------------------------------------------------------------------------------
# 1. clmfates_interfaceMod.F90
# --------------------------------------------------------------------------------
IFACE_USE_ANCHOR = "   use EDInitMod             , only : set_site_properties\n"
IFACE_USE = IFACE_USE_ANCHOR + (
    "   ! ---- FATES parity instrumentation (schema-v2 phase-tagged pdump) ----------\n"
    "   use FatesParityDumpMod    , only : fates_parity_dump, parity_nstep\n"
    "   use FatesParityDumpMod    , only : parity_site, parity_daily\n"
)

IFACE_PROC_ANCHOR = "      procedure, public :: init_coldstart\n"
IFACE_PROC = IFACE_PROC_ANCHOR + "      procedure, public :: parity_dump   ! FATES parity instrumentation\n"

# cold start (phase 0) — fires inside init_coldstart, right before its t_stopf.
IFACE_COLD_ANCHOR = "     call t_stopf('fates_initcoldstart')\n"
IFACE_COLD = (
    "     ! ---- FATES parity instrumentation: cold-start (phase 0) ------------------\n"
    "     do nc = 1, nclumps\n"
    "        call this%parity_dump(nc, 0)\n"
    "     end do\n"
    "     ! --------------------------------------------------------------------------\n"
) + IFACE_COLD_ANCHOR

# daily dynamics: tag the site ordinal + arm the EDMainMod sub-phase hooks.
IFACE_DYN_ANCHOR = """      do s = 1,this%fates(nc)%nsites

            call ed_ecosystem_dynamics(this%fates(nc)%sites(s),    &
                  this%fates(nc)%bc_in(s), &
                  this%fates(nc)%bc_out(s))

            call ed_update_site(this%fates(nc)%sites(s), &
                  this%fates(nc)%bc_in(s), &
                  this%fates(nc)%bc_out(s), &
                  is_restarting = .false.)
      enddo
"""
IFACE_DYN = """      do s = 1,this%fates(nc)%nsites

            ! ---- FATES parity instrumentation: arm the EDMainMod sub-phase hooks ----
            parity_site  = s
            parity_daily = .true.

            call ed_ecosystem_dynamics(this%fates(nc)%sites(s),    &
                  this%fates(nc)%bc_in(s), &
                  this%fates(nc)%bc_out(s))

            call ed_update_site(this%fates(nc)%sites(s), &
                  this%fates(nc)%bc_in(s), &
                  this%fates(nc)%bc_out(s), &
                  is_restarting = .false.)

            parity_daily = .false.
            ! ------------------------------------------------------------------------
      enddo
"""

# the new type-bound procedure body, appended just before `end module`.
IFACE_END_ANCHOR = "end module CLMFatesInterfaceMod\n"
IFACE_END = """
  ! ====================================================================================
  ! FATES parity instrumentation: dump every site of clump `nc` at (parity_nstep, phase).
  ! ====================================================================================
  subroutine parity_dump(this, nc, phase)
    class(hlm_fates_interface_type), intent(inout) :: this
    integer, intent(in) :: nc
    integer, intent(in) :: phase
    if (this%fates(nc)%nsites > 0) then
       call fates_parity_dump(this%fates(nc)%sites, this%fates(nc)%bc_in, &
            this%fates(nc)%nsites, parity_nstep, phase)
    end if
  end subroutine parity_dump

end module CLMFatesInterfaceMod
"""

patch("src/utils/clmfates_interfaceMod.F90", [
    (IFACE_USE_ANCHOR, IFACE_USE),
    (IFACE_PROC_ANCHOR, IFACE_PROC),
    (IFACE_COLD_ANCHOR, IFACE_COLD),
    (IFACE_DYN_ANCHOR, IFACE_DYN),
    (IFACE_END_ANCHOR, IFACE_END),
])

# --------------------------------------------------------------------------------
# 2. clm_driver.F90 — the per-timestep "fast" phase (phase 1)
# --------------------------------------------------------------------------------
DRV_USE_ANCHOR = "  use clm_time_manager       , only : get_nstep, is_beg_curr_day, is_beg_curr_year\n"
DRV_USE = DRV_USE_ANCHOR + \
    "  use FatesParityDumpMod     , only : parity_nstep   ! FATES parity instrumentation\n"

# Emitted at the top of the FATES dynamics branch: after WrapUpdateFatesRmean and
# wrap_update_hifrq_hist (the last per-timestep FATES updates), before the daily step.
DRV_FAST_ANCHOR = """          if( is_beg_curr_day() ) then

             ! --------------------------------------------------------------------------
             ! This is the main call to FATES dynamics
             ! --------------------------------------------------------------------------
"""
DRV_FAST = """          ! ---- FATES parity instrumentation: fast (sub-daily) phase = 1 -----------
          parity_nstep = get_nstep()
          call clm_fates%parity_dump(nc, 1)
          ! -------------------------------------------------------------------------

""" + DRV_FAST_ANCHOR

patch("src/main/clm_driver.F90", [
    (DRV_USE_ANCHOR, DRV_USE),
    (DRV_FAST_ANCHOR, DRV_FAST),
])

# --------------------------------------------------------------------------------
# 3. EDMainMod.F90 — the daily sub-phases (10..18)
# --------------------------------------------------------------------------------
def dump(phase, indent="    "):
    return f"{indent}call fates_parity_dump_site(currentSite, bc_in, parity_site, parity_nstep, {phase})\n"


ED_USE_ANCHOR = "  implicit none\n"
ED_USE = (
    "  ! ---- FATES parity instrumentation (schema-v2 phase-tagged pdump) ------------\n"
    "  use FatesParityDumpMod , only : fates_parity_dump_site, parity_nstep, parity_site\n"
    "\n"
) + ED_USE_ANCHOR

ED_EDITS = [
    (ED_USE_ANCHOR, ED_USE),
    # 10 dyn_in — the daily step's input state
    ("    call TotalBalanceCheck(currentSite, 0)\n",
     "    call TotalBalanceCheck(currentSite, 0)\n" + dump(10)),
    # 11 phenology
    ("      end if ! SP phenology\n    end if\n",
     "      end if ! SP phenology\n    end if\n" + dump(11)),
    # 12 disturbance rates (+ fire model)
    ("       call disturbance_rates(currentSite, bc_in)\n",
     "       call disturbance_rates(currentSite, bc_in)\n" + dump(12, "       ")),
    # 13 growth / allocation / PRT / mortality derivative
    ("       call ed_integrate_state_variables(currentSite, bc_in, bc_out )\n",
     "       call ed_integrate_state_variables(currentSite, bc_in, bc_out )\n" + dump(13, "       ")),
    # 14 recruitment
    ("       call TotalBalanceCheck(currentSite,1)\n",
     "       call TotalBalanceCheck(currentSite,1)\n" + dump(14, "       ")),
    # 15 cohort sort / terminate / fuse
    ("    call TotalBalanceCheck(currentSite,2)\n",
     "    call TotalBalanceCheck(currentSite,2)\n" + dump(15)),
    # 16 patch spawn (disturbance)
    ("       call TotalBalanceCheck(currentSite,3)\n",
     "       call TotalBalanceCheck(currentSite,3)\n" + dump(16, "       ")),
    # 17 patch fuse/terminate == end of ed_ecosystem_dynamics
    ("    call TotalBalanceCheck(currentSite,5)\n",
     "    call TotalBalanceCheck(currentSite,5)\n" + dump(17)),
    # 18 end of ed_update_site (canopy spread/structure/trim) == end of the daily step
    ("  end subroutine ed_update_site\n",
     dump(18) + "\n  end subroutine ed_update_site\n"),
]
patch("src/fates/main/EDMainMod.F90", ED_EDITS)

print("[ok] instrumentation complete")
