# `use_hydrstress` cannot be defaulted on yet — three blockers

**Status: the flip was attempted, measured, and REVERTED.** This is the honest
result of campaign #3 of the defaults programme (#256, row M3). The default
stays `false`; the reasons are below, with the numbers.

## What CTSM says

For `phys=clm5_0`, non-FATES, `configuration="clm"`, CTSM defaults
`use_hydrstress = .true.` — plant hydraulic stress is **on** in CLM5. CLM.jl
ships `false`, i.e. the BTRAN path, which gives different `gs`, transpiration
and GPP under water stress. So the port's default vegetation-water physics is
not CLM5's. That part of the audit stands.

## What happened when we flipped it

Full suite, `--check-bounds=yes`, against a baseline of
**26150 passed / 0 failed / 0 errored / 3 broken** (confirmed by two
independent runs):

| | baseline | `use_hydrstress = true` |
|---|---|---|
| passed | 26150 | 26088 |
| failed | 0 | 1 |
| **errored** | 0 | **18** |

Errors, not value mismatches — the code *throws*. Three distinct causes.

### Blocker 1 — PHS and FATES are mutually exclusive (3 errors)

`control.jl:102` raises `use_hydrstress and use_fates cannot both be true`,
which matches CTSM: the `.true.` default applies to the **non-FATES**
configuration only. The audit correctly labelled this row `COND-MISMATCH`.

*Fix when unblocked:* make the default conditional (`nothing` sentinel →
`!use_fates`), exactly as `use_bedrock` was handled in #252. This one is
straightforward and is not the reason the campaign stopped.

### Blocker 2 — the PHS path is not AD-differentiable (14 errors)

```
MethodError: no method matching CLM._Psn2Pft(
    ::Vector{Float64}, ::Vector{Float64}, ::Vector{Float64}, ::Vector{Float64},
    ::Vector{Float64}, ::Vector{Float64}, ::Vector{Dual{…}})
```

`_Psn2Pft{Vp}` (`photosynthesis.jl:1918`) declares **all seven fields with a
single shared type parameter**. Six arrive as `Vector{Float64}`; the seventh,
`lmr_intercept_atkin`, arrives as `Vector{Dual}` because it is the only field
built through the element type:

```julia
_lmratk = T.(params_inst.lmr_intercept_atkin)          # T = Dual under AD
lmr_intercept_atkin = slatop_pft isa Array ? _lmratk : convert(typeof(slatop_pft), _lmratk)
```

When `slatop_pft` is a plain `Array`, the ternary keeps the `Dual`-typed
`_lmratk`, so `Vp` cannot unify and construction fails. `lmr_intercept_atkin`
is a **constant parameter** — it has no derivative to carry, so promoting it to
`Dual` is wrong independently of this campaign.

This hits `test_ad_robustness` (all 6 scenarios), `test_parameter_recovery`,
`test_multisite_calibration` and `test_calibration` — i.e. **enabling PHS breaks
the entire AD/calibration surface**. That is a real port gap, not a
test-calibration artefact, and it is the blocker that stopped the flip.

*Fix:* build `lmr_intercept_atkin` to match the other pft vectors' type rather
than the AD element type (or give the struct per-field type parameters). Needs
its own PR + validation, since it touches the live photosynthesis path.

### Blocker 3 — unallocated PHS state on some fixtures (≈6 errors)

`BoundsError: attempt to access 0-element Vector{Float64} at index [2]` (and an
`Int64` variant) from `test_clm_driver` / `test_control`. The PHS branch reads
per-patch state that those fixtures never allocate, because they were built for
the BTRAN path. Needs a look at what `use_hydrstress=true` requires at init that
the non-PHS path does not.

## Why this is the right outcome

Nothing here was tuned to pass. The suite was allowed to fail honestly, and the
failures turned out to be **three real defects** — one config-constraint
mismatch and two genuine port gaps — rather than tests calibrated to the wrong
default. Making PHS the default *before* fixing blockers 2 and 3 would have
broken AD and calibration for every caller who trusts the default: exactly the
class of bug this whole programme exists to remove.

## Order of work to unblock

1. Fix `_Psn2Pft`'s type unification (blocker 2) — highest value on its own
   merit, since it makes the photosynthesis path AD-safe regardless of PHS.
2. Fix the PHS init/allocation gap (blocker 3).
3. Re-run this experiment; if clean, land the default as a **conditional**
   (`!use_fates`, blocker 1) with the same evidence discipline used for
   `use_bedrock` (#252) and `use_aquifer_layer` (#259).

`use_luna` (row M4) is untouched and remains queued behind this.
