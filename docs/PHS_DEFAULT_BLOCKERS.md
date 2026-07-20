# `use_hydrstress` — the three blockers, and how they closed

**Status: the conditional default is IN.** `use_hydrstress` now takes CTSM's
conditional namelist default (`nothing` → `!use_fates`), closing row M3 of
`docs/DRIVER_DEFAULTS_AUDIT.md` — the last live-physics row of the defaults
campaign.

The first attempt (#261) flipped it, measured **18 errors**, and REVERTED. That
was the right call, and the documented blocker report was the right deliverable:
the 18 errors turned out to be **three real defects**, not tests calibrated to
the wrong default. This document is that report, updated with how each closed.

## What CTSM says

For `phys=clm5_0`, non-FATES, `configuration="clm"`, CTSM defaults
`use_hydrstress = .true.` — plant hydraulic stress is **on** in CLM5. CLM.jl
shipped `false`, which is CTSM's **code fallback** (`clm_varctl.F90`), the value
the namelist overrides. Same root cause as #252 (`use_bedrock`), #259
(`use_aquifer_layer`) and #267 (`use_luna`): the audit found fourteen such
namelist-vs-code-fallback disagreements, and CLM.jl had copied the fallback in
every one.

## Blocker 1 — PHS and FATES are mutually exclusive — **CLOSED**

`control.jl:102` raises `use_hydrstress and use_fates cannot both be true`,
matching CTSM: the `.true.` default applies to the **non-FATES** configuration
only. So a bare `true` would have been just as wrong as `false`.

**Fix:** the `nothing` sentinel, resolving to `!use_fates`, at the three entry
points that carry the default — `clm_initialize!`, `clm_run!` and
`CLMDriverConfig`. `varctl.use_hydrstress` keeps its `false` **struct** default,
which mirrors CTSM's code fallback: it is the DERIVATION that carries the
namelist default, not the struct (the #265 principle, as for
`create_crop_landunit`).

Pinned in `test/test_driver_defaults_audit.jl` — the sentinel **and** the
resolved values (`CLMDriverConfig().use_hydrstress === true`, `use_fates=true` →
`false`, explicit `false` → `false`), so the pin cannot go vacuously green.

## Blocker 2 — the PHS path was not AD-differentiable — **CLOSED**

Originally 14 errors:

```
MethodError: no method matching CLM._Psn2Pft(
    ::Vector{Float64}, …, ::Vector{Dual{…}})
```

`_Psn2Pft{Vp}` (`photosynthesis.jl:1981`) declares **all eight pft fields with a
single shared type parameter**. If any one field arrives as `Vector{Dual}` while
the rest are `Vector{Float64}`, `Vp` cannot unify and construction throws.

#262 fixed the first offender, `lmr_intercept_atkin`. **This campaign found the
same bug one field over**, and it was still fatal:

```julia
# canopy_fluxes.jl — before
crop_pft_s = isempty(crop_pft) ?
    _to_backend_like(ivt_vec, FT, zeros(FT, length(medlynslope_pft))) : crop_pft
```

`FT` is the AD element type. When `crop_pft` is empty — which is the
**default, non-CN configuration**, i.e. exactly the AD/calibration path — the
substitute was built as `Vector{Dual}`, and PHS pass-2 construction threw. Fix:

```julia
crop_pft_s = isempty(crop_pft) ? zero(medlynslope_pft) : crop_pft
```

`zero(sibling)` keeps both the eltype (`Float64` — `crop_pft` is a CONSTANT
per-PFT parameter carrying no derivative) and the container (a device array
stays a device array, so the GPU path is unperturbed).

**Measured:** `test_ad_robustness` with PHS on goes from **6 errors / 0 passes**
to **96 passes / 0 errors**, with AD-vs-finite-difference relative error
**0.0%** on `d(LH)/d(T)`, `d(SH)/d(T)` and `d(T_grnd)/d(T)` in all six climate
scenarios. 96 is exactly the PHS-**off** baseline for that file, so this is a
RESTORATION, not an inflated count: each scenario previously threw at the first
PHS call, so its six assertions never ran. The suite total is unaffected by
this file.

A sweep for the same bug class found LUNA's `_LunaPft` promotes *all* its pft
fields to `FT` — wasteful, but it unifies, so it is not broken. `crop_pft` was
the unique unification breaker.

## Blocker 3 — unallocated PHS state on some fixtures — **CLOSED**

Reproduced exactly: **6 `BoundsError`s** in `test_clm_driver`, every one inside
`get_froot_carbon_patch`. Two distinct causes, and only one was the fixture's
fault.

**3a — a real port fragility (2 of the 6).** `clm_drv_core!` decided CN-vs-SP
*itself* and then called the facade, which decides **again** on its own
`veg.config.use_cn`:

```julia
config.use_cn ? get_froot_carbon_patch(veg, bc_patch)      # no PFT inputs
              : get_froot_carbon_patch(veg, bc_patch; tlai, slatop, froot_leaf, ivt)
```

When the two configs disagree, the driver's CN call lands in the facade's SP
branch with all four PFT inputs defaulted **empty** → `BoundsError` on
`ivtH[p]`. Fortran's `CNVegetationFacade` owns that choice; the port duplicated
it. Fixed by passing the SP inputs unconditionally and letting the facade
branch. The CN path is unchanged (the kwargs are read only by the SP branch).

**3b — fixture incoherence (4 of the 6), fixed the #271 way.** The
`make_driver_data` fixture hardcodes a BTRAN-era pftcon list; with PHS now the
default it runs the PHS path, which reads three parameters the fixture never
set. The fix is to give the fixture **real parameter values** — not to weaken
#263's guard, which is new and was doing its job:

| parameter | value | source |
|---|---|---|
| `froot_leaf` | 2.8143 | `clm5_params.nc`, all veg PFTs |
| `root_radius` | `ROOT_RADIUS_PARAM` = 0.29e-3 m | `pftcon.jl:107` (constant, not file-read) |
| `root_density` | `ROOT_DENSITY_PARAM` = 0.31e6 | `pftcon.jl:106` |

**A silent NaN found on the way.** `set_params_for_testing!` (the port of
Fortran's `setParamsForTesting`) sets only `ck` and `psi50`, while
`photo_params_init!` leaves `krmax`, `kmax` and `theta_cj` at **NaN**. Nothing
read them on the BTRAN default. The PHS path does — and NaN there is a *silent
wrong answer*, not an error: no `BoundsError`, no guard, just a quiet NaN
propagating into the PHS conductance chain. The fixture now seeds the real
`clm5_params.nc` values (`krmax` 7.943e-10, `kmax` 2e-8, `theta_cj` 0.98).

That is the `conservation-is-not-accuracy` lesson in miniature: the loud failure
was the easy one, and the thing worth finding was the quiet one next to it.

## Why this landed and #261 did not

Nothing was tuned to pass. Every one of the 18 original errors was traced to a
defect and fixed at the defect:

| blocker | original errors | disposition |
|---|---|---|
| 1 — PHS+FATES exclusivity | 3 | conditional default (`nothing` → `!use_fates`) |
| 2 — `_Psn2Pft` AD unification | 14 | **real port gap** — `crop_pft` built at the AD element type (#262, one field over) |
| 3 — unallocated PHS state | ~6 | **1 real port fragility** (duplicated CN branch) + fixture given real parameters |

Errors mean a port gap; value mismatches mean a test calibrated to the old
default. These were errors, and they were gaps.

## Full-suite measurement

See the PR body for the baseline-vs-after diff with every moved test named.
