# EDAccumulateFluxesMod.jl
# Julia port of FATES src/fates/biogeophys/EDAccumulateFluxesMod.F90
#
# Accumulates NPP, GPP and respiration of each cohort over the course of each
# 24-hour period. The per-timestep (_tstep) fluxes are computed in EDPhotosynthesis;
# this routine sums them into the daily accumulators (_acc). It cannot live in
# EDPhotosynthesis because that is a leaf-layer loop and would erroneously add
# these things up multiple times.  (Rosie Fisher, March 2014.)
#
# Translation notes:
#   * fates_r8 -> Float64. The Fortran multi-site arrays sites(nsites) /
#     bc_in(nsites) / bc_out(nsites) become Julia Vectors indexed 1:nsites.
#   * Patch / cohort linked-list pointer walks become the Union{...,Nothing}
#     traversals used by the merged FATES types (oldest_patch -> younger ;
#     shortest -> taller).
#   * `nocomp_bareground` is the FatesConstantsMod const; `filter_photo_pa` is on
#     bc_in (the per-patch photosynthesis filter flag — value 3 means run).
#   * The c13 GPP-weighted mean, the year_net_uptake 999-sentinel reset, and the
#     leaf-layer accumulation are reproduced verbatim.
#   * bc_out is in the signature for Fortran fidelity (unused here, as in Fortran).

"""
    AccumulateFluxes_ED(nsites, sites, bc_in, bc_out, dt_time)

Accumulate the per-cohort photosynthesis / respiration fluxes from per-timestep
(`*_tstep`, KgC/indiv/timestep) into the daily accumulators (`*_acc`,
KgC/indiv/day) for every site.

For each site, walk the patches (oldest -> younger), skip bare-ground patches,
increment the FATES-patch index `ifp`, and — only when `bc_in.filter_photo_pa[ifp] == 3`
— walk that patch's cohorts (shortest -> taller) accumulating:
  * `npp_acc  += npp_tstep`, `gpp_acc += gpp_tstep`, `resp_acc += resp_tstep`,
  * `sym_nfix_daily += sym_nfix_tstep`,
  * `c13disc_acc` as the GPP-weighted mean of `c13disc_clm` (zero if total GPP is 0),
  * `year_net_uptake[iv] += ts_net_uptake[iv]` over leaf layers `1:nv`, resetting
    the 999-sentinel ("there were leaves in this layer this year") to 0 first.

`dt_time` (the timestep interval) is accepted for Fortran-signature fidelity;
`bc_out` is likewise accepted and unused, matching the Fortran.
"""
function AccumulateFluxes_ED(nsites::Integer,
                             sites::Vector{ed_site_type},
                             bc_in::Vector{bc_in_type},
                             bc_out::Vector{bc_out_type},
                             dt_time::Real)

    for s in 1:nsites

        ifp = 0

        # Note: Do not attempt to accumulate or log any heterotrophic respiration
        # fluxes from the HLM here. It is likely this has not been calculated yet
        # (ELM/CLM).

        cpatch = sites[s].oldest_patch
        while cpatch !== nothing
            if cpatch.nocomp_pft_label != nocomp_bareground
                ifp = ifp + 1

                if bc_in[s].filter_photo_pa[ifp] == 3
                    ccohort = cpatch.shortest
                    while ccohort !== nothing

                        # Accumulate fluxes from hourly to daily values.
                        # _tstep fluxes are KgC/indiv/timestep, _acc are KgC/indiv/day.
                        # The daily accumulator starts from 0 and sums finite per-step
                        # fluxes; guard the sums so an unset (NaN) accumulator or a cohort
                        # whose instantaneous flux was not computed this step (NaN _tstep)
                        # contributes 0 rather than poisoning the whole-site mass balance.
                        _fz(x) = isfinite(x) ? x : 0.0
                        ccohort.npp_acc  = _fz(ccohort.npp_acc)  + _fz(ccohort.npp_tstep)
                        ccohort.gpp_acc  = _fz(ccohort.gpp_acc)  + _fz(ccohort.gpp_tstep)
                        ccohort.resp_acc = _fz(ccohort.resp_acc) + _fz(ccohort.resp_tstep)

                        ccohort.sym_nfix_daily = _fz(ccohort.sym_nfix_daily) + _fz(ccohort.sym_nfix_tstep)

                        # weighted mean of D13C by gpp
                        if (ccohort.gpp_acc + ccohort.gpp_tstep) == 0.0
                            ccohort.c13disc_acc = 0.0
                        else
                            ccohort.c13disc_acc =
                                ((ccohort.c13disc_acc * ccohort.gpp_acc) +
                                 (ccohort.c13disc_clm * ccohort.gpp_tstep)) /
                                (ccohort.gpp_acc + ccohort.gpp_tstep)
                        end

                        for iv in 1:ccohort.nv
                            # note that there were leaves in this layer this year.
                            if ccohort.year_net_uptake[iv] == 999.0
                                ccohort.year_net_uptake[iv] = 0.0
                            end
                            ccohort.year_net_uptake[iv] =
                                ccohort.year_net_uptake[iv] + ccohort.ts_net_uptake[iv]
                        end

                        ccohort = ccohort.taller
                    end  # while(associated(ccohort))
                end
            end  # not bare ground
            cpatch = cpatch.younger
        end  # while(associated(cpatch))
    end

    return nothing
end
