# test_fates_prtlossfluxes.jl
# Tests for FATES Batch 3: PRTLossFluxesMod — the PARTEH loss-flux routines that
# operate on the generic prt_vartypes state.
#
# Strategy: build a small CNP hypothesis descriptor (leaf, fineroot, sapwood,
# structure, storage, repro) x (carbon12, nitrogen, phosphorus), install it as
# the module global, allocate prt_params with known retranslocation/turnover/
# stoichiometry values, then drive each loss routine on a synthetic plant with
# known masses and assert:
#   * mass conservation (mass removed from a pool == leaves-plant + retranslocated)
#   * the leaves-plant part lands in the turnover/burned/damaged diagnostic
#   * the retranslocated part lands in the storage pool
#   * correct retranslocation partitioning (carbon never retranslocates; N/P use
#     turnover_{nitr,phos}_retrans).

using Test
using CLM

# --- helper: build + install a CNP prt_global descriptor -----------------------
# organ order in this test hypothesis: leaf, fnrt, sapw, store, repro, struct.
# Each organ holds 3 elements (C12, N, P). Leaf carbon uses 2 age positions.
function _build_cnp_global(; nleaf_pos::Int = 2)
    g = CLM.prt_global_type()
    CLM.ZeroGlobal!(g)
    g.hyp_name = "test_cnp"
    g.hyp_id   = CLM.prt_cnp_flex_allom_hyp   # == 2  (<= 2 path)

    organs   = [CLM.leaf_organ, CLM.fnrt_organ, CLM.sapw_organ,
                CLM.store_organ, CLM.repro_organ, CLM.struct_organ]
    elements = [CLM.carbon12_element, CLM.nitrogen_element, CLM.phosphorus_element]

    g.num_vars     = length(organs) * length(elements)
    g.num_bc_in    = 0
    g.num_bc_out   = 0
    g.num_bc_inout = 0
    g.state_descriptor = [CLM.state_descriptor_type() for _ in 1:g.num_vars]

    i_var = 0
    for org in organs
        for el in elements
            i_var += 1
            # only leaf-carbon gets multiple age positions; everything else is 1
            npos = (org == CLM.leaf_organ) ? nleaf_pos : 1
            CLM.RegisterVarInGlobal!(g, i_var, "v$(org)_$(el)", "s$(org)_$(el)",
                                     org, el, npos)
        end
    end
    return g
end

# total plant mass (sum of val over all variables/positions) for one element
function _plant_mass(prt, g, element_id)
    m = 0.0
    for i_var in 1:g.num_vars
        if g.state_descriptor[i_var].element_id == element_id
            m += sum(prt.variables[i_var].val)
        end
    end
    return m
end

@testset "FATES PRTLossFluxesMod" begin

    old_global = CLM.prt_global[]
    g = _build_cnp_global()
    CLM.prt_global[] = g

    # ---- allocate + fill prt_params for a single PFT -----------------------
    # norgan here is the number of param-file organs; we make organ_param_id an
    # identity map over the 6 PRT organs so organ_param_id[organ] == organ.
    npft = 1; norgan = CLM.num_organ_types; nleafage = 2
    CLM.allocate_prt_params!(CLM.prt_params, npft, norgan, nleafage;
                             nall_organs = CLM.num_organ_types)

    ipft = 1
    p = CLM.prt_params
    # identity organ_param_id so every organ has a (positive) stoichiometry slot
    for o in 1:CLM.num_organ_types
        p.organ_param_id[o] = o
    end
    p.woody[ipft]     = 0          # non-woody so all organs allowed in decid drop
    p.evergreen[ipft] = CLM.itrue  # so maintenance leaf turnover is active

    # retranslocation fractions (PFT x organ): C ignores these; N/P use them.
    retr_n = 0.40
    retr_p = 0.25
    for o in 1:CLM.num_organ_types
        p.turnover_nitr_retrans[ipft, o] = retr_n
        p.turnover_phos_retrans[ipft, o] = retr_p
    end
    # stoichiometry (used by PhenologyFlush nutrient targets)
    for o in 1:CLM.num_organ_types
        p.nitr_stoich_p1[ipft, o] = 0.05
        p.phos_stoich_p1[ipft, o] = 0.005
    end
    # turnover timescales [yr]
    p.branch_long[ipft] = 50.0
    p.root_long[ipft]   = 1.0
    p.leaf_long[ipft, :]        .= 1.5
    p.leaf_long_ustory[ipft, :] .= 2.0
    p.senleaf_long_fdrought[ipft] = 0.5

    # convenience: set every pool to a known starting value, then ZeroRates!
    function fresh_plant(; leaf_c=(2.0, 1.0), default=1.0)
        prt = CLM.prt_vartypes()
        CLM.InitPRTVartype!(prt)
        for i_var in 1:g.num_vars
            org = g.state_descriptor[i_var].organ_id
            el  = g.state_descriptor[i_var].element_id
            for ip in 1:g.state_descriptor[i_var].num_pos
                if org == CLM.leaf_organ && el == CLM.carbon12_element
                    prt.variables[i_var].val[ip] = leaf_c[ip]
                else
                    prt.variables[i_var].val[ip] = default
                end
            end
        end
        CLM.ZeroRates!(prt)
        return prt
    end

    try
        # ===============================================================
        @testset "PRTBurnLosses!: no retrans, all mass leaves the plant" begin
            prt = fresh_plant()
            organ = CLM.fnrt_organ
            frac = 0.3
            # snapshot fnrt masses by element
            before = Dict(el => CLM.GetState(prt, organ, el)
                          for el in (CLM.carbon12_element, CLM.nitrogen_element,
                                     CLM.phosphorus_element))
            CLM.PRTBurnLosses!(prt, organ, frac)
            for (el, b) in before
                @test CLM.GetState(prt, organ, el) ≈ (1 - frac) * b
                @test CLM.GetBurned(prt, organ, el)  ≈ frac * b
                # removed mass == burned mass exactly (no retranslocation)
                @test (b - CLM.GetState(prt, organ, el)) ≈ CLM.GetBurned(prt, organ, el)
            end
            # mass conservation invariant holds with burned as the loss term
            @test CLM.CheckMassConservation(prt, ipft, 1) === nothing
        end

        # ===============================================================
        @testset "PRTDamageLosses!: no retrans, tracked in damaged" begin
            prt = fresh_plant()
            organ = CLM.sapw_organ
            frac = 0.5
            b = CLM.GetState(prt, organ, CLM.nitrogen_element)
            CLM.PRTDamageLosses!(prt, organ, frac)
            @test CLM.GetState(prt, organ, CLM.nitrogen_element) ≈ (1 - frac) * b
            @test CLM.GetBurned(prt, organ, CLM.nitrogen_element) ≈ 0.0
            # the damaged diagnostic carries the loss; conservation holds
            @test CLM.CheckMassConservation(prt, ipft, 1) === nothing
        end

        # ===============================================================
        @testset "PRTReproRelease!: returns mass_out, sets val0 for balance" begin
            prt = fresh_plant()
            frac = 0.6
            el = CLM.carbon12_element
            b = CLM.GetState(prt, CLM.repro_organ, el)
            mass_out = CLM.PRTReproRelease!(prt, CLM.repro_organ, el, frac)
            @test mass_out ≈ frac * b
            @test CLM.GetState(prt, CLM.repro_organ, el) ≈ (1 - frac) * b
            # the val0 hack keeps val - val0 == net_alloc (mass balance preserved)
            @test CLM.CheckMassConservation(prt, ipft, 1) === nothing
            # calling on a non-repro organ errors
            @test_throws Exception CLM.PRTReproRelease!(prt, CLM.leaf_organ, el, 0.1)
        end

        # ===============================================================
        @testset "DeciduousTurnover!: retrans partition C/N/P -> storage" begin
            prt = fresh_plant()
            organ = CLM.fnrt_organ
            frac = 0.5

            # element -> expected retranslocation fraction
            exp_retr = Dict(CLM.carbon12_element => 0.0,
                            CLM.nitrogen_element  => retr_n,
                            CLM.phosphorus_element => retr_p)

            # snapshots before
            org_before = Dict(el => CLM.GetState(prt, organ, el) for el in keys(exp_retr))
            store_before = Dict(el => CLM.GetState(prt, CLM.store_organ, el)
                                for el in keys(exp_retr))
            plant_before = Dict(el => _plant_mass(prt, g, el) for el in keys(exp_retr))

            CLM.PRTDeciduousTurnover!(prt, ipft, organ, frac)

            for (el, retr) in exp_retr
                b = org_before[el]
                turnover_mass = (1 - retr) * frac * b
                retrans_mass  = retr * frac * b

                # organ pool drops by both pieces
                @test CLM.GetState(prt, organ, el) ≈ b - (turnover_mass + retrans_mass)
                # the leaves-plant piece is in the turnover diagnostic
                @test CLM.GetTurnover(prt, organ, el) ≈ turnover_mass
                # the retranslocated piece moved to storage (val grew by retrans_mass)
                @test CLM.GetState(prt, CLM.store_organ, el) ≈ store_before[el] + retrans_mass

                # ---- mass conservation: removed == turnover(out) + retrans(stays)
                removed = b - CLM.GetState(prt, organ, el)
                @test removed ≈ turnover_mass + retrans_mass

                # ---- whole-plant: only the turnover piece truly leaves the plant
                @test _plant_mass(prt, g, el) ≈ plant_before[el] - turnover_mass
            end

            # per-variable flux/state invariant still holds for the whole plant
            @test CLM.CheckMassConservation(prt, ipft, 1) === nothing
        end

        # ===============================================================
        @testset "MaintTurnover!: evergreen leaf + branch/root retrans" begin
            prt = fresh_plant()

            # expected base turnover rates [1/day]
            ypd = CLM.years_per_day
            bt_fnrt = ypd / p.root_long[ipft]
            bt_branch = ypd / p.branch_long[ipft]   # sapw/struct/store
            # canopy layer (icanlayer==1), not drought -> leaf uses last leaf_long
            bt_leaf = ypd / p.leaf_long[ipft, end]

            # snapshot fnrt-N before (retranslocated fraction = retr_n)
            fnrt_n0 = CLM.GetState(prt, CLM.fnrt_organ, CLM.nitrogen_element)
            # whole-plant N before: only the turnover term truly leaves the plant;
            # retranslocation just relocates N to the storage pool.
            plant_n0 = _plant_mass(prt, g, CLM.nitrogen_element)
            store_n0 = CLM.GetState(prt, CLM.store_organ, CLM.nitrogen_element)
            # leaf-carbon last age bin before (only senescing position turns over)
            ileaf_c = CLM.sp_organ_map_get(g, CLM.leaf_organ, CLM.carbon12_element)
            nlpos = g.state_descriptor[ileaf_c].num_pos
            leaf_c_sen0 = prt.variables[ileaf_c].val[nlpos]
            leaf_c_young0 = prt.variables[ileaf_c].val[1]

            CLM.PRTMaintTurnover!(prt, ipft, 1, false)

            # --- fine-root nitrogen: carbon-retrans=0 for C, retr_n for N ---
            # turnover removed first (no retrans on the turnover term), then retrans
            # on the *post-turnover* value (matches Fortran ordering exactly).
            t_mass = (1 - retr_n) * bt_fnrt * fnrt_n0
            after_turn = fnrt_n0 - t_mass
            r_mass = retr_n * bt_fnrt * after_turn
            @test CLM.GetTurnover(prt, CLM.fnrt_organ, CLM.nitrogen_element) ≈ t_mass
            @test CLM.GetState(prt, CLM.fnrt_organ, CLM.nitrogen_element) ≈ after_turn - r_mass
            # storage strictly gains N (retranslocation from every turning-over organ)
            @test CLM.GetState(prt, CLM.store_organ, CLM.nitrogen_element) > store_n0

            # whole-plant N drops by exactly the sum of all turnover terms; the
            # retranslocated N is conserved within the plant (it moved to storage).
            total_turn_n = 0.0
            for i_var in 1:g.num_vars
                if g.state_descriptor[i_var].element_id == CLM.nitrogen_element
                    total_turn_n += sum(prt.variables[i_var].turnover)
                end
            end
            @test _plant_mass(prt, g, CLM.nitrogen_element) ≈ plant_n0 - total_turn_n

            # --- leaf carbon: only the senescing (last) age position turns over;
            #     carbon never retranslocates so the whole turnover leaves ---
            @test prt.variables[ileaf_c].val[1] ≈ leaf_c_young0       # young bin untouched
            lc_t = bt_leaf * leaf_c_sen0
            @test prt.variables[ileaf_c].val[nlpos] ≈ leaf_c_sen0 - lc_t
            @test prt.variables[ileaf_c].turnover[nlpos] ≈ lc_t

            @test bt_branch > 0   # sanity: branch turnover active
            # whole-plant mass-balance invariant
            @test CLM.CheckMassConservation(prt, ipft, 1) === nothing
        end

        # ===============================================================
        @testset "PhenologyFlush!: storage -> leaf C then nutrient top-up" begin
            prt = fresh_plant()
            organ = CLM.leaf_organ
            cfrac = 0.5

            i_store_c = CLM.sp_organ_map_get(g, CLM.store_organ, CLM.carbon12_element)
            store_c0 = prt.variables[i_store_c].val[1]
            ileaf_c  = CLM.sp_organ_map_get(g, organ, CLM.carbon12_element)
            leaf_c0  = prt.variables[ileaf_c].val[1]

            CLM.PRTPhenologyFlush!(prt, ipft, organ, cfrac)

            # carbon transfer = store_c0 * cfrac into leaf bin 1
            xfer = store_c0 * cfrac
            @test prt.variables[ileaf_c].val[1] ≈ leaf_c0 + xfer
            @test prt.variables[i_store_c].val[1] ≈ store_c0 - xfer
            # net_alloc mirrors the move (+ on leaf, - on storage)
            @test prt.variables[ileaf_c].net_alloc[1] ≈ xfer

            # nutrient top-up toward target = leaf_c(new) * stoich, demand-limited
            ileaf_n = CLM.sp_organ_map_get(g, organ, CLM.nitrogen_element)
            i_store_n = CLM.sp_organ_map_get(g, CLM.store_organ, CLM.nitrogen_element)
            leaf_c_new = prt.variables[ileaf_c].val[1]
            target = leaf_c_new * p.nitr_stoich_p1[ipft, p.organ_param_id[organ]]
            leaf_n_pre = 1.0   # default fresh value before flush
            demand = max(0.0, target - leaf_n_pre)
            # storage N before nutrient phase is still the fresh default (1.0)
            n_xfer = min(demand, 1.0)
            @test prt.variables[ileaf_n].val[1] ≈ leaf_n_pre + n_xfer
            @test prt.variables[i_store_n].val[1] ≈ 1.0 - n_xfer
        end

    finally
        CLM.prt_global[] = old_global
    end
end
