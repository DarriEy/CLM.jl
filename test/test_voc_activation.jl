# =============================================================================
# test_voc_activation.jl
#
# Tests the LIVE VOC/MEGAN ACTIVATION PLUMBING (not the kernel physics, which is
# oracle-validated in scripts/validate_voc_megan_fortran_parity.jl):
#   * megan_factors_file_init!   — read megan_factors_file NetCDF -> table + mf
#   * megan_config_init          — shr_megan_readnl + megan_factors_init analogue
#   * read_efisop_from_surfdata! — VOCEmission::iniTimeConst analogue (EF1_* read)
#   * CLMDriverConfig activation constructor path
#   * DEFAULT (VOC off) stays byte-identical / inert
#
# Uses SYNTHETIC NetCDF files written to a tempdir (CTSM variable layout) so the
# real NCDatasets reader code path is exercised without any external inputdata —
# portable to CI. (scripts/validate_voc_megan_activation.jl runs the same path on
# the REAL megan_factors_file + surfdata when present.)
# =============================================================================
using NCDatasets

# ---- helpers to write synthetic CTSM-layout NetCDF fixtures ------------------

# Write a minimal MEGAN2.1 emission-factors file (exact CDL dim/var layout of
# megan21_emis_factors_*.nc): char Comp_Name(Comp_Num,char), Comp_EF(PFT_Num,
# Comp_Num), Class_EF(PFT_Num,Class_Num), per-class LDF/Agro/.../Ceo.
function _write_megan_factors_file(path; npft=20, nclass=3,
                                   names=["isoprene","myrcene","limonene"],
                                   class_nums=[1,2,3], mw=[68.12,136.23,136.23])
    ncomp = length(names)
    NCDataset(path, "c") do ds
        defDim(ds, "Comp_Num", ncomp)
        defDim(ds, "Class_Num", nclass)
        defDim(ds, "PFT_Num", npft)
        nchar = 40
        defDim(ds, "char", nchar)
        # Comp_Name as a char(Comp_Num, char) field -> stored (char, Comp_Num),
        # NUL-padded, to exercise the name cleaner.
        cn = fill('\0', nchar, ncomp)
        for (j, nm) in enumerate(names), (i, ch) in enumerate(collect(nm))
            cn[i, j] = ch
        end
        v = defVar(ds, "Comp_Name", Char, ("char", "Comp_Num")); v[:, :] = cn
        defVar(ds, "Comp_MW", Float64, ("Comp_Num",))[:] = mw
        defVar(ds, "Class_Num", Int16, ("Comp_Num",))[:] = Int16.(class_nums)
        # Comp_EF(PFT_Num,Comp_Num) in CDL -> NCDatasets writes with dims named in
        # that order; give distinctive values so Comp_EF*Class_EF is checkable.
        comp_ef = [Float64(2 + i) for i in 1:ncomp, k in 1:npft]     # (ncomp, npft)
        defVar(ds, "Comp_EF", Float64, ("Comp_Num", "PFT_Num"))[:, :] = comp_ef
        class_ef = [Float64(10 * c) for c in 1:nclass, k in 1:npft]  # (nclass, npft)
        defVar(ds, "Class_EF", Float64, ("Class_Num", "PFT_Num"))[:, :] = class_ef
        # per-class coefficients (read from file, NOT the all-ones defaults)
        defVar(ds, "LDF",   Float64, ("Class_Num",))[:] = [0.999, 0.4, 0.2][1:nclass]
        defVar(ds, "Agro",  Float64, ("Class_Num",))[:] = fill(0.6, nclass)
        defVar(ds, "Amat",  Float64, ("Class_Num",))[:] = fill(1.125, nclass)
        defVar(ds, "Anew",  Float64, ("Class_Num",))[:] = fill(0.05, nclass)
        defVar(ds, "Aold",  Float64, ("Class_Num",))[:] = fill(1.0, nclass)
        defVar(ds, "betaT", Float64, ("Class_Num",))[:] = [0.13, 0.10, 0.10][1:nclass]
        defVar(ds, "ct1",   Float64, ("Class_Num",))[:] = fill(95.0, nclass)
        defVar(ds, "ct2",   Float64, ("Class_Num",))[:] = fill(230.0, nclass)
        defVar(ds, "Ceo",   Float64, ("Class_Num",))[:] = fill(2.0, nclass)
    end
    return path
end

# Write a minimal surfdata-like file carrying the 6 EF1_* gridded isoprene maps.
function _write_surfdata_with_ef(path; efvals=(10000.0, 2000.0, 2000.0, 4000.0, 800.0, 1.0))
    vars = ("EF1_BTR","EF1_FET","EF1_FDT","EF1_SHR","EF1_GRS","EF1_CRP")
    NCDataset(path, "c") do ds
        defDim(ds, "lsmlon", 1); defDim(ds, "lsmlat", 1)
        for (k, vn) in enumerate(vars)
            defVar(ds, vn, Float64, ("lsmlon", "lsmlat"))[:, :] = fill(efvals[k], 1, 1)
        end
    end
    return path
end

@testset "VOC/MEGAN activation plumbing" begin
    tmp = mktempdir()
    mff = _write_megan_factors_file(joinpath(tmp, "megan_factors.nc"))
    sd  = _write_surfdata_with_ef(joinpath(tmp, "surfdata_ef.nc"))

    @testset "megan_factors_file_init! reads table + per-class coeffs from file" begin
        tbl = CLM.MEGANFactorsTable{Float64}()
        mf  = CLM.MEGANFactors{Float64}()
        CLM.megan_factors_file_init!(mf, tbl, mff)
        @test tbl.npfts == 20
        # stored factor for isoprene = Comp_EF(iso)*Class_EF(class1) = (2+1)*(10*1)=30
        (f, cnum, w) = CLM.megan_factors_get(tbl, "isoprene")
        @test cnum == 1
        @test w ≈ 68.12
        @test f[1] ≈ 30.0
        # myrcene = (2+2)*(10*2) = 80 ; limonene = (2+3)*(10*3) = 150
        @test CLM.megan_factors_get(tbl, "myrcene")[1][1]  ≈ 80.0
        @test CLM.megan_factors_get(tbl, "limonene")[1][1] ≈ 150.0
        # per-class coefficients came from the FILE (not the all-ones defaults)
        @test mf.n_classes == 3
        @test mf.LDF[1]   ≈ 0.999
        @test mf.betaT[1] ≈ 0.13
        @test mf.Anew[1]  ≈ 0.05
        @test mf.Ceo[1]   ≈ 2.0
    end

    @testset "megan_config_init builds populated MEGANConfig via the reader path" begin
        spec = ["ISOP = isoprene", "BIGENE = isoprene + 0.5*myrcene"]
        cfg = CLM.megan_config_init(; use_voc=true, megan_specifier=spec,
                                    megan_factors_file=mff)
        @test length(cfg.meg_compounds) == 2     # isoprene, myrcene (deduped)
        @test length(cfg.mech_comps)    == 2     # ISOP, BIGENE
        @test cfg.mapped_emisfctrs == false      # CTSM default
        @test cfg.meg_compounds[1].name == "isoprene"
        @test cfg.meg_compounds[1].emis_factors[1] ≈ 30.0
        @test cfg.meg_compounds[2].emis_factors[1] ≈ 80.0
        @test cfg.megan_factors.LDF[1] ≈ 0.999   # coeffs threaded from the file
        # BIGENE maps isoprene + 0.5*myrcene
        @test cfg.mech_comps[2].n_megan_comps == 2
        @test cfg.meg_compounds[2].coeff ≈ 0.5

        # use_voc=false OR no specifier => inert (empty) config
        @test isempty(CLM.megan_config_init(; use_voc=false, megan_specifier=spec,
                                            megan_factors_file=mff).meg_compounds)
        @test isempty(CLM.megan_config_init(; use_voc=true, megan_specifier=nothing,
                                            megan_factors_file=mff).mech_comps)
    end

    @testset "read_efisop_from_surfdata! (present + absent)" begin
        voc = CLM.VOCEmisData(); CLM.vocemis_init!(voc, 1, 1, 2, 2)
        @test all(isnan, voc.efisop_grc)              # starts NaN
        ok = CLM.read_efisop_from_surfdata!(voc, sd, [1])
        @test ok
        @test all(isfinite, voc.efisop_grc)
        @test voc.efisop_grc[1, 1] ≈ 10000.0
        @test voc.efisop_grc[6, 1] ≈ 1.0
        # absent EF1_* -> filled with default (0), returns false, warns
        voc2 = CLM.VOCEmisData(); CLM.vocemis_init!(voc2, 1, 1, 2, 2)
        ok2 = @test_logs (:warn,) match_mode=:any CLM.read_efisop_from_surfdata!(voc2, mff, [1])
        @test ok2 == false
        @test all(iszero, voc2.efisop_grc)
    end

    @testset "DEFAULT CLMDriverConfig is VOC-off / inert (byte-identical)" begin
        c = CLM.CLMDriverConfig()
        @test c.use_voc === false
        @test isempty(c.megan.meg_compounds)
        @test isempty(c.megan.mech_comps)
        @test c.megan.mapped_emisfctrs === false
        # driver gate short-circuits => voc_emission! never called
        gate = c.use_voc && !isempty(c.megan.mech_comps) && !isempty(c.megan.meg_compounds)
        @test gate === false
    end

    @testset "CLMDriverConfig activation constructor + kernel emits nonzero flux" begin
        spec = ["ISOP = isoprene", "BIGENE = isoprene + 0.5*myrcene"]
        cfg = CLM.CLMDriverConfig(use_voc=true, megan_specifier=spec,
                                  megan_factors_file=mff)
        @test cfg.use_voc === true
        @test length(cfg.megan.meg_compounds) == 2
        gate = cfg.use_voc && !isempty(cfg.megan.mech_comps) &&
               !isempty(cfg.megan.meg_compounds)
        @test gate === true

        # Drive the kernel through the init-built descriptors (mapped=false => table EF).
        np, nc, ng = 1, 1, 1
        voc = CLM.VOCEmisData()
        CLM.vocemis_init!(voc, np, ng, length(cfg.megan.meg_compounds),
                          length(cfg.megan.mech_comps))
        CLM.read_efisop_from_surfdata!(voc, sd, [1])
        patch = CLM.PatchData(); CLM.patch_init!(patch, np)
        patch.itype .= [6]; patch.gridcell .= [1]; patch.column .= [1]  # broadleaf decid trop
        CLM.voc_emission!(voc, cfg.megan.meg_compounds, cfg.megan.mech_comps,
            cfg.megan.megan_factors, patch, 1:np, trues(np),
            fill(220.0, nc, 2), fill(95.0, ng, 2), fill(101325.0, nc), fill(45.0, ng),
            fill(180.0, np), fill(170.0, np), fill(85.0, np), fill(80.0, np),
            fill(0.5, np), fill(0.5, np), fill(0.5, np), fill(2.0, np), fill(1.8, np),
            fill(0.7*45.0, np, 1), fill(0.6*45.0, np, 1),
            fill(300.0, np), fill(299.0, np), fill(297.0, np), fill(0.8, np);
            use_mapped_emisfctrs=cfg.megan.mapped_emisfctrs)
        @test isfinite(voc.vocflx_tot_patch[1])
        @test voc.vocflx_tot_patch[1] > 0.0
        @test all(isfinite, voc.efisop_grc)
    end
end
