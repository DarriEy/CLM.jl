# Build the crop-tails surfdata: the USplains 1-pt crop surfdata (gridcell 352 of
# the ne3np4 78-pft file, scripts/build_crop_cft_surfdata.py) with the foddergrass
# CFT (pft 45, 15.1% of the crop landunit) re-assigned to WINTER WHEAT (pft 21).
#
# Why a re-assignment rather than an extract: in
#   surfdata_ne3np4_hist_2000_78pfts_c251022.nc
# PCT_CFT is identically 0 for every winter cereal (winter wheat 21/22, winter
# barley 27/28, winter rye 31/32) in all 488 gridcells (checked; the LUH3 wheat
# area is mapped entirely to spring wheat), so no real winter-wheat cell exists.
# pft 45 is merged into spring wheat (mergetoclmpft(45) = 19) by CTSM anyway, so
# moving it to 21 leaves corn/soy/spring-wheat untouched and adds one live
# winter-wheat patch. clm5_params has a complete winter-wheat parameter set
# (plant 1 Sep - 30 Nov, mxmat 330, mergetoclmpft(21) = 21).
#
# Usage: julia --project=. make_winterwheat_surfdata.jl [src] [dst]
using NCDatasets

sd  = get(ENV, "SYMFLUENCE_DATA", "/Users/darri/repos/SYMFLUENCE_data")
src = length(ARGS) >= 1 ? ARGS[1] : joinpath(sd, "crop_cft_surfdata", "surfdata_cropCFT_USplains_1pt.nc")
dst = length(ARGS) >= 2 ? ARGS[2] : joinpath(sd, "runs", "crop_tails_inputs", "surfdata_cropCFT_USplains_1pt_wwheat.nc")
mkpath(dirname(dst))
cp(src, dst; force = true)
chmod(dst, 0o644)

NCDataset(dst, "a") do ds
    v = ds["PCT_CFT"]                       # (lsmlon, lsmlat, cft), cft 1 <-> pft 15
    dn = dimnames(v)
    @assert dn[end] == "cft" "unexpected PCT_CFT dims $dn"
    a = Array(v[:, :, :])
    i45 = 45 - 15 + 1
    i21 = 21 - 15 + 1
    @assert a[1, 1, i21] == 0 "winter wheat already present?"
    moved = a[1, 1, i45]
    a[1, 1, i21] += moved
    a[1, 1, i45] = 0
    @assert isapprox(sum(a[1, 1, :]), 100; atol = 1e-6)
    v[:, :, :] = a
    ds.attrib["crop_tails_note"] =
        "CLM.jl crop-tail parity: PCT_CFT of pft 45 ($(moved)%) moved to winter wheat (pft 21)"
    println("moved $(moved)% of the crop landunit from pft 45 to pft 21 -> $dst")
end
