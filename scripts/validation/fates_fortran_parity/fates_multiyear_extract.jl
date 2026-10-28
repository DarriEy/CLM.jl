# fates_multiyear_extract.jl <h0.nc>  (OUT=<csv> to write the daily CSV)
# Reduce the multi-year Fortran FATES daily history (stage_multiyear.py) to the
# per-day CSV that scripts/fates_multiyear_compare.jl scores against:
# ncoh, npatch, vegc/structc [gC/m2], NPP/GPP [gC/m2/day], recruitment [/m2/day],
# per-PFT nplant/recruitment, and per-PFT largest occupied FATES_NPLANT_SZPF class.
using NCDatasets, Printf, Statistics
f = ARGS[1]
ds = NCDataset(f)
edges = [0,5,10,15,20,30,40,50,60,70,80,90,100]
nsz = length(edges)
ncoh = ds["FATES_NCOHORTS"][1, :]
npat = ds["FATES_NPATCHES"][1, :]
rec  = ds["FATES_RECRUITMENT_PF"][1, :, :]   # (pft, time)  m-2 yr-1
npf  = ds["FATES_NPLANT_PF"][1, :, :]
szpf = ds["FATES_NPLANT_SZPF"][1, :, :]      # (sz*pft, time) sz fastest
vegc = ds["FATES_VEGC"][1, :]
strc = ds["FATES_STRUCTC"][1, :]
npp  = ds["FATES_NPP"][1, :]
gpp  = ds["FATES_GPP"][1, :]
stem = ds["FATES_STEM_ALLOC"][1, :]
ddbh = ds["FATES_DDBH_SZPF"][1, :, :]
np_ = size(rec, 1); nt = size(rec, 2)
println("nt=$nt npft=$np_ nszpf=$(size(szpf,1))")
out = get(ENV, "OUT", "")
io = isempty(out) ? nothing : open(out, "w")
io === nothing || println(io, "day,ncoh,npatch,vegc,structc,npp_d,gpp_d,recruit_d,maxsz_lower,", join(["nplant_$p" for p in 1:np_], ","), ",", join(["recruit_$p" for p in 1:np_], ","), ",", join(["maxsz_$p" for p in 1:np_], ","))
for t in 1:nt
    maxsz = 0; maxszp = zeros(Int, np_)
    for p in 1:np_, s in 1:nsz
        v = szpf[(p-1)*nsz + s, t]
        if !ismissing(v) && v > 0
            maxsz = max(maxsz, edges[s]); maxszp[p] = max(maxszp[p], edges[s])
        end
    end
    r = coalesce.(rec[:, t], 0.0) ./ 365   # per m2 per day
    io === nothing || @printf(io, "%d,%g,%g,%.8g,%.8g,%.8g,%.8g,%.8g,%d,%s,%s,%s\n", t, ncoh[t], npat[t], vegc[t]*1000, strc[t]*1000,
        npp[t]*86400*1000, gpp[t]*86400*1000, sum(r), maxsz,
        join((@sprintf("%.6g", coalesce(x,0.0)) for x in npf[:, t]), ","),
        join((@sprintf("%.6g", x) for x in r), ","), join(maxszp, ","))
    if t % 30 == 0 || t == 1
        @printf("day %4d ncoh=%4d npatch=%d vegC=%8.1f gC/m2 structC=%8.1f NPP=%7.3f gC/m2/d rec/d=%.4g maxsz>=%d\n",
            t, ncoh[t], npat[t], vegc[t]*1000, strc[t]*1000, npp[t]*86400*1000, sum(r), maxsz)
    end
end
io === nothing || close(io)
