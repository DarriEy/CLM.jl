# Convert date windows to CROPDUMP_WINDOWS nstep ranges for a run that starts at
# START (default 2016-01-01, tod 14400) on CTSM's NOLEAP calendar, dt = 3600 s.
# nstep n's model time (get_curr_date, END of the step) = start + n*dt.
#
#   julia nstep_windows.jl 2017-07-03:2017-08-15 2017-10-02:2017-10-31 ...
#   START=2017-07-25T04 julia nstep_windows.jl ...
using Dates

noleap_days(d::Date) = (year(d) - 2000) * 365 + dayofyear(d) - ((isleapyear(d) && month(d) > 2) ? 1 : 0)
st = DateTime(get(ENV, "START", "2016-01-01T04"))
h0 = noleap_days(Date(st)) * 24 + hour(st)
nst(d::Date, h::Int) = noleap_days(d) * 24 + h - h0
wins = String[]
for a in ARGS
    lo, hi = Date.(split(a, ":"))
    push!(wins, "$(nst(lo, 0) + 1):$(nst(hi, 0) + 24)")   # 01:00 of lo .. 24:00 of hi
end
println(join(wins, ","))
