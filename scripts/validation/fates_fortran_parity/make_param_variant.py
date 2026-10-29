"""make_param_variant.py <in.cdl> <out.cdl> name=value [name=value ...]

Write a copy of a FATES parameter CDL with the named parameters set to `value`
(broadcast over every entry the variable already has, e.g. every PFT). Used by the
FATES research-mode references (scripts/validation/fates_fortran_parity/README.md "RESEARCH MODES"), whose switches are inert
with the default parameter values:

  seed dispersal: fates_seed_dispersal_fraction / _max_dist / _pdf_scale / _pdf_shape
                  are unset (`_`) in the default file;
  tree damage:    fates_damage_event_code = 1 means "no damage events".

The Julia port reads the .cdl directly; build the Fortran .nc with
`ncgen -o out.nc out.cdl`.
"""
import re, sys

src, dst = sys.argv[1], sys.argv[2]
s = open(src).read()
for kv in sys.argv[3:]:
    k, v = kv.split("=", 1)
    m = re.search(rf"(?m)^ {k} = ([^;]*);", s)
    assert m, k
    n = len([x for x in m.group(1).split(",") if x.strip()])
    s = s[:m.start()] + f" {k} = " + ", ".join([v] * n) + " ;" + s[m.end():]
open(dst, "w").write(s)
print("ok", dst)
