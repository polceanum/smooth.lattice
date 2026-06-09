# Certified k=8 Sums-Only Unrank

- Timestamp: `2026-06-09T08:44:52.980804+00:00`
- Git commit: `503a9893a389bd52604bd5624ddb9ce76cda286c`
- Git dirty: `False`
- P: `2,3,5,7,11,13,17,19`
- N: `1000000000000`
- Solver return code: `0`
- Auditor return code: `0`
- Rank certified: `True`
- Certified count_le: `1000000000000`
- Exponents: `[75, 28, 9, 16, 3, 22, 5, 1]`

| solver wall s | audit wall s | solver reported s | cert reported s | groupA | groupB | ambiguous | digits |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2.760041 | 20.606740 | 2.028203 | 19.980800 | 35137278 | 1717623 | 1 | 91 |

The solver is the exploratory high-k sums-only MITM path. The certificate is
the independent interval-log auditor with exact big-integer resolution of
boundary ambiguities. This artifact supports correctness of this fixed
instance; it is not a broad best-known claim.
