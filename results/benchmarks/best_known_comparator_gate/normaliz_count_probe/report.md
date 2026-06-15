# Normaliz Count Probe

- Tool: `PyNormaliz`
- Toy cases passed: `2/2`
- Target cases completed: `0/3`
- Target timeouts: `3`
- Target errors: `0`
- Timeout per count: `5.0` seconds

The target counts use rationalized log simplexes and are comparator
measurements, not independent certificates for the original irrational-log
rank. Correctness certification still comes from the interval auditor.

| P | d | mode | status | seconds | count | delta from certified N |
|---|---:|---|---|---:|---:|---:|
| `2,3,5,7,11` | 1 | rounded | timeout | 5.006109 |  | None |
| `2,3,5,7,11,13` | 1 | rounded | timeout | 5.004431 |  | None |
| `2,3,5,7,11,13,17,19` | 1 | rounded | timeout | 5.004594 |  | None |
