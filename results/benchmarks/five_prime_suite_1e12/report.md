# Certified Five-Prime Suite

- Timestamp: `2026-06-08T21:16:09.606220+00:00`
- Git commit: `87da37090b939d217a5ee2a51e2c08101d5d13ac`
- Git dirty: `False`
- N: `1000000000000`
- Cases: `6`
- Certified cases: `6`
- Layer wall-time wins: `6`
- Mean wall ratio X+Y/layer: `1.4962629718211538`

The `xplusy_adaptive` rows measure adaptive Cartesian-sum value selection only.
The `layer_compressed` rows measure full exponent-vector unranking.

| P | layer wall s | X+Y wall s | ratio | layer RSS KB | X+Y RSS KB | certified | exps |
|---|---:|---:|---:|---:|---:|---|---|
| `2,3,5,7,11` | 1.367318 | 2.122471 | 1.552288 | 23256 | 153456 | True | `[1052, 26, 33, 53, 4]` |
| `2,3,5,7,13` | 0.772592 | 1.109113 | 1.435574 | 23700 | 149692 | True | `[205, 279, 119, 131, 16]` |
| `2,3,5,11,13` | 0.851199 | 1.010214 | 1.186813 | 31228 | 138540 | True | `[254, 220, 258, 4, 52]` |
| `2,3,7,11,13` | 0.797991 | 1.176589 | 1.474439 | 29068 | 129264 | True | `[291, 331, 90, 84, 28]` |
| `2,5,7,11,13` | 1.032364 | 1.461805 | 1.415978 | 26244 | 159688 | True | `[106, 306, 164, 53, 32]` |
| `3,5,7,11,13` | 0.861198 | 1.647030 | 1.912485 | 30840 | 207688 | True | `[273, 219, 5, 219, 5]` |
