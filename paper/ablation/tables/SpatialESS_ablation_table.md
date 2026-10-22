# SpatialESS ablation table

Validated: 170/170 runs. Pending and incomplete values are provisional.

| Dataset | Cells | Permutations | Configuration | Graph (s) | Observed (s) | Permutation (s) | Inference (s) | Peak RSS (GiB) | Time / Full | Validated / planned | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| cosmx_5000 | 5000 | 20 | Full SpatialESS | 0.011 [0.010, 0.011] | 0.125 [0.125, 0.127] | 0.537 [0.534, 0.540] | 0.738 [0.738, 0.739] | 0.198 [0.198, 0.199] | 1.00 | 5/5 | COMPLETE |
| cosmx_5000 | 5000 | 20 | Without spatial hashing | 0.170 [0.167, 0.177] | 0.124 [0.124, 0.125] | 0.531 [0.521, 0.533] | 0.894 [0.888, 0.902] | 0.198 [0.198, 0.199] | 1.21 | 5/5 | COMPLETE |
| cosmx_5000 | 5000 | 20 | COO instead of CSR | 0.011 [0.011, 0.011] | 0.237 [0.233, 0.238] | 0.530 [0.520, 0.532] | 0.845 [0.841, 0.848] | 0.197 [0.197, 0.197] | 1.14 | 5/5 | COMPLETE |
| cosmx_5000 | 5000 | 20 | Without edge-score reuse | 0.011 [0.011, 0.011] | 0.135 [0.135, 0.137] | 0.752 [0.737, 0.764] | 0.973 [0.962, 0.975] | 0.199 [0.197, 0.199] | 1.32 | 5/5 | COMPLETE |
| cosmx_5000 | 5000 | 20 | Dense neighborhood baseline | 0.516 [0.510, 0.536] | 1.580 [1.570, 1.610] | 0.518 [0.517, 0.529] | 2.813 [2.705, 2.831] | 0.369 [0.369, 0.369] | 3.81 | 5/5 | COMPLETE |
| cosmx_10000 | 10000 | 20 | Full SpatialESS | 0.022 [0.021, 0.028] | 0.257 [0.256, 0.260] | 1.115 [1.115, 1.128] | 1.476 [1.469, 1.485] | 0.203 [0.203, 0.204] | 1.00 | 5/5 | COMPLETE |
| cosmx_10000 | 10000 | 20 | Without spatial hashing | 0.667 [0.665, 0.672] | 0.255 [0.255, 0.257] | 1.104 [1.098, 1.106] | 2.104 [2.102, 2.114] | 0.203 [0.203, 0.204] | 1.43 | 5/5 | COMPLETE |
| cosmx_10000 | 10000 | 20 | COO instead of CSR | 0.022 [0.022, 0.022] | 0.482 [0.478, 0.484] | 1.103 [1.101, 1.112] | 1.696 [1.689, 1.696] | 0.201 [0.201, 0.202] | 1.15 | 5/5 | COMPLETE |
| cosmx_10000 | 10000 | 20 | Without edge-score reuse | 0.021 [0.021, 0.021] | 0.282 [0.281, 0.285] | 1.588 [1.569, 1.592] | 1.965 [1.943, 1.975] | 0.203 [0.203, 0.204] | 1.33 | 5/5 | COMPLETE |
| cosmx_10000 | 10000 | 20 | Dense neighborhood baseline | 1.021 [1.013, 1.037] | 6.073 [5.947, 6.096] | 1.105 [1.098, 1.119] | 8.762 [8.629, 8.868] | 0.943 [0.943, 0.943] | 5.94 | 5/5 | COMPLETE |
| cosmx_10000 | 10000 | 100 | Full SpatialESS | 0.020 [0.020, 0.020] | 0.253 [0.253, 0.253] | 5.368 [5.363, 5.371] | 5.706 [5.701, 5.710] | 0.225 [0.225, 0.225] | 1.00 | 5/5 | COMPLETE |
| cosmx_10000 | 10000 | 100 | Without edge-score reuse | 0.020 [0.020, 0.020] | 0.276 [0.276, 0.276] | 7.682 [7.664, 7.693] | 8.042 [8.024, 8.053] | 0.225 [0.224, 0.225] | 1.41 | 5/5 | COMPLETE |
| cosmx_10000 | 10000 | 500 | Full SpatialESS | 0.020 [0.020, 0.020] | 0.267 [0.267, 0.270] | 26.776 [26.730, 26.961] | 27.138 [27.097, 27.328] | 0.276 [0.276, 0.277] | 1.00 | 5/5 | COMPLETE |
| cosmx_10000 | 10000 | 500 | Without edge-score reuse | 0.020 [0.020, 0.020] | 0.293 [0.292, 0.300] | 38.620 [38.454, 38.632] | 39.010 [38.844, 39.029] | 0.276 [0.275, 0.277] | 1.44 | 5/5 | COMPLETE |
| cosmx_25000 | 25000 | 20 | Full SpatialESS | 0.051 [0.051, 0.053] | 0.697 [0.696, 0.706] | 2.902 [2.901, 2.905] | 4.035 [4.024, 4.158] | 0.241 [0.240, 0.241] | 1.00 | 5/5 | COMPLETE |
| cosmx_25000 | 25000 | 20 | Without spatial hashing | 4.079 [4.076, 4.121] | 0.693 [0.690, 0.701] | 2.908 [2.885, 2.949] | 8.114 [7.998, 8.210] | 0.241 [0.240, 0.241] | 2.01 | 5/5 | COMPLETE |
| cosmx_25000 | 25000 | 20 | COO instead of CSR | 0.056 [0.055, 0.057] | 1.262 [1.261, 1.289] | 2.936 [2.904, 2.966] | 4.786 [4.703, 4.800] | 0.241 [0.241, 0.242] | 1.19 | 5/5 | COMPLETE |
| cosmx_25000 | 25000 | 20 | Without edge-score reuse | 0.052 [0.051, 0.053] | 0.763 [0.759, 0.764] | 4.257 [4.227, 4.261] | 5.559 [5.538, 5.581] | 0.241 [0.240, 0.241] | 1.38 | 5/5 | COMPLETE |
| cosmx_25000 | 25000 | 20 | Dense neighborhood baseline | 5.109 [5.076, 5.405] | 35.074 [34.797, 35.506] | 2.902 [2.897, 2.934] | 46.978 [46.571, 47.361] | 4.886 [4.885, 4.886] | 11.64 | 5/5 | COMPLETE |
| cervical_center_5000 | 5000 | 100 | Full SpatialESS | 0.017 [0.017, 0.017] | 0.346 [0.339, 0.357] | 5.317 [5.296, 5.409] | 5.997 [5.995, 6.040] | 0.245 [0.245, 0.246] | 1.00 | 5/5 | COMPLETE |
| cervical_center_5000 | 5000 | 100 | Without spatial hashing | 0.186 [0.172, 0.192] | 0.372 [0.366, 0.375] | 5.368 [5.347, 5.398] | 6.190 [6.151, 6.206] | 0.245 [0.245, 0.245] | 1.03 | 5/5 | COMPLETE |
| cervical_center_5000 | 5000 | 100 | COO instead of CSR | 0.018 [0.018, 0.020] | 0.935 [0.932, 0.942] | 5.362 [5.314, 5.385] | 6.618 [6.566, 6.640] | 0.246 [0.245, 0.247] | 1.10 | 5/5 | COMPLETE |
| cervical_center_5000 | 5000 | 100 | Without edge-score reuse | 0.018 [0.016, 0.018] | 0.359 [0.353, 0.363] | 5.382 [5.370, 5.423] | 6.062 [6.022, 6.130] | 0.246 [0.245, 0.246] | 1.01 | 5/5 | COMPLETE |
| cervical_center_5000 | 5000 | 100 | Dense neighborhood baseline | 0.574 [0.560, 0.580] | 0.544 [0.535, 0.551] | 5.263 [5.250, 5.292] | 6.821 [6.815, 6.894] | 0.425 [0.425, 0.425] | 1.14 | 5/5 | COMPLETE |
| cervical_stratified_5000 | 5000 | 100 | Full SpatialESS | 0.009 [0.008, 0.009] | 0.245 [0.242, 0.245] | 6.555 [6.517, 6.555] | 7.064 [7.023, 7.093] | 0.238 [0.238, 0.238] | 1.00 | 5/5 | COMPLETE |
| cervical_stratified_5000 | 5000 | 100 | Without spatial hashing | 0.166 [0.162, 0.167] | 0.239 [0.239, 0.242] | 6.750 [6.569, 6.755] | 7.447 [7.249, 7.488] | 0.238 [0.237, 0.238] | 1.05 | 5/5 | COMPLETE |
| cervical_stratified_5000 | 5000 | 100 | COO instead of CSR | 0.010 [0.009, 0.012] | 0.226 [0.225, 0.228] | 6.612 [6.566, 6.648] | 7.101 [7.088, 7.146] | 0.237 [0.237, 0.238] | 1.01 | 5/5 | COMPLETE |
| cervical_stratified_5000 | 5000 | 100 | Without edge-score reuse | 0.009 [0.008, 0.011] | 0.248 [0.236, 0.250] | 6.725 [6.628, 6.726] | 7.297 [7.149, 7.307] | 0.238 [0.238, 0.238] | 1.03 | 5/5 | COMPLETE |
| cervical_stratified_5000 | 5000 | 100 | Dense neighborhood baseline | 0.555 [0.550, 0.556] | 0.603 [0.596, 0.606] | 6.506 [6.494, 6.546] | 8.118 [8.049, 8.191] | 0.420 [0.420, 0.421] | 1.15 | 5/5 | COMPLETE |
| cosmx_50000 | 50000 | 20 | Full SpatialESS | 0.199 [0.198, 0.199] | 1.403 [1.402, 1.407] | 5.993 [5.990, 6.036] | 7.830 [7.824, 7.867] | 0.321 [0.320, 0.321] | 1.00 | 5/5 | COMPLETE |
| cosmx_50000 | 50000 | 20 | Without spatial hashing | 16.140 [16.113, 16.157] | 1.393 [1.391, 1.396] | 6.021 [5.946, 6.086] | 23.814 [23.730, 23.829] | 0.320 [0.320, 0.321] | 3.04 | 5/5 | COMPLETE |
| cosmx_50000 | 50000 | 20 | COO instead of CSR | 0.203 [0.203, 0.203] | 2.610 [2.586, 2.611] | 6.010 [5.972, 6.045] | 9.053 [8.981, 9.090] | 0.323 [0.322, 0.323] | 1.16 | 5/5 | COMPLETE |
| cosmx_50000 | 50000 | 20 | Without edge-score reuse | 0.198 [0.198, 0.199] | 1.547 [1.544, 1.554] | 8.874 [8.845, 8.951] | 10.847 [10.818, 10.934] | 0.321 [0.321, 0.321] | 1.39 | 5/5 | COMPLETE |
| cosmx_50000 | 50000 | 20 | Dense neighborhood baseline | 14.648 [13.986, 16.037] | 107.003 [106.339, 115.844] | 6.141 [6.098, 6.148] | 142.172 [141.298, 152.957] | 18.928 [18.928, 18.929] | 18.16 | 5/5 | COMPLETE |

Times and RSS: median [25th percentile, 75th percentile] of numerical-validation PASS runs.
Time / Full = ratio of median inference times, compared within the same input and nboot.
Incomplete groups and ratios are provisional. Empty cells are not zero or inferred measurements.
Graph, observed and permutation times exclude input loading and output serialization.
Whole-process wall time and peak RSS include loading, frozen permutations and output writing.
Graph payload is measured numeric-array bytes; it is not process RSS.
CSR-to-COO tests graph format and traversal. CSR-to-dense tests sparse representation and traversal.
No-reuse recomputes LR edge probabilities from cached invariant molecular factors and edge support.
Canonical sender/receiver accumulation order, frozen label permutations and filtering are preserved.
Five repetitions are technical measurements, not independent biological replicates.
