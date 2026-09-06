# Stage 1A sparse exact triMean validation

Date: 2026-07-31

The observed and permuted triMean paths now share a zero-aware R type-7 kernel.
The production kernel uses repeated `std::nth_element` calls on one mutable
positive-value buffer and never expands structural zeros. A full-sort
implementation is retained only in an internal benchmark oracle.

## Kernel benchmark

The benchmark crossed five group sizes (100 to 1,000,000 cells) with six zero
fractions (0 to 0.99). All 30 sort-versus-multi-select comparisons had an
absolute numerical difference of zero.

For combinations that required positive order statistics, multi-select was
approximately 1.5- to 21.2-fold faster on medium and large groups. One
100-cell dense case was approximately 4% slower. At zero fractions high enough
to force all triMean quantiles to zero, both kernels return immediately.

Files:

- `benchmarks/results/sparse_type7_kernel/sort_vs_multiselect.tsv`
- `benchmarks/results/sparse_type7_kernel/time_verbose.txt`

## Million-cell regression

The new kernel recomputed 534 referenced genes across 2,357 groups in the
1,156,091-cell Xenium object. All 1,258,638 output values were bitwise
identical to the frozen sort result:

- maximum absolute difference: 0;
- differing values: 0;
- multi-select kernel time: 5.444 seconds;
- independent-process peak RSS: 4,250,680 kB;
- process exit status: 0.

The corresponding frozen sort run recorded 6.730 seconds for the same triMean
stage. The observed full-object speedup is therefore 1.24-fold; loading the
large CellChat object remains outside this kernel time.

Files:

- `benchmarks/results/multiselect_xenium1156091/summary.tsv`
- `benchmarks/results/multiselect_xenium1156091/time_verbose.txt`
- `benchmarks/results/multiselect_xenium1156091/group_trimean_multiselect.rds`

