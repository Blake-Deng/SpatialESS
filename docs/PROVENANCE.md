# Reference provenance

Reference snapshots were recorded on 2026-07-30.

- CellChatESS v0.2.1: commit
  `551d1d01e303ba7bc0b6870a1d422297ce918bee`
- kidney_autoimmune: commit
  `b90251686e25117f2db9104a06eefd7b17f591e9`
- FastCCC: commit
  `e63c386bfc406b1a7aa590411db64f6e8b506d9d`
- CytoSignal: commit
  `b2c1a3cfda258fa57f6c8c3a92635dfe4b0505a0`
- CytoSignal figure repository: commit
  `4e38ba9ca915be292b289f4b07d47ac839a7ef93`

Local read-only snapshots are stored under
`/home/dzf/cellchat_acceleration/SpatialESS_research_20260730/references`.


The public package does not include raw expression matrices, RDS objects,
H5/H5AD files or full benchmark outputs. The reproducible compact tables are
under results/tables/; their source roots and regeneration commands are in
docs/REPRODUCE_MAIN.md. CosMx queue and output reconciliation is recorded in
results/tables/cosmx_queue_reconciliation.tsv.

The release benchmark environment was Ubuntu 24.04.4 LTS on two AMD EPYC
7763 sockets (128 physical cores, 256 logical CPUs) with approximately 1 TiB
RAM. It used R 4.5.3, Matrix 1.7.5, Rcpp 1.1.1.1.1, testthat 3.3.2 and C++17
(x86_64-conda-linux-gnu-c++). Main comparisons used one native thread and
one BLAS/OpenMP thread, with matching cells, genes, LR table, coordinates,
seed and 20 permutations between methods. The package check completed with
Status: OK.
