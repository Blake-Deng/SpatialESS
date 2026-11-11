# Current-source integrity

`SHA256SUMS.txt` covers every regular file in this source tree except itself and
Git internals. `SYMLINKS.tsv` records link paths and literal targets; links are
not dereferenced while building the manifest. On systems with the source fully
extracted, verify regular files with `shasum -a 256 -c SHA256SUMS.txt`.

Earlier root checksum snapshots are preserved under
`docs/archive_20261111_checksums/`. Any `SUBMISSION_SHA256.tsv` and older nested
archive manifests describe historical submission material, not this development
update. GitHub-generated ZIPs and R source packages have separate archive hashes
in the delivery bundle; their contents need not match a complete repository tree.
