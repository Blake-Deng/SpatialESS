# Minimal example

After installing SpatialESS and this integration package:

```bash
./examples/minimal/run_demo.sh
```

The example creates a 24-cell SPARKLE-like corrected H5AD file, exports only
the two genes required by a synthetic LR pair, runs ten deterministic
permutations, and writes a sparse input bundle plus communication records to
`examples/minimal/output/`.
