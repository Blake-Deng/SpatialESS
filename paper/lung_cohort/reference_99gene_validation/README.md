# Lung reference validation: 99-gene input

This numerical-fidelity analysis uses a 99-gene signalling preparation across
45 inputs. It is separate from the 343-gene pulmonary-fibrosis cohort model.
Code and diagnostic logs are available here; download and mount the result
tables using the [reproduction guide](../../../docs/REPRODUCE.md).

There are 42 evaluable pairs: 41 use the unmodified reference, and one uses a
reference repaired only to preserve matrix dimensions. Three reference
attempts failed. Of the 42 evaluable pairs, 41 meet the strict criterion;
VUILD142MA retains one permutation-P difference of 0.05 without a changed
significance decision.

Tables are under `current_release/` after mounting. The patch, failed-attempt
logs and diagnostics document reference execution and are not successful
performance observations. These results do not establish full numerical
agreement for the distinct 343-gene cohort preparation.
