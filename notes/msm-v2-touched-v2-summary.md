# MSM v2 touched-v2 summary

## Scope

Compared four variants on GCP `c3-standard-22` with `taskset -c 0-15`, release build (`-O3 -DNDEBUG`, AVX512 IFMA), and 3 repetitions with median CPU time:

- `baseline` (commit `2a758b3` + current benchmark harness)
- `epoch` (commit `d501d02` logic)
- `touched` (first touched-bucket reset version)
- `touched_v2` (touched-bucket + cached conflict resolve + buffer reuse)

Source files used for comparison:

- `notes/baseline_median_runs.json`
- `notes/epoch_median_runs.json`
- `notes/touched_median_runs.json`
- `notes/touched_v2_median_runs.json`

## Key result highlights vs baseline

### BM_VROOM_MSM_V2

- `2^18`: `+5.19%`
- `2^20`: `+3.49%`
- Small sizes mostly near-flat (`-0.09%` to `+0.37%`)

### BM_VROOM_MSM_V2_Parallel

- `2^20`: `+14.69%`
- `2^18`: `+0.92%`
- Mid/smaller sizes mixed (some small regressions)

## Notes

- The largest improvements are at large sizes, where memory and bucket management costs dominate.
- Parallel results still show noise at some sizes; for final publication-quality numbers, run 5-7 repetitions and report median plus spread.
