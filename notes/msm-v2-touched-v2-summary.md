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

## Latest comparison benchmark table

Comparison columns:

- `2a758b3_ms`: baseline median CPU time
- `50722be_ms`: touched_v2 median CPU time
- `latest_ms`: latest parallel tuning attempt median CPU time
- `latest_vs_50722be_pct`: positive means latest is faster than `50722be`
- `latest_vs_2a758b3_pct`: positive means latest is faster than `2a758b3`

| benchmark | exp | 2a758b3_ms | 50722be_ms | latest_ms | latest_vs_50722be_pct | latest_vs_2a758b3_pct |
|---|---:|---:|---:|---:|---:|---:|
| BM_VROOM_MSM_V2 | 8  | 3.055 | 3.050 | 3.057 | -0.22% | -0.07% |
| BM_VROOM_MSM_V2 | 10 | 9.562 | 9.540 | 9.549 | -0.10% | +0.14% |
| BM_VROOM_MSM_V2 | 12 | 31.209 | 31.147 | 31.374 | -0.73% | -0.53% |
| BM_VROOM_MSM_V2 | 14 | 108.990 | 109.091 | 109.735 | -0.59% | -0.68% |
| BM_VROOM_MSM_V2 | 16 | 386.353 | 384.906 | 385.440 | -0.14% | +0.24% |
| BM_VROOM_MSM_V2 | 18 | 1530.436 | 1451.041 | 1452.697 | -0.11% | +5.08% |
| BM_VROOM_MSM_V2 | 20 | 6211.082 | 5994.271 | 6043.485 | -0.82% | +2.70% |
| BM_VROOM_MSM_V2_Parallel | 8  | 3.059 | 3.051 | 3.058 | -0.24% | +0.01% |
| BM_VROOM_MSM_V2_Parallel | 10 | 1.596 | 1.600 | 1.609 | -0.53% | -0.77% |
| BM_VROOM_MSM_V2_Parallel | 12 | 4.216 | 4.353 | 4.374 | -0.46% | -3.75% |
| BM_VROOM_MSM_V2_Parallel | 14 | 13.824 | 13.829 | 13.858 | -0.20% | -0.24% |
| BM_VROOM_MSM_V2_Parallel | 16 | 54.492 | 54.749 | 55.198 | -0.82% | -1.29% |
| BM_VROOM_MSM_V2_Parallel | 18 | 228.207 | 226.106 | 229.284 | -1.41% | -0.47% |
| BM_VROOM_MSM_V2_Parallel | 20 | 888.680 | 758.149 | 911.754 | -20.26% | -2.60% |
