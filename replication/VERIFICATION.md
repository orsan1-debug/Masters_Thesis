# Verification

Written by `replication/verification/make_verification.R` from the logs of the last run of
`replication/run_checks.R` (16 Oct 2026 01:59) on macOS Tahoe 26.6.2, R version 4.6.1, balnet 0.0.4 and glmnet 5.1.

## The checks

| Link | Check | Result | Source |
|---|---|---|---|
| Code = the code that ran | `tests/run_tests.R`: the replication code against the code under `R/`, `runs/` and `handoff/`, and against the stored numbers | 9 files, 766 expectations, 0 failures | `logs/tests.txt` |
| Code = the code that ran, without `results/` | `tests/run_tests.R` in a copy of the package without `results/`, the tests reading the rows they use from `tests/fixtures/` | 9 files, 449 expectations, 0 failures, 3 skipped (the recomputations in `test_measures.R`, which need `results/`) | `logs/tests_without_results.txt` |
| Raw results = a rerun of the code | `smoke_test.R`: 15 replications of the first and last cell of every grid, n up to 50 000, rerun through `runs/` and compared row by row. Not part of `run_checks.R` (`SMOKE=1` adds it) | 54 stored files: 48 identical, 2 last digits, 4 means within 0.5 MCSE | `smoke_report.csv` (15 Oct 2026) |
| Summaries = raw results | `verify_summaries.R`: the Study 1 section summaries and the coverage table recomputed from `results/` | 17 recomputed, largest difference 1.9e-13 | `logs/summaries.txt` |
| Summaries = raw results | `handoff/scripts/build_summaries.R`: the appendix entries, the Study 2 summaries and `wc_compare` by a second computation | 36 of 36 passed | `output/thesis_summaries/verification.md` |
| PDF = a render from the summaries | `verify_render.R`: the thesis rendered from the stored summaries in a clean copy. `verification/compare_render.R`: the numbers of every table and figure, every number in the text and the words of every page, against `output/rendered/thesis/AAAAA-Working-Thesis.pdf` | 52 of 52 displays with identical numbers (3 665 numbers). In the whole text 25 numbers in one PDF and not the other. The 56 pages by their words: 44 identical, 3 reordered, 3 URL split elsewhere, 6 page break elsewhere. Both listed below | `verification/render_compare.csv`, `verification/render_numbers.csv`, `verification/render_text.csv` |

The steps of the last run:

```
setup      ok                   1.8 secs
tests      ok                   5.1 mins
summaries  ok                   7.3 mins
render     failed, status 1     25 secs
manifest   ok                   16 secs
Error: see replication/logs/ for the failed steps
Execution halted
```

The render step stops where there is no TeX. How the render was then made is below.

## Each table and figure

Page is the page of the stored PDF. Summary = results gives the largest difference between the
recomputed and the stored numbers, relative to the stored number where that exceeds 1 (above
1e-10 stops `verify_summaries.R`). Rerun = stored counts the `smoke_test.R` verdicts of the results
files behind the display. Render = PDF compares the numbers of the display's rows, in order.

| Display | Page | qmd chunk | Summary = results | Rerun = stored | Render = PDF |
|---|---|---|---|---|---|
| Table 1 | 8 | `tbl-spec` (line 223) | typed in the qmd |  | identical (no numbers) |
| Table 2 | 9 | `tbl-grids` (line 240) | typed in the qmd |  | identical (89 numbers) |
| Table 3 | 10 | `tbl-grids-s2` (line 297) | typed in the qmd |  | identical (50 numbers) |
| Table 4 | 15 | `tbl-n-adv` (line 528) | n_both 1.1e-14 | 4 identical | identical (108 numbers) |
| Figure 1 | 16 | `fig-n-refit` (line 549) | misspec_npath 1.9e-13 | 1 identical. Not rerun: aipw_misspec_n_dgp2n1mill | identical (20 numbers) |
| Figure 2 | 17 | `fig-n-misspec` (line 566) | misspec_npath 1.9e-13 | 1 identical. Not rerun: aipw_misspec_n_dgp2n1mill | identical (28 numbers) |
| Table 5 | 19 | `tbl-ov-adv` (line 630) | overlap 1.2e-13 | 4 identical | identical (24 numbers) |
| Table 6 | 19 | `tbl-ov-adv-hi` (line 647) | overlap 1.2e-13 | 4 identical | identical (30 numbers) |
| Table 7 | 20 | `tbl-ov-adv-lo` (line 660) | overlap 1.2e-13 | 4 identical | identical (24 numbers) |
| Figure 3 | 20 | `fig-ov-adv` (line 689) | overlap 1.2e-13 | 4 identical | identical (164 numbers) |
| Table 8 | 21 | `tbl-pn-adv` (line 765) | p_axis 3.8e-15 | 2 means within 0.5 MCSE | identical (40 numbers) |
| Figure 4 | 22 | `fig-pn-rmse` (line 781) | p_axis 3.8e-15. app_p_axis passed (second computation) | 2 means within 0.5 MCSE | identical (24 numbers) |
| Figure 5 | 23 | `fig-sp-change` (line 824) | sparsity 8.9e-16 | 4 identical | identical (9 numbers) |
| Table 9 | 25 | `tbl-sig-adv` (line 860) | signal 1.9e-15 | 2 identical | identical (6 numbers) |
| Figure 6 | 26 | `fig-n-coverage` (line 987) | coverage_empse:n 1.9e-15 | 2 identical. Not rerun: aipw_misspec_n_dgp2n1mill | identical (12 numbers) |
| Table 10 | 27 | `tbl-s2-counts` (line 1486) | cv_overlap, cv_noise, cv_noise_x_overlap, cv_alpha_x_overlap, cv_n, cv_nbig, cv_p, cv_misspec, cv_sparsity, cv_sparsity_misspec passed (second computation) | 21 identical, 1 last digits | identical (132 numbers) |
| Figure 7 | 28 | `fig-s2-path` (line 1498) | cv_noise_x_overlap passed (second computation) | 2 identical | identical (27 numbers) |
| Figure 8 | 29 | `fig-s2-lead` (line 1514) | cv_noise_x_overlap passed (second computation) | 2 identical | identical (12 numbers) |
| Table 11 | 30 | `tbl-concl` (line 1567) | n_both 1.1e-14, overlap 1.2e-13, p_axis 3.8e-15, sparsity 8.9e-16, signal 1.9e-15 | 14 identical, 2 means within 0.5 MCSE | identical (25 numbers) |
| Appendix Table 1 | 38 | `atbl-q1-adv` (line 1650) | correct_n 1.6e-15 | 2 identical | identical (12 numbers) |
| Appendix Table 2 | 38 | `atbl-app-n-diag-m-lin` (line 1686) | app_misspec_n passed (second computation) | 2 identical | identical (50 numbers) |
| Appendix Table 3 | 39 | `atbl-n-path` (line 1700) | misspec_npath 1.9e-13 | 1 identical. Not rerun: aipw_misspec_n_dgp2n1mill | identical (675 numbers) |
| Appendix Table 4 | 40 | `atbl-app-avg-n-m` (line 1740) | misspec_n 1.1e-14 | 2 identical | identical (8 numbers) |
| Appendix Table 5 | 40 | `atbl-app-n-ref` (line 1758) | n_both 1.1e-14 | 4 identical | identical (36 numbers) |
| Appendix Figure 1 | 41 | `afig-ov-levels` (line 1787) | drawn in the chunk from R/dgp.R |  | identical (12 numbers) |
| Appendix Table 6 | 41 | `atbl-ov-levels` (line 1863) | closed form, tests/test_dgp.R |  | identical (38 numbers) |
| Appendix Table 7 | 41 | `atbl-app-ov-full-c` (line 1906) | app_overlap_full passed (second computation) | 4 identical | identical (30 numbers) |
| Appendix Table 8 | 42 | `atbl-app-ov-full-m` (line 1919) | app_overlap_full passed (second computation) | 4 identical | identical (52 numbers) |
| Appendix Figure 2 | 42 | `afig-ov-coverage` (line 1932) | coverage_empse:overlap 3.7e-15 | 4 identical | identical (12 numbers) |
| Appendix Figure 3 | 43 | `afig-app-nc` (line 1974) | nc 2e-15 | 2 identical | identical (12 numbers) |
| Appendix Table 9 | 43 | `atbl-app-nc-full-m` (line 1987) | app_nc passed (second computation) | 2 identical | identical (50 numbers) |
| Appendix Table 10 | 44 | `atbl-app-nc-trunc` (line 2001) | nc 2e-15 | 2 identical | identical (125 numbers) |
| Appendix Figure 4 | 45 | `afig-app-pn-bias` (line 2028) | p_axis 3.8e-15. app_p_axis passed (second computation) | 2 means within 0.5 MCSE | identical (24 numbers) |
| Appendix Table 11 | 45 | `atbl-app-pn-full-c` (line 2040) | app_p_axis passed (second computation) | 2 means within 0.5 MCSE | identical (305 numbers) |
| Appendix Table 12 | 47 | `atbl-app-pn-full-m` (line 2051) | app_p_axis passed (second computation) | 2 means within 0.5 MCSE | identical (245 numbers) |
| Appendix Table 13 | 48 | `atbl-app-pn-diag-c-lin` (line 2064) | app_p_axis passed (second computation) | 2 means within 0.5 MCSE | identical (65 numbers) |
| Appendix Table 14 | 48 | `atbl-app-pn-diag-m-lin` (line 2076) | app_p_axis passed (second computation) | 2 means within 0.5 MCSE | identical (65 numbers) |
| Appendix Table 15 | 49 | `atbl-app-sp-cells` (line 2109) | app_sparsity, app_sparsity_equi passed (second computation) | 2 identical | identical (151 numbers) |
| Appendix Table 16 | 49 | `atbl-app-sp-bound` (line 2120) | sparsity 8.9e-16 | 4 identical | identical (82 numbers) |
| Appendix Table 17 | 50 | `atbl-app-sp-exact` (line 2135) | sparsity 8.9e-16 | 4 identical | identical (105 numbers) |
| Appendix Table 18 | 50 | `atbl-app-sp-trunc` (line 2155) | sparsity 8.9e-16 | 4 identical | identical (17 numbers) |
| Appendix Table 19 | 51 | `atbl-app-sp-tan` (line 2175) | app_sparsity_tan passed (second computation) | 2 identical | identical (284 numbers) |
| Appendix Table 20 | 52 | `atbl-app-sp-share` (line 2187) | closed form, tests/test_dgp.R |  | identical (9 numbers) |
| Appendix Table 21 | 52 | `atbl-app-sig-full-c` (line 2223) | app_signal passed (second computation) | 2 identical | identical (152 numbers) |
| Appendix Table 22 | 53 | `atbl-app-sig-diag-c-lin` (line 2236) | app_signal passed (second computation) | 2 identical | identical (18 numbers) |
| Appendix Figure 5 | 53 | `afig-sig-coverage` (line 2252) | coverage_empse:signal 1.9e-15 | 2 identical | identical (18 numbers) |
| Appendix Figure 6 | 54 | `afig-app-wc` (line 2302) | wc_compare passed (second computation) | 1 identical, 1 last digits | identical (12 numbers) |
| Appendix Table 23 | 54 | `atbl-app-wc` (line 2313) | wc_compare passed (second computation) | 1 identical, 1 last digits | identical (40 numbers) |
| Appendix Table 24 | 55 | `atbl-s2-selectors` (line 2324) | cv_overlap, cv_noise, cv_noise_x_overlap, cv_alpha_x_overlap, cv_n, cv_nbig, cv_p, cv_misspec, cv_sparsity, cv_sparsity_misspec passed (second computation) | 21 identical, 1 last digits | identical (40 numbers) |
| Appendix Table 25 | 55 | `atbl-s2-alpha` (line 2335) | cv_alpha_x_overlap passed (second computation) | 2 identical | identical (45 numbers) |
| Appendix Figure 7 | 56 | `afig-s2-design` (line 2346) | cv_alpha_x_overlap, cv_n, cv_nbig, cv_p passed (second computation) | 9 identical, 1 last digits | identical (23 numbers) |
| Appendix Figure 8 | 56 | `afig-s2-spec` (line 2360) | cv_noise_x_overlap, cv_misspec, cv_sparsity, cv_sparsity_misspec passed (second computation) | 8 identical | identical (29 numbers) |

## How the render was made

`thesis_render.pdf` was made on 15 Oct 2026. `verify_render.R` ran on the Mac, which has no TeX,
so Quarto stopped after writing the LaTeX, kept in `render_tex/`. That LaTeX was compiled with
XeLaTeX (TeX Live 2022) in a Linux virtual machine on the same Mac, with the Mac's Times New
Roman. The stored PDF was compiled elsewhere, with a later TeX Live (its xdvipdfmx is dated
2026-03-17, the virtual machine's 2021-06-09). In a few places a line breaks at another point,
so the comparison goes by numbers and words rather than layout.

## Numbers in one PDF and not the other

Counted over the whole text with page numbers dropped. Pages are where the number occurs in that PDF.

| Number | Found in | Times | Pages |
|---|---|---|---|
| `https://doi.org/10.1111/rssb.12268` | stored PDF | 1 | 34 |
| `https://doi.org/10.1097/EDE` | stored PDF | 1 | 34 |
| `0000000000001644` | stored PDF | 1 | 34 |
| `https://doi.org/10.1111/` | stored PDF | 1 | 35 |
| `rssb.12027` | stored PDF | 1 | 35 |
| `https://doi.org/10.1002/sim` | stored PDF | 1 | 35 |
| `8086` | stored PDF | 1 | 35 |
| `10.1093/biomet/asaa020` | stored PDF | 1 | 35 |
| `https://doi.org/10.1093/` | stored PDF | 1 | 36 |
| `biomet/70.1.41` | stored PDF | 1 | 36 |
| `https://doi.org/10` | stored PDF | 1 | 36 |
| `1093/biomet/asx069` | stored PDF | 1 | 36 |
| `https://doi.org/10.1111/rssb` | render | 1 | 34 |
| `.12268` | render | 1 | 34 |
| `https://doi.org/10.1` | render | 1 | 34 |
| `097/EDE.0000000000001644` | render | 1 | 34 |
| `https://doi.org/10.111` | render | 1 | 35 |
| `1/rssb.12027` | render | 1 | 35 |
| `https://doi.org/10.1002/` | render | 1 | 35 |
| `sim.8086` | render | 1 | 35 |
| `//doi.org/10.1093/biomet/asaa020` | render | 1 | 35 |
| `https://doi.org/10.1093/bi` | render | 2 | 36 |
| `omet/70.1.41` | render | 1 | 36 |
| `omet/asx069` | render | 1 | 36 |

## Pages whose words are not identical

Each page's words read in order, with words split at the end of a line joined. Reordered: the
same words read in another order, as the labels of a figure drawn in another order. URL split
elsewhere: the same letters, a long URL or DOI broken at another point. Page break elsewhere:
the same letters once the next or previous page is added, a line on the other side of a page break.

| Page | Verdict | Only in the stored PDF | Only in the render |
|---|---|---|---|
| 9 | page break elsewhere | `model,` |  |
| 10 | page break elsewhere |  | `model,` |
| 16 | page break elsewhere |  | `part` |
| 17 | page break elsewhere | `part` |  |
| 20 | reordered |  |  |
| 21 | page break elsewhere |  | `over at = holds to AIPW p/n from with lead (Appendix Its 0.5 ...` |
| 22 | page break elsewhere | `p = Figure p/n at growing with from 0.5 (Appendix 4). Its le ...` |  |
| 34 | URL split elsewhere | `https://doi.org/10.1111/rssb.12268. https://doi.org/10.1097/ ...` | `https://doi.org/10.1111/rssb .12268. https://doi.org/10.1 09 ...` |
| 35 | URL split elsewhere | `https://doi.org/10.1111/ rssb.12027. https://doi.org/10.1002 ...` | `https: https://doi.org/10.111 1/rssb.12027. https://doi.org/ ...` |
| 36 | URL split elsewhere | `https://doi.org/10.1093/ biomet/70.1.41. https://web.stanfor ...` | `https://doi.org/10.1093/bi https://doi.org/10.1093/bi omet/7 ...` |
| 38 | reordered |  |  |
| 42 | reordered |  |  |

## Not rerun here

- The full simulations. The smoke test reruns the first and last cell of every grid.
- The n = 1 million cell (`aipw_misspec_n_dgp2n1mill`), too large for the smoke test.
  `tests/test_runs.R` regenerates its first draw, and the summaries that read it are
  recomputed from the stored file.

