[READ] Horvitz–Thompson and Hájek across every paper in the project — per-paper quote compendium, dot points, for me to write the section from.

Read first, in this order:
1. thesis_state.md — rules in "decided": deliverables are dot-point source notes (not prose), saved as a NEW file in the repo root `C:\Users\otisr\Documents\Thesis 2026\Masters_Thesis\` (request folder access; the repo is otherwise read-only) and mirrored as a project doc under `claude/`.
2. claude/methods_estimators_dotpoints_2026-09-16.md, sections 1–2 — what has already been found for HT and Hájek. Start from it, re-verify every marker against the paper file, and extend it; do not redo it from scratch.
3. paper_facts.md — the entry format and the page-marker conventions.

Task. For EVERY paper file in the project, read the whole file with project_read (RAG search misses passages — do not rely on project_search alone) and pull out every passage that defines, names, attributes, or characterises (a) the Horvitz–Thompson / unnormalised IPW estimator and (b) the Hájek / normalised / ratio / "stabilized-weight" IPW estimator — ATE and ATT forms, and the DR-HT / DR-Hájek variants where a paper gives them. For each passage record, as one dot point:
- a short verbatim quote (≤ 2 sentences, in quotation marks; never a paraphrase presented as a quote);
- the equation, transcribed in Quarto display maths in the thesis notation ($W_i$, $Y_i$, $Y_i(1)$, $e(x)$, $\hat e(X_i)$, $\tau$), with the paper's own symbols noted once (e.g. "paper writes $Z_i$, $T_i$, $\pi$");
- the name the paper uses for it (HT, IPW, IPW-POP, ratio IPW, rIPW, Hájek, "linear weighted difference", stabilized weights, $T_{\mathrm{lin}}$ …);
- whom the paper attributes it to (Horvitz & Thompson 1952, Hájek 1971, Hirano–Imbens–Ridder 2003, Rosenbaum 1987, …), tagged [UNSOURCED] when that origin paper is not in the project;
- any property or result claimed (unbiasedness, invariance to $Y + c$, stability, variance, sample boundedness, efficiency, simulation results with the table number);
- the page marker in that paper's convention (below).

Page-marker conventions (verified 2026-09-16; use these, and state the convention once per paper block):
- Ding 2023 "A First Course in Causal Inference" = the file `AIPW-IPTW-Covariate-Balancing-LLM.md` (ch. 11–14 excerpts): cite printed p.; printed = PDF − 26.
- Tan 2020 RCAL (`Tan_2020_Regularized_calibrated_...`): arXiv v1 PDF pages.
- Tan 2020 AoS (`Tan - 2020 - Model-assisted ...`): [CHECK] whether the file carries journal or PDF pages; state which you use.
- Wager 2025 (ch. 2–3–7 excerpt): printed p. (= PDF − 1).
- Wang & Zubizarreta 2020: printed p. (file gives "PDF page N = printed page M").
- Chattopadhyay 2020: printed p. (= PDF + 3226).
- Kang & Schafer 2007: printed p. (= PDF + 522).
- Imai & Ratkovic 2014: printed p. (= PDF + 242).
- Rosenbaum & Rubin 1983: printed p. (= PDF + 40).
- Sverdrup & Hastie 2026, Shang 2025, Zhao 2017, Ben-Michael 2021, Bruns-Smith 2025, Keele 2025, Wyss 2026, Hainmueller 2012, Sun & Tan 2020, Ertefaie 2023, Austin & Stuart 2015, Friedman 2010, Morris 2019: use the anchor the file provides (PDF page and printed/journal page where both are given); state which.

Papers, in this order (all 21 in the project): Rosenbaum & Rubin 1983; Kang & Schafer 2007; Friedman 2010; Hainmueller 2012; Imai & Ratkovic 2014; Austin & Stuart 2015; Zhao 2017; Morris 2019; Chattopadhyay 2020; Sun & Tan 2020; Tan 2020 RCAL; Tan 2020 AoS; Wang & Zubizarreta 2020; Ben-Michael 2021; Ding 2023 (`AIPW-IPTW-Covariate-Balancing-LLM.md`); Ertefaie 2023; Bruns-Smith 2025; Keele 2025; Shang 2025; Wager 2025; Wyss 2026; Sverdrup & Hastie 2026. Where a paper says nothing about HT or Hájek, write one line: "nothing — searched for: Horvitz, Thompson, Hájek/Hajek, ratio, normalized/normalised, stabilized, IPW-POP, $T_{\mathrm{lin}}$, weighted difference".

Already found (re-verify, then build on): Ding §11.2.2 p. 160 (both forms; "Horvitz and Thompson (1952) proposed it in survey sampling and Rosenbaum (1987a) used in causal inference"; "Hajek estimator due to Hájek (1971)"), Prop. 11.1 p. 160 (HT not invariant to $Y+c$), p. 161 (blow-up near 0/1, truncation), §13.2 p. 184 (ATT HT and Hájek with odds weights), Prop. 14.1 p. 197 (Hájek = WLS coefficient, Imbens 2004), Problem 12.2 p. 178 (Hájek-form DR, Robins et al. 2007); Chattopadhyay §2.4 (i)–(iv) pp. 3231–3232, eq. (8) p. 3234, §3.3 p. 3235 (sample boundedness), fn ** p. 3243 (DR-HT ≈ DR-Hájek), Tables 2–3 pp. 3244, 3246; Kang & Schafer §2.1 p. 528 (HT origin; "Precision is often enhanced if we use a denominator of $\sum_i t_i\hat\pi_i^{-1}$ rather than $n$"; IPW-POP eq. (3)), p. 529 and Table 1 ("bias and RMSE actually get worse as the sample size grows!"); Imai & Ratkovic §3.1 p. 251 (HT = Horvitz & Thompson 1952; IPW = Hirano et al. 2003; formulas), Table 1 p. 253; Tan RCAL PDF p. 6 ($\hat\mu^1_{\mathrm{IPW}}$, $\hat\mu^1_{\mathrm{rIPW}}$), p. 7 (the two coincide under calibration), p. 8 eq. (10) (weights sum to $n$), §4 p. 20 (ratio form used throughout), §6 pp. 31–32 (ATT $\hat\nu^0_{\mathrm{IPW}}$, $\hat\nu^0_{\mathrm{rIPW}}$); Wager eqs. (2.11)–(2.12) p. 22, Theorem 2.2 p. 23, eq. (2.16) p. 25, eq. (2.17) p. 25, ch. 2 notes p. 28, eq. (7.1) p. 84, eq. (7.11) p. 87, eq. (7.21) p. 95; Shang §3.1 eqs. (6)–(7) PDF p. 5 ("Hájek estimator (also called stabilized weights)", refs 46–47); Sverdrup & Hastie PDF p. 3 (identification $E[W_iY_i/e(X_i)]$), Remark 2 PDF p. 5, ATT weights PDF p. 7; Keele PDF pp. 3–4 ($\hat\Delta_w$ weighted-difference form; Horvitz & Thompson 1952 cited); Hainmueller §2.1 (odds weights $d_i = \hat p/(1-\hat p)$, HT survey analogy) — page marker still needed.

Deliverable (one .qmd, Quarto maths kept, dot points only, no prose paragraphs):
1. One block per paper in the order above, sub-headed "HT", "Hájek", "Names and attributions", "Properties / results", or the one-line "nothing" statement.
2. A cross-paper table: paper | name used for HT | name used for Hájek | attribution given | forms given (ATE / ATT / DR) | page(s).
3. An [UNSOURCED] register: every origin paper cited by the project papers but not in the project, with which project paper (and page) makes the attribution.
4. A [CHECK] list.
5. Proposed paper_facts.md entries in its format ("[Author Year §/eq., p. N] claim — verified against … — used for: Methods; …") — list them; do not write to paper_facts.md until I say.

Rules: every dot point carries a page marker; verify before flagging — read the file, never assume; no interpretation and no numbers from results/; notation as in thesis/Final draft.qmd; one term per concept (HT for the unnormalised form, Hájek for the normalised form — record the papers' own names beside them). Save as `HT_Hajek_paper_sweep_<today>.qmd` in the repo root and as `claude/HT_Hajek_paper_sweep_<today>.md` in the project; show the file in chat. Return the file, then a ≤ 100-word note on what could not be verified.
