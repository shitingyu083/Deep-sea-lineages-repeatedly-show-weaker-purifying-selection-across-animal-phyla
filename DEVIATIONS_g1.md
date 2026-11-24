# DEVIATIONS_g1.md — append-only

Deviations from `PREREG_G1_gene_family_convergence.md`
(md5 `8a19c101207c7016ee1affacb1ee8943`, frozen 2026-11-03).

PREREG_G1 itself is never edited. Every departure from it, every discovered defect in the
statistic or its null, and every decision taken after seeing a number is recorded here, in
the order it happened. Nothing below was written before the relevant number was seen; each
entry states what was seen and what was then done.

---

## D1 — 2026-11-03 — the permutation battery does not preserve K

**Seen:** `g1_kdiagnostic.py` scanned every tree on disk (obs, `null1..5`, `nullp1..10`)
per clade and recorded K = number of maximal deep subtrees, deep-branch count, deep-tip
count, R, `rate = R/(N_fam*(K-1))`, and a K-invariant binomial recurrence count.

| battery      | K matched to obs?        | deep branches matched? | deep tips matched? |
|--------------|--------------------------|------------------------|--------------------|
| `nullp1..10` | **no** in 4 of 6 clades  | no                     | yes                |
| `null1..5`   | yes in 4 of 6 clades     | no                     | no                 |

Measured: Annelida obs K=8 vs perm 11–20; Mollusca_rescue obs 18 vs perm 26–32;
Cnidaria obs 12 vs perm 13–20; Actinopteri_m_rescue obs 36 vs perm 25–40.
Deep branch counts move too: Annelida 34→22, Mollusca 58→44, Cnidaria 30→24,
Echinodermata 6→12.

**Why it matters:** R = Σ_f max(0, max(e_f,c_f) − 1) increases mechanically with K. A raw-R
comparison between trees of different K is invalid — it partly tests "how much deep-labelled
tree does this tree have", not "do directions repeat across origins".

**Done:** R is still reported unnormalised because the PREREG fixes it as the primary
statistic, but no p-value is taken from raw R against a K-unmatched battery. Every p below is
taken on the PREREG-normalised `rate = R/(N_fam*(K-1))`, or on the K-matched `null1..5`
battery where it exists. Both are reported side by side; neither is selected after the fact.

---

## D2 — 2026-11-03 — the structure battery does not preserve structure either

**Seen:** same diagnostic. `null1..5` preserve K but change the number of deep-labelled
branches (Echinodermata 6 → 12, i.e. doubled) and the number of deep-labelled tips.

**Done:** `null1..5` is labelled a **K-matched** null, not a structure-preserving null. It is
reported as a secondary check only, and its own mismatch is stated wherever it is used.

---

## D3 — 2026-11-03 — the count of BH-significant families is a power floor, not evidence

**Seen:** `n_binom_bh ≈ 0` in every clade.

**Diagnosis:** the smallest attainable two-sided exact binomial p for a family with K origins
is `2·0.5^K`. Against the BH threshold `0.05/N`:

| K  | p_min = 2·0.5^K | BH threshold (Annelida) | detectable? |
|----|-----------------|--------------------------|-------------|
| 6  | 3.1e-2          | 4.3e-6                   | impossible  |
| 8  | 7.8e-3          | 4.1e-6                   | impossible  |
| 12 | 4.9e-4          | 5.1e-6                   | impossible  |
| 36 | 2.9e-11         | 2.9e-6                   | possible    |

For K < 24 no family can be detected however strong its signal. Only Actinopteri_m_rescue
(K=36) is in the detectable regime, and there obs `n_binom_bh = 1` vs nulls 0–1.

**Done:** the near-zero count is reported as an arithmetic power floor and is explicitly **not**
offered as evidence of absence. No claim in the paper rests on it.

---

## D4 — 2026-11-03 — a sign-flip null was built and then REJECTED as miscalibrated

**Built:** `g1_signflip.py`, a structure-preserving paired sign-flip (holds component
structure, sizes and |Δ| exactly fixed; flips only the sign alignment across origins;
B=2000, seed 20261103).

**Seen:** p = 1/2001 = 0.0005 in **all six clades**.

**Why that result was not reported as convergence:** the Δ sign marginals are strongly
contraction-biased in every clade (Chondrichthyes +4322/−8855; Mollusca +25374/−64244;
Cnidaria +21606/−14768; Echinodermata +5297/−6549; Actinopteri +30139/−84513). A symmetric
sign-flip null destroys that global bias, so the test measures "contractions outnumber
expansions", not "directions agree across independent origins". Independent corroboration
that the null is displaced by construction: the permutation R range brackets the observation
(e.g. Chondrichthyes perm [2777, 5384] vs obs 3676) whereas the sign-flip null predicts 3175.8.

**Done:** this test is **discarded**. Its p-values are not reported anywhere and must not be
cited. The files remain on disk as a record of the diagnostic.

---

## D5 — 2026-11-03 — `_rescue` clades have no `prep_report.txt`; N0 borrowed from the sibling

**Seen:** `g1_family_convergence.py --bridge` died with `cannot locate N0 for Annelida`
(cause: `prep_report.txt` field 2 reads `"N0 /abs/path/N0.tsv"`, and the leading `N0 ` key
was not stripped — fixed) and then with `cannot locate N0 for Mollusca_rescue`.

**Why the rescue dirs have no report:** they were prepared as re-runs of a sibling clade.

**Verification performed before borrowing:**

- species name sets are **identical**: Mollusca vs Mollusca_rescue 79/79 names, diff empty;
  Actinopteri_m vs Actinopteri_m_rescue 75/75 names, diff empty;
- the family IDs are drawn verbatim from the sibling's N0:
  Mollusca_rescue **16712 / 16712** family IDs present in the Mollusca N0;
  Actinopteri_m_rescue **17157 / 17157** present in the Actinopteri_m N0.

**Done:** `n0_path()` gained an explicit `SIBLING_N0` fallback
(`Mollusca_rescue → Mollusca`, `Actinopteri_m_rescue → Actinopteri_m`) that logs every borrow.
This affects the eggNOG bridge only (an annotation lookup); no copy number or Δ value is read
from N0, so the within-clade statistics C1 are untouched.

---

## D6 — 2026-11-03 — C2's frozen null is miscalibrated; a marginal-matched null was added

**PREREG C2 (frozen, verbatim):** "The cross-clade count of eggNOG-bridged families that recur
in the same direction in ≥ 2 independent clades, against the same permutation null restricted
to bridged families." This was implemented literally (`g1_c2.py`).

**Seen (frozen null, run as specified):** S_obs = 5726; permutation null
[4843, 4871, 4894, 4916, 4957, 5022, 5076, 5195, 5235, 5266]; p = 1/11 = **0.091** — the
resolution floor of a 10-replicate battery. K-matched `null1..5` (4 clades): obs 2732 vs
[1771, 1775, 1795, 1828, 1849], p = 1/6 = **0.167**.

**Diagnosis — the excess is a tie-rate artefact.** Measured share of bridged families that
receive a non-neutral direction, obs vs 10-replicate mean:

| clade | obs | perm |
|---|---|---|
| Annelida | 91.6 % | 81.6 % |
| Chondrichthyes | 97.2 % | 86.4 % |
| Mollusca_rescue | 91.7 % | 86.4 % |
| Cnidaria | 90.9 % | 81.4 % |
| Echinodermata | 94.9 % | 67.9 % |
| Actinopteri_m_rescue | 94.4 % | 91.4 % |

A neutral family (e=f) cannot contribute to "two clades agree", so a null with more ties
mechanically produces a smaller S. This is D1's K-confound arriving by another route: larger
K pushes e and c toward ties.

**Added (`g1_c2_marginal.py`):** a marginal-matched independence null. Each clade's observed
(+,−,0) marginal over its bridged families and each OG's clade membership are held fixed;
only the alignment of directions *between* clades is destroyed. B=2000, seed 20261103.

**Seen (marginal-matched):** S_obs = 5689; null mean 5694.4, sd 33.4, range [5592, 5807];
p = **0.567**, z = −0.16. Pairwise concordance among doubly-called (OG, clade) pairs:
observed **0.597** vs conditional expectation **0.600**.

**Done:** the frozen permutation p (0.091 / 0.167) is reported as specified, with its tie-rate
bias stated next to it, and C2 is **not** claimed. The marginal-matched test is the one that
decides the claim, and it is null. This is recorded here rather than silently substituting one
null for another.

**Bug found and fixed inside this entry, disclosed for completeness:** the first
marginal-matched run assigned the `P(+)` interval to `−1` (cuts applied in the wrong order),
producing a spurious z = 164. Caught because the pairwise sub-statistic disagreed with the
aggregate (0.597 observed vs 0.590 expected — impossible pairing). After the fix the two
sub-statistics agree with each other (0.597 vs 0.600) and with the aggregate. The conditional
expectation `(p+a·p+b + p−a·p−b)/((p+a+p−a)(p+b+p−b))` was also corrected at the same time.

---

## D7 — 2026-11-03 — outcome, stated before the paper text is written

**C1:** no clade beats its null on any valid test. On the normalised rate, four of six clades
sit exactly at the 0.091 floor of a 10-replicate battery (which is the floor, not
significance). Echinodermata is the only clade whose permutation K matches exactly; it gives
rate p = 0.091, but its K-matched structure null gives **p = 1.000** with the null producing
*more* recurrence than the observation.

**C2:** null. Observed shared-direction count (5689) is indistinguishable from the
marginal-matched expectation (5694.4, p = 0.567); pairwise concordance 0.597 vs 0.600.

**C3 (the pre-committed dissociation statement):** PREREG_G1 §1 fixed in advance that if the
family-size level does not show convergence, the paper "reports a dissociation: no convergence
at either the rate level (R6-5) or the family-size level. That is a clean negative on the
'recurrent genomic convergence' hypothesis and must be stated as such, not buried."

That is the outcome. It is reported as a clean negative, not reframed.

---

## D8 — 2026-11-03 — E4 (timing) run twice; the first run's numbers are discarded, and the second run's control overturns its own signal

PREREG_G1 §8 registered E4 as exploratory, verbatim: "E4 Timing: whether recurring
changes concentrate on early versus late deep branches."  No frozen file was edited; the
pre-registration md5 is still `8a19c101207c7016ee1affacb1ee8943`.  Every number in this
entry is exploratory, uncorrected for multiple testing, hypothesis-generating, and cannot
become a headline.

**Why G1's own nulls could not be reused.** Measured on Annelida: the observed deep tree
carries 34 deep branches, the `lambda_nullp` battery carries 22–31 over its 40 replicates,
and the structure-matched battery carries 22 over its 5.  Any null that changes the
deep-branch set cannot compare ages on it.  E4 therefore uses a null that holds the
deep-branch set fixed and permutes branch age among the branches.

### D8.1 — first run (v0) used the DISTAL end; discarded in full

The script computed `a_b = dist(node)/height`, the root-to-distal-end distance, while its
own docstring said proximal.  On a strictly ultrametric tree every terminal branch then has
`a_b = 1.000` exactly, so the statistic collapses into "terminal branch or internal branch".
Echinodermata exposes it: 6 deep branches, all six terminal, ages all 1.000,
`A_obs = −0.00000`, N1 null sd `0.00000`, N2 p exactly `1.0000`.  Caught before any number
was reported.  The v0 outputs are kept only as
`work/addon/G1/E4_timing_v0_DISCARDED_distal.{tsv,log}` and
`E4_timing_perclade_v0_DISCARDED_distal.tsv`, and their numbers are never quoted.

**Fixed** to the proximal end, `a_b = dist(parent(node))/height`, i.e. the age at which the
branch began.  Terminal and internal branches are then treated identically.

### D8.2 — second null N2 is invalid, measured

N2 permutes ages only within an origin component.  A component is a chain — the deep
branches on one descending lineage — so its proximal ages are **monotonically increasing
along the chain by construction**.  Measured: 35 of 35 non-singleton chains have strictly
increasing ages (Annelida 6/6, Chondrichthyes 4/4, Mollusca_rescue 6/6, Cnidaria 5/5,
Actinopteri_m_rescue 8/8; Echinodermata has none).  Within-chain permutation therefore
destroys a constraint the tree guarantees, the N2 null mean sits at ≈ −0.08 instead of 0,
and every clade's N2 z turns spuriously positive (Stouffer Z = +3.99 against N1's −2.64).

N2 is retained in the tables as a diagnostic and is **not used for inference**.  This is
written down so a later reader cannot pick N2 because it gives the larger number.

### D8.3 — shallow-branch control is clean, and the topology control overturns the signal

Two controls were fixed before the v2 run.

**Shallow control** (same statistic over the clade's shallow branches, N1): flat in all six
clades, z from −1.26 to +0.90, Stouffer Z = −0.46 (ALL) / −0.29 (RECUR).  So a general
"CAFE puts turnover near tips" explanation is ruled out.

**Topology control**: on an ultrametric tree a deep branch that leads **directly to a tip**
can start late for a purely topological reason.  81–98% of the deep change mass sits on
exactly that class — Annelida 42,982/52,864; Chondrichthyes 18,118/20,281; Mollusca
142,034/164,243; Cnidaria 82,878/96,542; Echinodermata 21,403/21,403; Actinopteri
153,263/155,939.  Restricting the statistic to deep branches that lead to **internal**
nodes (where a late start cannot be produced by tip topology alone):

| clade | ALL A_obs (all deep) | N1 p | ALL A_obs (internal only) | N1 p |
|---|---|---|---|---|
| Annelida | +0.01529 | 0.6466 | +0.09563 | 0.2236 |
| Chondrichthyes | +0.03682 | 0.2015 | −0.01530 | 0.8363 |
| Mollusca_rescue | −0.06121 | 0.0173 | **+0.07902** | 0.0911 |
| Cnidaria | −0.06630 | 0.0223 | **+0.00949** | 0.6066 |
| Echinodermata | −0.05455 | 0.1856 | undefined (0 internal-leading deep branches) | — |
| Actinopteri_m_rescue | −0.01616 | 0.4234 | +0.05516 | 0.5215 |

Both clades that were nominally significant flip **sign** and lose significance.  Cross-clade
Stouffer on the topology-restricted arm: Z = +1.70 (p = 0.0887, ALL) and +1.59 (p = 0.1129,
RECUR) over 5 clades, 4 of 5 late-weighted — direction reversed from the unrestricted arm,
magnitude exploratory only.

### D8.4 — outcome

**E4 does not find a timing signal.**  The unrestricted statistic measures where the deep
node labels happen to fall relative to tips, not when deep-sea turnover happened; once that
class is removed there is nothing.  Combined with D7 this leaves G1 with a clean negative on
both its confirmatory tests and no supporting exploratory result.  It is reported as such.

Files: `scripts/g1_e4_timing.py`; outputs `work/addon/G1/E4_timing.tsv`,
`E4_timing_perclade.tsv`, `E4_timing.log`.
