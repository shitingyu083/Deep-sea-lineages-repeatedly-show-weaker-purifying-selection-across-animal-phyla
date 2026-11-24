# DEVIATIONS_g2.md — append-only

Deviations from `PREREG_G2_protein_physchem.md`
(md5 `b86f0b87cd1f17b97c62f241844d91b2`).

PREREG_G2 is never edited. Each entry records what was seen and what was then done.

---

## E1 — 2026-11-03 — the pre-registered shallow-partition null is INFEASIBLE for 3 of 6 clades

**Seen (feasibility probe, run before any physicochemical effect was computed):**

| clade | K | deep tips | shallow tips | component sizes | shallow pool ≥ deep tips? |
|---|---|---|---|---|---|
| Annelida | 8 | 21 | 17 | 2,5,4,3,1,3,1,2 | **no** |
| Chondrichthyes | 6 | 11 | 6 | 3,2,1,1,2,2 | **no** |
| Mollusca_rescue | 18 | 38 | 41 | 1,2,4,1,3,1,9,2,6,1×9 | yes |
| Cnidaria | 12 | 21 | 25 | 1,2,2,1,2,1,5,1,3,1,1,1 | yes |
| Echinodermata | 6 | 6 | 16 | 1,1,1,1,1,1 | yes |
| Actinopteri_m_rescue | 36 | 46 | 29 | sum 46 | **no** |

PREREG_G2 §4.1 drew K groups of the observed component sizes **from the shallow pool**. In
Annelida, Chondrichthyes and Actinopteri_m_rescue the shallow pool is smaller than the sum of
those sizes, so that partition cannot be formed at all — not "rarely", but never. The
preregistered null would have produced zero replicates for half the clades.

**Also a bias in the original design, found at the same moment:** every Δ_k shares the same
denominator x̄_shallow. If the shallow reference is a noisy draw from a wide distribution,
all K values shift together and produce spurious concordance. Re-partitioning only the
shallow pool cannot detect this, because it never re-randomises the reference.

**Done — amended null (fixed before any effect was computed): pooled partition null.**

Pool all species of the clade that carry data for this family (deep + shallow). For each
replicate, draw a uniform random permutation of the pooled species and assign, in order:

- K pseudo-origin groups whose sizes equal the observed deep component sizes, then
- one pseudo-reference group whose size equals the observed shallow species count.

Recompute Δ'_k, e', c', R' exactly as in §2. The pooled set is exactly the union of the
observed groups, so the assignment is always complete and feasible.

This preserves: K, every component size, the reference size, the family set, the trait values,
the metric distribution and the number of species. It destroys only **which species are the
deep ones** — the hypothesis under test — and it now also randomises the reference, removing
the shared-denominator artifact above.

Nothing else in PREREG_G2 changes. The §4.2 secondary label-permutation battery is still
reported beside it and still cannot decide anything (DEVIATIONS_g1.md D1).

---
## E2 — 2026-11-03 — the §1 coverage table was measured on the wrong trait table; no prereg change follows

**Seen.** Before any effect was computed, the first attempt to apply §1 used
`work/physchem/gene_metrics.tsv.gz`, the table named in §1's "Trait source" line. It
returned 688–881 passing families per clade instead of the 2,843–6,479 in §1's coverage
table, and only 286 families in Echinodermata. `g2_join_diag.py` gave the reason: that
table holds compleasm **BUSCO-Single** hits only — ~620 genes per species out of 30,000+
— covering 2.0–3.2 % of the genes in each clade's `N0.tsv`. It was built
(`30_physchem.py`) for *species-level* composition, where restricting to a single-copy
benchmark keeps annotation route from driving the comparison. G2 compares genes
**within** a homologous family, where that constraint does not apply, and it caps the
analysis at ~700 families instead of the thousands §1 assumed.

**Done — the trait table was rebuilt, and §1 itself was not touched.**

`34_physchem_full.py` recomputes the same metric set (formulae copied verbatim from
`30_physchem.py`) for **every** gene in the N0 union of the six clades: 282 species,
6,640,854 gene models, one gzipped TSV per species under `work/physchem_full/`.
Re-applying §1's rule 1+2 to those tables (`g2_incl_count2.py`) gives

| clade | K | families in N0 | passing | PREREG §1 |
|---|---|---|---|---|
| Annelida | 8 | 42,300 | **4,148** | 4,148 |
| Chondrichthyes | 6 | 21,193 | **6,477** | 6,479 |
| Echinodermata | 6 | 24,318 | **6,029** | 6,029 |
| Cnidaria | 12 | 54,323 | **3,927** | 3,927 |
| Mollusca_rescue | 18 | 82,651 | **2,843** | 2,843 |
| Actinopteri_m_rescue | 36 | 30,814 | **3,781** | 3,781 |

Five of six match exactly and Chondrichthyes differs by 2 of 6,479 — round-off, not a
coverage gap. So §1's table was correct as written; it describes full-proteome
coverage, and the "families in N0" column differs slightly from the current `N0.tsv`
counts because that column was read at freeze time. **§1 is not overstated. No E2 was
needed for the reason the earlier working note expected, and none is claimed here.**
The entry exists only to record that the first execution used a table the prereg's own
trait-source line did not intend for this analysis, and that nothing downstream
depends on it.

**Also fixed, and worth recording separately from E1.** `g2_run.py` as first executed
on the server still carried a no-op permutation: `G[perm, j*ng + lab[perm]] = 1.0`
assigns row `perm[s]` to its own label, so every replicate reproduced the observed
assignment and the null had exactly zero variance (R = null mean exactly, p = 1.000,
z = 0.00, in all 32 metrics). Corrected to `G[np.arange(nsp), j*ng + lab[perm]] = 1.0`.
After the fix the Echinodermata null has sd ≈ 1,600 on R ≈ 18,000 and range
16,351–26,712 at B = 200; all reported G2 numbers come from the corrected script
(md5 `c820bb9b3a60bb1e405b5f0925b09007`).

---
## E3 — 2026-11-03 — the pooled partition null scored dataless families as ties; null rebuilt (v3)

**Seen.** PREREG_G2 §1 keeps a family only if every deep component and the shallow
reference hold a species carrying that metric. The observed family set is therefore
built under an inclusion rule. The permuted null in `g2_run.py` was not: it re-used the
observed family set and re-partitioned species, so a group of size 1–2 could receive
only species with no data for that family. Its group mean is then NaN. Because
`(NaN > 0)` and `(NaN < 0)` are both False, the family is scored as `e == c` — counted
as a **tie**, silently and at full weight.

Measured directly (`g2_tie_diag.py`, Annelida × gravy, B = 4,000):

| quantity | value |
|---|---|
| NaN families per replicate | **427.7** (mean) |
| replicates containing ≥1 NaN family | **4,000 / 4,000** |
| observed tie rate | 0.097 |
| null tie rate | **0.179** |
| within-null R-vs-tie slope | −20,223 |
| within-null R-vs-tie corr | **−0.824** |

The artefact inflates the null tie rate and deflates null R in lockstep, so the v2 raw
`p = (1 + #{R_null ≥ R_obs}) / (1 + B)` was **anti-conservative**. Every v2 raw p is
withdrawn and must not be quoted.

**Separately void.** `g2_nullcheck2.py` reported "null rejection rate 0.000 in all six
clades". That script ran under the no-op permutation recorded at the end of E2, so it
measured nothing; the number carries no information and is not cited.

**Done — the null was rebuilt as v3** (`g2_run_v3.py`, md5
`a62e8914c6fc3dfcf6e27844b12c3a88`), with two changes and nothing else:

1. §1's inclusion rule is applied **inside each replicate**:
   `valid = (cnt[:,:,:K] > 0).all(axis=2) & (cnt[:,:,K] >= 2)`. Families failing it are
   **dropped as undefined**, not scored as ties. The permutation stays one draw per
   replicate applied to all families, because that is what preserves the between-family
   correlation of the observed design; permuting per family would collapse the null
   variance toward zero and manufacture significance.
2. The comparison is on the normalised `rate = R / (n_valid · (K−1))`, so numerator and
   denominator carry the identical inclusion rule on both sides.

**A consequence of the fix, recorded because it changes how the numbers are read.** Under
E1 the reference is re-randomised, so a real effect *lowers* the observed tie rate. The
tie gap therefore co-varies with the effect by construction: corr(tie_gap, z) = **+0.876**
confirmatory, +0.706 exploratory, **+0.702** control. A tie-matched conditioning would
remove part of the effect along with the nuisance. So the deciding statistic is `p` on
normalised rate. A ridge diagnostic (null rate regressed on null tie, observed residual
read off the same line) is reported as `p_ridge` — a **calibration readout, not a
deciding test**, and `rate@obs_tie` is reported beside it so the reader can see the price.

**Calibration check, which v2 failed.** Annelida `log10_len` — the §3.2 negative control —
has observed tie 0.152 against null tie 0.156 (gap 0.004) and ridge residual **−0.0012**,
`p_ridge` 0.5771. On a metric with no effect the observed tie rate lands on the null's, as
it must. Annelida gravy on the same readout: residual +0.0137.

**Corrected outcome, B = 2,000, 6 clades × 32 metrics × 2 aggregations, 12 clean jobs.**
Control: **0/6** at p < 0.05, z ∈ [−1.82, +0.26]. Confirmatory: **1/18**. Exploratory:
**4/168**, all four Annelida and all in the same direction (f_E 0.0195 z +2.63;
charged_DEKR 0.0205 z +2.65; aliphatic 0.0210 z +2.37; f_K 0.0360 z +2.19).

The one confirmatory cell is Annelida × gravy: rate 0.7290 against null 0.6399
(sd 0.0247), p = 0.0045, z = +3.61; Holm within clade over the three confirmatory metrics
= 0.0135; median aggregation reproduces it (rate 0.7284, p 0.0045, z +3.59). Every other
confirmatory cell has p ≥ 0.2029. **G2-C1 is met once, in one clade, for one metric** — and
E4 below then removes it.

---

## E4 — 2026-11-03 — the Annelida × gravy hit is a proteome-wide compositional shift, not family-level convergence

**Seen, before anything was written up.** The family-level table for Annelida × gravy
(`G2v3_fam_Annelida_gravy_mean.tsv`) is 86.5% positive: of 4,148 families, 3,588 have
mean Δ > 0, overall mean +0.0594, median +0.0546, and the (e, c) distribution is (7,1)
n=1,091, (6,2) n=983, (5,3) n=678, (8,0) n=530 with a thin contraction tail. The
registered statistic is direction-agnostic by design and cannot distinguish bidirectional
family-level convergence from a uniform offset applied to every family.

That distinction is decisive, because E1's permutation only shuffles *which species are
deep*. A clade-wide compositional difference between the deep and shallow proteomes
reproduces the same signal in every family and is not reproduced by the null.

**Test run** — within-species centring (`g2_ann_gravy_diag.py`): subtract each species'
own proteome-wide mean for the metric from all of its genes, then recompute rate and
null. If the recurrence is the compositional offset it disappears; if a subset of
families shifted beyond their species baseline it survives. Same seed, same null, B=1,000.

| | rate | R | tie | frac(meanΔ>0) | p | z |
|---|---|---|---|---|---|---|
| raw | 0.7290 | 21,166 | 0.097 | 0.865 | **0.0030** | **+3.84** |
| within-species centred | 0.6577 | 19,098 | 0.164 | 0.606 | **0.1618** | **+0.99** |
| null (raw) | 0.6386 ± 0.0235 | | | | | |
| null (centred) | 0.6379 ± 0.0200 | | | | | |

The whole effect is the offset. Deep annelid proteomes average gravy **−0.3357** (n=21)
against shallow **−0.3778** (n=17), a difference of **+0.0421** — the same sign and scale as
the family-level mean Δ of +0.0594. Centring drops the positive fraction from 0.865 to
0.606, close to the 0.597 expected under the null, while the dispersion of Δ across
families is unchanged (sd 0.0607 → 0.0606). Nothing concentrates in particular families;
every family moves by about the same amount.

**Not one origin carrying it, either.** Per-origin raw mean Δ is positive in all eight
origins (+0.015 to +0.116), and leave-one-origin-out leaves the rate at 0.7155–0.7544
against a null of 0.6386 — dropping any single origin changes nothing. That is the
signature of a clade-wide shift, not of independent origins converging.

**Done — G2-C1 does not stand.** The Annelida × gravy result is recorded as a
**proteome-wide amino-acid compositional difference between deep and shallow annelids**,
not as recurrent protein-property convergence in homologous families, and it is not
reported as a convergence claim. §4.3's reliability rule would in any case have flagged it:
the observed tie rate moves 0.097 → 0.164 on centring.

The centring test is a diagnostic, not a confirmatory test, and it was not pre-registered;
it is reported in full with its numbers so the negative result can be checked as easily as
a positive one would have been. It is applied to the other five clades' confirmatory cells
(E3: all p ≥ 0.2029) only for the record — those are already null.

**Net G2 outcome.** Neither the rate level (R6-5), nor the family-size level, nor the
protein-property level shows recurrent genomic convergence across independent deep-sea
origins. The negative control passes cleanly in all six clades, the 4/168 exploratory
cells are all Annelida and all explainable by the same compositional offset, and the single
confirmatory hit is an artefact of that offset.
