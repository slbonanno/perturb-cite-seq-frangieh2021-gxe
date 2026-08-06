# Analysis plan

## Framing

The Norman 2019 CRISPRa screen supports a perturbation × perturbation
interaction analysis: does the double-perturbation phenotype deviate from the
additive expectation of its two singles? The machinery for that is a log2FC
signature matrix plus a model of expected additivity.

Frangieh 2021 reuses that machinery against a different second axis. Instead of
a second gene, the modifier is the cell's environment — unstimulated, IFN-γ
stimulated, or under attack by autologous TILs. The question becomes: is the
effect of knocking out gene *X* the same in all three environments?

This is the same statistical shape and a genuinely different biological
question, which is what makes it a next step rather than a repeat.

---

## Notebook 01 — Ingest, QC, schema

**Goal:** establish what the data actually contains and which arms are powered.

1. Download both h5ads (`scripts/download_data.py`).
2. `describe_schema()` on RNA and ADT. **Do not write analysis code before
   this.** The `schema:` block in `config.yaml` is a guess; correct it now.
3. `check_schema()` — fail loudly if config and data disagree.
4. `align_modalities()` — the two files are distributed separately and are not
   guaranteed to share cells or ordering. Intersect on barcodes.
5. Write the observed ADT panel to `config/panels_observed.yaml`; reconcile
   with `panels.yaml` by hand and commit.
6. **The power table.** `group_sizes()` gives cells per
   perturbation × condition. Plot it as a heatmap. This single figure
   determines what the project can and cannot claim — arms with 12 cells will
   not support a context-dependence call no matter what the p-value says.
   Apply `min_cells_per_perturbation_per_condition`, per condition, not
   globally.
7. Standard QC: counts, genes, mito fraction, distribution by condition. Expect
   the co-culture arm to look different — those cells are being killed.
8. Guide multiplicity check: confirm low MOI.

**Deliverable:** `data/interim/frangieh_qc.h5mu` plus the power heatmap.

---

## Notebook 02 — Condition-stratified differential expression

**Goal:** a log2FC signature for every (perturbation, condition) pair.

**The central design constraint.** Every contrast is perturbation vs.
**condition-matched** control. Pooling controls across conditions would let the
enormous IFN-γ response leak into every perturbation's signature and produce a
screen where everything is a hit. Encode this in the code path, not just in a
comment.

1. HVG selection **on control cells only**. Selecting genes on the full matrix
   selects partly on the effect being measured — the same leakage failure mode
   as feature selection before cross-validation splits.
2. `assign_pseudoreplicates()` → `make_pseudobulk()`, summing raw counts.
3. `pydeseq2` per contrast, LFC shrinkage on.
4. Sanity checks that must pass before proceeding:
   - Each perturbation's own target gene should be down in its own arm
     (transcript-level effect of frameshift-induced NMD). Perturbations failing
     this are editing failures, and knowing which ones is worth a figure.
   - IFN-γ vs. control in *unperturbed* cells should recover JAK/STAT and
     antigen-presentation induction. If it doesn't, something upstream broke.
5. `signature_matrix()` → `data/processed/signatures.parquet`.

---

## Notebook 03 — Context dependence

**Goal:** the hit list.

Two independent readouts, deliberately not one:

- **Signature divergence.** `cross_condition_correlation()` — per
  perturbation, correlate its signature in IFN-γ and in co-culture against its
  signature in control. Also track magnitude ratio: an effect can be
  *amplified* in one environment without changing direction, which is a
  different phenomenon from redirection and should not be collapsed with it.
- **E-distance** (`pertpy`, permutation null). A model-free answer to "did this
  perturbation do anything at all in this condition", computed in PCA space.

`call_context_dependent()` requires **both**: a real effect somewhere, and
divergence between environments. Requiring only divergence would fill the hit
list with noise, since two null signatures are also uncorrelated. This
conjunction is the analytical core of the project.

Then:
- Cluster perturbations by their context-dependence profile (Leiden on the
  signature manifold) and ask whether the resulting modules are coherent —
  IFN-γ receptor/JAK-STAT components should group.
- Pathway enrichment on the context-dependent set (`decoupler` + MSigDB
  hallmark/reactome).
- Reuse the STRING PPI approach from the Norman project on the hit set. Cheap,
  already-written, and a good visual.

---

## Notebook 04 — The protein layer

**Goal:** the payoff. Does the transcriptional story hold at protein level?

1. CLR-normalise ADT (per cell, margin 0).
2. Recompute perturbation effects in ADT space, condition-stratified, same
   matched-control logic.
3. For each ADT feature with an RNA counterpart in `panels.yaml`, compare RNA
   and protein log2FC across all perturbations and conditions. Handle
   many-to-one mappings explicitly — an MHC-I antibody detecting HLA-A/B/C
   cannot be naively correlated against a single transcript.
4. **The four quadrants** are the interesting output:
   - concordant down — clean loss of function
   - RNA down, protein flat — buffering, protein half-life, or ADT floor
   - RNA flat, protein down — post-transcriptional regulation; mechanistically
     the most interesting quadrant
   - concordant up — induction
5. Foreground CD58, MHC-I, PD-L1. The published mechanism is that CD58 protein
   is not IFN-γ-inducible while MHC-I is, and that CD58 loss confers evasion
   without compromising MHC. Recovering that independently validates the
   pipeline and gives a concrete interview anecdote.

This is where the earlier PBMC CITE-seq work pays off: the same RNA/protein
discordance question, but now with a causal handle rather than a
cross-sectional one.

---

## Notebook 05 — Selection readout

**Goal:** a second, orthogonal phenotype from the same object.

In co-culture, cells are being killed. A perturbation's *representation* in
that arm relative to control is therefore a fitness phenotype — the classic
pooled-screen readout — and it costs almost nothing to compute.

1. `guide_enrichment()` — Fisher exact per perturbation per condition,
   BH-corrected.
2. Caveat honestly: this is not a true timecourse-controlled dropout screen,
   and differential representation can also reflect infection efficiency or
   proliferation differences unrelated to immune pressure. The control arm
   partially handles this; say so rather than overclaiming.
3. **Join the two readouts.** Scatter transcriptional effect size (E-distance
   in co-culture) against selection log-odds. The four quadrants again:
   - depleted + strong transcriptional response — sensitiser, responding and
     dying
   - enriched + strong response — candidate evasion mechanism
   - enriched + no transcriptional response — evasion without a transcriptional
     signature; the most surprising quadrant and worth dwelling on
   - neither — likely a non-functional guide

The disagreements between the two readouts are the point. Transcriptional
response and survival are not the same phenotype, and a screen that only
measures one is blind to the other.

---

## What "done" looks like

A README with five figures and three claims:

1. A ranked, defensible list of context-dependent perturbations, with the
   IFN-γ/JAK-STAT and antigen-presentation modules recovered as positive
   controls.
2. An RNA-vs-protein quadrant analysis that independently rediscovers the CD58
   discordance.
3. A comparison of transcriptional and selection readouts showing they
   disagree, with a worked example.

If those three exist and the limitations section is honest, the project is
finished. Resist adding a prediction model.
