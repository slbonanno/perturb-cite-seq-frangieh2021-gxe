# perturb-cite-seq-frangieh2021-gxe

Reanalysis of **Frangieh et al. 2021** Perturb-CITE-seq asking a
gene-by-environment question: *which CRISPR knockouts change what they do
depending on immune context, and do those changes show up at the protein
level?*

> Frangieh CJ, Melms JC, Thakore PI, et al. Multimodal pooled Perturb-CITE-seq
> screens in patient models define mechanisms of cancer immune evasion.
> *Nature Genetics* 53:332–341 (2021). doi:10.1038/s41588-021-00779-1

---

## The question

Classical Perturb-seq analysis asks about perturbation × perturbation
interaction — epistasis. This dataset supports a different second axis.
~218,000 patient-derived melanoma cells carry knockouts of 248 genes, and each
cell sits in one of three environments: unstimulated control, IFN-γ
stimulation, or co-culture with autologous tumour-infiltrating lymphocytes.

So the interaction of interest is perturbation × **environment**. Does knocking
out gene *X* mean the same thing to a cell sitting alone in a dish as it does
to a cell under T-cell attack? For some genes the answer is obviously no, and
those genes are where immune-evasion mechanism lives.

Because the assay is CITE-seq, every transcriptional answer has an independent
protein-level check across a 20-plex surface panel — including the paper's
headline finding, CD58, whose protein is *not* IFN-γ-induced and whose loss
spares MHC. That makes it a built-in positive control for RNA/protein
discordance rather than a curiosity.

## Analysis plan

| Notebook | Question |
|---|---|
| `01_ingest_qc_schema.ipynb` | What is actually in these files, and which perturbation × condition arms have enough cells to say anything? |
| `02_condition_stratified_de.ipynb` | What does each perturbation do, *within* each environment, against condition-matched controls? |
| `03_context_dependence.ipynb` | Which perturbations behave differently across environments — by signature divergence and by E-distance? |
| `04_protein_layer.ipynb` | Do the transcriptional calls survive at the protein level, and where do RNA and ADT disagree? |
| `05_selection_readout.ipynb` | In co-culture, cells are under T-cell selection. Which perturbations are enriched or depleted, and does that agree with the expression readout? |

Full reasoning, including the design traps, is in
[`docs/analysis_plan.md`](docs/analysis_plan.md).

## Deliberately out of scope

- **Perturbation-response prediction** (GEARS, CPA, scGen). A real and
  interesting problem, and a *different project*. Folding it in here would
  turn one clean question into two half-answered ones.
- Cross-dataset integration with Norman 2019.
- Causal / GRN inference from interventional data.

## Setup

```bash
conda env create -f environment.yml
conda activate frangieh_gxe
python scripts/download_data.py        # ~3-5 GB into data/raw/
jupyter lab
```

Data comes from the [scPerturb](https://scperturb.org) harmonised release
(Zenodo record `10044268`), which applies uniform QC and feature harmonisation
across 44 perturbation datasets. Original deposit is Broad Single Cell Portal
`SCP1064`.

## Configuration

Everything tunable lives in YAML; nothing is hard-coded in notebooks.

- `config/config.yaml` — paths, seed, data URLs, `.obs` schema, QC thresholds,
  DE parameters, E-distance settings, plotting.
- `config/panels.yaml` — ADT↔RNA feature mapping and pathway gene sets.

Notebooks read config via `src/config.py`, which resolves the repo root by
walking up to `config/config.yaml`, so they run from any working directory.

> **Before trusting anything:** the `schema:` block in `config.yaml` is an
> *assumption* about scPerturb's column names. `nb01` runs
> `describe_schema()` and `check_schema()` to verify it against the real files.
> Correct the YAML there, then everything downstream inherits the fix.

## Known limitations

State these plainly rather than burying them:

1. **No biological replicates.** DE uses random pseudo-replicate splits of a
   single sample, which gives the negative-binomial model a variance term but
   systematically understates biological variance. DE p-values are
   anti-conservative — treat the ranking as a ranking. The permutation-tested
   E-distance in `nb03` is the honest effect-size check.
2. **Knockout is incomplete.** CRISPR-Cas9 produces a mixture of frameshift,
   in-frame and unedited alleles. A perturbation with no phenotype may be a
   real negative or an editing failure; guide-level consistency in `nb01` is a
   partial check, not a solution.
3. **Selection confounds the co-culture arm.** Cells surviving T-cell killing
   are a non-random subset. This is a nuisance for the expression analysis and
   the *signal* for `nb05` — the same fact wearing two hats.
4. **One patient-derived line.** Nothing here establishes generality.

## Repository layout

```
config/     YAML configuration (the control surface)
src/        importable helpers: config, io, pseudobulk, stats, plotting
notebooks/  01-05, run in order
scripts/    download_data.py
data/       raw / interim / processed  (gitignored)
results/    figures / tables           (gitignored)
docs/       analysis_plan.md
```

## License

MIT (code). The underlying data is redistributed by scPerturb under CC-BY;
cite Frangieh et al. 2021 and Peidli et al. 2024 for any reuse.
