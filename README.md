# perturb-cite-seq-frangieh2021-gxe

Reanalysis of **Frangieh et al. 2021** Perturb-CITE-seq on a patient-derived melanoma cell line x stimuli.

> Frangieh CJ, Melms JC, Thakore PI, et al. Multimodal pooled Perturb-CITE-seq
> screens in patient models define mechanisms of cancer immune evasion.
> *Nature Genetics* 53:332–341 (2021). doi:10.1038/s41588-021-00779-1

#### Experimental Design:
- **248 CRISPR-KO targets** immune-checkpoint-resistance / ICR program
   - **~744 targeting guides** (~3 guides/gene)
   - **74 control guides** (37 non-targeting + 37 intergenic)
- **virus, MOI ~0.1** each infected cell receives one guide
- **puromycin selection, transduced cells**
- **3 environmental arms:** Control, IFNg-treated, autologous TIL co-culture
- **readout:** scRNA-seq, 20 ADT panel

#### Processed Dataset:
- **~218k cells~** immune-checkpoint-resistance / ICR program

#### Data Exploration:
The interaction of interest is KO target × environment.
How do expression profiles (and ADT signal) differ, from cells targeted by same guide
but in each context: control, IFNg, or TILs (leukocyte attack)

Note that leukocyte attack on cells likely cause big changes in gene expression that will compete with
and interact with the CRISPR-KO

This project is in python using scanpy.  Each AnnData object handles one sparse matrix (and transformations of it).
MuData will hold both the RNA and ADT layers.

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

### Currently out of scope

- **Perturbation-response prediction** (GEARS, CPA, scGen).
- Cross-dataset integration with Norman 2019.
- Causal / GRN inference from interventional data.

### Setup

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

### Configuration

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

### Known limitations

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

### Repository layout

```
config/     YAML configuration (the control surface)
src/        importable helpers: config, io, pseudobulk, stats, plotting
notebooks/  01-05, run in order
scripts/    download_data.py
data/       raw / interim / processed  (gitignored)
results/    figures / tables           (gitignored)
docs/       analysis_plan.md
```

### License

MIT (code). The underlying data is redistributed by scPerturb under CC-BY;
cite Frangieh et al. 2021 and Peidli et al. 2024 for any reuse.


## Part I: data exploration and QC
schema of .h5ad processed data was explored, naming of columns etc adjusted in config.yaml and panels.yaml
A power table was constructed to look at number of pseudoreplicates (cells) per {perturbation1 x perturbation2} condition

### Figure 1 — Screen design and cell filtering

![Figure 1](results/figures/01_qc_overview.png)

**(A)** MOI distribution: number of guides detected per cell across all 218k cells.
Although the screen was designed for low multiplicity, only 58% of cells carry exactly one guide;
31% carry two or more (up to 19), and 11% have no guide captured (which may be seq depth issue or no guide delivery - they survived puro).
**(B)** For cells receiving more than one guide: no {perturbation1 x perturbation2} has enough cells to proceed with.
Gene1 x gene2 comparisons are not feasible in this experimental design.
**(C)** dataset filtered for cells with MOI=1 (1 guide received). Sufficient power is defined as ~30 cells within
{perturbation1 x perturbation2} bin.
**(D)** cell counts for perturbations where n<=30 cells in at least one {perturbation2}.
**(E)** sgRNA representation within each target {perturbation1}: for a given target, are there multiple guides represented in the data?

Panels A and B describe all cells (218k), since they motivate the filter; C–E describe
the 126k single-guide cells retained for downstream analysis.

#### Figure 1 conclusions:

Wet lab experimental design caveat: 42% of cells captured are filtered out (MOI not equal to 1; insufficient cell number to do gene x gene comparisons)

sgRNA representation within-target is a small issue, over all.  Though a few targets are dominated by just one sgRNA.

### Figure 2 — Cell-level quality control

![Figure 2](results/figures/01_qc_violins_moi_filter.png)

UMI count, gene count, and mitochondrial and ribosomal fractions, split by
condition, before **(A)** and after **(B)** restricting to single-guide cells.
White points mark medians; group sizes are given beneath each violin. UMI and
gene counts are on log scales.

QC cutoffs (i.e. pct.mito < ~20) determined by scPerturb "harmonised" release.
These filters can be treated as quality-neutral, since they look the same before and after filtering for single-guide cells.

#### Figure 2 conclusions:
No further QC filtering on single-guide-cells dataset is warranted at thsi point.
We **do** expect to see differences in "QC metrics" like pct.mito when we get to comparing TILs to control - the cells are under attach and dying.
This will need to be handled with biology in mind, conducting analysis downstream.





