# Analysis code for: cross-endpoint genetic localization of hearing loss in mouse cochlear single-cell atlases

Zhengyi Tang¹, Ting Liao¹, Chengmeng Wei¹, Yibei Jiao¹, Ge Wang¹

¹ Department of Otorhinolaryngology and Head & Neck Surgery, 923rd Hospital of PLA, Nanning 530021, China

Corresponding author: Zhengyi Tang — doctang@hotmail.com (co-corresponding author: Ge Wang)

This repository contains the code used to produce every number, table and figure
reported in the accompanying manuscript submitted to *BMC Genomics*. It does not
contain the input data, which are all publicly available (see [Input data](#input-data)).

---

## 1. What is here

| Directory | Contents |
|---|---|
| `scripts/` | Data preparation, QC, single-cell atlas construction, scDRS scoring, cross-endpoint comparison, and sensitivity analyses |
| `figure_scripts/` | Code that renders Figures 1–4 from the analysis outputs. File names match the published figure numbers |
| `environment/` | Environment snapshot (`ENVIRONMENT.md`), the anndata/SIGILL workaround note, and `requirements.txt` |

Files whose names begin with `_` are ad-hoc probes and are deliberately excluded.

## 2. Environment

Python 3.13.14, on an isolated virtual environment (`--` paths in `environment/ENVIRONMENT.md`
are machine-specific and must be adapted). Key pins:

```
anndata==0.13.4   h5py==3.16.0      numpy==2.5.3     pandas==2.3.3
scanpy==1.12.4    scdrs==1.0.2      scipy==1.18.1    scikit-learn==1.9.1
matplotlib==3.11.2  seaborn==0.13.2  statsmodels==0.15.0  numba==0.67.0
```

See `environment/ENVIRONMENT.md` for the full snapshot and
`environment/ENVIRONMENT_FIX_anndata.md` for the anndata crash workaround that the
environment requires.

## 3. Input data

All inputs are public; no individual-level human data are used.

| Input | Source |
|---|---|
| FinnGen R13 GWAS summary statistics (three endpoints: `H8_HL_SEN_NAS`, `H8_HL_CON_NAS`, `H8_NOISEINNER`) | `https://storage.googleapis.com/finngen-public-data-r13/summary_stats/` |
| Endpoint definitions and case/control counts | `finngen_R13_pheno_n.tsv` (same bucket) |
| Miao 2024 mouse whole-cochlea snRNA-seq | GEO **GSE281324** (GSM8618786–GSM8618789) |
| Sun 2023 mouse cochlear atlas | NGDC GSA **CRA004814** (also via gEAR) |
| Milon 2021 mouse lateral wall / spiral ganglion | GEO **GSE168041** |
| Taukulis 2021 mouse stria vascularis (cisplatin) | GEO **GSE165662** |
| Korrapati 2019 mouse stria vascularis (naive) | GEO **GSE136196** |
| Mouse–human orthology | MGI `HOM_MouseHumanSequence.rpt`, <https://informatics.jax.org> |
| LD reference panel | 1000 Genomes Phase 3 EUR (QCed PLINK mirror) |

## 4. Reproduction order

The pipeline is linear. Each step writes an output that the next step reads.

**Step 1 — prepare FUMA/MAGMA inputs (GRCh38, three columns)**

```bash
python scripts/make_fuma_input_grch38.py      # writes *_grch38_3col.tsv.gz
python scripts/qc_fuma_input.py               # full, non-sampled QC of the prepared files
```

**Step 2 — gene-level test (MAGMA) via the FUMA web platform**

The gene-based test was executed on the **FUMA web platform**, not locally. The exact
job settings are recorded in `fuma_output/<endpoint>/<endpoint>.params.config`
(retained with the results; `magma_window = 0`, `leadP = 5e-8`, GRCh37 gene
coordinates, Ensembl v102, MHC excluded by annotation). Returns are
`magma.genes.out` per endpoint. Local alternative scaffolding, retained for
transparency: `scripts/prepare_magma_input.py`, `scripts/export_fuma_input.py`,
`scripts/make_fuma_input_rsid.py`, `scripts/qc_significant_snps.py`.

**Step 3 — bridge human gene-level statistics to mouse symbols**

```bash
python scripts/build_ortholog_map.py          # MGI HOM_MouseHumanSequence -> strict 1:1 subset
python scripts/prep_gene_z_from_fuma.py       # human Z -> mouse symbol matrix (GENE leading column)
```

**Step 4 — build the three integrated mouse cochlear atlases**

```bash
python scripts/qc_scrna.py
python scripts/build_scrna_atlas.py           # sv_atlas, whole_cochlea_atlas, sgn_atlas
python scripts/label_transfer.py              # cell-type label harmonisation
```

**Step 5 — score cells and compare endpoints**

```bash
bash scripts/run_main_scdrs.sh                # scDRS per-cell scores, n_ctrl = 1000
python scripts/compare_endpoints.py           # cross-endpoint D = q95(A) - q95(B), MC P, BH-FDR
python scripts/jackknife_endpoints.py         # block jackknife
```

**Step 6 — sensitivity analyses**

```bash
bash scripts/run_sensitivity_a.sh
bash scripts/run_sensitivity_gs2000.sh
python scripts/sensitivity_reaggregate.py     # cell-level re-aggregation, leave-one-dataset
bash scripts/run_gs2000_noise_cpu.sh          # (attempted; noise endpoint at 2,000 genes did not complete)
```

**Step 7 — render figures**

```bash
python figure_scripts/fig1_gene_level.py          # -> Figure 1
python figure_scripts/fig2_cross_endpoint.py      # -> Figure 2
python figure_scripts/fig3_celltype_enrichment.py # -> Figure 3
python figure_scripts/fig4_sensitivity.py         # -> Figure 4
```

## 5. Conventions used throughout the manuscript

- **Cross-endpoint statistic** `D = q95(A) − q95(B)`, where `A` is the endpoint
  listed first in the comparison. A negative `D` means endpoint `B` has the higher
  enrichment. Values are taken as they appear in the `diff` column of
  `compare_pairwise_tests.csv`; **no sign conversion is applied anywhere**.
- Monte-Carlo P values come from full enumeration of control-gene-set 95th-percentile
  differences (10⁶ pairs), two-sided.
- BH-FDR is applied within each atlas over all pairwise tests of that atlas.
- `q′` is a prespecified correction for the number of candidate compartments:
  `q′ = 2q` (correction A, taken as primary) and `q′ = 3q` (correction B).

## 6. Known limitations of the released code

- The 2,000-gene configuration for the noise-induced endpoint could not be
  completed (eleven consecutive crashes); only 28 of 84 comparisons exist at that
  gene-set size. This is reported in the manuscript.
- Paths inside the shell wrappers are machine-specific and need editing before a
  clean-machine run.
- The gene-level step depends on an external web service (FUMA) whose availability
  and version are outside this repository.

## 7. Licence

MIT — see `LICENSE`.
