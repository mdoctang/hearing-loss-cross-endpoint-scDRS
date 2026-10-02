# Figure scripts — script file → output file → final figure number

> Every script reads the analysis products listed below read-only and writes
> one PNG (400 dpi) plus one PDF (vector, fonts embedded) into
> `C:/Users/docta/WorkBuddy/纯生信/figures/`.
> Run with the project venv and **always with `-u`**:
> `C:/Users/docta/.workbuddy/binaries/python/envs/default/Scripts/python.exe -u <script>`

## Mapping (script file name → output file → final figure number → content)

| script file | output file(s) | final figure | content |
|---|---|---|---|
| `fig1_gene_level.py` | `Fig1_gene_level.png/.pdf` | **Figure 1** | Gene-level MAGMA results: three Manhattan panels (GRCh37) + gene z-score distributions |
| `fig3_cross_endpoint.py` | `Fig2_cross_endpoint.png/.pdf` | **Figure 2** | Cross-endpoint differences D = q95(A) − q95(B): top-12 forest + all-84 lollipop |
| `fig2_celltype_enrichment.py` | `Fig3_celltype_enrichment.png/.pdf` | **Figure 3** | Cell-type-level scDRS enrichment: q95 dot plot + CON-vs-SEN scatter |
| `fig4_sensitivity.py` | `Fig4_sensitivity.png/.pdf` | **Figure 4** | Sensitivity of the IPhC_IBC CON-vs-SEN difference (cell-level re-aggregation / gene-set size / gene-set construction) |

⚠ **Figure numbering follows the order of first mention in the main text**
(R1 → Fig 1, R5 → Fig 2, R7 → Fig 3, R11 → Fig 4), so the numbers in the
output file names are NOT alphabetical with the script file names. Two script
files still carry an older number in their own name:

- `fig3_cross_endpoint.py` → produces **Figure 2** (`Fig2_cross_endpoint.*`)
- `fig2_celltype_enrichment.py` → produces **Figure 3** (`Fig3_celltype_enrichment.*`)

Each script file's docstring states its final figure number. Caption text for
all four figures lives in `figures/figure_legends.md`.

## Shared modules and diagnostics

| file | role |
|---|---|
| `figstyle.py` | Shared rcParams (light theme, DejaVu Sans, Okabe-Ito palette), `save()` helper writing PNG+PDF at 400 dpi. Sets `mathtext.default: regular`; text that must survive PDF text extraction is written as plain Unicode (e.g. U+2032 PRIME) rather than mathtext. |
| `explore_data.py` | Read-only exploration used to confirm every number printed in the figures against the products. |
| `_legendprobe.py` | Diagnostic that rendered the prior-correction legend string to show why the mathtext prime became `q0` in the PDF text layer. Retained as the evidence trail for MI-22; not needed to rebuild figures. |

## BMC layout constraints applied to all four figures

- Width ≈ 170 mm (BMC full-page width) for every figure, so nominal font size
  equals effective font size; all text ≥ 6.8 pt.
- Height ≤ 225 mm (BMC figure-and-legend cap).
- Minimum line width 0.4 pt (> the 0.25 pt requirement); TrueType fonts
  embedded (subset type `XXXXXX+DejaVuSans`).
