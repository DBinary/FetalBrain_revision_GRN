# T577 GRN Reproducibility Package

This package contains the code, direct inputs, and figure outputs for the current six-panel analysis, together with the full and five half-count GRN results. The original working files were left intact. Large expression matrices, full result H5AD files, and large reference databases are not included.

The current figures use five TF-expression groups. The two expression panels use cividis; the Top-k panel uses a purple gradient. Older four-group, fixed full-signature AUCell, and statistical-annotation figures are outside this package version.

## Contents

| Path | Contents |
|---|---|
| `T577_SpaGRN_visualization.ipynb` | Code and saved outputs for the six final figures. |
| `figures/` | Six standalone PNG figures and six vector PDF figures. |
| `data/figure_inputs/` | Confirmed plotting inputs and analysis baselines; these are not overwritten by recomputation. |
| `data/analysis_recomputed/` | Outputs from lightweight recomputation, kept separate for comparison with the baselines. |
| `data/grn/full/`, `data/grn/rep01/`–`rep05/` | Native regulon JSON and motif results, plus weighted final-edge tables for all six runs. |
| `notebooks/analysis/` | Edge export, network and Top-k analysis, native-score summaries, expression features, and package validation. |
| `notebooks/downsampling/`, `notebooks/grn/` | Downsampling and six inference notebooks with their saved historical outputs. |
| `environment/` | Python/package version snapshots, the Arial font file, and a SpaGRN Python-source snapshot. Version snapshots do not establish that these environments were reinstalled or validated on another machine. |
| `reference/` | TF list, macaque reference description, and related build material. |
| `config.json` | Runtime paths and parameters for analysis and training. |
| `CAPTION.md`, `CAPTION.tex` | Figure caption in Markdown and LaTeX; replace “Fig. X” with the final figure number. |

## Configuration and running

Run notebooks from the package root. `config.json` separates `analysis` settings, including paths for the full run and replicate runs, from `training` settings, including input matrices, reference resources, downsampling parameters, and output paths. Update these values for the local machine before running workflows that use external data. The package includes small plotting inputs and final networks, but it does not include the large source matrices or reference databases needed to recreate those inputs and networks.

To run the lightweight plotting and package checks, install the notebook dependencies in a Jupyter environment and use the corresponding local kernel. For example:

```bash
git clone https://github.com/DBinary/FetalBrain_revision_GRN.git
cd FetalBrain_revision_GRN
python - <<'PY'
from pathlib import Path
import nbformat
from nbclient import NotebookClient

root = Path.cwd()
kernel = "hulei_spatrack"  # Replace with a kernel installed on this machine.
for name in [
    "notebooks/analysis/02_network_topk.ipynb",
    "notebooks/analysis/05_TF_expression_quintiles.ipynb",
    "T577_SpaGRN_visualization.ipynb",
    "notebooks/analysis/06_validate_package.ipynb",
]:
    path = root / name
    notebook = nbformat.read(path, as_version=4)
    NotebookClient(notebook, kernel_name=kernel, timeout=600).execute(cwd=str(root))
    nbformat.write(notebook, path)
PY
```

The root plotting notebook reads packaged tables and the bundled font; it does not require the external expression matrices. The plots use 16-point type, 4.5 × 4.5 inch panels, 600 dpi PNGs, and PDFs with embedded TrueType fonts. The font is explicitly loaded from the package.

Lightweight recomputation writes to `data/analysis_recomputed/` so its tables can be compared with `data/figure_inputs/`. `notebooks/analysis/01_export_final_edges.ipynb` can regenerate weighted edge tables from motif results and requires a compatible PySCENIC/SpaGRN environment; the edge tables are already included.

## Figure definitions and interpretation

| Figure | Packaged input and plotted quantity |
|---|---|
| Relative edge count | Total or shared edges in each run divided by the full network's 68,492 edges. |
| Relative regulon count | Total or shared regulons in each run divided by the full network's 260 regulons. |
| Shared-regulon Pearson | Native AUCell Pearson scores from `native_shared_regulon_pearson.tsv`. |
| TF-expression module recovery | Modules retained per expression group and run, divided by 52 TFs. |
| TF-expression target recovery | Retained full-network TF–target edges divided by full-network edges in that group. |
| Top-k Jaccard | Median regulon-level Jaccard overlap of Top-k target sets within each repeat's fixed eligible regulon subset. |

Pearson values compare native scores for each run's own inferred regulons, not scores from fixed full-network signatures. Each correlation pairs all 117,867 bins. The five runs have 183, 181, 181, 178, and 175 shared regulons, respectively, for 898 correlations. Of the 1,300 full-by-run regulon pairs, 402 modules are absent from a run; these have zero recovery and undefined Pearson values. Negative correlations are retained.

Expression groups are defined from mean full-run CP10K, including zero-expression bins. Each quintile contains 52 TFs and is shown from Q5 to Q1. Target recovery counts TF–target edges, so the same target linked to different TFs contributes once per TF; this differs from averaging per-TF recovery equally.

Top-k targets are sorted by native weight descending, then target name ascending; weights are not converted to absolute values. A shared regulon is eligible only if it has at least 100 final targets in both the full and downsampled networks. The five repeats contain 67, 60, 61, 64, and 62 eligible regulons. Each repeat uses the same eligible regulon list for k = 20, 40, 60, 80, and 100; the lists can differ between repeats. This conditional comparison does not by itself represent stability across all modules.

## GRN results and workflow scope

For each run, `regulons.json` and `motifs.csv` are copied native outputs. `final_edges.tsv.gz` is derived from motif `gene2weight` values and restricted to members in the final JSON regulons. Its columns are `run, regulon, TF, target, weight`. Self-edges and native negative weights are retained. This is the final network, not the original spatial candidate-adjacency table; weights are not edge-level P values or FDR values. The full run has 260 regulons and 68,492 edges.

The six inference notebooks preserve historical outputs and were path-adapted; training was not repeated during package cleanup. The original difference in constant-gene handling between full and downsampled runs is retained. Downsampling keeps all 117,867 bins and 15,541 gene columns, drawing without replacement per bin to `floor(total_counts / 2)`. Seeds are 20261016–20261020. The workflow does not require detected gene counts to be halved; the rep01 historical consistency check is retained.

Counts are reconstructed integers and have not been confirmed as original UMI counts. The five repeats are technical downsamplings of one T577 input, not five independent individuals. File integrity and numerical rechecks do not establish cross-individual scientific validity.

## External data and reference requirements

The following operations need data not included in the package:

- Regenerating the five downsampled matrices requires the full input matrix and the historical rep01 matrix for the included consistency check. Rerunning GRN inference requires the selected full or downsampled input matrix and reference resources.
- Recomputing native Pearson scores requires the six original result H5AD files and their AUCell matrices.
- Recomputing full-run gene-expression features requires the full expression matrix.
- Rebuilding the macaque cisTarget references requires the genome, gene annotation, motif archive and annotation, and ranking-generation tools. The two ranking files are approximately 549–550 MiB each; the motif annotation table is approximately 57 MiB. They are not bundled.

The packaged reference material describes *Macaca fascicularis*, assembly Macaca_fascicularis_5.0, Ensembl release 91, and the v10nr public motif collection. Historical source-file ordering for the two ranking windows is not established here; do not infer which file was consumed first from the listed order. Configure the required file paths in `config.json` before rebuilding or training.

After the configuration and English-language updates, `01_export_final_edges.ipynb`, `02_network_topk.ipynb`, and `06_validate_package.ipynb` were executed successfully, with their actual outputs saved. The checks confirmed native edge membership and weights, agreement with the packaged network/Top-k baselines, and the six figure files. Other notebook outputs were retained. Training, downsampling, native Pearson recomputation from H5AD, and full-expression recomputation were not rerun during cleanup. Environment snapshots are records of prior versions, not proof of a successful installation on another machine.
