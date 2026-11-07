# T577 GRN Package Scope and Acceptance

## Purpose

Preserve a portable research package for the current six figures and the associated full and five downsampled GRN runs. Keep the original working files unchanged. This plan documents the archival scope; it does not claim that workflows not actually run have been rerun or validated.

## Package scope

The package retains:

- Six final figures as PNG and vector PDF, with the plotting notebook and direct plotting inputs.
- Native regulon and motif outputs for the full run and five repeats, plus weighted final-edge tables and their export code.
- Downsampling, inference, and analysis notebooks with their historical outputs.
- `config.json` for runtime analysis paths and training inputs, references, downsampling parameters, and output paths.
- Reference notes, TF list, environment snapshots, and the figure caption.

The current figure set uses purple for Top-k and cividis for both TF-expression panels. Older four-group, fixed full-signature AUCell, and statistical-annotation figures are excluded.

## Acceptance criteria

A package review should confirm that:

1. The six PNG/PDF pairs and root plotting notebook correspond to the current figure definitions and use the included plotting inputs.
2. Packaged `final_edges.tsv.gz` tables correspond to the native final regulon members and retain the documented weight, self-edge, and negative-weight conventions.
3. Lightweight recomputation and package validation use packaged inputs and keep recomputed outputs separate from the confirmed figure-input baselines.
4. `config.json` clearly separates analysis settings from training settings and makes external input, reference, downsampling, and output paths explicit.
5. Documentation distinguishes packaged historical outputs from workflows run during the current cleanup and identifies the external files required for full retraining and recomputation.

## Execution boundary

Following the configuration and English-language updates, the final-edge export, network/Top-k comparison, and package-validation notebooks were rerun successfully and their actual outputs saved. The data and figures were retained. Training, downsampling, native Pearson recomputation from H5AD, and full-expression recomputation were not rerun during cleanup. No scientific result should be described as newly reproduced unless its corresponding workflow was actually executed.
