# Helper functions

[← Main README](../../README.md)

Post-processing of RZooRoH results, used after any of the analyses in [`DAMONA`](../DAMONA), [`Labradors`](../Labradors) or [`NordicRed`](../NordicRed).

For each animal, both scripts compute the proportion of the genome in HBD classes up to each of the 14 thresholds (T = 2, 4, 8 … 16384).

| Script | Based on | Output | Usage |
|---|---|---|---|
| `pHBD.R` | HBD probabilities | `pHBD.RData` | `Rscript pHBD.R zr_6K_Array.RData` |
| `pvit.R` | HBD segments (Viterbi) | `pvit.RData` | `Rscript pvit.R zr_6K_Array.RData` |

`pvit.R` needs HBD segments, so it cannot be used on runs with `vit = FALSE`, such as the Labrador runs. The output is written to the current working directory.

