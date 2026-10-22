# RZooRoH post-processing helpers

Both scripts take an RZooRoH result `.RData` file as their first command-line argument.

- `pvit.R` uses the Viterbi HBD segments and saves `pvit.RData`.
- `pHBD.R` uses HBD probabilities and saves `pHBD.RData`.

Example: `Rscript pvit.R 'zr_6K_Array.RData'` or `Rscript pHBD.R 'zr_6K_Array.RData'`.

Run these after the relevant RZooRoH analysis. Each script writes its output to the current working directory.
