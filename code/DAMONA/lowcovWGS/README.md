# DAMONA – low-coverage WGS

[← DAMONA](../README.md) · [Main README](../../../README.md)

Simulates low-coverage sequencing by randomly downsampling the reads of the DAMONA WGS data. At low coverage the genotype calls are unreliable, so RZooRoH uses the read counts directly (`ad` format) instead of called genotypes.

## Order of steps

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | [`../WGS/get_snp_panel_format.sh`](../WGS) | `DAMONA132_CEPE_unrefined.vcf.gz` | `Damona132_snp_AD` |
| 2a | `makelowcov.sh` (compiles `makelowcov.f90`) | `Damona132_snp_AD` | `out_vcf_lcWGS_<cov>x_gen.txt`, all animals at the same coverage |
| 2b | `makelowcov_combined.sh` (compiles `makelowcov_combined.f90`) | `Damona132_snp_AD` | `out_vcf_lcWGS_2x5x10x_gen.txt`, one third of the animals at each coverage |
| 3 | `run_zooroh_lcWGS.R` | Step 2 output | `zr_lcwgs_<cov>x_<ini>-<fin>.RData` |

Then summarise the results with [`helper_functions`](../../helper_functions).

## Choosing the coverage (step 2a)

Set `pcov` and `seed` in `makelowcov.f90`, update the file names in the `.f90` and `.sh` files, then recompile:

| Coverage | `pcov` | `seed` |
|---|---|---|
| 0.1x | 0.005 | 541986 |
| 0.2x | 0.01 | 541987 |
| 0.5x | 0.025 | 541985 |
| 1x | 0.05 | 541984 |
| 2x | 0.10 | 541983 |
| 5x | 0.25 | 541988 |

## Running RZooRoH (step 3)

`Rscript run_zooroh_lcWGS.R <genotype_file> <coverage> <subset>`, for example `Rscript run_zooroh_lcWGS.R out_vcf_lcWGS_0.5x_gen.txt 0.5 2`.

Each subset has 12 animals (subsets 1–11 = 132 animals). To submit it on SLURM, [`../WGS/run_zooroh_wgs.sh`](../WGS) can be used as a template.
