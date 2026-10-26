# CEPE (Cost Effective Panel Evaluation) scripts

This repository contains the analysis scripts associated with the manuscript:

Forneris N, Wragg D, Fischer D, Schoenebeck J, Clements D, Iso-Touru T, Gautier M, Druet T. (in press). Comparing Whole-Genome Sequencing and Cost-Effective Genotyping Strategies for Estimating Homozygosity-by-Descent. *Heredity*.

## Overview

Homozygosity-by-descent (HBD) and inbreeding are estimated with [RZooRoH](https://cran.r-project.org/web/packages/RZooRoH) on several genotyping strategies, and the results are compared with whole-genome sequencing (WGS). Every analysis follows the same three steps:

1. **Prepare the genotypes** for a given strategy (WGS, SNP array, RRS/GBS, low-coverage WGS).
2. **Run RZooRoH** to estimate HBD segments and inbreeding coefficients.
3. **Summarise the results** with the scripts in [`helper_functions`](code/helper_functions).

## Repository map

Each folder has its own README that describes its files and the order in which to run them.

| Folder | Dataset | Content |
|---|---|---|
| [`code/DAMONA`](code/DAMONA) | Holstein cattle, 132 animals, high-coverage WGS | WGS reference and the strategies simulated from it |
| ├─ [`WGS`](code/DAMONA/WGS) | | WGS reference analysis; **start here**, it creates the inputs for the other DAMONA folders |
| ├─ [`Arrays`](code/DAMONA/Arrays) | | Simulated 6K and 30K SNP arrays |
| ├─ [`RRS`](code/DAMONA/RRS) | | Simulated reduced-representation sequencing (PstI), at full and low coverage |
| └─ [`lowcovWGS`](code/DAMONA/lowcovWGS) | | Simulated low-coverage WGS (0.1x to 5x) |
| [`code/Labradors`](code/Labradors) | Labrador Retrievers, real WGS at 0.9x, 3.8x and 43.5x | Low-coverage sequencing with genotype likelihoods |
| [`code/NordicRed`](code/NordicRed) | Nordic Red cattle, real GBS and WGS on the same animals | GBS compared with WGS |
| [`code/helper_functions`](code/helper_functions) | All | Post-processing of RZooRoH results |

## Software

- R with [RZooRoH](https://cran.r-project.org/web/packages/RZooRoH) v0.4.1 (all analyses use a 14-class HBD model)
- BCFtools, VCFtools
- Beagle 5.4 (RRS only)
- Intel Fortran (`ifort`) for the `.f90` programs
- GBSX for the in silico digestion (RRS only, run beforehand)

## General notes

- Most `.sh` files are example **SLURM** submission scripts. Paths, modules and resources refer to the cluster we used and need to be adapted.
- Some scripts handle several panels or coverages: you choose which one to run by commenting/uncommenting lines or editing parameters. Each folder README says which ones.
