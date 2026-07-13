# CEPE (Cost Effective Panel Evaluation) scripts

This repository contains the analysis scripts associated with a manuscript currently under review:

*Comparing Whole-Genome Sequencing and Cost-Effective Genotyping Strategies for Estimating Homozygosity-by-Descent*

All the code underlying the manuscript is already shared here. Some parts may still be refined before the review process concludes.

## Guide to the `code` folder

This section summarizes the contents of the [`code/`](https://github.com/nsforneris/CEPE/tree/main/code) folder.

The folder is organized into five subfolders, one per dataset/analysis (`DAMONA`, `Labradors`, `NordicRed`) plus one for helper functions (`helper_functions`). Inside `DAMONA` there are, in turn, four subfolders (`WGS`, `Arrays`, `RRS`, `lowcovWGS`).

Broadly, the repository's workflow is: **(1)** prepare/filter genotypes according to the genotyping strategy (WGS, array, RRS/GBS, low-coverage WGS), **(2)** run `RZooRoH` to estimate HBD (Homozygosity-By-Descent) segments and inbreeding coefficients, and **(3)** post-process/summarize those results.

---

## `code/DAMONA/WGS/`
Whole-genome-sequencing (WGS) Holstein cattle data (DAMONA).

| File | What it does |
|---|---|
| `get_snp_panel_format.sh` | SLURM script that extracts/converts the WGS genotype file in the desired format from a VCF file. |
| `run_beagle.sh` | SLURM script (array job, one task per chromosome 1–29) that phases genotypes with Beagle; reads a VCF file already filtered by chromosome. |
| `run_zooroh_wgs.R` | Runs `RZooRoH` (14-class HBD model) on WGS data for a subset of 12 of the 132 Damona samples (useful when `localhbd = TRUE`) and saves the result to a `.RData` file. |
| `run_zooroh_wgs.sh` | SLURM script (array job, one task per subset of samples 1–11) that runs `run_zooroh_wgs.R` on the cluster. |
| `het_in_window.sh` | SLURM script (array job, one task per chromosome 1–29) that compiles and runs `het_in_window_25.f90` / `het_in_window_50.f90` on the WGS genotypes, computing, for each SNP, the heterozygosity of each individual within a physical window around it. |
| `het_in_window_50.f90` | Fortran program that, for each SNP, computes each individual's heterozygosity within a **50 kb** window (±25 kb) around it (minimum 100 SNPs in the window). |
| `het_in_window_25.f90` | Same as above but with **25 kb** windows (±12.5 kb, minimum 50 SNPs). |


## `code/DAMONA/Arrays/`
6K and 30K SNP genotyping panels derived from the DAMONA's WGS data.

| File | What it does |
|---|---|
| `run_zooroh_arrays.R` | Runs `RZooRoH` (14-class HBD model) on the array genotypes (`6K_Array_gen.txt` or `30K_Array_gen.txt`) and saves the result to a `.RData` file. |


## `code/DAMONA/RRS/`
Generation and analysis of simulated RRS (*Reduced Representation Sequencing*) panels, built from DAMONA's WGS data after an in silico digestion with the PstI enzyme.

| File | What it does |
|---|---|
| `RRS_panels_generation.R` | Main script: starting from the PstI digestion fragments and the sequence SNPs, it detects SNPs at the enzyme's cut site (which cause *allelic dropout*), corrects the affected haplotypes, computes call rate and MAF, filters (call rate ≥0.95, MAF ≥0.01), and exports the final VCF files for the **RRS50K**, **RRS30K**, and **RRS15K** panels (the latter obtained by subsampling 50% of the RRS30K fragments). |
| `RRS_panels_generation.sh` | SLURM script that runs `RRS_panels_generation.R` on the cluster, copying inputs to local scratch, and then converts the output VCFs into the `gt` format (0/1/2/9) used by RZooRoH for the three panels (RRS15K/30K/50K). |
| `RRS_panels_generation_lowcov.R` | Adapted version of the previous script for combining RRS with low-coverage sequencing: instead of exporting an already-corrected VCF, it generates two text files (positions to downsample and read-count corrections for dropout) that are later used by `makeRRSlowcov.f90`. |
| `RRS_panels_generation_lowcov.sh` | SLURM script that runs `RRS_panels_generation_lowcov.R` on the cluster. |
| `makeRRSlowcov.f90` | Fortran program that reads the read counts (AD) from the full sequence, filters down to the positions of the RRS panel of interest, applies the dropout corrections computed in R, and downsamples the reads to simulate low coverage (e.g. RRS15K@2x or @5x). |
| `makeRRSlowcov.sh` | SLURM script that compiles and runs `makeRRSlowcov.f90`, and formats the resulting VCF as a final genotype file. |
| `makeRRSlowcov.f90_seeds` | List of the random seeds used for each panel × coverage combination (RRS15K@2x, RRS15K@5x, RRS30K@2x, RRS30K@5x), to be set manually in `makeRRSlowcov.f90` before compiling. |
| `run_zooroh_rrs.R` | Runs `RZooRoH` on the full-coverage RRS panels (RRS15K/30K/50K, genotype format). |
| `run_zooroh_rrs_lowcov.R` | Runs `RZooRoH` on the low-coverage RRS panels (`ad` format). |


## `code/DAMONA/lowcovWGS/`
Simulation of low-coverage sequencing from DAMONA's WGS data.

| File | What it does |
|---|---|
| `makelowcov.f90` | Fortran program that, from the read counts (AD) of the full WGS data, downsamples each animal's reads to a uniform target coverage (0.1x, 0.2x, 0.5x, 1x, 2x, or 5x, depending on the `pcov`/`seed` values noted in the comments) and writes the resulting low-coverage VCF. |
| `makelowcov.sh` | SLURM script that compiles and runs `makelowcov.f90` and formats the output file as a genotype file for RZooRoH. |
| `makelowcov_combined.f90` | Variant of `makelowcov.f90` that generates a single dataset in which animals are split into three equal groups with different, simultaneous coverage levels (2x, 5x, and 10x), randomly assigned. |
| `makelowcov_combined.sh` | SLURM script that compiles and runs `makelowcov_combined.f90`. |
| `run_zooroh_lcWGS.R` | Runs `RZooRoH` on the low-coverage WGS data (`ad` format), splitting the analysis into subsets of individuals (useful when `localhbd = TRUE`). It is called from the command line with the genotype file, coverage, and subset number as arguments. |


## `code/Labradors/`
Preparation of Labrador Retriever data sequenced at different coverage levels (0.9x, 3.8x, 43.5x) and their analysis with RZooRoH using genotype likelihoods (PL).

| File | What it does |
|---|---|
| `prep40x.sh` | SLURM script for the 43.5x data: concatenates, filters by depth (`--min-meanDP`/`--max-meanDP`), selects biallelic SNPs, and generates both the PL file and the list of selected positions (`selected_positions.txt`) later used by `prep1x.sh`/`prep4x.sh`. |
| `prep4x.sh` | SLURM script: concatenates the per-chromosome VCFs of the 3.8x sequence data, filters biallelic SNPs with `vcftools`/`bcftools`, keeps only the previously selected positions (`selected_positions.txt`), and formats the phred-scaled genotype likelihoods (PL) into a text file (`Labradors4x_SNPs4_PL.txt`) ready for RZooRoH. |
| `prep1x.sh` | Same as `prep4x.sh` but for the 0.9x data. |
| `Run_40X.R` | Runs `RZooRoH` (genotype-likelihood format `"gl"`) on the 43.5x data. |
| `Run_1X.R` | Same as above but on `Labradors1x_SNPs4_PL.txt`. |
| `Run_4X.R` | Same as above but on `Labradors4x_SNPs4_PL.txt`. |


## `code/NordicRed/`
Filtering of Nordic Red cattle data genotyped by GBS (18K/20K) and by WGS (with two filtering stringency levels), and their subsequent analysis with RZooRoH.

| File | What it does |
|---|---|
| `wgs_more_stringent_filtering.sh` | Filters the WGS VCF (`bcftools`), keeping only positions also present in the DAMONA set (stricter filtering), retains biallelic autosomal SNPs, and exports genotypes in gt format (0/1/2/9) ready for RZooRoH (`filtered_DAMONA_biallelic_snps_auto_gt.gen`). |
| `wgs_less_stringent_filtering.sh` | Filters the same WGS VCF, keeping only variants flagged `PASS` (less strict filtering, without restricting to the DAMONA positions), and exports the analogous genotype file (`filtered_PASS_biallelic_snps_auto_gt.gen`). |
| `GBS20K_filtering.sh` | Filters the original GBS VCF (`PASS`, biallelic, autosomal) and exports the **GBS20K** genotype panel. |
| `GBS18K_filtering.sh` | Starting from the filtered GBS20K panel, keeps only the SNPs also present in `filtered_WGS18K_gt.gen`, producing the **GBS18K** panel (GBS ∩ WGS intersection). |
| `WGS18K_filtering.sh` | Starting from the filtered WGS data, keeps only the SNPs present in the filtered GBS20K panel, producing the **WGS18K** panel (same intersection, from the sequence side). |
| `run_zooroh_nordicred.R` | Runs `RZooRoH` (`"gt"` format) on any of the panels generated above (strictly/loosely filtered WGS, WGS18K, GBS20K, or GBS18K — the active file is selected via comments). |


## `code/helper_functions/`

| File | What it does |
|---|---|
| `pvit.R` | Post-processing helper code: from a `RZooRoH` results `.RData`, it uses the HBD segments from the Viterbi approach to compute, for each individual, the proportion of the genome in HBD classes cumulated up to each of the 14 HBD thresholds (2, 4, 8 ... 16384), and saves the result to `pvit.RData`. Run as `Rscript pvit.R 'file.RData'`. |
| `pHBD.R` | Post-processing helper code: from a `RZooRoH` results `.RData`, it uses the HBD probabilities to compute, for each individual, the proportion of the genome in HBD classes cumulated up to each of the 14 HBD thresholds (2, 4, 8 ... 16384), and saves the result to `pHBD.RData`. Run as `Rscript pHBD.R 'file.RData'`. |

---

### General notes
- Almost all `.sh` files are examples of **SLURM** submission scripts that stage data in a local scratch directory, run the corresponding step (R, compiled Fortran, or `bcftools`/`vcftools`), and copy the results back to shared scratch.
- The `.f90` programs were be compiled with `ifort` (Intel Fortran) before running; the comments in each file indicate the compilation command and, where relevant, which parameters (`pcov`, `seed`, file names) need to be adjusted for the desired panel/coverage.
- All `run_zooroh_*.R` scripts use the **RZooRoH** package v0.4.1 with a **14-class HBD model** (`K = 14`), varying the genotype input format (`gt` = called genotype, `gl`/`PL` = genotype likelihood, `ad` = read counts) depending on the genotyping strategy.
- RZooRoH is available at https://cran.r-project.org/web/packages/RZooRoH
