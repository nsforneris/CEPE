# DAMONA – WGS

[← DAMONA](../README.md) · [Main README](../../../README.md)

WGS reference analysis of the 132 Holstein animals. This folder also creates the input files for [`Arrays`](../Arrays), [`lowcovWGS`](../lowcovWGS) and the low-coverage part of [`RRS`](../RRS), so **run step 1 first**.

## Order of steps

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | `get_snp_panel_format.sh` | `DAMONA132_CEPE_refined.vcf.gz` | `seq_evalset_gen.txt` (WGS genotypes, RZooRoH `gt` format) |
| | | `DAMONA132_CEPE_unrefined.vcf.gz` + [`../Arrays/*_ids.txt`](../Arrays) | `Damona132_snp_AD` (read counts, used by `lowcovWGS` and `RRS`), `seq_unrefined_gen.txt`, `6K_Array_gen.txt`, `30K_Array_gen.txt` |
| 2 | `run_zooroh_wgs.sh` (calls `run_zooroh_wgs.R`) | `seq_evalset_gen.txt.gz` | `zr_wgs_<ini>-<fin>.RData` |
| 3 (optional) | `het_in_window.sh` (compiles `het_in_window_25.f90` or `het_in_window_50.f90`) | `seq_evalset_gen.txt.gz` | Per-SNP heterozygosity of each animal in a window around the SNP |

Then summarise the RZooRoH results with [`helper_functions`](../../helper_functions).

## Files

| File | What it does |
|---|---|
| `get_snp_panel_format.sh` | Converts the VCFs into the WGS genotype file, the read-count (AD) file and the 6K/30K array genotype files. |
| `run_zooroh_wgs.R` | Runs RZooRoH on one subset of 12 animals. Usage: `Rscript run_zooroh_wgs.R <genotype_file> <subset>`. |
| `run_zooroh_wgs.sh` | SLURM array job (subsets 1–11 = 132 animals) that runs `run_zooroh_wgs.R`. The animals are split into subsets because `localhbd = TRUE` uses a lot of memory. |
| `het_in_window.sh` | SLURM array job (one task per chromosome) that compiles and runs the heterozygosity program. Keep only the `ifort` line for the window size you want: both lines create the same executable. |
| `het_in_window_50.f90` | Heterozygosity in a 50 kb window (±25 kb, minimum 100 SNPs). |
| `het_in_window_25.f90` | Heterozygosity in a 25 kb window (±12.5 kb, minimum 50 SNPs). |
