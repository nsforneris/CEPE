# DAMONA – SNP arrays

[← DAMONA](../README.md) · [Main README](../../../README.md)

Simulated low-density (6K) and medium-density (30K) SNP arrays, obtained by keeping only the array SNPs from the DAMONA sequence data.

## Order of steps

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | [`../WGS/get_snp_panel_format.sh`](../WGS) | `DAMONA132_CEPE_unrefined.vcf.gz` + `6K_ids.txt` / `30K_ids.txt` | `6K_Array_gen.txt`, `30K_Array_gen.txt` |
| 2 | `run_zooroh_arrays.R` | `6K_Array_gen.txt` or `30K_Array_gen.txt` | `zr_6K_Array.RData` or `zr_30K_Array.RData` |

Then summarise the results with [`helper_functions`](../../helper_functions).

## Files

| File | What it does |
|---|---|
| `6K_ids.txt`, `30K_ids.txt` | SNPs of each array (`chr_pos`). |
| `run_zooroh_arrays.R` | Runs RZooRoH on one array. The 6K panel is active; to run the 30K panel, uncomment the 30K input and output lines. |
