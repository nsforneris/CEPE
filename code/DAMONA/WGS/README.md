# DAMONA WGS

## Workflow

1. Run `get_snp_panel_format.sh` with `DAMONA132_CEPE_refined.vcf.gz` and `DAMONA132_CEPE_unrefined.vcf.gz`. It creates `seq_evalset_gen.txt`, `Damona132_snp_AD`, `seq_unrefined_gen.txt`, `6K_Array_gen.txt`, and `30K_Array_gen.txt`.
2. For the WGS RZooRoH analysis, provide `seq_evalset_gen.txt.gz` at the location configured in `run_zooroh_wgs.sh`, then run `run_zooroh_wgs.sh`. It calls `run_zooroh_wgs.R` for sample subsets.
3. For window heterozygosity, provide `seq_evalset_gen.txt.gz` at the location configured in `het_in_window.sh`, then run that script. `het_in_window_25.f90` and `het_in_window_50.f90` contain the two window-size implementations.

`Damona132_snp_AD` is also an input to the low-coverage WGS and RRS simulations.
