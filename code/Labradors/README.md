# Labrador Retriever WGS

## Workflow

1. Run `prep40x.sh` first. It creates `selected_positions.txt` and `Labradors50x_SNPs4_PL.txt`.
2. Run `prep4x.sh` and/or `prep1x.sh`. Both use `selected_positions.txt` from Step 1 and create `Labradors4x_SNPs4_PL.txt` and `Labradors1x_SNPs4_PL.txt`, respectively.
3. Run the matching analysis script: `Run_40X.R`, `Run_4X.R`, or `Run_1X.R`. Each script reads its corresponding PL file and saves an RData result.

The preparation scripts use per-chromosome VCFs from the paths configured inside each script.
