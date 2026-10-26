# Labrador Retrievers

[← Main README](../../README.md)

Real WGS data from Labrador Retrievers sequenced at three coverages: 0.9x, 3.8x and 43.5x. RZooRoH is run on genotype likelihoods (PL, `gl` format) so that genotype uncertainty at low coverage is taken into account.

Files named `50x` contain the 43.5x data.

## Order of steps

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | `prep40x.sh` | Per-chromosome VCFs (43.5x), `biallelic_positions1.txt` | `Labradors50x_SNPs4_PL.txt`, `selected_positions.txt` |
| 2 | `prep4x.sh`, `prep1x.sh` | Per-chromosome VCFs (3.8x / 0.9x), `selected_positions.txt` | `Labradors4x_SNPs4_PL.txt`, `Labradors1x_SNPs4_PL.txt` |
| 3 | `Run_40X.R`, `Run_4X.R`, `Run_1X.R` | Matching PL file | `Lab14L.RData`, `Lab14L_4x.RData`, `Lab14L_1x.RData` |

Step 1 must run first, because it defines the SNPs (`selected_positions.txt`) used for all three coverages.

Then summarise the results with [`helper_functions`](../helper_functions). These runs use `vit = FALSE` and produce no HBD segments, so use `pHBD.R`, not `pvit.R`.

## Files

| File | What it does |
|---|---|
| `prep40x.sh` | Concatenates the VCFs, keeps biallelic SNPs with mean depth between 33 and 48, writes the PL file and the list of selected positions. |
| `prep4x.sh`, `prep1x.sh` | Concatenate the VCFs, keep biallelic SNPs among the selected positions, write the PL file. |
| `Run_40X.R`, `Run_4X.R`, `Run_1X.R` | Run RZooRoH (`gl` format) on each coverage. |
