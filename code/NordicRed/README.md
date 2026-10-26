# Nordic Red cattle

[← Main README](../../README.md)

The same Nordic Red animals were genotyped by GBS and by WGS. Unlike [`DAMONA/RRS`](../DAMONA/RRS), these GBS data are real, so no simulation or phasing is needed. All scripts keep biallelic autosomal SNPs and write genotypes in RZooRoH `gt` format.

## Order of steps

**A. Full WGS panels.** The two scripts are independent of each other.

| Script | Input | Output |
|---|---|---|
| `wgs_more_stringent_filtering.sh` | `fine_tuning_GBS_data-Cow_WGS.vcf.gz` + DAMONA SNP list | `filtered_DAMONA_biallelic_snps_auto_gt.gen` (only positions also in the DAMONA set) |
| `wgs_less_stringent_filtering.sh` | `fine_tuning_GBS_data-Cow_WGS.vcf.gz` | `filtered_PASS_biallelic_snps_auto_gt.gen` (all `PASS` variants) |

**B. GBS panels and GBS–WGS intersection.** Run these in this order.

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | `GBS20K_filtering.sh` | `fine_tuning_GBS_data-Cow_GBS-250.vcf.gz` | `filtered_GBS20K_gt.gen`, `filtered_PASS_GBS_biallelic_snps_auto.vcf.gz` |
| 2 | `WGS18K_filtering.sh` | `fine_tuning_GBS_data-Cow_WGS.vcf.gz` + `filtered_GBS20K_gt.gen` | `filtered_WGS18K_gt.gen` (WGS SNPs also in GBS20K) |
| 3 | `GBS18K_filtering.sh` | Step 1 VCF + `filtered_WGS18K_gt.gen` | `filtered_GBS18K_gt.gen` (GBS SNPs also in WGS18K) |

**C. RZooRoH.** In `run_zooroh_nordicred.R`, uncomment the genotype file, its sample-ID file and the matching output name:
- WGS panels use `sample_ids.txt`.
- GBS panels use `sample_ids_gbs.txt`.

Several filtering scripts write a sample-ID file with the same name, so keep the matching file when you switch between panels.

Then summarise the results with [`helper_functions`](../helper_functions).
