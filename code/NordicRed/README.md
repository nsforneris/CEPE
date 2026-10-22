# Nordic Red

## Workflow and dependencies

1. Run `GBS20K_filtering.sh` on `fine_tuning_GBS_data-Cow_GBS-250.vcf.gz`. It creates `filtered_GBS20K_gt.gen` and `filtered_PASS_GBS_biallelic_snps_auto.vcf.gz`.
2. Run `WGS18K_filtering.sh` on `fine_tuning_GBS_data-Cow_WGS.vcf.gz`. It uses `filtered_GBS20K_gt.gen` from Step 1 and creates `filtered_WGS18K_gt.gen`.
3. Run `GBS18K_filtering.sh`. It uses `filtered_PASS_GBS_biallelic_snps_auto.vcf.gz` from Step 1 and `filtered_WGS18K_gt.gen` from Step 2, creating `filtered_GBS18K_gt.gen`.
4. The alternative WGS panels are generated independently by `wgs_more_stringent_filtering.sh` and `wgs_less_stringent_filtering.sh`; both use `fine_tuning_GBS_data-Cow_WGS.vcf.gz`.
5. Run `run_zooroh_nordicred.R` after selecting the desired genotype file and its matching sample-ID file in the script.

The filtering scripts write some sample-ID files under the same name (`sample_ids.txt` or `sample_ids_gbs.txt`); preserve the matching file when switching between panels.
