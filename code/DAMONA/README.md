# DAMONA

[← Back to main README](../../README.md)

Whole-genome sequencing data from 132 Holstein cattle (8,417,679 biallelic autosomal SNPs, chromosomes 1–29). These data are the reference, and the other genotyping strategies are simulated from them.

| Folder | Strategy | Needs input from |
|---|---|---|
| [`WGS`](WGS) | WGS reference: **run first** | — |
| [`Arrays`](Arrays) | 6K and 30K SNP arrays | `WGS` |
| [`RRS`](RRS) | Reduced-representation sequencing (RRS15K/30K/50K), full and low coverage | `WGS` (low coverage only) |
| [`lowcovWGS`](lowcovWGS) | Low-coverage WGS (0.1x–5x, and a 2x/5x/10x mix) | `WGS` |

The two input VCFs are `DAMONA132_CEPE_refined.vcf.gz` (genotypes) and `DAMONA132_CEPE_unrefined.vcf.gz` (read counts and array genotypes).
