# DAMONA – RRS panels

[← DAMONA](../README.md) · [Main README](../../../README.md)

Simulated reduced-representation sequencing (RRS/GBS) panels built from the DAMONA WGS data after an in silico digestion with the *PstI* enzyme. The simulation accounts for **allelic dropout**, which happens when a SNP sits in the restriction site.

## Principle

1. **In silico digestion.** GBSX takes the reference genome (FASTA) and *PstI*, finds the restriction sites and returns the fragments (BED). This is done beforehand and is not part of this folder.
2. **Candidate RRS SNPs.** These are the WGS SNPs that fall inside the selected fragments.
3. **Phasing with Beagle.** Allelic dropout acts on a whole haplotype: if one chromosome carries a SNP that destroys the restriction site, that chromosome's fragment is not sequenced, so every SNP in the fragment loses the alleles on that chromosome. Phased haplotypes are needed to know which alleles sit on the same chromosome as the cut-site allele.
4. **Dropout correction.** This applies to fragments with a SNP within 3 bp of the cut site:
   - Heterozygous at the cut site: only the haplotype with the REF allele is kept, so the animal looks homozygous across the fragment.
   - Homozygous ALT: the fragment is not produced, so the genotypes are missing.
   - Homozygous REF: nothing changes.
5. **Quality control.** SNPs are kept if call rate ≥ 0.95 and MAF ≥ 0.01.
6. **Panels by fragment length.**
   - RRS50K: fragments of 200–300 bp.
   - RRS30K: fragments of 250–300 bp.
   - RRS15K: a random 50 % of the RRS30K fragments.
7. **HBD estimation** with RZooRoH.

For **RRS at low coverage** (RRS15K/30K at 2x and 5x), the dropout detection is the same, but the corrections are applied to the WGS read counts, which are then downsampled.

## Inputs prepared beforehand

| File | Content |
|---|---|
| `genome.PstI.1000nt.digest.bed` | *PstI* fragments from the GBSX in silico digestion |
| `RRS_PstI_snp_all.bed` | WGS SNPs inside the candidate RRS fragments |
| `SNP_seq.bed` | All WGS SNPs, used to find variants at the restriction sites |
| `DAMONA132_CEPE_refined.vcf.gz` (+ index) | WGS genotypes, input for Beagle |
| `Damona132_snp_AD` | WGS read counts, created by [`../WGS/get_snp_panel_format.sh`](../WGS) (low coverage only) |

## Order of steps

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | `run_beagle.sh` | `DAMONA132_CEPE_refined.vcf.gz` | `beagle_chr{1..29}.vcf.gz` (phased) |
| 2 | `RRS_panels_generation.sh` (calls `RRS_panels_generation.R`) | BED files + phased VCFs | `RRS15K_gen.txt`, `RRS30K_gen.txt`, `RRS50K_gen.txt` |
| 3 | `run_zooroh_rrs.R` | `RRS*K_gen.txt` | `zr_RRS*K_Array.RData` |
| *Low coverage only:* | | | |
| 4 | `RRS_panels_generation_lowcov.sh` (calls `RRS_panels_generation_lowcov.R`) | BED files + phased VCFs | `RRS15K_positions_to_downsample.txt`, `RRS30K_positions_to_downsample.txt`, `reads_to_correct.txt` |
| 5 | `makeRRSlowcov.sh` (compiles `makeRRSlowcov.f90`) | Step 4 files + `Damona132_snp_AD` | `out_vcf_RRS*K_*x_gen.txt` |
| 6 | `run_zooroh_rrs_lowcov.R` | `out_vcf_RRS*K_*x_gen.txt` | `zr_RRS*K_*x.RData` |

Then summarise the results with [`helper_functions`](../../helper_functions). The file `order_of_steps` lists the same order in short form.

## Files

| File | What it does |
|---|---|
| `run_beagle.sh` | SLURM array job (one task per chromosome) that extracts each chromosome and phases it with Beagle. |
| `RRS_panels_generation.R` | Detects cut-site SNPs, applies the dropout correction to the phased haplotypes, filters on call rate and MAF, and writes the RRS50K, RRS30K and RRS15K VCFs. |
| `RRS_panels_generation.sh` | Runs the R script and converts the VCFs to RZooRoH `gt` format. |
| `RRS_panels_generation_lowcov.R` | Same detection as above, but it writes the positions of each panel and the read-count corrections instead of corrected genotypes. |
| `RRS_panels_generation_lowcov.sh` | Runs the R script above. |
| `makeRRSlowcov.f90` | Keeps the read counts at the panel positions, applies the dropout corrections and downsamples the reads to the target coverage. |
| `makeRRSlowcov.sh` | Compiles and runs the Fortran program and formats the output for RZooRoH. |
| `makeRRSlowcov.f90_seeds` | Random seed used for each panel × coverage combination. |
| `run_zooroh_rrs.R` | RZooRoH on the full-coverage panels (`gt` format). |
| `run_zooroh_rrs_lowcov.R` | RZooRoH on the low-coverage panels (read counts, `ad` format). |

## Choosing the panel and coverage

- In `run_zooroh_rrs.R` and `run_zooroh_rrs_lowcov.R`, uncomment the input file and the matching output name.
- In `makeRRSlowcov.f90`, set the following, then recompile once for each panel × coverage combination:
  - the positions file (RRS15K or RRS30K);
  - `pcov` (0.10 for 2x, 0.25 for 5x);
  - `seed` (see `makeRRSlowcov.f90_seeds`);
  - the output file names.
- In `makeRRSlowcov.sh`, change the file names to match.
