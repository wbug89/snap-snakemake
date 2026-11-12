# Changelog

## v1.0.0 — first formal release

This is the first formally versioned release of **SNaP — Snakemake for NGS and Polyploids**.

### Core workflow

- Supports `SE` and `PE` reads from one Snakefile.
- Supports `preprocess`, `gvcf`, `joint`, `joint_shard`, `joint_finalize`, and `full` stages.
- Supports multi-batch sequencing units and biological-sample grouping.
- Preserves RAD-seq/DArT-seq duplicate-flagged reads for variant calling:
  - `MarkDuplicates --REMOVE_DUPLICATES false`
  - `--REMOVE_SEQUENCING_DUPLICATES false`
  - HaplotypeCaller `--disable-read-filter NotDuplicateReadFilter`
- Supports configurable ploidy, including hexaploid calling.
- Uses Trimmomatic bundled adapters first and user-supplied adapter paths as fallback.
- Supports optional Telegram start/success/failure notifications.

### Joint genotyping

- Uses per-contig `GenomicsDBImport -> GenotypeGVCFs` instead of cohort-wide `CombineGVCFs`.
- Supports distributed joint calling across multiple servers by assigning reference contigs/LGs to partitions.
- Supports central `joint_finalize` concatenation and filtering in reference `.fai` order.
- Explicitly retains `--max-alternate-alleles 6` and `--max-genotype-count 1024` for polyploid memory safety.

### Variant identifiers

- SNP IDs are assigned before final SNP filtering as `CHROM_POS` using:

  ```bash
  bcftools annotate --set-id '%CHROM_%POS'
  ```

- Example IDs: `LG1_123456`, `LG5_9876543`.
- The final `final_combined_filtered_snps.vcf.gz` retains these IDs.

### Versioning

- Formal release version reset from the internal development series (`v3.x`) to semantic version `v1.0.0`.
- Future public releases should increment from `v1.0.0` using semantic versioning.
