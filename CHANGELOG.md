# Changelog

## v1.0.1 — documentation and configuration templates
- Standardized the default `sample_group_regex` example as `^((?:SM)[0-9]{3})` and documented how to adapt capture group 1 to custom naming schemes, including `IB001`, `IBC1`/`TN71`, and `ABC_001`-style IDs.

- Expanded English and Traditional Chinese README files with reference-genome indexing commands, stage-by-stage usage, config guidance, Snakemake commands, resource controls, and Telegram setup.
- Expanded every `configs/*.example.yaml` template to expose all user-adjustable settings while retaining stage-specific example paths.
- Standardized sample naming examples to `SM001`–`SM999` and the corresponding prefix-capture regex.
- Removed `docs/GIT_ON_MAC.md` from the public project.
- No SNP-calling algorithm or core workflow logic changed from v1.0.0.

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
