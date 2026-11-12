# SNaP v1.0.0

### Snakemake for NGS and Polyploids

**SNaP makes SNP calling as simple as a snap.**

SNaP is a reproducible Snakemake workflow for SNP discovery from RAD-seq-style NGS datasets, including ddRAD-seq and DArT-seq. It supports SE/PE reads, configurable ploidy, GVCF-based variant calling, per-contig GenomicsDB joint genotyping, and distributed execution across multiple Linux servers.

The public project is maintained as **SNaP**, with the GitHub repository name **`snap-snakemake`**. Version `v1.0.0` is the first formal release.

Key design principles:

- one workflow for SE and PE data;
- configurable ploidy, including hexaploid calling;
- RAD/DArT-aware duplicate handling;
- GVCF calling by biological sample;
- joint genotyping by reference contig/LG, enabling PC1/PC3 distributed execution;
- final SNP IDs standardized as `CHROM_POS`, for example `LG1_123456`;
- resumable execution through Snakemake and per-contig GenomicsDB workspaces.

## Stages

| stage | purpose |
|---|---|
| `preprocess` | FASTQ → trim → BWA → sort → MarkDuplicates |
| `gvcf` | merge all units for each biological sample → HaplotypeCaller GVCF |
| `joint` | single-host all-contig GenomicsDBImport → GenotypeGVCFs → concat/filter |
| `joint_shard` | process only reference contigs assigned to this server |
| `joint_finalize` | collect all per-contig VCFs → reference-order concat → final filtering |
| `full` | single-host end-to-end workflow |

The distributed rule is important: gVCF calling may be divided by **sample**, but joint calling is divided by **contig**, not by sample. Every PC1/PC3 `joint_shard` run must see the complete cohort of single-sample gVCFs.


## Final SNP ID convention

After joint genotyping and SNP selection, the workflow assigns every SNP VCF `ID` as:

```text
CHROM_POS
```

For example:

```text
LG1_123456
LG5_9876543
```

The rule uses:

```bash
bcftools annotate --set-id '%CHROM_%POS'
```

The SNP path is:

```text
allsample.vcf.gz
  → SelectVariants (SNP)
allsample_snps.vcf.gz
  → set SNP ID = CHROM_POS
allsample_snps.ID.vcf.gz
  → hard filtering
final_combined_filtered_snps.vcf.gz
```

Therefore, `final_combined_filtered_snps.vcf.gz` retains the `CHROM_POS` SNP IDs for downstream PLINK, GWAS, annotation, and cross-analysis matching.

## Distributed joint calling

Copy the same `joint_contig_assignment.tsv` to both servers:

```text
contig  server
LG1     PC1
...
LG9     PC3
...
LG15    PC3
```

PC1:

```bash
snakemake -s Snakefile \
  --configfile config.joint_shard_PC1.yaml \
  --cores 32 \
  --resources joint_jobs=1 mem_mb=200000 \
  --printshellcmds --rerun-incomplete
```

PC3:

```bash
snakemake -s Snakefile \
  --configfile config.joint_shard_PC3.yaml \
  --cores 30 \
  --resources joint_jobs=1 mem_mb=120000 \
  --printshellcmds --rerun-incomplete
```

`joint_jobs` must be supplied as a **global Snakemake resource limit**. The per-rule `resources:` declaration alone does not cap whole-machine parallelism.

Each shard outputs:

```text
results/per_contig/000001.vcf.gz
results/per_contig/000001.vcf.gz.tbi
```

The shard-to-contig mapping is written to `metadata/joint_contigs.tsv` using reference `.fai` order.

After PC1 and PC3 finish, copy only their per-contig VCFs and indexes to the final server. Do not copy GenomicsDB workspaces. Configure `per_contig_dirs` and run:

```bash
snakemake -s Snakefile \
  --configfile config.joint_finalize.yaml \
  --cores 8 \
  --printshellcmds --rerun-incomplete
```

`joint_finalize` validates that every reference contig has exactly one per-contig VCF, concatenates in reference order, then runs the existing SNP/INDEL split, ID assignment, hard filtering, indexing, and SNP count.

## Polyploid genotype-state control

v1.0.0 explicitly keeps GATK's conservative defaults:

```text
--max-alternate-alleles 6
--max-genotype-count 1024
```

Warnings such as `too many genotypes ... current limit 1024` can occur for hexaploid multiallelic sites. v1.0.0 does not raise this limit because genotype-state growth is CPU- and memory-intensive.

## RAD/DArT duplicate handling

The existing behavior is retained:

```text
MarkDuplicates --REMOVE_DUPLICATES false --REMOVE_SEQUENCING_DUPLICATES false
HaplotypeCaller --disable-read-filter NotDuplicateReadFilter
```

## Resume

Use `--rerun-incomplete`. A rerun of `GenomicsDBImport` recreates only that shard's workspace, leaving completed shards untouched.

```bash
snakemake -s Snakefile --configfile CONFIG --cores 30 \
  --resources joint_jobs=1 mem_mb=120000 \
  --rerun-incomplete -n -p
```

For details, see the Traditional Chinese README and `CHANGELOG.md`.
