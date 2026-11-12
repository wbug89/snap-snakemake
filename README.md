# SNaP v1.0.1

### Snakemake for NGS and Polyploids

**SNaP makes SNP calling as simple as a snap.**

SNaP is a reproducible Snakemake workflow for SNP discovery from RAD-seq-style NGS datasets, including ddRAD-seq and DArT-seq. It supports SE/PE reads, configurable ploidy, GVCF-based variant calling, per-contig GenomicsDB joint genotyping, and distributed joint calling across multiple Linux servers.

Repository: `snap-snakemake`

## Workflow stages

| Stage | Main input | Main work | Main output |
|---|---|---|---|
| `preprocess` | FASTQ | Trimmomatic → BWA-MEM → sort/index → MarkDuplicates | `marked/*.marked.bam` |
| `gvcf` | marked BAMs | merge units by biological sample → HaplotypeCaller | `gvcf_file/*.g.vcf.gz` |
| `joint` | complete cohort gVCFs | per-contig GenomicsDBImport → GenotypeGVCFs → concat/filter | final cohort VCF |
| `joint_shard` | complete cohort gVCFs + contig assignment | process only contigs assigned to one server | `results/per_contig/*.vcf.gz` |
| `joint_finalize` | all per-contig VCFs | concat in reference order → SNP/INDEL split → filtering | final cohort VCF |
| `full` | FASTQ | single-host end-to-end execution | final cohort VCF |

For distributed work, gVCF calling may be divided by **sample**, but joint calling is divided by **reference contig**, not by sample. Each `joint_shard` server must see the complete cohort of single-sample gVCFs.

---

## 1. Installation

Create the Conda environment:

```bash
conda env create -f environment.yaml
conda activate snap-snakemake-v1.0.1
```

The workflow expects the configured executables for Java, BWA, samtools, bcftools, GATK, and Trimmomatic to be available.

---

## 2. Prepare and index the reference genome

Before running SNaP, index the same reference FASTA on every server that will use it.

Assume:

```bash
REF=/path/to/reference.fa
```

Create the BWA index (required by `preprocess` and `full`):

```bash
bwa index "$REF"
```

Create the FASTA index used by samtools/GATK:

```bash
samtools faidx "$REF"
```

Create the GATK sequence dictionary:

```bash
gatk CreateSequenceDictionary \
  -R "$REF" \
  -O "${REF%.*}.dict"
```

For `reference.fa`, the expected files include:

```text
reference.fa
reference.fa.fai
reference.dict
reference.fa.amb
reference.fa.ann
reference.fa.bwt
reference.fa.pac
reference.fa.sa
```

All servers must use the same reference sequence. When distributing work, verify it if needed:

```bash
sha256sum "$REF"
```

---

## 3. Sample naming and `sample_group_regex`

The example configuration uses biological sample IDs `SM001`–`SM999`.

Default example:

```yaml
sample_group_regex: '^((?:SM)[0-9]{3})'
```

Examples:

```text
SM001_R1.fastq.gz                  → SM001
SM001_R2.fastq.gz                  → SM001
SM001_laneA_R1.fastq.gz            → SM001
SM125_sequencingBatch2_R1.fastq.gz → SM125
```

The first capture group is used as the biological sample name. In the default expression, `((?:SM)[0-9]{3})` captures IDs such as `SM001`, `SM125`, and `SM999` from the beginning of the raw-unit filename.

### Adapting `sample_group_regex` to your own sample names

Set `sample_group_regex` so that **capture group 1** contains exactly the biological sample ID that should be shared across sequencing lanes, batches, or read files. The regex is applied from the beginning of each raw-unit filename, so the leading `^` should normally be kept. Information such as lane, sequencing batch, `_R1`, `_R2`, or other suffixes should remain outside capture group 1.

The default example is:

```yaml
sample_group_regex: '^((?:SM)[0-9]{3})'
```

This means:

- `^` = start matching at the beginning of the filename/unit name.
- `((?:SM)[0-9]{3})` = capture the biological sample ID as group 1.
- `(?:SM)` = match the literal prefix `SM` without creating an extra capture group.
- `[0-9]{3}` = match exactly three digits.

Examples of adapting the pattern:

```yaml
# SM001, SM125, SM999
sample_group_regex: '^((?:SM)[0-9]{3})'

# IB001, IB125, IB999
sample_group_regex: '^((?:IB)[0-9]{3})'

# Either IBC1-IBC9 or TN followed by digits, e.g. IBC1 and TN71
sample_group_regex: '^((?:IBC[0-9]+|TN[0-9]+))'

# ABC_001, ABC_125, XYZ_999
sample_group_regex: '^([A-Z]{3}_[0-9]{3})'
```

For example, with:

```yaml
sample_group_regex: '^((?:SM)[0-9]{3})'
```

all of these raw units are grouped as biological sample `SM001`:

```text
SM001_R1.fastq.gz
SM001_R2.fastq.gz
SM001_lane1_R1.fastq.gz
SM001_batchB_R1.fastq.gz
```

The important rule is that **the first capture group must return only the biological sample ID**. Do not include lane, batch, or read suffixes inside that group. For example, if `SM001_lane1_R1.fastq.gz` and `SM001_lane2_R1.fastq.gz` belong to the same biological sample, both must resolve to `SM001`, not to `SM001_lane1` and `SM001_lane2`.

Keep the outer capturing parentheses `(...)`. Non-capturing groups such as `(?:SM)` may be used inside the sample pattern without changing which value is returned as group 1. If the regex does not capture a biological sample ID, SNaP stops with an error rather than silently assigning an incorrect sample.

---

## 4. Configuration model

`config.yaml` contains the complete parameter set. Every file under `configs/*.example.yaml` also exposes all user-adjustable parameters, while providing stage-specific example paths.

Recommended workflow:

```bash
cp configs/config.preprocess.example.yaml config.preprocess.yaml
```

Edit the copied local config. Local run configs are ignored by Git.

### Main top-level parameters

| Key | Meaning |
|---|---|
| `workflow_stage` | `preprocess`, `gvcf`, `joint`, `joint_shard`, `joint_finalize`, or `full` |
| `read_layout` | `SE` or `PE`; used by `preprocess/full` |
| `batch_id` | unique ID for an independently preprocessed sequencing batch |
| `ploidy` | HaplotypeCaller ploidy, e.g. `2` or `6` |
| `project_dir` | output root for the current run |
| `fastq_dir` | FASTQ directory for `preprocess/full` |
| `reference` | reference FASTA |
| `sample_group_regex` | regex capturing biological sample ID in group 1 |
| `exclude_units` | exact raw sequencing units to skip |
| `include_biological_samples` | optional biological-sample allowlist |
| `exclude_biological_samples` | optional biological-sample blocklist |
| `include_biological_samples_file` | optional one-ID-per-line include file |
| `exclude_biological_samples_file` | optional one-ID-per-line exclude file |
| `marked_dirs` | marked BAM directories for `gvcf` |
| `sample_assignment_file` | optional sample-to-server TSV for distributed `gvcf` |
| `calling_partition` | partition selected from `sample_assignment_file` |
| `gvcf_dirs` | complete-cohort gVCF directories for `joint/joint_shard` |
| `joint_contig_assignment_file` | contig-to-server TSV for `joint_shard` |
| `joint_partition` | partition selected from the contig assignment TSV |
| `per_contig_dirs` | completed per-contig VCF directories for `joint_finalize` |

### Tool paths

```yaml
tools:
  java: java
  bwa: bwa
  samtools: samtools
  bcftools: bcftools
  gatk: gatk
  trimmomatic: trimmomatic
```

Replace any value with an absolute executable path if required.

### Trimmomatic

```yaml
trimmomatic:
  phred: -phred33
  threads: 8
  leading: 3
  trailing: 3
  slidingwindow: '4:30'
  minlen: 36
  builtin_adapter_se: NexteraPE-PE.fa
  builtin_adapter_pe: NexteraPE-PE.fa
  user_adapter_se: null
  user_adapter_pe: null
  illuminaclip_se: '2:30:10'
  illuminaclip_pe: '2:30:10'
```

SNaP first searches adapters packaged with the active Trimmomatic installation, then uses `user_adapter_*` as fallback if configured.

### CPU / Java / GenomicsDB resources

```yaml
resources:
  bwa_threads: 12
  samtools_view_threads: 4
  samtools_sort_threads: 8
  haplotypecaller_threads: 2

java_options:
  markduplicates: '-Xmx20g'
  haplotypecaller: '-Xmx20g'
  genomicsdbimport: '-Xmx16g -XX:ParallelGCThreads=2'
  genotypegvcfs: '-Xmx32g -XX:ParallelGCThreads=2'

genomicsdb:
  batch_size: 30
  reader_threads: 1
  import_mem_mb: 24000
  genotype_mem_mb: 40000
```

`import_mem_mb` and `genotype_mem_mb` are **per-job Snakemake resource requests**. They do not limit total machine-wide usage by themselves. For joint calling, also set a global runtime ceiling, for example:

```bash
--resources joint_jobs=1 mem_mb=120000
```

### Calling/filtering

```yaml
calling:
  minimum_mapping_quality: 30
  max_alternate_alleles: 6
  max_genotype_count: 1024
  hard_filter_expression: 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0'
```

The genotype-state limits are intentionally conservative for polyploid multiallelic sites.

---

## 5. Stage: `preprocess`

Purpose:

```text
FASTQ
→ Trimmomatic
→ BWA-MEM
→ samtools sort/index
→ GATK MarkDuplicates
→ marked BAM
```

For RAD/DArT data, duplicate-flagged reads are retained:

```text
--REMOVE_DUPLICATES false
--REMOVE_SEQUENCING_DUPLICATES false
```

Prepare config:

```bash
cp configs/config.preprocess.example.yaml config.preprocess.yaml
```

Edit at minimum:

```yaml
workflow_stage: preprocess
read_layout: PE
batch_id: batch1
project_dir: /path/to/preprocess_output
fastq_dir: /path/to/rawData
reference: /path/to/reference.fa
ploidy: 6
sample_group_regex: '^((?:SM)[0-9]{3})'
```

Dry run:

```bash
snakemake \
  -s Snakefile \
  --configfile config.preprocess.yaml \
  --cores 30 \
  -n -p
```

Run:

```bash
snakemake \
  -s Snakefile \
  --configfile config.preprocess.yaml \
  --cores 30 \
  --printshellcmds \
  --rerun-incomplete
```

---

## 6. Stage: `gvcf`

This stage discovers marked BAMs from `marked_dirs`, groups all units/batches belonging to the same biological sample, merges them, then runs HaplotypeCaller once per sample.

HaplotypeCaller retains duplicate-flagged RAD/DArT reads with:

```text
--disable-read-filter NotDuplicateReadFilter
```

Prepare config:

```bash
cp configs/config.gvcf.example.yaml config.gvcf.yaml
```

Typical settings:

```yaml
workflow_stage: gvcf
project_dir: /path/to/gvcf_output
reference: /path/to/reference.fa
ploidy: 6
marked_dirs:
  - /path/to/batch01/marked
  - /path/to/batch02/marked
```

### Optional sample distribution across servers

Create:

```text
biological_sample    server
SM001                PC1
SM002                PC1
SM003                PC3
SM004                PC3
```

Then configure PC1:

```yaml
sample_assignment_file: sample_assignment.tsv
calling_partition: PC1
```

and PC3:

```yaml
sample_assignment_file: sample_assignment.tsv
calling_partition: PC3
```

Run:

```bash
snakemake \
  -s Snakefile \
  --configfile config.gvcf.yaml \
  --cores 30 \
  --printshellcmds \
  --rerun-incomplete
```

---

## 7. Stage: `joint` — single-server joint calling

Use this when one server will process all reference contigs.

Every `gvcf_dirs` entry may contain part of the cohort, but together they must provide exactly one indexed gVCF for every biological sample.

```yaml
workflow_stage: joint
project_dir: /path/to/joint_output
reference: /path/to/reference.fa
gvcf_dirs:
  - /path/to/gvcf_PC1/gvcf_file
  - /path/to/gvcf_PC3/gvcf_file
```

Dry run:

```bash
snakemake \
  -s Snakefile \
  --configfile config.joint.yaml \
  --cores 30 \
  --resources joint_jobs=1 mem_mb=120000 \
  -n -p
```

Run:

```bash
snakemake \
  -s Snakefile \
  --configfile config.joint.yaml \
  --cores 30 \
  --resources joint_jobs=1 mem_mb=120000 \
  --printshellcmds \
  --rerun-incomplete
```

Increase `joint_jobs` only after monitoring actual RAM usage.

---

## 8. Stage: `joint_shard` — distributed joint calling

Joint calling is split by reference contig, not by sample. Every shard server must see the **complete cohort gVCFs**.

Create a local assignment file from `joint_contig_assignment.example.tsv`:

```text
contig    server
LG1       PC1
LG2       PC1
LG3       PC3
...
```

Each reference contig must appear exactly once.

PC1 config:

```yaml
workflow_stage: joint_shard
project_dir: /path/to/joint_PC1
reference: /path/to/reference.fa
gvcf_dirs:
  - /path/to/gvcf_PC1/gvcf_file
  - /path/to/gvcf_PC3/gvcf_file
joint_contig_assignment_file: joint_contig_assignment.tsv
joint_partition: PC1
```

PC3 uses the same assignment file but:

```yaml
joint_partition: PC3
```

Run PC1, for example:

```bash
snakemake \
  -s Snakefile \
  --configfile config.joint_shard_PC1.yaml \
  --cores 32 \
  --resources joint_jobs=1 mem_mb=200000 \
  --printshellcmds \
  --rerun-incomplete
```

Run PC3, for example:

```bash
snakemake \
  -s Snakefile \
  --configfile config.joint_shard_PC3.yaml \
  --cores 30 \
  --resources joint_jobs=1 mem_mb=120000 \
  --printshellcmds \
  --rerun-incomplete
```

Each server produces:

```text
results/per_contig/000001.vcf.gz
results/per_contig/000001.vcf.gz.tbi
...
```

`metadata/joint_contigs.tsv` records the actual shard-to-contig mapping.

---

## 9. Stage: `joint_finalize`

After all `joint_shard` runs complete, collect the per-contig VCFs and `.tbi` files on one server. GenomicsDB workspaces do not need to be transferred.

Example:

```yaml
workflow_stage: joint_finalize
project_dir: /path/to/joint_final
reference: /path/to/reference.fa
per_contig_dirs:
  - /path/to/joint_PC1/results/per_contig
  - /path/to/joint_PC3/results/per_contig
```

Run:

```bash
snakemake \
  -s Snakefile \
  --configfile config.joint_finalize.yaml \
  --cores 8 \
  --printshellcmds \
  --rerun-incomplete
```

The workflow validates complete/non-duplicated contig coverage, concatenates in reference `.fai` order, then performs SNP/INDEL selection and filtering.

---

## 10. Stage: `full`

`full` executes the complete workflow on one host from FASTQ to final filtered variants.

```bash
cp configs/config.full.example.yaml config.full.yaml
```

Edit `fastq_dir`, `reference`, `project_dir`, `read_layout`, `batch_id`, `ploidy`, and sample naming as needed.

Run:

```bash
snakemake \
  -s Snakefile \
  --configfile config.full.yaml \
  --cores 30 \
  --resources joint_jobs=1 mem_mb=120000 \
  --printshellcmds \
  --rerun-incomplete
```

---

## 11. Telegram notifications

All example configs include:

```yaml
telegram:
  enabled: false
  bot_token: ''
  chat_id: ''
```

### Recommended method: environment variables

Keep secrets out of Git-tracked YAML files:

```bash
export TG_BOT_TOKEN='YOUR_BOT_TOKEN'
export TG_CHAT_ID='YOUR_CHAT_ID'
```

Then enable notifications in the local config:

```yaml
telegram:
  enabled: true
  bot_token: ''
  chat_id: ''
```

SNaP uses config values first; empty values fall back to `TG_BOT_TOKEN` and `TG_CHAT_ID`.

Alternatively, credentials can be placed directly in a **local, Git-ignored config**:

```yaml
telegram:
  enabled: true
  bot_token: 'YOUR_BOT_TOKEN'
  chat_id: 'YOUR_CHAT_ID'
```

Do not commit real credentials.

If Telegram is enabled but credentials are missing, SNaP prints a warning and continues without notifications.

---

## 12. Final SNP ID convention

After joint genotyping and SNP selection, SNaP sets the VCF `ID` field with:

```bash
bcftools annotate --set-id '%CHROM_%POS'
```

Example IDs:

```text
LG1_123456
LG5_9876543
```

The final:

```text
final_combined_filtered_snps.vcf.gz
```

retains these `CHROM_POS` IDs.

---

## 13. Resume and troubleshooting

Dry-run before large jobs:

```bash
snakemake -s Snakefile --configfile CONFIG.yaml --cores 30 -n -p
```

Resume incomplete jobs:

```bash
snakemake -s Snakefile --configfile CONFIG.yaml --cores 30 --rerun-incomplete
```

For joint calling, always pair this with an explicit whole-machine resource ceiling appropriate for the server:

```bash
--resources joint_jobs=1 mem_mb=120000
```

For long jobs, run Snakemake inside a persistent terminal such as `tmux`.

## Version

`v1.0.1` is a documentation/configuration maintenance release. Core SNP-calling logic is unchanged from v1.0.0.
