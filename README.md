# SNaP

**Snakemake for NGS and Polyploids**

> **SNaP makes SNP calling as simple as a snap.**

SNaP is a Snakemake workflow for RAD-seq/ddRAD-seq/DArT-seq SNP discovery. It supports SE/PE data, configurable ploidy, GVCF calling, GenomicsDB joint genotyping, reference interval selection, distributed joint calling, and memory-aware execution profiles.

Current version: **v1.0.2**

## Highlights in v1.0.2

- Raw FASTQ: `.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`.
- Robust PE R1/R2 validation and exclusion of trimmed/paired/unpaired intermediate FASTQs.
- Preflight validation of `sample_group_regex`, reference indexes, inputs, assignments, and sharding.
- Reference sequence selection: `all`, `regex`, `min_length`, `list`, `file`.
- Joint sharding: one-sequence-per-shard (`contig`) or deterministic length-balanced shards (`balanced`).
- Per-rule memory declarations and dynamic Java heap sizing.
- Built-in `low-memory`, `workstation`, and `server` workflow profiles.
- Linux, macOS, and Windows via WSL2.
- Concise preflight/dry-run workflow.
- Run-level history/logging through `scripts/run_snap.sh`.
- Final SNP IDs remain `CHROM_POS` (e.g. `LG1_123456`).

## 1. Installation

Create the stable environment name:

```bash
conda env create -f environment.yaml
conda activate snap-snakemake
```

Update an existing environment:

```bash
conda env update -n snap-snakemake -f environment.yaml --prune
```

SNaP versions are tracked by `VERSION` and Git tags, not by the Conda environment name.

Supported platforms: Linux, macOS, and Windows 10/11 through WSL2. Native Windows CMD/PowerShell execution is not a supported bioinformatics runtime. Apple Silicon support depends on availability of the corresponding Bioconda builds for each dependency.

## 2. Prepare the reference genome

```bash
REF=/path/to/reference.fa
bwa index "$REF"
samtools faidx "$REF"
gatk CreateSequenceDictionary -R "$REF" -O /path/to/reference.dict
```

`preprocess/full` require the BWA indexes. Other stages require the FASTA, `.fai`, and GATK `.dict` but do not require BWA mapping indexes.

## 3. FASTQ inputs

Supported extensions:

```text
.fastq.gz  .fq.gz  .fastq  .fq
```

PE filenames must contain `_R1` and `_R2` immediately before the optional `_001` and FASTQ extension, for example:

```text
SM001_R1.fastq.gz / SM001_R2.fastq.gz
SM002_L001_R1_001.fq.gz / SM002_L001_R2_001.fq.gz
```

Files containing workflow-output tokens such as `trim`, `trimmed`, `paired`, or `unpaired` are excluded from raw FASTQ discovery.

## 4. Biological sample parsing

Default:

```yaml
sample_group_regex: '^((?:SM)[0-9]{3})'
```

Capture group 1 must contain only the biological sample ID. Examples:

```yaml
# SM001
sample_group_regex: '^((?:SM)[0-9]{3})'
# IB001
sample_group_regex: '^((?:IB)[0-9]{3})'
# IBC1 / TN71
sample_group_regex: '^((?:IBC[0-9]+|TN[0-9]+))'
# ABC_001
sample_group_regex: '^([A-Z]{3}_[0-9]{3})'
```

A sample with multiple lanes/runs should resolve to the same biological ID so its BAM units can be merged before HaplotypeCaller.

## 5. Reference interval selection

```yaml
reference_intervals:
  mode: all
  regex: null
  min_length: null
  sequences: []
  file: null
```

Modes:

- `all`: all sequences in the FASTA `.fai`.
- `regex`: e.g. `'^chr[0-9]+$'` or `'^LG[0-9]+$'`.
- `min_length`: retain sequences at or above a bp threshold.
- `list`: explicit YAML list of sequence names.
- `file`: one sequence name per line.

SNaP v1.0.2 deliberately does not guess whether a sequence is a chromosome or contig. Use the naming scheme/metadata of the chosen reference.

## 6. Joint sharding

Default compatibility mode:

```yaml
joint_sharding:
  strategy: contig
  shards: null
```

One selected reference sequence becomes one joint shard.

For references with many scaffolds/contigs:

```yaml
joint_sharding:
  strategy: balanced
  shards: 16
```

SNaP uses deterministic longest-first greedy balancing by `.fai` sequence length. `metadata/reference_sequences.tsv` records each sequence, its length, selection status, exclusion reason, and shard.

For distributed `joint_shard`, all sequences inside a balanced shard must be assigned to the same partition in `joint_contig_assignment.tsv`; SNaP fails fast if a shard is split across servers.

## 7. Resource profiles

Profiles are under:

```text
profiles/low-memory/   # ~16 GiB RAM
profiles/workstation/  # ~32-64 GiB RAM
profiles/server/       # >=128 GiB RAM
```

Heavy rules declare `mem_mb`. SNaP marks `mem_mb` and `joint_jobs` as global scheduling resources so profile-level limits constrain concurrent high-memory jobs. GATK `-Xmx` is calculated from approximately 75% of the final per-rule `mem_mb`, leaving headroom for native/off-heap memory.

Profiles reduce the risk of whole-machine OOM, but they cannot guarantee that a single unusually large job fits into a given machine.

## 8. Preflight, dry-run, and execution

The recommended entry point is the run wrapper.

Preflight:

```bash
bash scripts/run_snap.sh --config config.full.yaml --profile workstation --mode preflight
```

Concise dry-run:

```bash
bash scripts/run_snap.sh --config config.full.yaml --profile workstation --mode dry-run
```

The dry-run uses Snakemake `--dry-run --quiet`, avoiding thousands of per-job shell commands.

Normal run:

```bash
bash scripts/run_snap.sh --config config.full.yaml --profile workstation --mode run
```

Only use `--printshellcmds` for debugging:

```bash
snakemake -s Snakefile --configfile config.full.yaml \
  --workflow-profile profiles/workstation \
  --dry-run --printshellcmds
```

## 9. Run history and logs

Each wrapper invocation creates:

```text
<project_dir>/logs/runs/YYYY-MM-DD_HHMMSS/
├── run.log
├── command.txt
├── config.snapshot.yaml
├── preflight_summary.txt
├── environment.txt
├── git_version.txt
└── run_metadata.tsv
```

`logs/runs/latest` points to the newest run when symbolic links are supported.

This records the command, configuration snapshot, Git commit/tag, Conda package list, preflight summary, Snakemake output, final status, and exit code.

## 10. Stages

### preprocess

```text
FASTQ → Trimmomatic → BWA-MEM → samtools sort/index → MarkDuplicates
```

```yaml
workflow_stage: preprocess
read_layout: PE
batch_id: batch1
project_dir: /path/to/preprocess
fastq_dir: /path/to/raw
reference: /path/to/reference.fa
```

```bash
bash scripts/run_snap.sh --config config.preprocess.yaml --profile workstation --mode run
```

RAD/DArT duplicate reads are marked for QC but are not removed.

### gvcf

```text
marked BAM units → merge by biological sample → HaplotypeCaller GVCF
```

```yaml
workflow_stage: gvcf
marked_dirs:
  - /path/to/batch1/marked
  - /path/to/batch2/marked
```

Optional distributed sample calling:

```yaml
sample_assignment_file: sample_assignment.tsv
calling_partition: PC1
```

HaplotypeCaller uses configurable ploidy and disables `NotDuplicateReadFilter` so duplicate-flagged RAD/DArT reads remain available.

### joint

Single-machine joint genotyping over the complete sample cohort:

```yaml
workflow_stage: joint
gvcf_dirs:
  - /path/to/gvcf_PC1
  - /path/to/gvcf_PC3
```

```bash
bash scripts/run_snap.sh --config config.joint.yaml --profile server --mode run
```

### joint_shard

Distribute genomic shards (not sample cohorts) across servers:

```yaml
workflow_stage: joint_shard
joint_contig_assignment_file: joint_contig_assignment.tsv
joint_partition: PC1
```

Every shard server must see the complete single-sample gVCF cohort.

### joint_finalize

Collect all shard VCFs:

```yaml
workflow_stage: joint_finalize
per_contig_dirs:
  - /path/to/PC1/results/per_contig
  - /path/to/PC3/results/per_contig
```

Balanced shards can contain non-adjacent reference sequences, so finalization performs `bcftools concat` followed by `bcftools sort`.

### full

End-to-end FASTQ to final VCF on one machine:

```bash
bash scripts/run_snap.sh --config config.full.yaml --profile low-memory --mode run
```

Choose `workstation` or `server` when appropriate.

## 11. Main configuration keys

Core:

```yaml
workflow_stage: full
read_layout: PE
batch_id: batch1
ploidy: 6
project_dir: /path/to/project
fastq_dir: /path/to/raw
reference: /path/to/reference.fa
sample_group_regex: '^((?:SM)[0-9]{3})'
```

Sample selection:

```yaml
exclude_units: []
include_biological_samples: []
exclude_biological_samples: []
include_biological_samples_file: null
exclude_biological_samples_file: null
```

gVCF distribution:

```yaml
marked_dirs: []
sample_assignment_file: null
calling_partition: null
```

Joint/finalize:

```yaml
gvcf_dirs: []
joint_contig_assignment_file: null
joint_partition: null
per_contig_dirs: []
```

Per-rule defaults:

```yaml
resources:
  trimmomatic_threads: 2
  trimmomatic_mem_mb: 2048
  bwa_threads: 4
  bwa_mem_mb: 4096
  samtools_view_threads: 1
  samtools_sort_threads: 2
  samtools_sort_mem_mb: 4096
  markduplicates_mem_mb: 6144
  merge_bam_mem_mb: 4096
  haplotypecaller_threads: 2
  haplotypecaller_mem_mb: 8192
  bcftools_mem_mb: 2048
```

GenomicsDB:

```yaml
genomicsdb:
  batch_size: 30
  reader_threads: 1
  import_mem_mb: 8192
  genotype_mem_mb: 12288
```

Calling/filtering:

```yaml
calling:
  minimum_mapping_quality: 30
  max_alternate_alleles: 6
  max_genotype_count: 1024
  hard_filter_expression: 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0'
```

## 12. Telegram

Recommended configuration:

```yaml
telegram:
  enabled: true
  bot_token: ''
  chat_id: ''
```

Keep credentials out of Git:

```bash
export TG_BOT_TOKEN='YOUR_BOT_TOKEN'
export TG_CHAT_ID='YOUR_CHAT_ID'
```

Missing Telegram credentials generate a warning and disable notifications; they do not abort the analysis.

## 13. SNP IDs and outputs

Final SNP IDs are assigned with:

```text
%CHROM_%POS
```

Example: `LG1_123456`.

Primary output:

```text
results/final_combined_filtered_snps.vcf.gz
results/final_combined_filtered_snps.vcf.gz.tbi
results/final_combined_filtered_indels.vcf.gz
results/final_combined_filtered_indels.vcf.gz.tbi
```

## 14. Lock/recovery

Confirm no other Snakemake process is running:

```bash
pgrep -af snakemake
```

Remove a stale lock after an interrupted run:

```bash
snakemake -s Snakefile --configfile config.yaml --unlock
```

`--unlock` removes workflow locking metadata; it does not delete analysis outputs. Resume with the wrapper; normal runs use `--rerun-incomplete`.

## 15. Core analysis logic retained in v1.0.2

```text
Trimmomatic
→ BWA-MEM
→ samtools sort/index
→ MarkDuplicates (mark only; do not remove RAD/DArT duplicates)
→ merge technical units by biological sample
→ HaplotypeCaller (configurable ploidy, GVCF, duplicate reads retained)
→ GenomicsDBImport
→ GenotypeGVCFs
→ SNP/INDEL selection
→ CHROM_POS SNP IDs
→ hard filtering
```

v1.0.2 is primarily a portability, resource-safety, input-flexibility, reference-control, and usability release; it does not intentionally change the core SNP-calling model.
