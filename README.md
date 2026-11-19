# SNaP

**Snakemake for NGS and Polyploids**

> **SNaP makes SNP calling as simple as a snap.**

SNaP is a Snakemake workflow for RAD-seq/ddRAD-seq/DArT-seq SNP discovery. It supports SE/PE data, configurable ploidy, GVCF calling, GenomicsDB joint genotyping, reference interval selection, distributed joint calling, and memory-aware execution profiles.

Current version: **v1.0.3**

## Highlights in v1.0.3

- Fixes the Full-stage default-target bug discovered after the v1.0.2 GitHub release.
- `rule all` is explicitly marked as the Snakemake default target.
- `scripts/run_snap.sh` explicitly targets `all` for dry-run and production runs.
- Corrects success reporting so a successful preflight alone cannot be reported as a completed Full-stage workflow.
- Adds the `minimal` resource profile for systems with about 8 GiB of visible memory.
- Adds profile-aware RAM warnings, `PASS WITH WARNINGS`, and WSL2 memory guidance.
- Core SNP-calling methodology and filtering remain unchanged from v1.0.2.

## Highlights introduced in v1.0.2

- Raw FASTQ: `.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`.
- Robust PE R1/R2 validation and exclusion of trimmed/paired/unpaired intermediate FASTQs.
- Preflight validation of `sample_group_regex`, reference indexes, inputs, assignments, and sharding.
- Reference sequence selection: `all`, `regex`, `min_length`, `list`, `file`.
- Joint sharding: one-sequence-per-shard (`contig`) or deterministic length-balanced shards (`balanced`).
- Per-rule memory declarations and dynamic Java heap sizing.
- Built-in `minimal`, `low-memory`, `workstation`, and `server` workflow profiles.
- Linux, macOS, and Windows via WSL2.
- Concise preflight/dry-run workflow.
- Run-level history/logging through `scripts/run_snap.sh`.
- Final SNP IDs remain `CHROM_POS` (e.g. `LG1_123456`).


## Quick navigation

- [Installation](#1-installation)
- [Prepare the reference genome](#2-prepare-the-reference-genome)
- [FASTQ inputs](#3-fastq-inputs)
- [Biological sample parsing](#4-biological-sample-parsing)
- [Reference interval selection](#5-reference-interval-selection)
- [Joint sharding](#6-joint-sharding)
- [Resource profiles](#7-resource-profiles)
- [Preflight, dry-run, and execution](#8-preflight-dry-run-and-execution)
- [Run history and logs](#9-run-history-and-logs)
- [Workflow stages](#10-stages)
- [Main configuration keys](#11-main-configuration-keys)
- [Telegram](#12-telegram)
- [SNP IDs and outputs](#13-snp-ids-and-outputs)
- [Lock and recovery](#14-lockrecovery)

---

## Quick Start: Full stage on a desktop or laptop

Use this guide to run SNaP from raw FASTQ files to the final SNP/INDEL VCF on a single computer.

The Full stage performs:

```text
FASTQ
→ Trimmomatic
→ BWA-MEM
→ samtools sort/index
→ MarkDuplicates
→ merge technical units
→ HaplotypeCaller GVCF
→ GenomicsDBImport
→ GenotypeGVCFs
→ SNP / INDEL selection
→ CHROM_POS SNP IDs
→ hard filtering
→ final VCF
```

### Step 1 — Enter the SNaP repository and create the environment

For the first installation:

```bash
cd /path/to/snap-snakemake

conda env create -f environment.yaml
conda activate snap-snakemake
```

To update an existing environment:

```bash
conda env update \
  -n snap-snakemake \
  -f environment.yaml \
  --prune

conda activate snap-snakemake
```

Check the main tools:

```bash
snakemake --version
bwa 2>&1 | head
samtools --version
bcftools --version
gatk --version
```

### Step 2 — Prepare the reference genome

```bash
REF=/absolute/path/to/reference.fa
```

Create the required indexes:

```bash
bwa index "$REF"
samtools faidx "$REF"
gatk CreateSequenceDictionary \
  -R "$REF" \
  -O "${REF%.*}.dict"
```

The Full stage requires the FASTA, `.fai`, GATK sequence dictionary, and BWA indexes.

### Step 3 — Prepare FASTQ files

Supported extensions:

```text
.fastq.gz
.fq.gz
.fastq
.fq
```

PE example:

```text
SM001_R1.fastq.gz
SM001_R2.fastq.gz
SM002_R1.fastq.gz
SM002_R2.fastq.gz
```

The default biological-sample parser is:

```yaml
sample_group_regex: '^((?:SM)[0-9]{3})'
```

For example:

```text
SM001_L001_R1.fastq.gz → SM001
SM001_L002_R1.fastq.gz → SM001
```

Modify `sample_group_regex` if your sample naming scheme differs.

### Step 4 — Create a Full-stage config

Copy the complete example:

```bash
cp configs/config.full.example.yaml config.full.yaml
```

Edit it:

```bash
nano config.full.yaml
```

At minimum, review:

```yaml
workflow_stage: full

read_layout: PE
batch_id: batch1
ploidy: 2

project_dir: /absolute/path/to/snap_output
fastq_dir: /absolute/path/to/fastq
reference: /absolute/path/to/reference.fa

sample_group_regex: '^((?:SM)[0-9]{3})'

reference_intervals:
  mode: all

joint_sharding:
  strategy: contig
  shards: null

telegram:
  enabled: false
  bot_token: ''
  chat_id: ''
```

Other advanced options can initially remain at the values supplied in `config.full.example.yaml`.

### Step 4.1 — Optional: Enable Telegram notifications

To receive Telegram notifications when the workflow starts, succeeds, or fails, create a Telegram Bot and obtain its `bot_token` and your target `chat_id`.

You can place them directly in the local config:

```yaml
telegram:
  enabled: true
  bot_token: 'YOUR_BOT_TOKEN'
  chat_id: 'YOUR_CHAT_ID'
```

To disable Telegram:

```yaml
telegram:
  enabled: false
  bot_token: ''
  chat_id: ''
```

Do not commit local configuration files containing real credentials to GitHub.

A safer alternative is to use environment variables:

```bash
export TG_BOT_TOKEN='YOUR_BOT_TOKEN'
export TG_CHAT_ID='YOUR_CHAT_ID'
```

and keep the config as:

```yaml
telegram:
  enabled: true
  bot_token: ''
  chat_id: ''
```

SNaP uses values supplied in the config first and falls back to `TG_BOT_TOKEN` and `TG_CHAT_ID` when the config values are empty. If Telegram is enabled but credentials are unavailable, SNaP prints a warning and skips notifications rather than terminating the workflow.

The next preflight step also verifies whether the Telegram configuration can be resolved.

### Step 5 — Choose a resource profile

| Profile | Recommended hardware | Typical use |
|---|---|---|
| `minimal` | ~8 GiB RAM, 2–4 CPU threads | Small datasets, teaching, low-memory laptops/WSL2 tests |
| `low-memory` | ~16 GiB RAM, 4–8 CPU threads | Ordinary laptop or desktop |
| `workstation` | ~32–64 GiB RAM, 8–16 threads | Research workstation |
| `server` | ≥128 GiB RAM | High-memory server |

If SNaP sees only ~7–8 GiB RAM, use `minimal`. For ~16 GiB, use `low-memory`. Large cohorts, polyploid datasets, or highly multiallelic regions may require substantially more memory than the `minimal` profile can provide.

Avoid starting with:

```bash
snakemake --cores all
```

for a memory-intensive variant-calling workflow.

### Step 6 — Run preflight

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile minimal \
  --mode preflight
```

Fix any reported input, reference, sample-parsing, interval, sharding, resource, or Telegram configuration errors before continuing.

### Step 7 — Run a concise dry-run

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile minimal \
  --mode dry-run
```

This validates the complete `all` DAG without executing the analysis. The Job stats table remains visible; detailed per-sample job output is retained in the run log instead of flooding the terminal.

For debugging only:

```bash
snakemake \
  -s Snakefile \
  --configfile config.full.yaml \
  --workflow-profile profiles/minimal \
  all \
  --dry-run \
  --printshellcmds
```

### Step 8 — Run the Full stage

For a ~8 GiB computer or constrained WSL2 environment:

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile minimal \
  --mode run
```

For a ~16 GiB laptop or desktop:

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile low-memory \
  --mode run
```

For a 32–64 GiB workstation:

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile workstation \
  --mode run
```

For a high-memory server:

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile server \
  --mode run
```

### Step 9 — Inspect run history

Each wrapper invocation creates:

```text
<project_dir>/logs/runs/YYYY-MM-DD_HHMMSS/
```

Inspect the latest run:

```bash
less /absolute/path/to/snap_output/logs/runs/latest/run.log
```

The run directory records the command, config snapshot, environment, Git version, preflight summary, Snakemake output, and exit status.

### Step 10 — Check the final outputs

Primary outputs:

```text
<project_dir>/results/final_combined_filtered_snps.vcf.gz
<project_dir>/results/final_combined_filtered_snps.vcf.gz.tbi
<project_dir>/results/final_combined_filtered_indels.vcf.gz
<project_dir>/results/final_combined_filtered_indels.vcf.gz.tbi
```

SNP IDs follow `CHROM_POS`, for example:

```text
LG1_123456
chr05_987654
```

Inspect the final SNP VCF:

```bash
bcftools view \
  -H \
  /absolute/path/to/snap_output/results/final_combined_filtered_snps.vcf.gz \
  | head
```

### Minimal Full-stage SOP

Once the environment, reference indexes, and `config.full.yaml` are ready:

```bash
conda activate snap-snakemake
```

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile minimal \
  --mode preflight
```

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile minimal \
  --mode dry-run
```

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile minimal \
  --mode run
```

For ~16 GiB RAM, replace `minimal` with `low-memory`. For 32–64 GiB, use `workstation`; for ≥128 GiB, use `server`.

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
profiles/minimal/      # ~8 GiB RAM
profiles/low-memory/   # ~16 GiB RAM
profiles/workstation/  # ~32-64 GiB RAM
profiles/server/       # >=128 GiB RAM
```

The `minimal` profile uses 2 cores, a 6000 MiB global memory budget, and `joint_jobs: 1`. It is intended for small datasets, teaching, smoke tests, and machines where SNaP sees roughly 7–8 GiB RAM. Large cohorts, polyploid data, or highly multiallelic regions may still require more memory.

Heavy rules declare `mem_mb`. SNaP marks `mem_mb` and `joint_jobs` as global scheduling resources so profile-level limits constrain concurrent high-memory jobs. GATK `-Xmx` is calculated from approximately 75% of the final per-rule `mem_mb`, leaving headroom for native/off-heap memory.

Profiles reduce the risk of whole-machine OOM, but they cannot guarantee that a single unusually large job fits into a given machine.

## 8. Preflight, dry-run, and execution

The recommended entry point is the run wrapper. The wrapper explicitly targets `preflight` for validation and explicitly targets `all` for both dry-run and real execution. This prevents Snakemake rule order from accidentally changing the requested workflow target.

Preflight reports the selected resource profile, its global memory budget, detected CPU threads, and detected memory. If the profile memory budget exceeds detected memory, validation becomes `PASS WITH WARNINGS` and SNaP recommends a smaller profile. Under WSL2, the detected memory is the memory visible to WSL2, which may be lower than the RAM installed in Windows.

Preflight:

```bash
bash scripts/run_snap.sh --config config.full.yaml --profile workstation --mode preflight
```

Concise dry-run:

```bash
bash scripts/run_snap.sh --config config.full.yaml --profile workstation --mode dry-run
```

The wrapper explicitly requests target `all` and runs a full Snakemake dry-run. The complete native dry-run output is saved to `run.log`, while the console is filtered to retain the DAG/Job stats summary and important status/error lines without flooding the terminal with per-sample job details. This avoids Snakemake 9.x `--quiet` behavior that can also suppress the Job stats summary. Use `--printshellcmds` only when debugging.

Normal run:

```bash
bash scripts/run_snap.sh --config config.full.yaml --profile workstation --mode run
```

Only use `--printshellcmds` for debugging:

```bash
snakemake -s Snakefile --configfile config.full.yaml \
  --workflow-profile profiles/workstation \
  all --dry-run --printshellcmds
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
