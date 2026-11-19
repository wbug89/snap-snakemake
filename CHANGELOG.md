# Changelog

## v1.0.3

### Fixed

- Fixed the Full-stage execution bug where the newly added `preflight` rule could be treated as Snakemake's implicit default target after the v1.0.2 release.
- Marked `rule all` with `default_target: True`.
- `scripts/run_snap.sh` now explicitly requests target `all` for both dry-run and production runs.
- Corrected run-mode success reporting so `SNaP SUCCESS` is only reported after the explicit `all` target exits successfully.

### Added

- Added the `minimal` resource profile for systems with about 8 GiB of visible memory.
- Added profile-aware RAM validation and `PASS WITH WARNINGS` when the selected profile exceeds detected memory.
- Added WSL2-specific memory guidance and profile recommendations.

### Unchanged

- Core SNP-calling methodology and filtering remain unchanged: Trimmomatic → BWA-MEM → samtools → MarkDuplicates → merged BAM → HaplotypeCaller GVCF → GenomicsDBImport → GenotypeGVCFs → SNP/INDEL selection → `CHROM_POS` IDs → hard filtering.
- RAD/DArT duplicate retention, configurable ploidy, FASTQ format support, reference interval control, balanced sharding, Telegram notifications, and run-level logging remain unchanged.

## v1.0.2

### Release-blocking corrections

- Marked `rule all` with `default_target: True` so adding `preflight` cannot silently change the default workflow target.
- `run_snap.sh` now explicitly requests `preflight` for validation and `all` for both dry-run and real execution.
- Dry-run now builds the full `all` DAG while retaining concise Job stats output.
- Added a `minimal` profile for ~8 GiB RAM systems.
- Preflight now reports the selected profile and memory budget, warns when the profile exceeds detected RAM, and reports `PASS WITH WARNINGS` in that case.
- Added WSL2-specific memory guidance.
- Run-mode success is only reported after Snakemake exits successfully for explicit target `all`.

### Portability, resource safety, and usability

- Added raw FASTQ support for `.fastq.gz`, `.fq.gz`, `.fastq`, and `.fq`.
- Added strict PE mate validation and exclusion of trimmed/paired/unpaired intermediate FASTQs from raw discovery.
- Kept the default biological-sample parser as `'^((?:SM)[0-9]{3})'` and added preflight validation/preview.
- Added reference prerequisite validation with stage-aware BWA index requirements.
- Added `reference_intervals` modes: `all`, `regex`, `min_length`, `list`, and `file`.
- Added deterministic `contig` and length-balanced joint sharding.
- Added `metadata/reference_sequences.tsv`, selected interval metadata, and a concise preflight report.
- Added per-rule memory resources for heavy steps and global `mem_mb`/`joint_jobs` scheduling scopes.
- Java `-Xmx` is now derived from the final rule memory request, allowing workflow profiles to safely scale GATK jobs.
- Added `minimal`, `low-memory`, `workstation`, and `server` workflow profiles.
- Added Linux/macOS/Windows-via-WSL2 documentation and Windows `Zone.Identifier` ignore rule.
- Added the `preflight` target and concise dry-run console filtering while preserving full dry-run output in `run.log`.
- Added `scripts/run_snap.sh` with timestamped run-level logs, config snapshots, Git/version metadata, Conda environment snapshots, and exit status history.
- Stabilized the Conda environment name as `snap-snakemake`.
- Expanded lock/recovery documentation.
- Core RAD/DArT duplicate handling, configurable ploidy, GVCF/GenomicsDB calling, hard filtering, and `CHROM_POS` SNP IDs remain unchanged.

## v1.0.1

- Expanded user documentation and example configuration coverage.
- Documented reference indexing, stage-specific execution, Telegram setup, and `sample_group_regex` customization.

## v1.0.0

- First formal SNaP release.
