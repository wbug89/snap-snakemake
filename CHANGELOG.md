# Changelog

## v1.0.2 — Portability, resource safety, and usability

- Added raw FASTQ support for `.fastq.gz`, `.fq.gz`, `.fastq`, and `.fq`.
- Added strict PE mate validation and exclusion of trimmed/paired/unpaired intermediate FASTQs from raw discovery.
- Kept the default biological-sample parser as `'^((?:SM)[0-9]{3})'` and added preflight validation/preview.
- Added reference prerequisite validation with stage-aware BWA index requirements.
- Added `reference_intervals` modes: `all`, `regex`, `min_length`, `list`, and `file`.
- Added deterministic `contig` and length-balanced joint sharding.
- Added `metadata/reference_sequences.tsv`, selected interval metadata, and a concise preflight report.
- Added per-rule memory resources for heavy steps and global `mem_mb`/`joint_jobs` scheduling scopes.
- Java `-Xmx` is now derived from the final rule memory request, allowing workflow profiles to safely scale GATK jobs.
- Added `low-memory`, `workstation`, and `server` workflow profiles.
- Added Linux/macOS/Windows-via-WSL2 documentation and Windows `Zone.Identifier` ignore rule.
- Added the `preflight` target and concise dry-run guidance (`--dry-run --quiet`).
- Added `scripts/run_snap.sh` with timestamped run-level logs, config snapshots, Git/version metadata, Conda environment snapshots, and exit status history.
- Stabilized the Conda environment name as `snap-snakemake`.
- Expanded lock/recovery documentation.
- Core RAD/DArT duplicate handling, configurable ploidy, GVCF/GenomicsDB calling, hard filtering, and `CHROM_POS` SNP IDs remain unchanged.

## v1.0.1

- Expanded user documentation and example configuration coverage.
- Documented reference indexing, stage-specific execution, Telegram setup, and `sample_group_regex` customization.

## v1.0.0

- First formal SNaP release.
