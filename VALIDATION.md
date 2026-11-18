# SNaP v1.0.2 validation

Validation performed for the v1.0.2 portability/resource-safety release.

## Passed static checks

- `VERSION` is `1.0.2`.
- Conda environment name is the stable `snap-snakemake`.
- `config.yaml`, every `configs/*.example.yaml`, `environment.yaml`, and all workflow profiles parse as valid YAML.
- All example configs expose the same configurable key tree as `config.yaml`.
- `workflow_inputs.py` compiles successfully with Python.
- `scripts/run_snap.sh` passes `bash -n` and avoids Bash-4-only syntax so it is compatible with the Bash 3.2 shipped by macOS.
- Run-wrapper smoke test with a stub Snakemake executable successfully created a timestamped run-history directory, captured the command/config/environment/Git metadata, streamed a run log, and retained the exit status.
- Default `sample_group_regex` is consistently `^((?:SM)[0-9]{3})`.
- Reference interval helper tests passed for `all`, `regex`, `min_length`, and `list`; `file` uses the same validated selection path.
- Deterministic `contig` and length-balanced joint sharding helper tests passed.
- `CombineGVCFs` is absent from executable workflow code; joint genotyping remains GenomicsDBImport -> GenotypeGVCFs.
- Final SNP IDs are still assigned with `bcftools annotate --set-id '%CHROM_%POS'` before hard filtering.
- `docs/GIT_ON_MAC.md` is not included.
- `.gitignore` includes macOS/Windows metadata exclusions, including `*:Zone.Identifier`.

## Implemented v1.0.2 checks/features

- FASTQ discovery supports `.fastq.gz`, `.fq.gz`, `.fastq`, and `.fq`.
- PE discovery validates R1/R2 pairs and rejects ambiguous/missing pairs.
- FASTQ names containing workflow-intermediate tokens (`trim`, `trimmed`, `paired`, `unpaired`) are excluded from raw-input discovery.
- Reference preflight checks FASTA, `.fai`, and GATK `.dict`; BWA indexes are required only for `preprocess` and `full`.
- Reference sequence selection supports `all`, `regex`, `min_length`, `list`, and `file`.
- Joint sharding supports one-sequence-per-shard (`contig`) and deterministic length-balanced shards (`balanced`).
- Preflight writes `metadata/preflight_summary.txt`, `metadata/reference_sequences.tsv`, `metadata/selected_intervals.list`, and `metadata/fastq_units.tsv`.
- Heavy rules declare `mem_mb`; `mem_mb` and `joint_jobs` are global scheduling resources. Java heap is derived from the final per-rule memory allocation with headroom for native/off-heap memory.
- Low-memory, workstation, and server workflow profiles are included.
- The recommended wrapper uses concise `--dry-run --quiet` output and does not enable `--printshellcmds` unless the user explicitly passes it.

## Limitation of this build-time validation

The artifact-build environment does not contain the Snakemake executable or the bioinformatics toolchain (BWA, samtools, bcftools, GATK, Trimmomatic). Therefore a true `snakemake --lint` / DAG dry-run and end-to-end biological-data run could not be executed here.

Before production use, run on a machine with the SNaP Conda environment:

```bash
conda activate snap-snakemake
bash scripts/run_snap.sh --config <your-config.yaml> --profile low-memory --mode preflight
bash scripts/run_snap.sh --config <your-config.yaml> --profile low-memory --mode dry-run
```

Then choose `low-memory`, `workstation`, or `server` for the production run according to the host.
