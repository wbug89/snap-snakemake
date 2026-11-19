# SNaP v1.0.3 validation

This file records static and wrapper-level validation performed for the v1.0.3 patch release.

## Release-blocking target correction

Validated statically:

- `rule all` contains `default_target: True`.
- wrapper `preflight` mode explicitly targets `preflight`.
- wrapper `dry-run` mode explicitly targets `all`.
- wrapper `run` mode explicitly targets `all`.
- run mode reports `SNaP SUCCESS` only when the explicit target is `all` and Snakemake returns exit code 0.

A wrapper smoke test with a fake Snakemake executable confirmed the actual argument sequences:

```text
preflight -> ... preflight --forcerun preflight --quiet
dry-run  -> ... all --dry-run
run      -> ... all --rerun-incomplete
```

The dry-run smoke test also confirmed that the console retains a `Job stats:` table while per-job details are suppressed, and that the complete native dry-run output remains in `run.log`.

## Resource profiles

All profile YAML files parse successfully. Current global limits are:

| Profile | Cores | Global mem_mb | joint_jobs |
|---|---:|---:|---:|
| minimal | 2 | 6000 | 1 |
| low-memory | 4 | 12000 | 1 |
| workstation | 12 | 48000 | 1 |
| server | 30 | 150000 | 2 |

For an illustrative WSL2-visible memory size of 7.3 GiB:

- `minimal` (6000 MiB) -> expected `PASS`.
- `low-memory` (12000 MiB) -> expected `PASS WITH WARNINGS`.
- larger profiles -> expected `PASS WITH WARNINGS`.

The wrapper passes both the selected profile name and profile memory budget into the preflight configuration as `snap_resource_profile` and `snap_profile_mem_mb`.

## Preflight warning behavior

Static validation confirms that preflight now reports:

- detected CPU threads;
- detected physical/visible memory;
- selected profile;
- profile memory budget;
- WSL2-specific memory guidance when applicable;
- `Validation: PASS WITH WARNINGS` when the selected profile memory budget exceeds detected memory.

## Syntax and configuration checks

Passed in the build environment:

- `bash -n scripts/run_snap.sh`;
- `python -m py_compile workflow_inputs.py`;
- YAML parsing for `config.yaml`, `environment.yaml`, all example configs, and all four resource profiles;
- static checks for `default_target: True`, explicit wrapper targets, WSL2 warning text, and `PASS WITH WARNINGS`;
- wrapper run-history smoke test with timestamped logs and metadata.

## Core workflow invariants

No intended changes were made to:

- Trimmomatic -> BWA-MEM -> samtools -> MarkDuplicates;
- RAD/DArT duplicate retention behavior;
- configurable ploidy;
- HaplotypeCaller GVCF calling;
- GenomicsDBImport / GenotypeGVCFs joint genotyping;
- `max_alternate_alleles` / `max_genotype_count`;
- final SNP ID convention `CHROM_POS`;
- FASTQ extension support;
- reference interval selection;
- balanced sharding;
- Telegram fallback behavior;
- timestamped run-level logging.

## Remaining end-to-end validation

The artifact build environment does not contain Snakemake/GATK/BWA/samtools/bcftools/Trimmomatic, so a real full DAG parse and tool execution could not be performed here. Before publishing the GitHub release, rerun the following on the WSL2 test machine using the 4-sample / 13-sequence test dataset:

1. `minimal` preflight: expect 4 biological samples, 13 selected sequences, 13 contig shards, and `PASS`.
2. `low-memory` preflight on ~7.3 GiB visible RAM: expect `PASS WITH WARNINGS`.
3. full dry-run: expect the complete `all` Job stats including trimming, mapping, BAM processing, HaplotypeCaller, GenomicsDBImport, GenotypeGVCFs, final filtering, and `all`.
4. full run: confirm execution begins with Trimmomatic/BWA rather than stopping after preflight.
5. rerun after all final outputs exist: `Nothing to be done` is valid and should exit 0.
