# SNaP v1.0.4-rc1 validation

This file records static and wrapper-level validation for the v1.0.4-rc1 release candidate. Real Snakemake DAG regression tests are also bundled for execution inside the `snap-snakemake` Conda environment.

## Release-blocking corrections

Validated statically:

- `rule all` still contains `default_target: True`.
- wrapper preflight target is separated as `-- preflight`.
- wrapper dry-run target is separated as `-- all`.
- wrapper production target is separated as `-- all`.
- wrapper `EXTRA_ARGS` remain before the target delimiter, so options such as `--keep-going` are not misinterpreted as targets.
- `--config` internal metadata is placed before the delimiter.

The delimiter is required because Snakemake's `--config` accepts multiple values and otherwise can consume a following target as a malformed config item.

## SE / PE wildcard regression

`UNIT_REGEX` is constructed from the retained `local_units` with `re.escape()` and applied as a global `unit` wildcard constraint. This prevents the SE output pattern from matching PE outputs by treating `_R1_paired` / `_R2_paired` as part of `{unit}`.

The output namespace itself is unchanged.

## Workflow config isolation

The Snakefile no longer contains `configfile: "config.yaml"`. Workflow configuration must be supplied explicitly with `--configfile` (or through `scripts/run_snap.sh --config`). This prevents implicit merging of the repository-root config with a user-selected config.

The workflow profiles remain execution/resource-only profiles; none contains a `configfile:` entry.

## Profile metadata

The wrapper reads `resources.mem_mb` from the selected `profiles/<name>/profile.yaml` and passes:

- `snap_resource_profile=<name>`
- `snap_profile_mem_mb=<MiB>`

to Snakemake before the target delimiter. Preflight therefore reports the selected profile and memory budget without a duplicated hard-coded profile table.

## Bundled regression dry-runs

Run:

```bash
bash tests/regression/run_regression_tests.sh
```

Expected coverage:

1. PE one unit: one `trim_fastq_pe` job and no scheduled `trim_fastq_se` job.
2. SE one unit: one `trim_fastq_se` job and no scheduled `trim_fastq_pe` job.
3. PE two technical units / one biological sample: two trim/mapping unit jobs but one merge/HaplotypeCaller biological-sample job.
4. No `InputFunctionException in rule trim_fastq_se` from PE paired-output suffixes.
5. Dry-run reaches `Job stats:` and includes target `all`.
6. Preflight metadata contains the selected profile and memory budget.

## Build-environment static checks

The release build checks:

- `bash -n scripts/run_snap.sh`;
- `bash -n tests/regression/run_regression_tests.sh`;
- `python -m py_compile workflow_inputs.py`;
- YAML parsing for `config.yaml`, `environment.yaml`, all example configs, and all four resource profiles;
- absence of `configfile:` in workflow profiles;
- explicit `-- preflight` / `-- all` wrapper targets;
- retained-unit wildcard constraint;
- `default_target: True`;
- final SNP ID convention `%CHROM_%POS`;
- duplicate-retention settings.

## Core workflow invariants

No intended changes were made to:

- Trimmomatic → BWA-MEM → samtools → MarkDuplicates;
- RAD/DArT duplicate retention behavior;
- technical-unit BAM merge;
- configurable ploidy;
- HaplotypeCaller GVCF calling;
- GenomicsDBImport / GenotypeGVCFs joint genotyping;
- `max_alternate_alleles` / `max_genotype_count`;
- final SNP ID convention `CHROM_POS`;
- FASTQ extension support and intermediate-FASTQ exclusion;
- reference interval selection;
- balanced sharding;
- Telegram fallback behavior;
- timestamped run-level logging.

## External validation still required before release

On the clean Snakemake 9.26.1 test server, rerun:

```bash
bash scripts/run_snap.sh --config configs/test.full.yaml --profile minimal --mode preflight
bash scripts/run_snap.sh --config configs/test.full.yaml --profile minimal --mode dry-run
```

For PE, the dry-run must show `trim_fastq_pe` and must not schedule `trim_fastq_se`. After the dry-run passes, proceed to a real Full-stage run and verify Trimmomatic → BWA → BAM processing → HaplotypeCaller → GenomicsDB → final VCF.
