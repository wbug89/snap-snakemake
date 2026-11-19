# SNaP v1.0.4-rc1 Release Candidate Notes

SNaP v1.0.4-rc1 is a release candidate for the v1.0.4 patch release, prepared after real-world testing of v1.0.3 on a clean Linux installation with Snakemake 9.26.1. It fixes two release-blocking DAG/CLI regressions and removes an implicit workflow-config merge. The core SNP-calling methodology is unchanged.

## Release-blocking fixes

### 1. Explicit target delimiter

Snakemake's `--config` option accepts multiple `name=value` arguments. In v1.0.3, a following target such as `preflight` could therefore be parsed as a malformed config item. The wrapper now places every workflow target after a literal `--`:

```text
preflight -> ... --config ... -- preflight
dry-run   -> ... --dry-run ... -- all
run       -> ... --rerun-incomplete ... -- all
```

`rule all` remains `default_target: True`, so explicit wrapper targets and the Snakefile default target provide independent protection.

### 2. SE / PE wildcard separation

The generic SE output pattern could previously interpret a PE output such as `<unit>_R1_paired.trim.fastq.gz` as an SE `{unit}` value. This caused `trim_fastq_se` to enter a PE DAG and fail in the SE input function with a `KeyError`.

SNaP v1.0.4-rc1 constrains `{unit}` to the actual retained unit keys discovered for the current run. PE paired/unpaired suffixes therefore cannot become synthetic unit IDs. Output directories and downstream paths remain unchanged.

## Config isolation

The Snakefile no longer auto-loads the repository-root `config.yaml`. Users must explicitly supply a workflow config with either:

```bash
bash scripts/run_snap.sh --config config.full.yaml ...
```

or direct Snakemake usage:

```bash
snakemake -s Snakefile --configfile config.full.yaml ...
```

This prevents a CLI-selected config from being implicitly extended by unrelated keys from the repository-root `config.yaml`.

## Profile metadata

The wrapper reads the selected workflow profile directly and passes:

- `snap_resource_profile`
- `snap_profile_mem_mb`

into preflight. This keeps `Selected profile` and `Profile memory budget` synchronized with the actual selected profile without maintaining a second hard-coded lookup table.

## Dry-run Job stats compatibility

Snakemake 9.x may suppress the `Job stats:` table when `--quiet` is used. SNaP v1.0.4-rc1 therefore does not pass `--quiet` for Full-stage dry-runs. The complete native dry-run output is saved to the run log, while the terminal stream is filtered to retain the DAG / Job stats summary and important errors without printing every routine per-sample detail.

## Regression checks

`tests/regression/run_regression_tests.sh` provides synthetic dry-run checks for:

- one-unit PE;
- one-unit SE;
- two PE technical units mapping to one biological sample;
- target delimiter behavior;
- profile metadata;
- Job stats generation;
- the `_R1_paired` / `_R2_paired` suffix collision.

Run it inside the SNaP Conda environment:

```bash
bash tests/regression/run_regression_tests.sh
```

The regression harness prints each case as it runs. On failure it reports the failed assertion, prints the captured output, and retains the temporary `/tmp/snap-regression.*` directory for inspection. Successful runs remove temporary fixtures automatically.

## Unchanged core analysis

```text
Trimmomatic
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

Also unchanged:

- duplicate reads are marked but not removed;
- HaplotypeCaller retains duplicate-flagged reads;
- ploidy remains configurable;
- GenomicsDB joint calling and genotype-complexity controls are unchanged;
- FASTQ support remains `.fastq.gz`, `.fq.gz`, `.fastq`, and `.fq`;
- reference interval selection and balanced sharding are unchanged;
- Telegram notifications and timestamped run-level logging are unchanged.

## Upgrade recommendation

Testing users of v1.0.3 are encouraged to evaluate v1.0.4-rc1 before the final v1.0.4 release. Keep existing v1.0.3 tags/releases available for reproducibility rather than moving or overwriting the old tag.
