# SNaP regression checks

## 1. Wrapper target parsing (no bioinformatics tools required)

```bash
bash tests/regression/test_wrapper_targets.sh
```

This uses a fake `snakemake` executable to verify the wrapper's argument construction:

- preflight ends with `-- preflight`;
- dry-run ends with `-- all`;
- run ends with `-- all`;
- selected profile metadata is supplied before the delimiter;
- extra Snakemake options such as `--keep-going` stay before `-- all`.

## 2. Real Snakemake DAG regression dry-runs

Run inside the `snap-snakemake` Conda environment:

```bash
bash tests/regression/run_regression_tests.sh
```

The script creates temporary synthetic FASTQ/reference fixtures and checks:

- PE one-unit DAG: `trim_fastq_pe = 1`, no scheduled `trim_fastq_se` job.
- SE one-unit DAG: `trim_fastq_se = 1`, no scheduled `trim_fastq_pe` job.
- PE technical-unit case: two raw units under one biological sample schedule two PE trim/mapping/MarkDuplicates jobs but one merge and one HaplotypeCaller job.
- wrapper target delimiter: `command.txt` contains `-- preflight` / `-- all`.
- preflight profile metadata: `Selected profile: minimal` and a memory budget are present.
- dry-run reaches the complete DAG summary (`Job stats:` and `all = 1`).
- PE paired-output suffixes do not reproduce the `trim_fastq_se` `KeyError` collision.

These are DAG/regression checks only; they do not execute Trimmomatic, BWA, GATK, samtools, or bcftools.

### Failure behavior

Each case prints its start and PASS/FAIL status. If a regression check fails, the script prints the failed assertion and captured output, then retains the temporary test directory (for example `/tmp/snap-regression.XXXXXX`) for inspection. Successful runs remove the temporary fixtures automatically.
