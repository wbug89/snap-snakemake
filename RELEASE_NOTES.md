# SNaP v1.0.2 Release Notes

SNaP v1.0.2 is a portability, resource-safety, input-flexibility, reference-control, and usability release. It intentionally preserves the core SNP-calling model while making the workflow safer on laptops/workstations and more predictable across Linux, macOS, and Windows via WSL2.

## Major additions

1. Four FASTQ extensions: `.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`.
2. Reference sequence selection (`all`, `regex`, `min_length`, `list`, `file`).
3. Deterministic balanced joint sharding for references with many contigs/scaffolds.
4. Stage-aware reference preflight and a user-facing `preflight` target.
5. Memory-aware heavy rules and three hardware profiles.
6. Dynamic Java heap sizing based on the scheduled rule memory.
7. Stable Conda environment name: `snap-snakemake`.
8. Concise dry-run workflow.
9. Timestamped run-level provenance/logging via `scripts/run_snap.sh`.

## Compatibility notes

- Existing `contig` joint sharding remains the default.
- `reference_intervals.mode: all` is the default and preserves whole-reference behavior.
- Existing final SNP IDs remain `CHROM_POS`.
- For distributed balanced sharding, all sequences within a generated balanced shard must be assigned to the same partition.
- `joint_finalize` must use the same reference-selection and sharding settings as the shard-producing runs.
- Java heap is no longer intended to be configured as a fixed `-Xmx` in user configs; it is calculated from `mem_mb`. Optional JVM flags can be added through `java_options.*_extra`.
