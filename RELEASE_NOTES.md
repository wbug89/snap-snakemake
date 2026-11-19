# SNaP v1.0.3 Release Notes

SNaP v1.0.3 is a patch release for the already published v1.0.2 release. It fixes a release-blocking Full-stage target-selection bug discovered during Windows WSL2 user testing. The core SNP-calling methodology is unchanged.

## Fixed

- `rule all` is explicitly marked with `default_target: True`, so the presence or position of the `preflight` rule can no longer change the default workflow target.
- `scripts/run_snap.sh` explicitly targets `all` during both dry-run and production execution.
- Run-mode success reporting now requires successful completion of the explicit `all` target. A successful preflight alone is not reported as a completed Full-stage workflow.
- Dry-run now builds the complete `all` DAG and retains the Job stats summary while suppressing routine per-job console noise. Full native dry-run output is still saved in `run.log`.

## Added

- `profiles/minimal/profile.yaml` for systems with about 8 GiB of visible memory.
- Profile-aware preflight memory checks. If a profile's global memory budget exceeds detected memory, preflight reports `PASS WITH WARNINGS` and recommends a safer profile or additional memory.
- WSL2-specific guidance explaining that SNaP sees the memory assigned to WSL2, which may be lower than the Windows host's installed RAM.

## Resource profiles

| Profile | Recommended RAM | Intended use |
|---|---:|---|
| `minimal` | ~8 GiB | Small datasets, teaching, low-memory laptops, WSL2 testing |
| `low-memory` | ~16 GiB | General laptops/desktops |
| `workstation` | 32–64 GiB | Research workstations |
| `server` | ≥128 GiB | High-memory servers |

The `minimal` profile is intentionally conservative. Large cohorts, polyploid datasets, and highly multiallelic loci can still require substantially more memory.

## Unchanged core analysis

The following analysis chain is unchanged from v1.0.2:

```text
Trimmomatic
→ BWA-MEM
→ samtools
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

- RAD/DArT duplicate reads are marked but not removed.
- HaplotypeCaller retains duplicate-flagged reads.
- Ploidy remains configurable.
- GVCF and GenomicsDB joint genotyping remain unchanged.
- `max_alternate_alleles` and `max_genotype_count` behavior remains unchanged.
- FASTQ support remains `.fastq.gz`, `.fq.gz`, `.fastq`, and `.fq`.
- Reference interval control and balanced sharding remain unchanged.
- Telegram notification fallback behavior remains unchanged.
- Timestamped run-level logging remains unchanged.

## Upgrade recommendation

Users of v1.0.2 should upgrade to v1.0.3, especially when using `workflow_stage: full` or the `scripts/run_snap.sh` wrapper. Existing v1.0.2 Git tags/releases should remain available for reproducibility; do not move or overwrite the v1.0.2 tag.
