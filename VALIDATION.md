# SNaP v1.0.1 validation notes

- YAML parsing validated for `config.yaml` and all `configs/*.example.yaml`.
- Every example config exposes the same complete adjustable key set as the master config.
- `workflow_inputs.py` compiles successfully.
- Generic sample regex validated for `SM001`–`SM999` prefixes and rejects `SM000`/`SM1000`.
- Reference-indexing instructions added for BWA, samtools, and GATK.
- Final SNP IDs remain `CHROM_POS` via `bcftools annotate --set-id '%CHROM_%POS'`.
- No executable `CombineGVCFs` rule is present.
- Core SNP-calling / GenomicsDB logic is unchanged from v1.0.0.
- `docs/GIT_ON_MAC.md` is not included in the public project.
