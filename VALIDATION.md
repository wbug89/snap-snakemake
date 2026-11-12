# v1.0.0 validation notes

Performed in the artifact environment:

- `workflow_inputs.py` compiled with Python 3.13.
- Smoke-tested reference `.fai` parsing and multi-directory gVCF discovery with synthetic inputs.
- Verified example TSVs contain real tab characters.
- Verified the v1.0.0 `Snakefile` contains `joint_shard`, `joint_finalize`, per-contig `GenomicsDBImport`, per-contig `GenotypeGVCFs`, and `concat_joint_vcfs` rules.
- Verified executable workflow code contains no `CombineGVCFs` rule/reference.
- Parsed all YAML configuration examples with PyYAML.

Limitation: Snakemake/GATK/bcftools executables are not installed in the artifact container, so an actual Snakemake DAG dry-run and real GATK execution could not be performed here. Run `snakemake -n -p` on PC1/PC3 before production execution.
