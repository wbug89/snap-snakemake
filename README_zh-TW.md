# SNaP v1.0.0

### Snakemake for NGS and Polyploids

**SNaP makes SNP calling as simple as a snap.**

SNaP 是一套以 Snakemake 建立、可重現的 SNP calling workflow，主要面向 RAD-seq 類型的 NGS 資料，包括 ddRAD-seq 與 DArT-seq。它支援 SE/PE reads、可調整 ploidy、GVCF-based variant calling、per-contig GenomicsDB joint genotyping，以及跨多台 Linux server 的分散式運算。

正式專案名稱為 **SNaP**，GitHub repository 名稱為 **`snap-snakemake`**；`v1.0.0` 為第一個正式版本。名稱中的 SNaP 取自 **Snakemake for NGS and Polyploids**，並以「像 snap 一樣簡單完成 SNP calling」作為專案概念。

核心原則：

- 同一套 workflow 支援 SE 與 PE；
- ploidy 可設定，包含 hexaploid SNP calling；
- 保留適合 RAD-seq / DArT-seq 的 duplicate handling；
- gVCF stage 依 biological sample 分工；
- joint calling 依 reference contig/LG 分工，可供 PC1/PC3 分散式運算；
- 最終 SNP ID 統一為 `CHROM_POS`，例如 `LG1_123456`；
- 透過 Snakemake 與 per-contig GenomicsDB workspace 支援中斷後續跑。

## Stage 一覽

| stage | 主要輸入 | 主要工作 | 主要輸出 |
|---|---|---|---|
| `preprocess` | FASTQ | trim → BWA → sort → MarkDuplicates | `marked/*.marked.bam` |
| `gvcf` | `marked_dirs` | biological-sample merge → HaplotypeCaller | `gvcf_file/*.g.vcf.gz` |
| `joint` | 完整 cohort gVCFs | 單機所有 LG：GenomicsDBImport → GenotypeGVCFs → concat/filter | final cohort VCF |
| `joint_shard` | 完整 cohort gVCFs + LG assignment | 只處理本 server 被分配的 LG | `results/per_contig/*.vcf.gz` |
| `joint_finalize` | PC1 + PC3 per-contig VCFs | concat → SNP/INDEL → filter/index | final cohort VCF |
| `full` | FASTQ | 單機 end-to-end | final cohort VCF |


## 最終 SNP ID 命名規則

在 joint calling 完成並分離 SNP 後，workflow 會使用 `bcftools annotate` 將每個 SNP 的 VCF `ID` 欄位統一設定為：

```text
CHROM_POS
```

例如：

```text
LG1_123456
LG5_9876543
```

實際規則為：

```bash
bcftools annotate --set-id '%CHROM_%POS'
```

處理順序為：

```text
allsample.vcf.gz
  ↓ SelectVariants (SNP)
allsample_snps.vcf.gz
  ↓ set SNP ID = CHROM_POS
allsample_snps.ID.vcf.gz
  ↓ hard filtering
final_combined_filtered_snps.vcf.gz
```

因此最終 `final_combined_filtered_snps.vcf.gz` 會保留 `CHROM_POS` 格式的 SNP ID。這個命名方式便於後續 PLINK、GWAS、annotation 與跨分析結果比對。

## v1.0.0 正式版重點

本正式版沿用已驗證的 per-contig GenomicsDB 架構，並將「跨兩台 server 依 LG 分工」納入正式 stage，以避免單機一次啟動過多大型 `GenotypeGVCFs` jobs。

v1.0.0 不提高 GATK polyploid multiallelic genotype 上限；預設明確保留：

```text
--max-alternate-alleles 6
--max-genotype-count 1024
```

高倍體、多 allele 位點仍可能看到 `too many genotypes ... current limit 1024` warning；這不等同於 workflow fatal error。

## 安裝與檔案

```text
snap-snakemake/
├── Snakefile
├── workflow_inputs.py
├── config.yaml
├── environment.yaml
├── joint_contig_assignment.example.tsv
├── sample_assignment.example.tsv
├── configs/
│   ├── config.joint.example.yaml
│   ├── config.joint_shard_PC1.example.yaml
│   ├── config.joint_shard_PC3.example.yaml
│   └── config.joint_finalize.example.yaml
├── README.md
├── README_zh-TW.md
├── CHANGELOG.md
└── VERSION
```

Conda：

```bash
conda env create -f environment.yaml
conda activate snap-snakemake-v1.0.0
```

每台 server 使用的 reference 必須相同，並具有：

```text
reference.fa
reference.fa.fai
reference.dict
```

`preprocess/full` 額外需要 BWA index。

建議 PC1 / PC3 比對：

```bash
sha256sum /path/to/reference.fa
```

## 1. preprocess

保留 RAD-seq/DArT-seq duplicate 策略：

```text
MarkDuplicates --REMOVE_DUPLICATES false --REMOVE_SEQUENCING_DUPLICATES false
HaplotypeCaller --disable-read-filter NotDuplicateReadFilter
```

同一 biological sample 不同 sequencing batch 會保留不同 read group，之後在 gVCF stage 合併。

## 2. gVCF stage：按 biological sample 分 PC1 / PC3

兩台可共用：

```text
sample_assignment.tsv
```

格式：

```text
biological_sample    server
IB001                PC1
IB002                PC1
IB098                PC3
IBC1                 PC3
TN71                 PC3
```

PC1 設：

```yaml
workflow_stage: gvcf
sample_assignment_file: sample_assignment.tsv
calling_partition: PC1
```

PC3 改成 `PC3`。

重要：同一 sample 的所有 batch/unit marked BAM 必須放到被分配的 calling server，再只做一次 HaplotypeCaller。

## 3. 準備 distributed joint calling

### 3.1 兩台都需要完整 cohort gVCFs

joint stage 不能再按 sample 分 cohort。PC1 / PC3 都必須能讀到同一套完整 biological samples：

```text
IB001.g.vcf.gz + .tbi
...
TN71.g.vcf.gz + .tbi
```

建議使用本機磁碟；若用 NAS，需留意 I/O。

### 3.2 建立 LG assignment

複製：

```bash
cp joint_contig_assignment.example.tsv joint_contig_assignment.tsv
```

格式：

```text
contig    server
LG1       PC1
LG2       PC1
...
LG9       PC3
...
LG15      PC3
```

每一個 reference `.fai` contig 必須恰好出現一次。v1.0.0 啟動時會檢查缺少、重複、或不存在於 reference 的 contig。

範例 8/7 分配只是一個模板；正式分工可依 `.fai` 長度與實際運算負荷調整。

## 4. PC1：joint_shard

```bash
cp configs/config.joint_shard_PC1.example.yaml config.joint_shard_PC1.yaml
```

確認：

```yaml
workflow_stage: joint_shard
joint_contig_assignment_file: joint_contig_assignment.tsv
joint_partition: PC1
```

先 dry-run：

```bash
snakemake \
  -s Snakefile \
  --configfile config.joint_shard_PC1.yaml \
  --cores 32 \
  --resources joint_jobs=1 mem_mb=200000 \
  -n -p
```

正式：

```bash
mkdir -p logs
snakemake \
  -s Snakefile \
  --configfile config.joint_shard_PC1.yaml \
  --cores 32 \
  --resources joint_jobs=1 mem_mb=200000 \
  --printshellcmds \
  --rerun-incomplete \
  2>&1 | tee -a logs/PC1_joint_shard.log
```

PC1 RAM 較多時，在確認實際 RSS 後可測試 `joint_jobs=2`；不要僅靠 `--cores` 控制 Java 大型 jobs。

## 5. PC3：joint_shard

```bash
cp configs/config.joint_shard_PC3.example.yaml config.joint_shard_PC3.yaml
```

建議先用最穩定設定：

```bash
snakemake \
  -s Snakefile \
  --configfile config.joint_shard_PC3.yaml \
  --cores 30 \
  --resources joint_jobs=1 mem_mb=120000 \
  --printshellcmds \
  --rerun-incomplete \
  2>&1 | tee -a logs/PC3_joint_shard.log
```

### 為什麼仍要 `--resources`？

Snakefile 中每個 GenomicsDBImport/GenotypeGVCFs job 都宣告：

```text
joint_jobs = 1
mem_mb = ...
```

這只描述「一個 job 需要多少資源」。命令列的：

```text
--resources joint_jobs=1 mem_mb=120000
```

才是整個 Snakemake run 的 machine-wide resource ceiling。若省略，Snakemake 仍可能因 15 個 LG 互相獨立而同時啟動很多 Java jobs。

## 6. joint_shard 輸出

每台會得到：

```text
<project_dir>/genomicsdb/000001/
<project_dir>/results/per_contig/000001.vcf.gz
<project_dir>/results/per_contig/000001.vcf.gz.tbi
...
```

`000001` 是 reference `.fai` 中第一條 contig；mapping 會寫入：

```text
metadata/joint_contigs.tsv
```

不要用 shard 編號猜 LG，應以這份 manifest 為準。

## 7. 把 PC1 per-contig VCF 集中到 final server

只需要搬：

```text
results/per_contig/*.vcf.gz
results/per_contig/*.vcf.gz.tbi
```

不需要搬 GenomicsDB workspace。

例如 PC1 → PC3：

```bash
rsync -avhP \
  /path/on/PC1/joint_PC1/results/per_contig/ \
  user@PC3:/mnt/22T/CWL_pc3/20260907_rerun_snpcalling/joint_from_PC1/results/per_contig/
```

## 8. joint_finalize

PC3 config 例：

```yaml
workflow_stage: joint_finalize
project_dir: /mnt/22T/CWL_pc3/20260907_rerun_snpcalling/joint_final
reference: /mnt/22T/CWL_pc3/refgenome/pasi3/pasi3.fa

per_contig_dirs:
  - /mnt/.../joint_from_PC1/results/per_contig
  - /mnt/.../joint_PC3/results/per_contig
```

SNaP v1.0.0 會檢查：

- 每一個 reference contig 都有對應 shard VCF。
- 沒有 shard 在兩個目錄重複。
- 每個 `.vcf.gz` 都有 `.tbi`。

Dry-run：

```bash
snakemake \
  -s Snakefile \
  --configfile config.joint_finalize.yaml \
  --cores 8 \
  -n -p
```

正式：

```bash
snakemake \
  -s Snakefile \
  --configfile config.joint_finalize.yaml \
  --cores 8 \
  --printshellcmds \
  --rerun-incomplete \
  2>&1 | tee -a logs/joint_finalize.log
```

它會按照 reference `.fai` 順序執行：

```text
per-contig VCFs
  ↓
bcftools concat
  ↓
allsample.vcf.gz
  ↓
SelectVariants SNP / INDEL
  ↓
set ID
  ↓
hard filter
  ↓
final_combined_filtered_snps.vcf.gz
final_combined_filtered_indels.vcf.gz
```

## 9. 保留單機 joint

如果某次分析只想用一台 server，仍可：

```yaml
workflow_stage: joint
```

這會對所有 `.fai` contigs 建 GenomicsDB / GenotypeGVCFs，最後 concat/filter。

務必加 machine-wide limit，例如 PC3：

```bash
snakemake \
  -s Snakefile \
  --configfile configs/config.joint.example.yaml \
  --cores 30 \
  --resources joint_jobs=1 mem_mb=120000 \
  --rerun-incomplete -p
```

## 10. 恢復與完整性

Snakemake lock：

```bash
snakemake -s Snakefile --configfile CONFIG --unlock
```

恢復：

```bash
snakemake -s Snakefile --configfile CONFIG --cores 30 \
  --resources joint_jobs=1 mem_mb=120000 \
  --rerun-incomplete -n -p
```

GenomicsDBImport rule 在重跑該 shard 時會先刪除**該 shard 自己的 workspace**，避免 GATK 因目錄已存在而失敗；其他已完成 shard 不受影響。

檢查 per-contig VCF：

```bash
for f in results/per_contig/*.vcf.gz; do
  bcftools index --stats "$f" >/dev/null || echo "BAD: $f"
done
```

## 11. tmux

建議先進 session，再執行 Snakemake：

```bash
tmux new -s joint_PC3
```

而不是把 Snakemake 當成 tmux session 的唯一 command。也可：

```bash
tmux set-option -t joint_PC3 remain-on-exit on
```

讓 pane 中 command 結束後保留畫面，方便看 exit 狀態。

## 12. 從內部 v3.x 開發版遷移至 SNaP v1.0.0

v3.1 原本：

```text
all cohort gVCFs
→ 15 × GenomicsDBImport
→ 15 × GenotypeGVCFs
→ concat
→ filter
```

SNaP v1.0.0 可改成：

```text
                         ┌→ PC1 joint_shard → assigned LG VCFs ─┐
all cohort gVCFs ────────┤                                      ├→ joint_finalize → final VCF
                         └→ PC3 joint_shard → assigned LG VCFs ─┘
```

這不是把 cohort 分成兩半。PC1/PC3 的每個 LG joint call 都仍使用**完整 cohort 的所有 biological samples**；分割的是 genomic intervals/LGs。
