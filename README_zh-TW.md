# SNaP v1.0.1

### Snakemake for NGS and Polyploids

**SNaP makes SNP calling as simple as a snap.**

SNaP 是一套以 Snakemake 建立、可重現的 SNP calling workflow，主要用於 RAD-seq 類型的 NGS 資料，包括 ddRAD-seq 與 DArT-seq。它支援 SE/PE reads、可調整 ploidy、GVCF-based variant calling、per-contig GenomicsDB joint genotyping，以及跨多台 Linux server 的 distributed joint calling。

GitHub repository：`snap-snakemake`

## Workflow stages

| Stage | 主要輸入 | 主要工作 | 主要輸出 |
|---|---|---|---|
| `preprocess` | FASTQ | Trimmomatic → BWA-MEM → sort/index → MarkDuplicates | `marked/*.marked.bam` |
| `gvcf` | marked BAM | 依 biological sample 合併 units → HaplotypeCaller | `gvcf_file/*.g.vcf.gz` |
| `joint` | 完整 cohort gVCFs | per-contig GenomicsDBImport → GenotypeGVCFs → concat/filter | final cohort VCF |
| `joint_shard` | 完整 cohort gVCFs + contig assignment | 只分析分配給該 server 的 contigs | `results/per_contig/*.vcf.gz` |
| `joint_finalize` | 所有 per-contig VCF | 依 reference 順序 concat → SNP/INDEL → filtering | final cohort VCF |
| `full` | FASTQ | 單機 end-to-end | final cohort VCF |

分散式運算時需特別注意：**gVCF stage 可以按 sample 分機器；joint calling 必須按 reference contig 分機器，不能把 cohort 按 sample 拆成兩群。** 每一台執行 `joint_shard` 的 server 都必須能讀取完整 cohort 的 single-sample gVCFs。

---

## 1. 安裝環境

建立 Conda environment：

```bash
conda env create -f environment.yaml
conda activate snap-snakemake-v1.0.1
```

workflow 需要可呼叫 Java、BWA、samtools、bcftools、GATK 與 Trimmomatic；若不在 PATH，可在 config 的 `tools:` 指定絕對路徑。

---

## 2. Reference genome 必須先建立 index

在執行 SNaP 前，每台會使用 reference 的 server 都必須準備相同 reference FASTA 與 index。

假設：

```bash
REF=/path/to/reference.fa
```

### BWA index

`preprocess` 與 `full` 需要：

```bash
bwa index "$REF"
```

### samtools FASTA index

```bash
samtools faidx "$REF"
```

會產生：

```text
reference.fa.fai
```

### GATK sequence dictionary

```bash
gatk CreateSequenceDictionary \
  -R "$REF" \
  -O "${REF%.*}.dict"
```

若檔名為 `reference.fa`，最後應包含：

```text
reference.fa
reference.fa.fai
reference.dict
reference.fa.amb
reference.fa.ann
reference.fa.bwt
reference.fa.pac
reference.fa.sa
```

若 PC1 / PC3 都參與分析，必須確認兩台使用完全相同 reference。必要時：

```bash
sha256sum "$REF"
```

---

## 3. Sample 命名與 `sample_group_regex`

v1.0.1 的 generic example 使用 `SM001`–`SM999`。

預設範例：

```yaml
sample_group_regex: '^((?:SM)[0-9]{3})'
```

例如：

```text
SM001_R1.fastq.gz                  → SM001
SM001_R2.fastq.gz                  → SM001
SM001_laneA_R1.fastq.gz            → SM001
SM125_sequencingBatch2_R1.fastq.gz → SM125
```

SNaP 會使用 regex 的**第一個 capture group** 作為 biological sample ID。預設的 `((?:SM)[0-9]{3})` 會從 raw-unit 檔名開頭擷取 `SM001`、`SM125`、`SM999` 這類樣本名稱。

### 如何依自己的 sample 命名方式修改 `sample_group_regex`

原則是讓 **第一個 capture group** 精確擷取 biological sample ID。這個 regex 會從 raw-unit 檔名開頭開始比對，因此通常應保留最前面的 `^`。lane、sequencing batch、`_R1`、`_R2` 或其他 suffix 都應留在 capture group 之外。

預設範例為：

```yaml
sample_group_regex: '^((?:SM)[0-9]{3})'
```

各部分意義：

- `^`：從檔名或 unit name 的開頭開始比對。
- `((?:SM)[0-9]{3})`：第一個 capture group，作為 biological sample ID。
- `(?:SM)`：匹配固定字首 `SM`，但不建立額外 capture group。
- `[0-9]{3}`：匹配三位數字。

如果 sample 命名方式不同，可以依格式修改，例如：

```yaml
# SM001、SM125、SM999
sample_group_regex: '^((?:SM)[0-9]{3})'

# IB001、IB125、IB999
sample_group_regex: '^((?:IB)[0-9]{3})'

# IBC1-IBC9 或 TN + 數字，例如 IBC1、TN71
sample_group_regex: '^((?:IBC[0-9]+|TN[0-9]+))'

# ABC_001、ABC_125、XYZ_999
sample_group_regex: '^([A-Z]{3}_[0-9]{3})'
```

以預設設定：

```yaml
sample_group_regex: '^((?:SM)[0-9]{3})'
```

下列 raw units 都會被歸為 biological sample `SM001`：

```text
SM001_R1.fastq.gz
SM001_R2.fastq.gz
SM001_lane1_R1.fastq.gz
SM001_batchB_R1.fastq.gz
```

最重要的是：**第一個 capture group 必須只回傳 biological sample ID**，不要把 lane、batch 或 read suffix 一起包進去。例如 `SM001_lane1_R1.fastq.gz` 與 `SM001_lane2_R1.fastq.gz` 如果屬於同一 biological sample，就都應解析成 `SM001`，而不是 `SM001_lane1` 與 `SM001_lane2`。

請保留最外層用來擷取 sample ID 的括號 `(...)`。像 `(?:SM)` 這種 non-capturing group 可以用來組織 regex，但不會取代第一個 capture group。若 regex 無法擷取 biological sample ID，SNaP 會直接報錯停止，而不是默默指定錯誤樣本。

---

## 4. Config 使用方式

`config.yaml` 包含完整可調參數；`configs/*.example.yaml` 也都列出所有 user-adjustable fields，再依 stage 提供不同 input/output path 範例。

建議不要直接修改 `.example.yaml`，而是複製：

```bash
cp configs/config.preprocess.example.yaml config.preprocess.yaml
```

再修改自己的 local config。這些 local config 已由 `.gitignore` 排除。

### 主要 top-level 參數

| 參數 | 說明 |
|---|---|
| `workflow_stage` | `preprocess` / `gvcf` / `joint` / `joint_shard` / `joint_finalize` / `full` |
| `read_layout` | `SE` 或 `PE`，供 `preprocess/full` 使用 |
| `batch_id` | 一批獨立 preprocessing sequencing batch 的唯一名稱 |
| `ploidy` | HaplotypeCaller ploidy，例如 diploid=`2`、hexaploid=`6` |
| `project_dir` | 此次 stage 的 output root |
| `fastq_dir` | `preprocess/full` 的 FASTQ folder |
| `reference` | reference FASTA |
| `sample_group_regex` | biological sample capture regex，group 1 為 sample ID |
| `exclude_units` | preprocess 時要排除的 exact raw sequencing units |
| `include_biological_samples` | biological sample allowlist；空 list 表示全部 |
| `exclude_biological_samples` | biological sample blocklist |
| `include_biological_samples_file` | 每行一個 sample ID 的 include file |
| `exclude_biological_samples_file` | 每行一個 sample ID 的 exclude file |
| `marked_dirs` | `gvcf` stage 要讀取的 marked BAM folders |
| `sample_assignment_file` | 分散式 gVCF calling 的 sample-to-server TSV |
| `calling_partition` | 選擇 assignment TSV 中哪個 server partition |
| `gvcf_dirs` | `joint/joint_shard` 的完整 cohort gVCF folders |
| `joint_contig_assignment_file` | `joint_shard` 的 contig-to-server TSV |
| `joint_partition` | 該 server 要處理的 contig partition |
| `per_contig_dirs` | `joint_finalize` 要收集的 per-contig VCF folders |

### Tools

```yaml
tools:
  java: java
  bwa: bwa
  samtools: samtools
  bcftools: bcftools
  gatk: gatk
  trimmomatic: trimmomatic
```

若工具不在 PATH，改成絕對路徑即可。

### Trimmomatic

```yaml
trimmomatic:
  phred: -phred33
  threads: 8
  leading: 3
  trailing: 3
  slidingwindow: '4:30'
  minlen: 36
  builtin_adapter_se: NexteraPE-PE.fa
  builtin_adapter_pe: NexteraPE-PE.fa
  user_adapter_se: null
  user_adapter_pe: null
  illuminaclip_se: '2:30:10'
  illuminaclip_pe: '2:30:10'
```

SNaP 會先找 active Trimmomatic installation 內附的 adapter；若找不到且設定 `user_adapter_*`，才改用 fallback adapter FASTA。

### CPU / Java / GenomicsDB

```yaml
resources:
  bwa_threads: 12
  samtools_view_threads: 4
  samtools_sort_threads: 8
  haplotypecaller_threads: 2

java_options:
  markduplicates: '-Xmx20g'
  haplotypecaller: '-Xmx20g'
  genomicsdbimport: '-Xmx16g -XX:ParallelGCThreads=2'
  genotypegvcfs: '-Xmx32g -XX:ParallelGCThreads=2'

genomicsdb:
  batch_size: 30
  reader_threads: 1
  import_mem_mb: 24000
  genotype_mem_mb: 40000
```

`import_mem_mb` / `genotype_mem_mb` 是**單一 job**的 Snakemake resource declaration，不代表整台機器的總上限。joint calling 還應在 command line 加上例如：

```bash
--resources joint_jobs=1 mem_mb=120000
```

### Variant calling / filtering

```yaml
calling:
  minimum_mapping_quality: 30
  max_alternate_alleles: 6
  max_genotype_count: 1024
  hard_filter_expression: 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0'
```

polyploid multiallelic site 的 genotype-state 數量會快速增加，因此預設不提高 `max_genotype_count`。

---

## 5. `preprocess`

流程：

```text
FASTQ
→ Trimmomatic
→ BWA-MEM
→ samtools sort/index
→ GATK MarkDuplicates
→ marked BAM
```

RAD-seq / DArT-seq 的相同 alignment coordinates 不一定是 PCR duplicate，因此只標記、不移除：

```text
--REMOVE_DUPLICATES false
--REMOVE_SEQUENCING_DUPLICATES false
```

建立 local config：

```bash
cp configs/config.preprocess.example.yaml config.preprocess.yaml
```

至少調整：

```yaml
workflow_stage: preprocess
read_layout: PE
batch_id: batch1
project_dir: /path/to/preprocess_output
fastq_dir: /path/to/rawData
reference: /path/to/reference.fa
ploidy: 6
sample_group_regex: '^((?:SM)[0-9]{3})'
```

Dry-run：

```bash
snakemake \
  -s Snakefile \
  --configfile config.preprocess.yaml \
  --cores 30 \
  -n -p
```

正式執行：

```bash
snakemake \
  -s Snakefile \
  --configfile config.preprocess.yaml \
  --cores 30 \
  --printshellcmds \
  --rerun-incomplete
```

---

## 6. `gvcf`

此 stage 從 `marked_dirs` 尋找 BAM，將同一 biological sample 的不同 unit/batch 合併，再對每個 sample 執行一次 HaplotypeCaller。

為保留 RAD/DArT duplicate-flagged reads，HaplotypeCaller 使用：

```text
--disable-read-filter NotDuplicateReadFilter
```

建立：

```bash
cp configs/config.gvcf.example.yaml config.gvcf.yaml
```

設定例如：

```yaml
workflow_stage: gvcf
project_dir: /path/to/gvcf_output
reference: /path/to/reference.fa
ploidy: 6
marked_dirs:
  - /path/to/batch01/marked
  - /path/to/batch02/marked
```

### PC1 / PC3 按 sample 分工

`sample_assignment.tsv`：

```text
biological_sample    server
SM001                PC1
SM002                PC1
SM003                PC3
SM004                PC3
```

PC1：

```yaml
sample_assignment_file: sample_assignment.tsv
calling_partition: PC1
```

PC3 改為：

```yaml
calling_partition: PC3
```

正式執行：

```bash
snakemake \
  -s Snakefile \
  --configfile config.gvcf.yaml \
  --cores 30 \
  --printshellcmds \
  --rerun-incomplete
```

---

## 7. `joint`：單台 server joint calling

適合所有 contigs 都在同一台 server 做 joint calling。

`gvcf_dirs` 合併後必須剛好提供完整 cohort，每個 biological sample 只能有一份 indexed gVCF。

```yaml
workflow_stage: joint
project_dir: /path/to/joint_output
reference: /path/to/reference.fa
gvcf_dirs:
  - /path/to/gvcf_PC1/gvcf_file
  - /path/to/gvcf_PC3/gvcf_file
```

Dry-run：

```bash
snakemake \
  -s Snakefile \
  --configfile config.joint.yaml \
  --cores 30 \
  --resources joint_jobs=1 mem_mb=120000 \
  -n -p
```

正式：

```bash
snakemake \
  -s Snakefile \
  --configfile config.joint.yaml \
  --cores 30 \
  --resources joint_jobs=1 mem_mb=120000 \
  --printshellcmds \
  --rerun-incomplete
```

只有在實際監測 RAM 後確認安全，才建議提高 `joint_jobs`。

---

## 8. `joint_shard`：兩台 server 依 contig 分工

joint calling 必須按 contig 分工，每台都要看到**完整 cohort gVCFs**。

從範例建立：

```text
contig    server
LG1       PC1
LG2       PC1
LG3       PC3
...
```

reference `.fai` 中每個 contig 必須恰好出現一次。

PC1：

```yaml
workflow_stage: joint_shard
project_dir: /path/to/joint_PC1
reference: /path/to/reference.fa
gvcf_dirs:
  - /path/to/gvcf_PC1/gvcf_file
  - /path/to/gvcf_PC3/gvcf_file
joint_contig_assignment_file: joint_contig_assignment.tsv
joint_partition: PC1
```

PC3 使用同一 assignment file，但：

```yaml
joint_partition: PC3
```

PC1 範例：

```bash
snakemake \
  -s Snakefile \
  --configfile config.joint_shard_PC1.yaml \
  --cores 32 \
  --resources joint_jobs=1 mem_mb=200000 \
  --printshellcmds \
  --rerun-incomplete
```

PC3 範例：

```bash
snakemake \
  -s Snakefile \
  --configfile config.joint_shard_PC3.yaml \
  --cores 30 \
  --resources joint_jobs=1 mem_mb=120000 \
  --printshellcmds \
  --rerun-incomplete
```

輸出：

```text
results/per_contig/000001.vcf.gz
results/per_contig/000001.vcf.gz.tbi
...
```

`metadata/joint_contigs.tsv` 會記錄 shard 對應哪一條 contig。

---

## 9. `joint_finalize`

所有 `joint_shard` 完成後，只需要把各 server 的：

```text
results/per_contig/*.vcf.gz
results/per_contig/*.vcf.gz.tbi
```

集中到 final server；**不需要搬 GenomicsDB workspace**。

設定：

```yaml
workflow_stage: joint_finalize
project_dir: /path/to/joint_final
reference: /path/to/reference.fa
per_contig_dirs:
  - /path/to/joint_PC1/results/per_contig
  - /path/to/joint_PC3/results/per_contig
```

執行：

```bash
snakemake \
  -s Snakefile \
  --configfile config.joint_finalize.yaml \
  --cores 8 \
  --printshellcmds \
  --rerun-incomplete
```

workflow 會先檢查 contig 是否完整且沒有重複，再依 reference `.fai` 順序 concat，最後執行 SNP/INDEL selection、ID assignment、hard filtering 與 index。

---

## 10. `full`

`full` 適合單台 server 從 FASTQ 一路跑到 final VCF。

```bash
cp configs/config.full.example.yaml config.full.yaml
```

依資料修改 `fastq_dir`、`reference`、`project_dir`、`read_layout`、`batch_id`、`ploidy`、sample naming 等設定。

執行：

```bash
snakemake \
  -s Snakefile \
  --configfile config.full.yaml \
  --cores 30 \
  --resources joint_jobs=1 mem_mb=120000 \
  --printshellcmds \
  --rerun-incomplete
```

---

## 11. Telegram 通知

所有 example config 都包含：

```yaml
telegram:
  enabled: false
  bot_token: ''
  chat_id: ''
```

### 建議方式：用 environment variables 保存 credential

```bash
export TG_BOT_TOKEN='YOUR_BOT_TOKEN'
export TG_CHAT_ID='YOUR_CHAT_ID'
```

local config 只需要：

```yaml
telegram:
  enabled: true
  bot_token: ''
  chat_id: ''
```

SNaP 會優先讀 config；若 config 為空，才 fallback 到：

```text
TG_BOT_TOKEN
TG_CHAT_ID
```

也可以把 credential 寫在**不會 commit 的 local config**：

```yaml
telegram:
  enabled: true
  bot_token: 'YOUR_BOT_TOKEN'
  chat_id: 'YOUR_CHAT_ID'
```

不要把真實 Telegram credential push 到 GitHub。

若 `enabled: true` 但缺少 token/chat ID，workflow 只會警告並跳過 Telegram，不會因此停止分析。

---

## 12. 最終 SNP ID

joint genotyping 完成並分離 SNP 後，SNaP 使用：

```bash
bcftools annotate --set-id '%CHROM_%POS'
```

將 VCF `ID` 設成：

```text
CHROM_POS
```

例如：

```text
LG1_123456
LG5_9876543
```

最終：

```text
final_combined_filtered_snps.vcf.gz
```

會保留這些 `CHROM_POS` SNP IDs。

---

## 13. Dry-run、resume 與長時間分析

正式分析前建議先：

```bash
snakemake -s Snakefile --configfile CONFIG.yaml --cores 30 -n -p
```

中斷後續跑：

```bash
snakemake -s Snakefile --configfile CONFIG.yaml --cores 30 --rerun-incomplete
```

joint calling 請同時設定適合該 server 的全機資源上限，例如：

```bash
--resources joint_jobs=1 mem_mb=120000
```

長時間執行建議放在 `tmux` 中。

## Version

`v1.0.1` 為 documentation / example-config maintenance release；核心 SNP calling 與 GenomicsDB 邏輯沒有改變。
