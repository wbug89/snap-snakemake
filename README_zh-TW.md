# SNaP

**Snakemake for NGS and Polyploids**

> **SNaP makes SNP calling as simple as a snap.**

SNaP 是以 Snakemake 建立的 RAD-seq / ddRAD-seq / DArT-seq SNP calling workflow，支援 SE/PE、可調 ploidy、GVCF calling、GenomicsDB joint genotyping，以及依 reference sequence 分散式 joint calling。

目前版本：**v1.0.3**

## v1.0.3 修正重點

- 修正 v1.0.2 GitHub release 後實測發現的 Full-stage default target 問題。
- `rule all` 明確設定為 Snakemake default target。
- `scripts/run_snap.sh` 在 dry-run 與正式 run 時明確指定 `all`。
- 修正成功判斷，避免只有 preflight 完成卻誤報 Full-stage 成功。
- 新增適合約 8 GiB 可見記憶體環境的 `minimal` resource profile。
- 新增 profile-aware RAM warning、`PASS WITH WARNINGS` 與 WSL2 記憶體提示。
- 核心 SNP calling 與 filtering 邏輯維持 v1.0.2 不變。

## v1.0.2 導入的主要功能

- FASTQ 支援 `.fastq.gz`、`.fq.gz`、`.fastq`、`.fq`。
- PE 自動驗證 R1/R2 配對，避免 trimmed/paired/unpaired output 被重新當成 raw input。
- `sample_group_regex` 啟動前驗證與解析 preview。
- reference preflight：檢查 FASTA、`.fai`、GATK `.dict`；`preprocess/full` 額外檢查 BWA index。
- 可選擇全部 reference sequences、regex、最短長度、指定 list 或檔案。
- joint sharding 支援 `contig` 與 length-balanced 模式。
- 所有 heavy rules 宣告 `mem_mb`；Java `-Xmx` 由 rule memory 自動計算。
- 內建 `minimal`、`low-memory`、`workstation`、`server` 四種 workflow profiles。
- Linux、macOS，以及 Windows 10/11 via WSL2。
- `preflight` target 提供簡潔的執行前摘要。
- `scripts/run_snap.sh` 每次執行建立完整 run history/log。
- 最終 SNP ID 維持 `CHROM_POS`，例如 `LG1_123456`。


## 快速索引

- [平台支援](#1-平台支援)
- [建立 Conda 環境](#2-建立-conda-環境)
- [Reference genome 準備](#3-reference-genome-準備)
- [FASTQ input](#4-fastq-input)
- [`sample_group_regex`](#5-sample_group_regex)
- [Reference interval selection](#6-reference-interval-selection)
- [Joint sharding](#7-joint-sharding)
- [Resource profiles](#8-resource-profiles)
- [Preflight](#9-preflight)
- [簡潔 dry-run](#10-簡潔-dry-run)
- [執行歷史與 log](#11-每次執行自動保存歷史-log)
- [Workflow stages](#121-preprocess)
  - [`preprocess`](#121-preprocess)
  - [`gvcf`](#122-gvcf)
  - [`joint`](#123-joint)
  - [`joint_shard`](#124-joint_shard)
  - [`joint_finalize`](#125-joint_finalize)
  - [`full`](#126-full)
- [Config 主要參數](#13-config-主要參數)
- [Telegram 通知](#14-telegram-通知)
- [Final SNP ID](#15-final-snp-id)
- [中斷、lock 與 recovery](#16-中斷lock-與-recovery)
- [Debugging](#17-debugging)

---

## 快速開始：在一般桌機或筆電執行 Full stage

這個 Quick Start 適合希望直接從 FASTQ 執行到最終 SNP/INDEL VCF 的使用者。

Full stage 會執行：

```text
FASTQ
→ Trimmomatic
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

### Step 1 — 進入 SNaP 專案並建立環境

第一次使用：

```bash
cd /path/to/snap-snakemake

conda env create -f environment.yaml
conda activate snap-snakemake
```

如果 `snap-snakemake` Conda environment 已存在：

```bash
conda env update \
  -n snap-snakemake \
  -f environment.yaml \
  --prune

conda activate snap-snakemake
```

確認主要工具：

```bash
snakemake --version
bwa 2>&1 | head
samtools --version
bcftools --version
gatk --version
```

### Step 2 — 準備 reference genome

假設 reference genome 為：

```bash
REF=/absolute/path/to/reference.fa
```

建立必要 index：

```bash
bwa index "$REF"
samtools faidx "$REF"
gatk CreateSequenceDictionary \
  -R "$REF" \
  -O "${REF%.*}.dict"
```

`full` stage 需要 FASTA、`.fai`、GATK sequence dictionary 與 BWA indexes。

### Step 3 — 準備 FASTQ

SNaP 支援：

```text
.fastq.gz
.fq.gz
.fastq
.fq
```

PE 範例：

```text
SM001_R1.fastq.gz
SM001_R2.fastq.gz
SM002_R1.fastq.gz
SM002_R2.fastq.gz
```

預設：

```yaml
sample_group_regex: '^((?:SM)[0-9]{3})'
```

例如：

```text
SM001_L001_R1.fastq.gz → SM001
SM001_L002_R1.fastq.gz → SM001
```

如果 sample naming 不同，請參考 [`sample_group_regex`](#5-sample_group_regex) 章節修改。

### Step 4 — 建立 Full-stage config

先複製完整範例：

```bash
cp configs/config.full.example.yaml config.full.yaml
```

再編輯：

```bash
nano config.full.yaml
```

一般使用者至少需要確認：

```yaml
workflow_stage: full

read_layout: PE
batch_id: batch1
ploidy: 2

project_dir: /absolute/path/to/snap_output
fastq_dir: /absolute/path/to/fastq
reference: /absolute/path/to/reference.fa

sample_group_regex: '^((?:SM)[0-9]{3})'

reference_intervals:
  mode: all

joint_sharding:
  strategy: contig
  shards: null

telegram:
  enabled: false
  bot_token: ''
  chat_id: ''
```

其中：

- `read_layout`: `SE` 或 `PE`
- `ploidy`: 依物種設定，例如 diploid = `2`、hexaploid = `6`
- `project_dir`: 分析結果輸出位置
- `fastq_dir`: raw FASTQ 所在資料夾
- `reference`: reference FASTA
- `sample_group_regex`: biological sample ID 解析規則
- `reference_intervals.mode: all`: 預設分析全部 reference sequences
- `joint_sharding.strategy: contig`: 預設一條 selected sequence 一個 joint shard

其他進階參數可以先保留 `config.full.example.yaml` 的預設值。

### Step 4.1 — 可選：快速設定 Telegram 通知

如果希望 workflow 開始、成功或失敗時收到 Telegram 通知，先建立 Telegram Bot，取得：

- `bot_token`
- `chat_id`

可以直接寫入本機使用的 config：

```yaml
telegram:
  enabled: true
  bot_token: 'YOUR_BOT_TOKEN'
  chat_id: 'YOUR_CHAT_ID'
```

如果不使用 Telegram：

```yaml
telegram:
  enabled: false
  bot_token: ''
  chat_id: ''
```

不要把含真實 credentials 的 local config commit 到 GitHub。

更建議使用環境變數，避免把 token 寫入 config：

```bash
export TG_BOT_TOKEN='YOUR_BOT_TOKEN'
export TG_CHAT_ID='YOUR_CHAT_ID'
```

config 保持：

```yaml
telegram:
  enabled: true
  bot_token: ''
  chat_id: ''
```

SNaP 會優先使用 config 中的值；config 為空時，再使用 `TG_BOT_TOKEN` 與 `TG_CHAT_ID`。若 `enabled: true` 但 credentials 不完整，SNaP 只會顯示警告並略過通知，不會因此中止 workflow。

下一步的 preflight 也會協助確認 Telegram 設定是否可被正確解析。

### Step 5 — 選擇硬體 profile

SNaP 提供四種 profile：

| Profile | 建議硬體 | 用途 |
|---|---|---|
| `minimal` | 約 8 GiB RAM、2–4 CPU threads | 小型資料、教學、低記憶體筆電或 WSL2 測試 |
| `low-memory` | 約 16 GiB RAM、4–8 CPU threads | 一般筆電 / 桌機 |
| `workstation` | 約 32–64 GiB RAM、8–16 threads | 一般研究工作站 |
| `server` | ≥128 GiB RAM | 高記憶體 server |

若 SNaP 偵測到只有約 7–8 GiB RAM，建議使用 **`minimal`**；約 16 GiB 使用 `low-memory`。大型 cohort、多倍體或高 multiallelic complexity 資料仍可能需要比 `minimal` 更多的記憶體。

不要因為電腦有很多 CPU cores 就直接使用：

```bash
snakemake --cores all
```

記憶體通常比 CPU core 數更容易成為 SNP calling workflow 的限制因素。

### Step 6 — 執行 Preflight

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile minimal \
  --mode preflight
```

Preflight 會檢查 config、FASTQ、PE R1/R2 pairing、`sample_group_regex`、reference indexes、reference interval selection、joint sharding、resource configuration 與 Telegram 設定。

如果 preflight 失敗，應先修正錯誤，不要直接進行正式分析。

### Step 7 — 執行簡潔 Dry-run

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile minimal \
  --mode dry-run
```

Dry-run 明確指定 target `all`，只建立與檢查完整 workflow DAG，不會真正執行分析。wrapper 會把 Snakemake 原始 dry-run 完整輸出保存到 `run.log`，terminal 則過濾成 DAG / Job stats 與重要狀態或錯誤訊息，因此可以保留任務統計表，又不會被逐 sample job details 洗版。這也避免 Snakemake 9.x 的 bare `--quiet` 可能連 Job stats 一起隱藏的情況。

只有除錯時才需要完整 shell commands：

```bash
snakemake \
  -s Snakefile \
  --configfile config.full.yaml \
  --workflow-profile profiles/minimal \
  all \
  --dry-run \
  --printshellcmds
```

### Step 8 — 正式執行 Full stage

約 7–8 GiB RAM 的電腦或低記憶體 WSL2：

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile minimal \
  --mode run
```

約 16 GiB RAM 的一般筆電 / 桌機：

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile low-memory \
  --mode run
```

32–64 GiB RAM 工作站：

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile workstation \
  --mode run
```

高記憶體 server：

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile server \
  --mode run
```

### Step 9 — 查看執行紀錄

每次透過 `run_snap.sh` 執行都會建立：

```text
<project_dir>/logs/runs/YYYY-MM-DD_HHMMSS/
```

最新一次執行：

```bash
less /absolute/path/to/snap_output/logs/runs/latest/run.log
```

每個 run 會保存：

```text
run.log
command.txt
config.snapshot.yaml
preflight_summary.txt
environment.txt
git_version.txt
run_metadata.tsv
```

因此分析結束或失敗後，可以回頭確認使用的 config、SNaP Git commit/tag、Conda environment、執行 command、Snakemake progress、error 與 exit status。

### Step 10 — 確認最終輸出

主要輸出：

```text
<project_dir>/results/final_combined_filtered_snps.vcf.gz
<project_dir>/results/final_combined_filtered_snps.vcf.gz.tbi
<project_dir>/results/final_combined_filtered_indels.vcf.gz
<project_dir>/results/final_combined_filtered_indels.vcf.gz.tbi
```

SNP ID 格式為 `CHROM_POS`，例如：

```text
LG1_123456
chr05_987654
```

快速查看：

```bash
bcftools view \
  -H \
  /absolute/path/to/snap_output/results/final_combined_filtered_snps.vcf.gz \
  | head
```

### 最短版 Full-stage SOP

已完成 Conda environment、reference indexing 與 `config.full.yaml` 設定後，一般使用者只需要依序執行：

```bash
conda activate snap-snakemake
```

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile minimal \
  --mode preflight
```

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile minimal \
  --mode dry-run
```

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile minimal \
  --mode run
```

約 16 GiB RAM 可把 `minimal` 改成 `low-memory`；32–64 GiB 使用 `workstation`；≥128 GiB 使用 `server`。

---

## 1. 平台支援

建議環境：

- Linux：Ubuntu / Debian 等常見 distribution。
- macOS：Intel 或 Apple Silicon；實際套件可用性依 Bioconda build 而定。
- Windows 10/11：**透過 WSL2 執行**。
- 不建議直接由 Windows CMD / PowerShell native 執行 bioinformatics workflow。

Windows/WSL 可能產生 `*:Zone.Identifier` metadata，已加入 `.gitignore`。

---

## 2. 建立 Conda 環境

環境名稱固定為：

```text
snap-snakemake
```

建立：

```bash
conda env create -f environment.yaml
conda activate snap-snakemake
```

既有環境更新：

```bash
conda env update \
  -n snap-snakemake \
  -f environment.yaml \
  --prune
```

軟體版本由 `VERSION` 與 Git tag 管理，不再把版本號放入 Conda environment name。

---

## 3. Reference genome 準備

假設 reference：

```bash
REF=/path/to/reference.fa
```

### 3.1 BWA index

`preprocess` 與 `full` 需要：

```bash
bwa index "$REF"
```

會建立：

```text
reference.fa.amb
reference.fa.ann
reference.fa.bwt
reference.fa.pac
reference.fa.sa
```

### 3.2 samtools FASTA index

所有需要 reference sequence metadata 的 stage 都需要：

```bash
samtools faidx "$REF"
```

建立：

```text
reference.fa.fai
```

### 3.3 GATK sequence dictionary

```bash
gatk CreateSequenceDictionary \
  -R "$REF" \
  -O /path/to/reference.dict
```

例如 `reference.fa` 對應 `reference.dict`。

### Stage requirements

| Stage | FASTA | `.fai` | `.dict` | BWA index |
|---|---:|---:|---:|---:|
| preprocess | ✓ | ✓ | ✓ | ✓ |
| gvcf | ✓ | ✓ | ✓ | — |
| joint | ✓ | ✓ | ✓ | — |
| joint_shard | ✓ | ✓ | ✓ | — |
| joint_finalize | ✓ | ✓ | ✓ | — |
| full | ✓ | ✓ | ✓ | ✓ |

若缺少 prerequisite，SNaP 會在建立 DAG 前停止並顯示建立 index 的指令。

---

## 4. FASTQ input

支援：

```text
.fastq.gz
.fq.gz
.fastq
.fq
```

### PE 命名

必須能辨識 `_R1` / `_R2`：

```text
SM001_R1.fastq.gz
SM001_R2.fastq.gz

SM002_L001_R1_001.fq.gz
SM002_L001_R2_001.fq.gz
```

R1/R2 可以使用上述任一支援副檔名。缺少 mate 時會在啟動階段直接報錯。

### SE 命名

可接受：

```text
SM001.fastq.gz
SM001_R1.fastq.gz
SM001.fq
```

### 避免重新讀入 workflow output

檔名含有下列 output token 的 FASTQ 會從 raw discovery 排除：

```text
trim
trimmed
paired
unpaired
```

因此例如：

```text
SM001_R1_paired.trim.fastq.gz
SM001_R2_unpaired.trim.fastq.gz
```

不會被重新當成 raw sequencing unit。

---

## 5. `sample_group_regex`

預設：

```yaml
sample_group_regex: '^((?:SM)[0-9]{3})'
```

第一個 capture group `(...)` **必須只包含 biological sample ID**。

例如：

```text
SM001_L001_R1.fastq.gz  -> SM001
SM001_run2_R1.fastq.gz  -> SM001
```

不同命名方式可自行調整：

```yaml
# SM001
sample_group_regex: '^((?:SM)[0-9]{3})'

# IB001
sample_group_regex: '^((?:IB)[0-9]{3})'

# IBC1 / TN71
sample_group_regex: '^((?:IBC[0-9]+|TN[0-9]+))'

# ABC_001
sample_group_regex: '^([A-Z]{3}_[0-9]{3})'
```

不要把 lane/run suffix 包進 capture group，例如若同一 biological sample 有：

```text
SM001_lane1_R1.fastq.gz
SM001_lane2_R1.fastq.gz
```

應解析成同一 `SM001`，讓 BAM 在 HaplotypeCaller 前合併。

preflight 會顯示前幾筆解析 preview。

---

## 6. Reference interval selection

設定：

```yaml
reference_intervals:
  mode: all
  regex: null
  min_length: null
  sequences: []
  file: null
```

### `all`

使用 `.fai` 中全部 sequence：

```yaml
reference_intervals:
  mode: all
```

### `regex`

例如只使用 chromosome：

```yaml
reference_intervals:
  mode: regex
  regex: '^chr[0-9]+$'
```

或 linkage groups：

```yaml
reference_intervals:
  mode: regex
  regex: '^LG[0-9]+$'
```

### `min_length`

例如只保留 ≥1 Mb sequences：

```yaml
reference_intervals:
  mode: min_length
  min_length: 1000000
```

### `list`

```yaml
reference_intervals:
  mode: list
  sequences:
    - chr1
    - chr2
    - chr3
```

### `file`

```yaml
reference_intervals:
  mode: file
  file: /path/to/chromosomes.txt
```

檔案每行一個 sequence name。

> SNaP v1.0.2 不會自行猜測「chr / contig」。應依 reference 命名使用 regex/list/file。以 NCBI assembly report 自動辨識 sequence role 留待後續版本。

HaplotypeCaller 與 joint genotyping 只分析 selected sequences。

---

## 7. Joint sharding

### `contig`

預設、與 v1.0.1 相容：

```yaml
joint_sharding:
  strategy: contig
  shards: null
```

一條 selected reference sequence = 一個 GenomicsDBImport + GenotypeGVCFs shard。

### `balanced`

大量 contigs/scaffolds 時建議：

```yaml
joint_sharding:
  strategy: balanced
  shards: 16
```

SNaP 依 `.fai` sequence length 使用 deterministic longest-first greedy balancing，把 selected sequences 分到指定數量 shards。

例如 5000 scaffolds 可以壓縮成 16 或 32 個 joint jobs，而不是建立 5000 個 jobs。

輸出：

```text
metadata/reference_sequences.tsv
```

包含：

```text
sequence
length
selected
exclusion_reason
shard
```

### Distributed `joint_shard` + balanced

`joint_contig_assignment.tsv` 仍以 sequence 指定 server/partition：

```text
sequence    server
LG1         PC1
LG2         PC1
LG3         PC3
```

如果使用 `balanced`，同一個 balanced shard 內的所有 sequences 必須分配到同一 partition。若同一 shard 被拆到 PC1/PC3，SNaP 會直接拒絕啟動。

`joint_finalize` 必須使用與 shard runs 完全一致的：

```yaml
reference_intervals:
joint_sharding:
```

設定。

---

## 8. Resource profiles

SNaP 提供：

```text
profiles/minimal/
profiles/low-memory/
profiles/workstation/
profiles/server/
```

### minimal

目標：約 8 GiB RAM、2–4 CPU threads。

```text
cores: 2
global mem_mb: 6000
joint_jobs: 1
```

適合小型資料、教學、smoke test 與低記憶體筆電 / WSL2。大型 cohort、多倍體或高 multiallelic complexity 資料仍可能需要更多 RAM。

### low-memory

目標：約 16 GiB RAM、4–8 CPU threads。

```text
cores: 4
global mem_mb: 12000
joint_jobs: 1
```

### workstation

目標：約 32–64 GiB RAM、8–16 threads。

```text
cores: 12
global mem_mb: 48000
joint_jobs: 1
```

### server

目標：128 GiB+ RAM。

```text
cores: 30
global mem_mb: 150000
joint_jobs: 2
```

所有 heavy rules 都具有 per-job `mem_mb`：

- Trimmomatic
- BWA
- samtools sort
- MarkDuplicates
- BAM merge
- HaplotypeCaller
- GenomicsDBImport
- GenotypeGVCFs
- joint VCF concat/sort

Snakefile 將 `mem_mb` 與 `joint_jobs` 設為 global scheduling resources。profile 的 global memory limit 因此會限制多個 high-memory jobs 同時執行。

GATK Java heap 不再固定為 20/32 GB，而是自動使用約：

```text
0.75 × rule mem_mb
```

保留約 25% 給 Java native memory、GenomicsDB/TileDB、compression buffers 與 OS。

> Resource declarations 是 scheduler limits，不是即時 RAM monitor；若單一資料集本身超過指定 memory，該 job 仍可能失敗。此設計主要避免多個 heavy jobs 同時把整台電腦 RAM 用盡。

---

## 9. Preflight

建議每次正式執行前先做。wrapper 會把所選 profile 名稱與 global memory budget 傳給 preflight；如果 profile 的 `mem_mb` 高於 SNaP 偵測到的實體記憶體，Validation 會顯示 `PASS WITH WARNINGS`。在 WSL2 中，SNaP 偵測到的是 WSL2 可見的記憶體，不一定等於 Windows 主機實際安裝的 RAM。

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile workstation \
  --mode preflight
```

會檢查：

- workflow stage
- SE/PE
- ploidy
- input directories
- FASTQ discovery / R1-R2 pairing
- `sample_group_regex`
- reference FASTA/indexes
- reference interval selection
- sample/gVCF discovery
- distributed assignment
- joint shard consistency

並產生：

```text
metadata/preflight_summary.txt
metadata/reference_sequences.tsv
metadata/selected_intervals.list
metadata/fastq_units.tsv
```

---

## 10. 簡潔 dry-run

大型 workflow 不建議預設加 `-p/--printshellcmds`。

使用 wrapper：

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile workstation \
  --mode dry-run
```

wrapper 會明確指定 `all --dry-run`。完整 Snakemake dry-run 輸出會保存到該次執行的 `run.log`，terminal 則只保留 DAG / Job stats 與重要狀態或錯誤訊息，因此不會被逐 sample 的 job details 洗版。

需要 debug command 時才直接使用：

```bash
snakemake \
  -s Snakefile \
  --configfile config.full.yaml \
  --workflow-profile profiles/workstation \
  all \
  --dry-run \
  --printshellcmds
```

---

## 11. 每次執行自動保存歷史 log

建議正式執行一律透過：

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile workstation \
  --mode run
```

每一次 wrapper invocation 會建立：

```text
<project_dir>/logs/runs/YYYY-MM-DD_HHMMSS/
├── run.log
├── command.txt
├── config.snapshot.yaml
├── preflight_summary.txt
├── environment.txt
├── git_version.txt
└── run_metadata.tsv
```

內容：

- `run.log`：Snakemake terminal output、job progress、錯誤。
- `command.txt`：實際執行 command。
- `config.snapshot.yaml`：該次分析 config 快照。
- `preflight_summary.txt`：輸入/reference/sample 摘要。
- `environment.txt`：OS、Python、Snakemake、Conda packages。
- `git_version.txt`：SNaP VERSION、Git commit/tag/branch。
- `run_metadata.tsv`：run ID、開始/結束時間、profile、status、exit code。

最近一次 run 可由：

```text
<project_dir>/logs/runs/latest
```

查看，例如：

```bash
less <project_dir>/logs/runs/latest/run.log
```

---

# 12. Workflow stages

## 12.1 `preprocess`

用途：

```text
FASTQ
→ Trimmomatic
→ BWA-MEM
→ samtools sort/index
→ GATK MarkDuplicates
```

RAD/DArT duplicate handling：duplicates **只標記、不移除**。

設定：

```yaml
workflow_stage: preprocess
read_layout: PE
batch_id: batch1
project_dir: /path/to/preprocess
fastq_dir: /path/to/raw
reference: /path/to/reference.fa
```

執行：

```bash
bash scripts/run_snap.sh \
  --config config.preprocess.yaml \
  --profile workstation \
  --mode run
```

主要 output：

```text
marked/*.marked.bam
marked/*.marked.bam.bai
marked/*.marked_dup_metrics.txt
metadata/<batch_id>.units.tsv
```

---

## 12.2 `gvcf`

用途：

```text
marked BAM units
→ group by biological sample
→ merge technical batches
→ HaplotypeCaller -ERC GVCF
```

設定：

```yaml
workflow_stage: gvcf
project_dir: /path/to/gvcf_output
reference: /path/to/reference.fa
ploidy: 6
marked_dirs:
  - /path/to/batch1/marked
  - /path/to/batch2/marked
```

單機：

```bash
bash scripts/run_snap.sh \
  --config config.gvcf.yaml \
  --profile workstation \
  --mode run
```

多 server 可使用：

```yaml
sample_assignment_file: sample_assignment.tsv
calling_partition: PC1
```

另一台使用相同 assignment file、改為 `PC3`。

HaplotypeCaller 保留 duplicate-flagged reads：

```text
--disable-read-filter NotDuplicateReadFilter
```

主要 output：

```text
gvcf_file/SM001.g.vcf.gz
gvcf_file/SM001.g.vcf.gz.tbi
```

---

## 12.3 `joint`

單機完整 joint genotyping：

```text
all single-sample gVCFs
→ GenomicsDBImport by shard
→ GenotypeGVCFs
→ concatenate + sort
→ SNP / INDEL split
→ ID assignment
→ hard filtering
```

設定：

```yaml
workflow_stage: joint
project_dir: /path/to/joint
gvcf_dirs:
  - /path/to/gvcf_PC1
  - /path/to/gvcf_PC3
reference: /path/to/reference.fa
```

每個 joint run 必須能看到**完整 cohort**。

執行：

```bash
bash scripts/run_snap.sh \
  --config config.joint.yaml \
  --profile server \
  --mode run
```

---

## 12.4 `joint_shard`

把 genomic shards 分給不同 server；不是把 samples 分 cohort。

例如：

```yaml
workflow_stage: joint_shard
joint_contig_assignment_file: joint_contig_assignment.tsv
joint_partition: PC1
```

PC1：

```bash
bash scripts/run_snap.sh \
  --config config.joint_shard_PC1.yaml \
  --profile server \
  --mode run
```

PC3 使用對應 config。

兩台都必須看到相同完整 single-sample gVCF cohort。

主要 output：

```text
results/per_contig/<shard>.vcf.gz
results/per_contig/<shard>.vcf.gz.tbi
```

即使目錄仍叫 `per_contig`（維持相容性），balanced mode 的一個 shard 可以包含多個 reference sequences。

---

## 12.5 `joint_finalize`

收集所有 shard VCF：

```yaml
workflow_stage: joint_finalize
per_contig_dirs:
  - /path/to/PC1/results/per_contig
  - /path/to/PC3/results/per_contig
```

執行：

```bash
bash scripts/run_snap.sh \
  --config config.joint_finalize.yaml \
  --profile workstation \
  --mode run
```

balanced shards 可能包含不連續 sequences，因此 finalize 使用：

```text
bcftools concat → bcftools sort
```

確保 final VCF 回到 reference header order。

---

## 12.6 `full`

單機 end-to-end：

```text
FASTQ → final SNP/INDEL VCF
```

設定：

```yaml
workflow_stage: full
read_layout: PE
project_dir: /path/to/project
fastq_dir: /path/to/raw
reference: /path/to/reference.fa
```

一般電腦：

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile low-memory \
  --mode run
```

workstation：

```bash
bash scripts/run_snap.sh \
  --config config.full.yaml \
  --profile workstation \
  --mode run
```

---

# 13. Config 主要參數

## Core

| Key | 用途 |
|---|---|
| `workflow_stage` | preprocess/gvcf/joint/joint_shard/joint_finalize/full |
| `read_layout` | SE / PE |
| `batch_id` | preprocessing batch ID |
| `ploidy` | HaplotypeCaller ploidy |
| `project_dir` | output root |
| `fastq_dir` | raw FASTQ directory |
| `reference` | reference FASTA |
| `sample_group_regex` | biological sample extraction |

## Sample filters

```yaml
exclude_units: []
include_biological_samples: []
exclude_biological_samples: []
include_biological_samples_file: null
exclude_biological_samples_file: null
```

## gVCF distribution

```yaml
marked_dirs: []
sample_assignment_file: null
calling_partition: null
```

## Joint

```yaml
gvcf_dirs: []
joint_contig_assignment_file: null
joint_partition: null
per_contig_dirs: []
```

## Resource defaults

```yaml
resources:
  trimmomatic_threads: 2
  trimmomatic_mem_mb: 2048
  bwa_threads: 4
  bwa_mem_mb: 4096
  samtools_view_threads: 1
  samtools_sort_threads: 2
  samtools_sort_mem_mb: 4096
  markduplicates_mem_mb: 6144
  merge_bam_mem_mb: 4096
  haplotypecaller_threads: 2
  haplotypecaller_mem_mb: 8192
  bcftools_mem_mb: 2048
```

一般情況優先選 `profiles/*`，不必逐項修改這些值。特殊硬體才調整 config 或 profile。

## GenomicsDB

```yaml
genomicsdb:
  batch_size: 30
  reader_threads: 1
  import_mem_mb: 8192
  genotype_mem_mb: 12288
```

## Calling

```yaml
calling:
  minimum_mapping_quality: 30
  max_alternate_alleles: 6
  max_genotype_count: 1024
  hard_filter_expression: 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0'
```

---

# 14. Telegram 通知

建議不要把 credential 寫進 Git-tracked config。

```yaml
telegram:
  enabled: true
  bot_token: ''
  chat_id: ''
```

在 shell：

```bash
export TG_BOT_TOKEN='YOUR_BOT_TOKEN'
export TG_CHAT_ID='YOUR_CHAT_ID'
```

若 `enabled: true` 但 credential 不完整，SNaP 只會警告並略過 Telegram，不會讓 workflow 因通知失敗而中止。

---

# 15. Final SNP ID

SNP 經：

```bash
bcftools annotate --set-id '%CHROM_%POS'
```

因此 final：

```text
results/final_combined_filtered_snps.vcf.gz
```

ID 例如：

```text
LG1_123456
chr03_987654
```

---

# 16. 中斷、lock 與 recovery

先確認沒有另一個 Snakemake：

```bash
pgrep -af snakemake
```

若前一次因斷電、kill、SSH/terminal 異常而留下 stale lock：

```bash
snakemake \
  -s Snakefile \
  --configfile config.yaml \
  --unlock
```

`--unlock` 只解除 Snakemake lock，**不會刪除 BAM/GVCF/VCF 結果**。

恢復：

```bash
bash scripts/run_snap.sh \
  --config config.yaml \
  --profile workstation \
  --mode run
```

profile 與 wrapper 已預設/加入 `--rerun-incomplete`。

---

# 17. Debugging

最新 run：

```bash
less /path/to/project/logs/runs/latest/run.log
```

單一 rule：

```bash
less /path/to/project/logs/SM001.haplotypecaller.log
```

需要完整 shell command 的 dry-run：

```bash
snakemake \
  -s Snakefile \
  --configfile config.yaml \
  --workflow-profile profiles/workstation \
  all \
  --dry-run \
  --printshellcmds
```

正常執行不建議 `--printshellcmds`，避免 tmux/terminal 被大量 command 洗掉。

---

# 18. 核心分析邏輯（v1.0.2 未改變）

```text
FASTQ
  ↓
Trimmomatic
  ↓
BWA-MEM
  ↓
samtools sort/index
  ↓
GATK MarkDuplicates
  └─ duplicate flags retained; reads are NOT removed
  ↓
merge technical units by biological sample
  ↓
GATK HaplotypeCaller
  ├─ configurable ploidy
  ├─ -ERC GVCF
  └─ NotDuplicateReadFilter disabled
  ↓
GenomicsDBImport
  ↓
GenotypeGVCFs
  ↓
SNP / INDEL selection
  ↓
CHROM_POS SNP IDs
  ↓
hard filtering
```

SNaP v1.0.2 的主要更新是 portability、resource safety、input flexibility、reference control 與 usability，而不是改變上述 SNP calling 科學流程。
