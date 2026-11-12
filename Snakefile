# SNaP — Snakemake for NGS and Polyploids — v1.0.1
# First formal release. Final SNP IDs use CHROM_POS (e.g. LG1_123456).
import glob
import os
import re
import shutil
import shlex
import sys
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime
from pathlib import Path

from snakemake.shell import shell

configfile: "config.yaml"
sys.path.insert(0, workflow.basedir)
from workflow_inputs import (reference_inputs, discover_marked, discover_gvcfs,
                             selected_samples, safe_id, directories)


############################################
# Configuration helpers
############################################

def cfg(section, key, default=None):
    return config.get(section, {}).get(key, default)


def as_abs(path):
    return os.path.abspath(os.path.expanduser(os.path.expandvars(str(path))))


WORKFLOW_STAGE = str(config.get("workflow_stage", "preprocess")).lower()
VALID_STAGES = {"preprocess", "gvcf", "joint", "joint_shard", "joint_finalize", "full"}
if WORKFLOW_STAGE not in VALID_STAGES:
    raise ValueError(
        "config workflow_stage must be one of: preprocess, gvcf, joint, "
        "joint_shard, joint_finalize, full"
    )

READ_LAYOUT = str(config.get("read_layout", "PE")).upper()
if READ_LAYOUT not in {"SE", "PE"}:
    raise ValueError("config read_layout must be either 'SE' or 'PE'.")

PLOIDY = int(config.get("ploidy", 2))
if PLOIDY < 1:
    raise ValueError("config ploidy must be >= 1.")

BATCH_ID = str(config.get("batch_id", "batch1"))
if not re.fullmatch(r"[A-Za-z0-9._-]+", BATCH_ID):
    raise ValueError(
        "batch_id may contain only letters, numbers, '.', '_' and '-'. "
        "Use a short unique ID for each independently processed sequencing batch."
    )

UNIT_DELIMITER = "__UNIT__"

############################################
# Tool paths
############################################

JAVA        = cfg("tools", "java", "java")
BWA         = cfg("tools", "bwa", "bwa")
SAMTOOLS    = cfg("tools", "samtools", "samtools")
BCFTOOLS    = cfg("tools", "bcftools", "bcftools")
GATK        = cfg("tools", "gatk", "gatk")
TRIMMOMATIC = cfg("tools", "trimmomatic", "trimmomatic")

############################################
# Data paths
############################################

path_start    = as_abs(config["project_dir"])
reference_dir = as_abs(config["reference"])
fastq_dir     = as_abs(config.get("fastq_dir", os.path.join(path_start, "rawData")))

############################################
# Output directories
############################################

trim_dir      = os.path.join(path_start, "trim")
mapped_dir    = os.path.join(path_start, "mapped")
sorted_dir    = os.path.join(path_start, "sorted")
marked_dir    = os.path.join(path_start, "marked")
metadata_dir  = os.path.join(path_start, "metadata")
merged_dir    = os.path.join(path_start, "merged")
gvcf_dir      = os.path.join(path_start, "gvcf_file")
results_dir   = os.path.join(path_start, "results")
log_dir       = os.path.join(path_start, "logs")
benchmark_dir = os.path.join(path_start, "benchmarks")

# Validate stage inputs before adapter discovery or output-directory creation.
if WORKFLOW_STAGE == "gvcf":
    directories(config, "marked_dirs")
elif WORKFLOW_STAGE in {"joint", "joint_shard"}:
    if cfg("joint", "marked_dirs", []):
        raise ValueError("Legacy joint.marked_dirs is no longer accepted; use gvcf + marked_dirs, then joint/joint_shard + gvcf_dirs")
    directories(config, "gvcf_dirs")
elif WORKFLOW_STAGE == "joint_finalize":
    directories(config, "per_contig_dirs")
REFERENCE_INPUTS, REFERENCE_CONTIGS = reference_inputs(reference_dir, WORKFLOW_STAGE)

############################################
# Optional notifications
############################################

ENABLE_TELEGRAM = bool(cfg("telegram", "enabled", False))

# Values set directly in config take precedence. Empty config values fall back
# to environment variables so secrets do not need to be stored in the workflow.
TELEGRAM_BOT_TOKEN = (
    str(cfg("telegram", "bot_token", "") or "").strip()
    or str(os.environ.get("TG_BOT_TOKEN", "") or "").strip()
)
TELEGRAM_CHAT_ID = (
    str(cfg("telegram", "chat_id", "") or "").strip()
    or str(os.environ.get("TG_CHAT_ID", "") or "").strip()
)

TELEGRAM_MISSING_SETTINGS = []
if not TELEGRAM_BOT_TOKEN:
    TELEGRAM_MISSING_SETTINGS.append("telegram.bot_token (or TG_BOT_TOKEN)")
if not TELEGRAM_CHAT_ID:
    TELEGRAM_MISSING_SETTINGS.append("telegram.chat_id (or TG_CHAT_ID)")

TELEGRAM_READY = ENABLE_TELEGRAM and not TELEGRAM_MISSING_SETTINGS
if ENABLE_TELEGRAM and TELEGRAM_MISSING_SETTINGS:
    print(
        "WARNING: Telegram notifications are enabled, but the following "
        "setting(s) are missing: "
        + ", ".join(TELEGRAM_MISSING_SETTINGS)
        + ". Telegram notifications will be skipped; the workflow will continue."
    )

############################################
# Shared parameters (SE/PE/ploidy use same values)
############################################

TRIMMOMATIC_PHRED   = cfg("trimmomatic", "phred", "-phred33")
TRIMMOMATIC_THREADS = int(cfg("trimmomatic", "threads", 8))
LEADING             = int(cfg("trimmomatic", "leading", 3))
TRAILING            = int(cfg("trimmomatic", "trailing", 3))
SLIDINGWINDOW       = str(cfg("trimmomatic", "slidingwindow", "4:30"))
MINLEN              = int(cfg("trimmomatic", "minlen", 36))
ILLUMINACLIP_SE     = str(cfg("trimmomatic", "illuminaclip_se", "2:30:10"))
ILLUMINACLIP_PE     = str(cfg("trimmomatic", "illuminaclip_pe", "2:30:10"))

BWA_THREADS           = int(cfg("resources", "bwa_threads", 12))
SAMTOOLS_VIEW_THREADS = int(cfg("resources", "samtools_view_threads", 4))
SAMTOOLS_SORT_THREADS = int(cfg("resources", "samtools_sort_threads", 8))
HC_THREADS            = int(cfg("resources", "haplotypecaller_threads", 2))
MIN_MAPPING_QUALITY   = int(cfg("calling", "minimum_mapping_quality", 30))

MARKDUP_JAVA_OPTS  = str(cfg("java_options", "markduplicates", "-Xmx20g"))
HC_JAVA_OPTS       = str(cfg("java_options", "haplotypecaller", "-Xmx20g"))
GENOMICSDB_JAVA_OPTS = str(cfg("java_options", "genomicsdbimport", "-Xmx16g -XX:ParallelGCThreads=2"))
GENOTYPE_JAVA_OPTS = str(cfg("java_options", "genotypegvcfs", "-Xmx32g -XX:ParallelGCThreads=2"))
GENOMICSDB_BATCH_SIZE = int(cfg("genomicsdb", "batch_size", 30))
GENOMICSDB_READER_THREADS = int(cfg("genomicsdb", "reader_threads", 1))
GENOMICSDB_MEM_MB = int(cfg("genomicsdb", "import_mem_mb", 24000))
GENOTYPE_MEM_MB = int(cfg("genomicsdb", "genotype_mem_mb", 40000))
MAX_ALTERNATE_ALLELES = int(cfg("calling", "max_alternate_alleles", 6))
MAX_GENOTYPE_COUNT = int(cfg("calling", "max_genotype_count", 1024))

VARIANT_FILTER = str(cfg(
    "calling",
    "hard_filter_expression",
    "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
))

############################################
# Trimmomatic adapter resolution
############################################

BUILTIN_ADAPTER_SE = str(cfg("trimmomatic", "builtin_adapter_se", "NexteraPE-PE.fa"))
BUILTIN_ADAPTER_PE = str(cfg("trimmomatic", "builtin_adapter_pe", "NexteraPE-PE.fa"))
USER_ADAPTER_SE = cfg("trimmomatic", "user_adapter_se", None)
USER_ADAPTER_PE = cfg("trimmomatic", "user_adapter_pe", None)


def trimmomatic_adapter_dirs():
    """Return likely adapter directories, preferring the active installation."""
    candidates = []

    env_adapter_dir = os.environ.get("TRIMMOMATIC_ADAPTER_DIR")
    if env_adapter_dir:
        candidates.append(as_abs(env_adapter_dir))

    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        candidates.extend(glob.glob(os.path.join(conda_prefix, "share", "trimmomatic*", "adapters")))
        candidates.append(os.path.join(conda_prefix, "share", "trimmomatic", "adapters"))

    resolved = shutil.which(TRIMMOMATIC) if os.path.sep not in str(TRIMMOMATIC) else as_abs(TRIMMOMATIC)
    if resolved:
        exe = Path(resolved).resolve()
        prefixes = [exe.parent, exe.parent.parent, exe.parent.parent.parent]
        for prefix in prefixes:
            candidates.extend(glob.glob(str(prefix / "share" / "trimmomatic*" / "adapters")))
            candidates.append(str(prefix / "share" / "trimmomatic" / "adapters"))

    candidates.extend([
        "/usr/share/trimmomatic/adapters",
        "/usr/local/share/trimmomatic/adapters",
        "/opt/trimmomatic/adapters",
    ])

    seen = set()
    unique = []
    for d in candidates:
        d = as_abs(d)
        if d not in seen:
            seen.add(d)
            unique.append(d)
    return unique


def resolve_adapter(layout):
    builtin_name = BUILTIN_ADAPTER_SE if layout == "SE" else BUILTIN_ADAPTER_PE
    user_path = USER_ADAPTER_SE if layout == "SE" else USER_ADAPTER_PE

    builtin_direct = as_abs(builtin_name)
    if os.path.isfile(builtin_direct):
        return builtin_direct

    searched = []
    for adapter_dir in trimmomatic_adapter_dirs():
        candidate = os.path.join(adapter_dir, builtin_name)
        searched.append(candidate)
        if os.path.isfile(candidate):
            return candidate

    if user_path:
        user_path = as_abs(user_path)
        if os.path.isfile(user_path):
            return user_path
        raise FileNotFoundError(
            f"Trimmomatic built-in adapter '{builtin_name}' was not found, and "
            f"configured user_adapter_{layout.lower()} does not exist: {user_path}\n"
            f"Built-in locations checked:\n  " + "\n  ".join(searched)
        )

    raise FileNotFoundError(
        f"Cannot find Trimmomatic built-in adapter '{builtin_name}' for {layout}, "
        f"and no user_adapter_{layout.lower()} fallback was configured.\n"
        f"Locations checked:\n  " + "\n  ".join(searched)
    )


# Gvcf and joint servers do not need Trimmomatic.
TRIMMOMATIC_ADAPTER = None
if WORKFLOW_STAGE in {"preprocess", "full"}:
    TRIMMOMATIC_ADAPTER = resolve_adapter(READ_LAYOUT)

############################################
# FASTQ/unit discovery and biological-sample grouping
############################################

FASTQ_RE = re.compile(r"\.(fastq|fq)\.gz$", flags=re.IGNORECASE)
R1_RE = re.compile(r"^(?P<unit>.+)_R1(?P<lane>_001)?(?P<ext>\.(?:fastq|fq)\.gz)$", flags=re.IGNORECASE)
R2_RE = re.compile(r"^(?P<unit>.+)_R2(?P<lane>_001)?(?P<ext>\.(?:fastq|fq)\.gz)$", flags=re.IGNORECASE)

SAMPLE_GROUP_REGEX = config.get("sample_group_regex", None)
BATCH_TOKEN_RE = re.compile(r"^(?=.*[A-Za-z])(?=.*\d)[A-Za-z0-9]{9,16}$")


def strip_fastq_suffix(filename):
    return FASTQ_RE.sub("", filename)


def get_sm_name(raw_unit):
    """Convert one FASTQ unit name to the biological sample (GATK SM)."""
    if SAMPLE_GROUP_REGEX:
        m = re.search(str(SAMPLE_GROUP_REGEX), raw_unit)
        if not m or m.lastindex is None:
            raise ValueError(
                f"sample_group_regex did not capture a biological sample for '{raw_unit}': "
                f"{SAMPLE_GROUP_REGEX}"
            )
        return safe_id(m.group(1))

    # Existing naming convention: S_TARI0166__2_H9 -> S_TARI0166
    if "__" in raw_unit:
        return raw_unit.split("__", 1)[0]

    # Standard Illumina lane suffixes.
    sm = re.sub(r"_S\d+_L\d{3}$", "", raw_unit, flags=re.IGNORECASE)
    sm = re.sub(r"_L\d{3}$", "", sm, flags=re.IGNORECASE)

    # Sequencing-run suffixes such as:
    # SM001_22HNHHLT4 -> SM001
    # SM001_22HNTCLT4 -> SM001
    # Only 9-16 character alphanumeric tokens containing BOTH letters and digits
    # are stripped automatically, reducing the chance of collapsing ordinary IDs.
    if "_" in sm:
        prefix, suffix = sm.rsplit("_", 1)
        if BATCH_TOKEN_RE.fullmatch(suffix):
            sm = prefix

    return safe_id(sm)


def discover_fastqs(fq_dir, layout):
    """
    Return raw_unit -> FASTQ path(s).

    SE accepts either:
      SM001.fastq.gz
      SM001_R1.fastq.gz
      SM001_R1_001.fastq.gz
    PE requires matching R1/R2 pairs.
    """
    if not os.path.isdir(fq_dir):
        raise FileNotFoundError(f"FASTQ directory does not exist: {fq_dir}")

    entries = sorted(
        f for f in os.listdir(fq_dir)
        if os.path.isfile(os.path.join(fq_dir, f)) and FASTQ_RE.search(f)
    )

    if layout == "SE":
        found = {}
        for f in entries:
            # Never use an R2 file as an SE input.
            if R2_RE.match(f):
                continue

            m = R1_RE.match(f)
            raw_unit = m.group("unit") if m else strip_fastq_suffix(f)
            if raw_unit in found:
                raise ValueError(
                    f"Multiple SE FASTQs resolve to the same unit '{raw_unit}' in {fq_dir}: "
                    f"{found[raw_unit]} and {os.path.join(fq_dir, f)}"
                )
            found[raw_unit] = os.path.join(fq_dir, f)
        return found

    found = {}
    for f in entries:
        m = R1_RE.match(f)
        if not m:
            continue

        raw_unit = m.group("unit")
        lane = m.group("lane") or ""
        ext = m.group("ext")
        r1 = os.path.join(fq_dir, f)
        r2 = os.path.join(fq_dir, f"{raw_unit}_R2{lane}{ext}")

        if not os.path.isfile(r2):
            raise FileNotFoundError(f"R1 exists but matching R2 is missing: {r1} -> {r2}")
        if raw_unit in found:
            raise ValueError(f"Duplicate PE FASTQ unit '{raw_unit}' in {fq_dir}")
        found[raw_unit] = (r1, r2)
    return found


def make_unit_key(sample, raw_unit):
    """
    Produce a globally unique output/read-group key.

    batch_id prevents collisions when two independent sequencing batches both use
    a generic filename such as SM001_R1.fastq.gz.
    """
    if UNIT_DELIMITER in sample or UNIT_DELIMITER in raw_unit:
        raise ValueError(
            f"Sample/unit names may not contain reserved delimiter '{UNIT_DELIMITER}': "
            f"sample={sample}, unit={raw_unit}"
        )
    safe_id(sample)
    safe_id(raw_unit)
    return f"{sample}{UNIT_DELIMITER}{BATCH_ID}__{raw_unit}"


RAW_FASTQ_MAP = {}
UNIT_FASTQ_MAP = {}
UNIT_META = {}
local_units = []
skipped_units = []

if WORKFLOW_STAGE in {"preprocess", "full"}:
    RAW_FASTQ_MAP = discover_fastqs(fastq_dir, READ_LAYOUT)

    retained_samples = set(selected_samples(config, {get_sm_name(u) for u in RAW_FASTQ_MAP}))

    # Backward-compatible name: exclude_samples is interpreted as exact raw-unit keys.
    EXCLUDE_UNITS = set(str(x) for x in config.get("exclude_units", config.get("exclude_samples", [])))
    EXCLUDE_BIOLOGICAL_SAMPLES = set(str(x) for x in config.get("exclude_biological_samples", []))

    for raw_unit in sorted(RAW_FASTQ_MAP):
        sm = get_sm_name(raw_unit)
        unit_key = make_unit_key(sm, raw_unit)

        if raw_unit in EXCLUDE_UNITS or unit_key in EXCLUDE_UNITS or sm not in retained_samples:
            skipped_units.append(raw_unit)
            continue

        if unit_key in UNIT_META:
            raise ValueError(f"Duplicate generated unit key: {unit_key}")

        UNIT_META[unit_key] = {
            "raw_unit": raw_unit,
            "sample": sm,
            "batch": BATCH_ID,
        }
        UNIT_FASTQ_MAP[unit_key] = RAW_FASTQ_MAP[raw_unit]
        local_units.append(unit_key)

    if not local_units:
        raise ValueError(f"No {READ_LAYOUT} FASTQ units were detected/retained in: {fastq_dir}")

############################################
# Stage-isolated discovery
############################################
joint_bams_by_sample = {}
joint_bais_by_sample = {}
external_gvcfs = {}
if WORKFLOW_STAGE == "gvcf":
    joint_bams_by_sample, joint_bais_by_sample = discover_marked(config, get_sm_name, REFERENCE_CONTIGS)
elif WORKFLOW_STAGE in {"joint", "joint_shard"}:
    external_gvcfs = discover_gvcfs(config, REFERENCE_CONTIGS)
elif WORKFLOW_STAGE == "full":
    for unit in local_units:
        sm = UNIT_META[unit]["sample"]
        bam = os.path.join(marked_dir, unit + ".marked.bam")
        joint_bams_by_sample.setdefault(sm, []).append(bam)
        joint_bais_by_sample.setdefault(sm, []).append(bam + ".bai")
    external_gvcfs = {sm: os.path.join(gvcf_dir, sm + ".g.vcf.gz")
                      for sm in sorted(joint_bams_by_sample)}

sample_names = sorted(
    joint_bams_by_sample if WORKFLOW_STAGE in {"gvcf", "full"} else external_gvcfs
)

# Reference-order shard IDs are stable across machines and are used only as filenames.
CONTIG_RECORDS = [
    {"shard": f"{i:06d}", "contig": name, "length": length}
    for i, (name, length) in enumerate(REFERENCE_CONTIGS, 1)
]
CONTIG_BY_SHARD = {x["shard"]: x for x in CONTIG_RECORDS}
SHARD_BY_CONTIG = {x["contig"]: x["shard"] for x in CONTIG_RECORDS}


def read_contig_assignment(path):
    if not path:
        return {}
    ap = as_abs(path)
    if not os.path.isfile(ap):
        raise ValueError(f"joint_contig_assignment_file does not exist: {ap}")
    assignment = {}
    with open(ap) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"\t+|\s+", line)
            if parts[0].lower() in {"contig", "chromosome", "lg"}:
                continue
            if len(parts) < 2:
                raise ValueError(f"Malformed contig assignment line {lineno}: {line}")
            contig, partition = parts[0], parts[1]
            if contig not in SHARD_BY_CONTIG:
                raise ValueError(f"Contig in assignment is absent from reference .fai: {contig}")
            if contig in assignment:
                raise ValueError(f"Duplicate contig in joint_contig_assignment_file: {contig}")
            assignment[contig] = partition
    missing = set(SHARD_BY_CONTIG) - set(assignment)
    if missing:
        raise ValueError("Contigs missing from joint_contig_assignment_file: " + ", ".join(sorted(missing)))
    return assignment


JOINT_PARTITION = str(config.get("joint_partition", "") or "").strip()
CONTIG_ASSIGNMENT = {}
if WORKFLOW_STAGE == "joint_shard":
    assignment_file = config.get("joint_contig_assignment_file")
    if not assignment_file or not JOINT_PARTITION:
        raise ValueError("joint_shard requires joint_contig_assignment_file and joint_partition")
    CONTIG_ASSIGNMENT = read_contig_assignment(assignment_file)
    ACTIVE_CONTIG_RECORDS = [
        x for x in CONTIG_RECORDS if CONTIG_ASSIGNMENT[x["contig"]] == JOINT_PARTITION
    ]
    if not ACTIVE_CONTIG_RECORDS:
        raise ValueError(f"No contigs are assigned to joint_partition={JOINT_PARTITION}")
elif WORKFLOW_STAGE in {"joint", "full"}:
    ACTIVE_CONTIG_RECORDS = list(CONTIG_RECORDS)
else:
    ACTIVE_CONTIG_RECORDS = []
ACTIVE_SHARDS = [x["shard"] for x in ACTIVE_CONTIG_RECORDS]

# joint_finalize receives completed per-contig VCFs from one or more joint_shard servers.
external_per_contig = {}
if WORKFLOW_STAGE == "joint_finalize":
    for d in directories(config, "per_contig_dirs"):
        for p in sorted(Path(d).glob("*.vcf.gz")):
            if p.name.endswith(".vcf.gz.tbi"):
                continue
            shard = p.name[:-len(".vcf.gz")]
            if shard not in CONTIG_BY_SHARD:
                continue
            ap = str(p.resolve())
            if shard in external_per_contig:
                raise ValueError(f"Duplicate per-contig VCF for shard {shard}: {external_per_contig[shard]} and {ap}")
            if not os.path.isfile(ap + ".tbi"):
                raise ValueError(f"Missing index for per-contig VCF: {ap}.tbi")
            external_per_contig[shard] = ap
    missing = [x["shard"] for x in CONTIG_RECORDS if x["shard"] not in external_per_contig]
    if missing:
        raise ValueError("joint_finalize is missing per-contig VCF shard(s): " + ", ".join(missing))

stage_dirs = {
    "preprocess": [trim_dir, mapped_dir, sorted_dir, marked_dir],
    "full": [trim_dir, mapped_dir, sorted_dir, marked_dir, merged_dir, gvcf_dir, results_dir],
    "gvcf": [merged_dir, gvcf_dir],
    "joint": [results_dir],
    "joint_shard": [results_dir],
    "joint_finalize": [results_dir],
}
for d in stage_dirs[WORKFLOW_STAGE] + [metadata_dir, log_dir, benchmark_dir]:
    os.makedirs(d, exist_ok=True)

GENOMICSDB_DIR = os.path.join(path_start, "genomicsdb")
PER_CONTIG_DIR = os.path.join(results_dir, "per_contig")
if WORKFLOW_STAGE in {"joint", "joint_shard", "full"}:
    os.makedirs(GENOMICSDB_DIR, exist_ok=True)
    os.makedirs(PER_CONTIG_DIR, exist_ok=True)

############################################
# Console summary
############################################

print(f"Workflow stage: {WORKFLOW_STAGE}")
print(f"Read layout: {READ_LAYOUT}")
print(f"Ploidy: {PLOIDY}")

if WORKFLOW_STAGE in {"preprocess", "full"}:
    print(f"Batch ID: {BATCH_ID}")
    print(f"Trimmomatic adapter: {TRIMMOMATIC_ADAPTER}")
    print(f"Detected FASTQ units: {len(RAW_FASTQ_MAP)}")
    print(f"Kept FASTQ units: {len(local_units)}")
    print(f"Skipped FASTQ units: {len(skipped_units)}")
    print("Biological samples represented locally:", len(set(UNIT_META[u]["sample"] for u in local_units)))

if WORKFLOW_STAGE == "gvcf":
    print(f"Biological samples entering HaplotypeCaller: {len(sample_names)}")
    print(f"Marked BAM units entering sample merge: {sum(len(v) for v in joint_bams_by_sample.values())}")
    print("Selected biological samples:", ", ".join(sample_names))
if WORKFLOW_STAGE in {"joint", "joint_shard"}:
    print(f"Single-sample gVCFs entering joint genotyping: {len(sample_names)}")
    print(f"Active joint contigs: {len(ACTIVE_CONTIG_RECORDS)}")
    if WORKFLOW_STAGE == "joint_shard":
        print(f"Joint partition: {JOINT_PARTITION}")
        print("Assigned contigs:", ", ".join(x["contig"] for x in ACTIVE_CONTIG_RECORDS))
if WORKFLOW_STAGE == "joint_finalize":
    print(f"Per-contig VCFs entering finalize: {len(external_per_contig)}")

############################################
# FASTQ path helpers
############################################


def se_fastq(wc):
    return UNIT_FASTQ_MAP[wc.unit]


def pe_r1(wc):
    return UNIT_FASTQ_MAP[wc.unit][0]


def pe_r2(wc):
    return UNIT_FASTQ_MAP[wc.unit][1]


def get_unit_sm(unit_key):
    return UNIT_META[unit_key]["sample"]

############################################
# Telegram helper
############################################


def send_telegram(message):
    if not TELEGRAM_READY:
        return

    url = f"https://api.telegram.org/bot{TELEGRAM_BOT_TOKEN}/sendMessage"
    data = urllib.parse.urlencode({
        "chat_id": TELEGRAM_CHAT_ID,
        "text": message
    }).encode("utf-8")

    try:
        req = urllib.request.Request(url, data=data, method="POST")
        with urllib.request.urlopen(req, timeout=20) as resp:
            _ = resp.read()
    except Exception as e:
        print("WARNING: Telegram send failed (details suppressed to protect credentials).")


onstart:
    msg = (
        f"[Snakemake] SNP pipeline ({WORKFLOW_STAGE}) started at "
        f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
    )
    print(msg)
    send_telegram(msg)


onsuccess:
    if WORKFLOW_STAGE == "preprocess":
        result_note = os.path.join(metadata_dir, f"{BATCH_ID}.units.tsv")
    elif WORKFLOW_STAGE == "gvcf":
        result_note = gvcf_dir
    elif WORKFLOW_STAGE == "joint_shard":
        result_note = PER_CONTIG_DIR
    else:
        result_note = os.path.join(results_dir, "final_combined_filtered_snps.vcf.gz")
    msg = (
        f"[Snakemake] SNP pipeline ({WORKFLOW_STAGE}) completed successfully at "
        f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        f"Output: {result_note}"
    )
    print(msg)
    send_telegram(msg)


onerror:
    msg = (
        f"[Snakemake] SNP pipeline ({WORKFLOW_STAGE}) failed at "
        f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        f"Please check logs under: {log_dir}"
    )
    print(msg)
    send_telegram(msg)

############################################
# Stage-aware targets
############################################


def trim_targets(wc):
    if READ_LAYOUT == "SE":
        return expand(os.path.join(trim_dir, "{unit}.trim.fastq.gz"), unit=local_units)
    return (
        expand(os.path.join(trim_dir, "{unit}_R1_paired.trim.fastq.gz"), unit=local_units)
        + expand(os.path.join(trim_dir, "{unit}_R2_paired.trim.fastq.gz"), unit=local_units)
    )


def trimmed_reads(wc):
    if READ_LAYOUT == "SE":
        return [os.path.join(trim_dir, f"{wc.unit}.trim.fastq.gz")]
    return [
        os.path.join(trim_dir, f"{wc.unit}_R1_paired.trim.fastq.gz"),
        os.path.join(trim_dir, f"{wc.unit}_R2_paired.trim.fastq.gz"),
    ]


def preprocess_targets(wc):
    targets = []
    targets += expand(os.path.join(marked_dir, "{unit}.marked.bam"), unit=local_units)
    targets += expand(os.path.join(marked_dir, "{unit}.marked.bam.bai"), unit=local_units)
    targets += expand(os.path.join(marked_dir, "{unit}.marked_dup_metrics.txt"), unit=local_units)
    if local_units:
        targets.append(os.path.join(metadata_dir, f"{BATCH_ID}.units.tsv"))
    return targets


def final_variant_targets(wc):
    return [
        os.path.join(results_dir, "final_combined_filtered_snps.vcf.gz"),
        os.path.join(results_dir, "final_combined_filtered_snps.vcf.gz.tbi"),
        os.path.join(results_dir, "final_combined_filtered_indels.vcf.gz"),
        os.path.join(results_dir, "final_combined_filtered_indels.vcf.gz.tbi"),
        os.path.join(results_dir, "final_combined_filtered_snps.count.txt"),
    ]


def per_contig_targets(records=None):
    records = ACTIVE_CONTIG_RECORDS if records is None else records
    shards = [x["shard"] for x in records]
    return (
        expand(os.path.join(PER_CONTIG_DIR, "{shard}.vcf.gz"), shard=shards)
        + expand(os.path.join(PER_CONTIG_DIR, "{shard}.vcf.gz.tbi"), shard=shards)
    )


def all_targets(wc):
    if WORKFLOW_STAGE == "full":
        return preprocess_targets(wc) + final_variant_targets(wc) + [os.path.join(metadata_dir, "full_inputs.tsv")]
    if WORKFLOW_STAGE == "preprocess":
        return preprocess_targets(wc)
    if WORKFLOW_STAGE == "gvcf":
        return (expand(os.path.join(gvcf_dir, "{sm}.g.vcf.gz"), sm=sample_names)
                + expand(os.path.join(gvcf_dir, "{sm}.g.vcf.gz.tbi"), sm=sample_names)
                + [os.path.join(metadata_dir, "gvcf_inputs.tsv")])
    if WORKFLOW_STAGE == "joint_shard":
        return per_contig_targets() + [
            os.path.join(metadata_dir, "joint_shard_inputs.tsv"),
            os.path.join(metadata_dir, "joint_contigs.tsv"),
            os.path.join(metadata_dir, "genomicsdb.sample_map.tsv"),
        ]
    if WORKFLOW_STAGE == "joint_finalize":
        return final_variant_targets(wc) + [os.path.join(metadata_dir, "joint_finalize_inputs.tsv")]
    return final_variant_targets(wc) + [os.path.join(metadata_dir, f"{WORKFLOW_STAGE}_inputs.tsv")]


rule all:
    input:
        all_targets


if WORKFLOW_STAGE in {"preprocess", "full"}:
    rule preprocess_all:
        input:
            preprocess_targets


    rule trim_all:
        input:
            trim_targets


    rule map_all:
        input:
            lambda wc: expand(os.path.join(mapped_dir, "{unit}.bam"), unit=local_units)


    rule sort_all:
        input:
            lambda wc: (
                expand(os.path.join(sorted_dir, "{unit}.sorted.bam"), unit=local_units)
                + expand(os.path.join(sorted_dir, "{unit}.sorted.bam.bai"), unit=local_units)
            )


    rule mark_all:
        input:
            preprocess_targets


if WORKFLOW_STAGE in {"gvcf", "full"}:
    rule merge_all:
        input:
            lambda wc: (
                expand(os.path.join(merged_dir, "{sm}.merged.bam"), sm=sample_names)
                + expand(os.path.join(merged_dir, "{sm}.merged.bam.bai"), sm=sample_names)
            )


    rule call_all:
        input:
            lambda wc: (
                expand(os.path.join(gvcf_dir, "{sm}.g.vcf.gz"), sm=sample_names)
                + expand(os.path.join(gvcf_dir, "{sm}.g.vcf.gz.tbi"), sm=sample_names)
            )


if WORKFLOW_STAGE in {"joint", "full", "joint_finalize"}:
    rule joint_all:
        input:
            final_variant_targets


if WORKFLOW_STAGE == "joint_shard":
    rule joint_shard_all:
        input:
            per_contig_targets

############################################
# Trimming (layout-specific) + per-unit mapping
############################################

if WORKFLOW_STAGE in {"preprocess", "full"}:
    rule trim_fastq_se:
        input:
            fq = se_fastq
        output:
            fq = os.path.join(trim_dir, "{unit}.trim.fastq.gz")
        params:
            adapter = TRIMMOMATIC_ADAPTER,
            clip = ILLUMINACLIP_SE
        log:
            os.path.join(log_dir, "{unit}.trimmomatic_se.log")
        benchmark:
            os.path.join(benchmark_dir, "{unit}.trimmomatic_se.txt")
        threads: TRIMMOMATIC_THREADS
        shell:
            r"""
            {TRIMMOMATIC:q} SE {TRIMMOMATIC_PHRED} -threads {threads} \
                {input.fq:q} \
                {output.fq:q} \
                ILLUMINACLIP:{params.adapter:q}:{params.clip} \
                LEADING:{LEADING} TRAILING:{TRAILING} SLIDINGWINDOW:{SLIDINGWINDOW} MINLEN:{MINLEN} \
                > {log:q} 2>&1
            """


    rule trim_fastq_pe:
        input:
            R1 = pe_r1,
            R2 = pe_r2
        output:
            R1P = os.path.join(trim_dir, "{unit}_R1_paired.trim.fastq.gz"),
            R1U = os.path.join(trim_dir, "{unit}_R1_unpaired.trim.fastq.gz"),
            R2P = os.path.join(trim_dir, "{unit}_R2_paired.trim.fastq.gz"),
            R2U = os.path.join(trim_dir, "{unit}_R2_unpaired.trim.fastq.gz")
        params:
            adapter = TRIMMOMATIC_ADAPTER,
            clip = ILLUMINACLIP_PE
        log:
            os.path.join(log_dir, "{unit}.trimmomatic_pe.log")
        benchmark:
            os.path.join(benchmark_dir, "{unit}.trimmomatic_pe.txt")
        threads: TRIMMOMATIC_THREADS
        shell:
            r"""
            {TRIMMOMATIC:q} PE {TRIMMOMATIC_PHRED} -threads {threads} \
                {input.R1:q} {input.R2:q} \
                {output.R1P:q} {output.R1U:q} \
                {output.R2P:q} {output.R2U:q} \
                ILLUMINACLIP:{params.adapter:q}:{params.clip} \
                LEADING:{LEADING} TRAILING:{TRAILING} SLIDINGWINDOW:{SLIDINGWINDOW} MINLEN:{MINLEN} \
                > {log:q} 2>&1
            """


    rule bwa_mapping:
        input:
            reads = trimmed_reads,
            reference_files = REFERENCE_INPUTS
        params:
            reference = reference_dir,
            rg = lambda wc: "'@RG\\tID:{unit}\\tSM:{sm}\\tLB:{sm}_{batch}\\tPL:ILLUMINA\\tPU:{batch}'".format(
                unit=wc.unit,
                sm=get_unit_sm(wc.unit),
                batch=BATCH_ID
            )
        output:
            bam = os.path.join(mapped_dir, "{unit}.bam")
        log:
            os.path.join(log_dir, "{unit}.bwa.log")
        benchmark:
            os.path.join(benchmark_dir, "{unit}.bwa.txt")
        threads: BWA_THREADS
        shell:
            r"""
            {BWA:q} mem -M -t {threads} -R {params.rg} {params.reference:q} {input.reads:q} \
            2> {log:q} \
            | {SAMTOOLS:q} view -@ {SAMTOOLS_VIEW_THREADS} -Sb - \
            > {output.bam:q}
            """

    ############################################
    # Per-unit BAM processing
    ############################################

    rule samtools_sort_index:
        input:
            os.path.join(mapped_dir, "{unit}.bam")
        output:
            bam = os.path.join(sorted_dir, "{unit}.sorted.bam"),
            bai = os.path.join(sorted_dir, "{unit}.sorted.bam.bai")
        log:
            os.path.join(log_dir, "{unit}.sort_index.log")
        benchmark:
            os.path.join(benchmark_dir, "{unit}.sort_index.txt")
        threads: SAMTOOLS_SORT_THREADS
        shell:
            r"""
            {SAMTOOLS:q} sort -@ {threads} {input:q} -o {output.bam:q} > {log:q} 2>&1
            {SAMTOOLS:q} quickcheck -v {output.bam:q} >> {log:q} 2>&1
            {SAMTOOLS:q} index {output.bam:q} {output.bai:q} >> {log:q} 2>&1
            """


    rule mark_duplicates:
        input:
            bam = os.path.join(sorted_dir, "{unit}.sorted.bam"),
            bai = os.path.join(sorted_dir, "{unit}.sorted.bam.bai")
        params:
            java_option = MARKDUP_JAVA_OPTS
        output:
            bam = os.path.join(marked_dir, "{unit}.marked.bam"),
            bai = os.path.join(marked_dir, "{unit}.marked.bam.bai"),
            metrics = os.path.join(marked_dir, "{unit}.marked_dup_metrics.txt")
        log:
            os.path.join(log_dir, "{unit}.markdup.log")
        benchmark:
            os.path.join(benchmark_dir, "{unit}.markdup.txt")
        shell:
            r"""
            # RAD-seq/DArT-seq: identical coordinates can represent genuine
            # restriction-fragment coverage. Mark for QC, but NEVER remove reads.
            {GATK:q} --java-options "{params.java_option}" MarkDuplicates \
                -I {input.bam:q} \
                -O {output.bam:q} \
                -M {output.metrics:q} \
                --CREATE_INDEX false \
                --REMOVE_DUPLICATES false \
                --REMOVE_SEQUENCING_DUPLICATES false \
                > {log:q} 2>&1

            {SAMTOOLS:q} quickcheck -v {output.bam:q} >> {log:q} 2>&1
            {SAMTOOLS:q} index {output.bam:q} {output.bai:q} >> {log:q} 2>&1
            """


    rule write_unit_manifest:
        input:
            bams = lambda wc: expand(os.path.join(marked_dir, "{unit}.marked.bam"), unit=local_units),
            bais = lambda wc: expand(os.path.join(marked_dir, "{unit}.marked.bam.bai"), unit=local_units)
        output:
            os.path.join(metadata_dir, f"{BATCH_ID}.units.tsv")
        run:
            with open(output[0], "w") as out:
                out.write("unit_key\tbiological_sample\tbatch_id\traw_unit\tread_layout\tmarked_bam\n")
                for unit_key in local_units:
                    meta = UNIT_META[unit_key]
                    out.write(
                        f"{unit_key}\t{meta['sample']}\t{meta['batch']}\t{meta['raw_unit']}\t"
                        f"{READ_LAYOUT}\t{os.path.join(marked_dir, unit_key + '.marked.bam')}\n"
                    )

############################################
# Central sample merge + HaplotypeCaller
############################################

if WORKFLOW_STAGE in {"gvcf", "full"}:
    rule merge_marked_bam:
        input:
            bams = lambda wc: joint_bams_by_sample[wc.sm],
            bais = lambda wc: joint_bais_by_sample[wc.sm]
        output:
            bam = os.path.join(merged_dir, "{sm}.merged.bam"),
            bai = os.path.join(merged_dir, "{sm}.merged.bam.bai")
        log:
            os.path.join(log_dir, "{sm}.merge_marked_bam.log")
        benchmark:
            os.path.join(benchmark_dir, "{sm}.merge_marked_bam.txt")
        run:
            bam_list = list(input.bams)
            shell("{SAMTOOLS:q} quickcheck -v {input.bams:q} > {log:q} 2>&1")

            if len(bam_list) == 1:
                src = bam_list[0]
                shell(r"""
                    cp {src:q} {output.bam:q}
                    {SAMTOOLS:q} index {output.bam:q} {output.bai:q} > {log:q} 2>&1
                """)
            else:
                shell(r"""
                    {SAMTOOLS:q} merge -f {output.bam:q} {input.bams:q} > {log:q} 2>&1
                    {SAMTOOLS:q} quickcheck -v {output.bam:q} >> {log:q} 2>&1
                    {SAMTOOLS:q} index {output.bam:q} {output.bai:q} >> {log:q} 2>&1
                """)


    rule gatk_haplotypecaller:
        input:
            bam = os.path.join(merged_dir, "{sm}.merged.bam"),
            bai = os.path.join(merged_dir, "{sm}.merged.bam.bai"),
            reference_files = REFERENCE_INPUTS
        params:
            java_option = HC_JAVA_OPTS,
            reference = reference_dir,
            ploidy = PLOIDY,
            min_mq = MIN_MAPPING_QUALITY
        output:
            gvcf = os.path.join(gvcf_dir, "{sm}.g.vcf.gz")
        log:
            os.path.join(log_dir, "{sm}.haplotypecaller.log")
        benchmark:
            os.path.join(benchmark_dir, "{sm}.haplotypecaller.txt")
        threads: HC_THREADS
        shell:
            r"""
            # HaplotypeCaller normally applies NotDuplicateReadFilter. Disable that
            # one default filter so SAM duplicate-flagged reads remain in analysis.
            {GATK:q} --java-options "{params.java_option}" HaplotypeCaller \
                -R {params.reference:q} \
                -I {input.bam:q} \
                --sample-ploidy {params.ploidy} \
                --disable-read-filter NotDuplicateReadFilter \
                --minimum-mapping-quality {params.min_mq} \
                --native-pair-hmm-threads {threads} \
                -ERC GVCF \
                --create-output-variant-index false \
                -O {output.gvcf:q} \
                > {log:q} 2>&1
            """


    rule index_gvcf:
        input:
            os.path.join(gvcf_dir, "{sm}.g.vcf.gz")
        output:
            os.path.join(gvcf_dir, "{sm}.g.vcf.gz.tbi")
        log:
            os.path.join(log_dir, "{sm}.index_gvcf.log")
        shell:
            r"""
            {BCFTOOLS:q} index --tbi --force {input:q} > {log:q} 2>&1
            {BCFTOOLS:q} index --stats {input:q} >> {log:q} 2>&1
            """

############################################
# Joint genotyping: GenomicsDB per reference contig
############################################


def contig_for_shard(wc):
    return CONTIG_BY_SHARD[wc.shard]["contig"]


def source_per_contig_vcf(wc):
    if WORKFLOW_STAGE == "joint_finalize":
        return external_per_contig[wc.shard]
    return os.path.join(PER_CONTIG_DIR, wc.shard + ".vcf.gz")


if WORKFLOW_STAGE in {"joint", "joint_shard", "full"}:
    rule write_joint_metadata:
        input:
            gvcfs = list(external_gvcfs.values()),
            tbis = [p + ".tbi" for p in external_gvcfs.values()],
            ref = reference_dir,
            reference_files = REFERENCE_INPUTS
        output:
            sample_map = os.path.join(metadata_dir, "genomicsdb.sample_map.tsv"),
            contigs = os.path.join(metadata_dir, "joint_contigs.tsv")
        run:
            with open(output.sample_map, "w") as out:
                for sm in sorted(external_gvcfs):
                    out.write(f"{sm}\t{external_gvcfs[sm]}\n")
            with open(output.contigs, "w") as out:
                out.write("shard\tcontig\tlength\tpartition\n")
                for rec in ACTIVE_CONTIG_RECORDS:
                    partition = CONTIG_ASSIGNMENT.get(rec["contig"], "ALL")
                    out.write(f"{rec['shard']}\t{rec['contig']}\t{rec['length']}\t{partition}\n")


    rule genomicsdb_import:
        input:
            sample_map = os.path.join(metadata_dir, "genomicsdb.sample_map.tsv"),
            contig_manifest = os.path.join(metadata_dir, "joint_contigs.tsv"),
            gvcfs = list(external_gvcfs.values()),
            tbis = [p + ".tbi" for p in external_gvcfs.values()],
            ref = reference_dir,
            reference_files = REFERENCE_INPUTS
        output:
            directory(os.path.join(GENOMICSDB_DIR, "{shard}"))
        params:
            java_option = GENOMICSDB_JAVA_OPTS,
            contig = contig_for_shard,
            batch_size = GENOMICSDB_BATCH_SIZE,
            reader_threads = GENOMICSDB_READER_THREADS
        resources:
            joint_jobs = 1,
            mem_mb = GENOMICSDB_MEM_MB
        log:
            os.path.join(log_dir, "{shard}.genomicsdbimport.log")
        benchmark:
            os.path.join(benchmark_dir, "{shard}.genomicsdbimport.txt")
        shell:
            r"""
            # GenomicsDBImport requires a fresh workspace on retry. Snakemake tracks
            # this directory as the rule output, so removing only this shard is safe.
            rm -rf {output:q}
            {GATK:q} --java-options "{params.java_option}" GenomicsDBImport \
                -R {input.ref:q} \
                -L {params.contig:q} \
                --sample-name-map {input.sample_map:q} \
                --validate-sample-name-map true \
                --batch-size {params.batch_size} \
                --reader-threads {params.reader_threads} \
                --genomicsdb-workspace-path {output:q} \
                > {log:q} 2>&1
            """


    rule genotype_gvcfs:
        input:
            db = os.path.join(GENOMICSDB_DIR, "{shard}"),
            ref = reference_dir,
            reference_files = REFERENCE_INPUTS
        output:
            vcf = os.path.join(PER_CONTIG_DIR, "{shard}.vcf.gz"),
            tbi = os.path.join(PER_CONTIG_DIR, "{shard}.vcf.gz.tbi")
        params:
            java_option = GENOTYPE_JAVA_OPTS,
            contig = contig_for_shard,
            db_uri = lambda wc: "gendb://" + os.path.join(GENOMICSDB_DIR, wc.shard)
        resources:
            joint_jobs = 1,
            mem_mb = GENOTYPE_MEM_MB
        log:
            os.path.join(log_dir, "{shard}.genotype_gvcfs.log")
        benchmark:
            os.path.join(benchmark_dir, "{shard}.genotype_gvcfs.txt")
        shell:
            r"""
            {GATK:q} --java-options "{params.java_option}" GenotypeGVCFs \
                -R {input.ref:q} \
                -V {params.db_uri:q} \
                -L {params.contig:q} \
                --max-alternate-alleles {MAX_ALTERNATE_ALLELES} \
                --max-genotype-count {MAX_GENOTYPE_COUNT} \
                --create-output-variant-index true \
                -O {output.vcf:q} \
                > {log:q} 2>&1
            """


if WORKFLOW_STAGE == "joint_finalize":
    rule write_finalize_manifest:
        input:
            vcfs = lambda wc: [external_per_contig[x["shard"]] for x in CONTIG_RECORDS],
            tbis = lambda wc: [external_per_contig[x["shard"]] + ".tbi" for x in CONTIG_RECORDS]
        output:
            os.path.join(metadata_dir, "joint_finalize_inputs.tsv")
        run:
            with open(output[0], "w") as out:
                out.write("shard\tcontig\tlength\tvcf\n")
                for rec in CONTIG_RECORDS:
                    out.write(f"{rec['shard']}\t{rec['contig']}\t{rec['length']}\t{external_per_contig[rec['shard']]}\n")


if WORKFLOW_STAGE in {"joint", "full", "joint_finalize"}:
    rule concat_joint_vcfs:
        input:
            vcfs = lambda wc: (
                [external_per_contig[x["shard"]] for x in CONTIG_RECORDS]
                if WORKFLOW_STAGE == "joint_finalize"
                else [os.path.join(PER_CONTIG_DIR, x["shard"] + ".vcf.gz") for x in CONTIG_RECORDS]
            ),
            tbis = lambda wc: (
                [external_per_contig[x["shard"]] + ".tbi" for x in CONTIG_RECORDS]
                if WORKFLOW_STAGE == "joint_finalize"
                else [os.path.join(PER_CONTIG_DIR, x["shard"] + ".vcf.gz.tbi") for x in CONTIG_RECORDS]
            )
        output:
            vcf = os.path.join(results_dir, "allsample.vcf.gz"),
            tbi = os.path.join(results_dir, "allsample.vcf.gz.tbi")
        log:
            os.path.join(log_dir, "concat_joint_vcfs.log")
        shell:
            r"""
            {BCFTOOLS:q} concat -a -Oz -o {output.vcf:q} {input.vcfs:q} > {log:q} 2>&1
            {BCFTOOLS:q} index -t --force {output.vcf:q} >> {log:q} 2>&1
            """


    rule select_snps:
        input:
            joint_vcf = os.path.join(results_dir, "allsample.vcf.gz"),
            tbi = os.path.join(results_dir, "allsample.vcf.gz.tbi"),
            ref = reference_dir,
            reference_files = REFERENCE_INPUTS
        output:
            os.path.join(results_dir, "allsample_snps.vcf.gz")
        log:
            os.path.join(log_dir, "select_snps.log")
        shell:
            r"""
            {GATK:q} SelectVariants \
                -R {input.ref:q} \
                -V {input.joint_vcf:q} \
                --select-type-to-include SNP \
                -O {output:q} \
                > {log:q} 2>&1
            """


    rule select_indels:
        input:
            joint_vcf = os.path.join(results_dir, "allsample.vcf.gz"),
            tbi = os.path.join(results_dir, "allsample.vcf.gz.tbi"),
            ref = reference_dir,
            reference_files = REFERENCE_INPUTS
        output:
            os.path.join(results_dir, "allsample_indels.vcf.gz")
        log:
            os.path.join(log_dir, "select_indels.log")
        shell:
            r"""
            {GATK:q} SelectVariants \
                -R {input.ref:q} \
                -V {input.joint_vcf:q} \
                --select-type-to-include INDEL \
                -O {output:q} \
                > {log:q} 2>&1
            """


    rule set_snp_id:
        input:
            os.path.join(results_dir, "allsample_snps.vcf.gz")
        output:
            vcf = os.path.join(results_dir, "allsample_snps.ID.vcf.gz"),
            tbi = os.path.join(results_dir, "allsample_snps.ID.vcf.gz.tbi")
        log:
            os.path.join(log_dir, "set_snp_id.log")
        shell:
            r"""
            {BCFTOOLS:q} annotate --set-id '%CHROM\\_%POS' {input:q} -Oz -o {output.vcf:q} > {log:q} 2>&1
            {BCFTOOLS:q} index -t --force {output.vcf:q} >> {log:q} 2>&1
            """


    rule set_indel_id:
        input:
            os.path.join(results_dir, "allsample_indels.vcf.gz")
        output:
            vcf = os.path.join(results_dir, "allsample_indels.ID.vcf.gz"),
            tbi = os.path.join(results_dir, "allsample_indels.ID.vcf.gz.tbi")
        log:
            os.path.join(log_dir, "set_indel_id.log")
        shell:
            r"""
            {BCFTOOLS:q} annotate --set-id '%CHROM\\_%POS' {input:q} -Oz -o {output.vcf:q} > {log:q} 2>&1
            {BCFTOOLS:q} index -t --force {output.vcf:q} >> {log:q} 2>&1
            """


    rule bcftools_filter:
        input:
            raw_snps = os.path.join(results_dir, "allsample_snps.ID.vcf.gz"),
            raw_indels = os.path.join(results_dir, "allsample_indels.ID.vcf.gz")
        params:
            filter = VARIANT_FILTER
        output:
            snps = os.path.join(results_dir, "final_combined_filtered_snps.vcf.gz"),
            indels = os.path.join(results_dir, "final_combined_filtered_indels.vcf.gz")
        log:
            os.path.join(log_dir, "bcftools_filter.log")
        shell:
            r"""
            {BCFTOOLS:q} filter -e {params.filter:q} {input.raw_snps:q} -Oz -o {output.snps:q} > {log:q} 2>&1
            {BCFTOOLS:q} filter -e {params.filter:q} {input.raw_indels:q} -Oz -o {output.indels:q} >> {log:q} 2>&1
            """


    rule index_filtered_variants:
        input:
            snps = os.path.join(results_dir, "final_combined_filtered_snps.vcf.gz"),
            indels = os.path.join(results_dir, "final_combined_filtered_indels.vcf.gz")
        output:
            snps_tbi = os.path.join(results_dir, "final_combined_filtered_snps.vcf.gz.tbi"),
            indels_tbi = os.path.join(results_dir, "final_combined_filtered_indels.vcf.gz.tbi")
        log:
            os.path.join(log_dir, "index_filtered_variants.log")
        shell:
            r"""
            {BCFTOOLS:q} index -t --force {input.snps:q} > {log:q} 2>&1
            {BCFTOOLS:q} index -t --force {input.indels:q} >> {log:q} 2>&1
            {BCFTOOLS:q} index --stats {input.snps:q} >> {log:q} 2>&1
            {BCFTOOLS:q} index --stats {input.indels:q} >> {log:q} 2>&1
            """


    rule count_filtered_snps:
        input:
            os.path.join(results_dir, "final_combined_filtered_snps.vcf.gz")
        output:
            os.path.join(results_dir, "final_combined_filtered_snps.count.txt")
        shell:
            r"""
            {BCFTOOLS:q} view -H {input:q} | wc -l > {output:q}
            """

if WORKFLOW_STAGE in {"gvcf", "joint", "joint_shard", "full"}:
    rule write_stage_manifest:
        input:
            files = ([p for paths in joint_bams_by_sample.values() for p in paths]
                     if WORKFLOW_STAGE in {"gvcf", "full"} else list(external_gvcfs.values()))
        output:
            os.path.join(metadata_dir, f"{WORKFLOW_STAGE}_inputs.tsv")
        params:
            partition = str(
                config.get("calling_partition", "")
                if WORKFLOW_STAGE == "gvcf"
                else (JOINT_PARTITION if WORKFLOW_STAGE == "joint_shard" else "ALL")
            )
        run:
            with open(output[0], "w") as out:
                out.write("biological_sample\tinput_path\tpartition\n")
                for sm in sample_names:
                    paths = joint_bams_by_sample[sm] if WORKFLOW_STAGE in {"gvcf", "full"} else [external_gvcfs[sm]]
                    for path in paths:
                        out.write(f"{sm}\t{path}\t{params.partition}\n")
