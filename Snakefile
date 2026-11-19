# SNaP — Snakemake for NGS and Polyploids — v1.0.4-rc1
# Portability/resource-safety release. Final SNP IDs use CHROM_POS (e.g. LG1_123456).
import glob
import os
import platform
import re
import subprocess
import shutil
import shlex
import sys
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime
from pathlib import Path

from snakemake.shell import shell

# Workflow configuration must be supplied explicitly with --configfile.
# This prevents a user-selected config from being implicitly merged with the
# repository-root config.yaml. The wrapper always supplies --configfile.
if not config:
    raise ValueError(
        "No workflow config was loaded. Use --configfile <config.yaml> "
        "or run SNaP via scripts/run_snap.sh --config <config.yaml>."
    )
sys.path.insert(0, workflow.basedir)
from workflow_inputs import (reference_inputs, select_reference_sequences, build_joint_shards,
                             discover_marked, discover_gvcfs, selected_samples, safe_id, directories)

# Treat memory as a whole-machine scheduling constraint for local/workstation use.
# Profiles or CLI --resources values then prevent concurrent heavy jobs from overcommitting RAM.
resource_scopes:
    mem_mb="global",
    joint_jobs="global"


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
preflight_summary_path = os.path.join(metadata_dir, "preflight_summary.txt")
reference_metadata_path = os.path.join(metadata_dir, "reference_sequences.tsv")
selected_intervals_path = os.path.join(metadata_dir, "selected_intervals.list")
fastq_metadata_path = os.path.join(metadata_dir, "fastq_units.tsv")

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
SELECTED_REFERENCE_CONTIGS, REFERENCE_AUDIT, REFERENCE_INTERVAL_MODE = select_reference_sequences(
    REFERENCE_CONTIGS, config
)
JOINT_SHARDS, JOINT_SHARD_STRATEGY = build_joint_shards(SELECTED_REFERENCE_CONTIGS, config)
SELECTED_SEQUENCE_SET = {name for name, _ in SELECTED_REFERENCE_CONTIGS}
REFERENCE_TOTAL_LENGTH = sum(length for _, length in REFERENCE_CONTIGS)
SELECTED_REFERENCE_LENGTH = sum(length for _, length in SELECTED_REFERENCE_CONTIGS)

def detected_memory_mb():
    try:
        pages = os.sysconf("SC_PHYS_PAGES")
        page_size = os.sysconf("SC_PAGE_SIZE")
        return int(pages * page_size / 1024 / 1024)
    except (AttributeError, ValueError, OSError):
        if platform.system() == "Darwin":
            try:
                value = subprocess.check_output(["sysctl", "-n", "hw.memsize"], text=True).strip()
                return int(int(value) / 1024 / 1024)
            except Exception:
                return None
        return None


DETECTED_CPU_THREADS = os.cpu_count() or 1
DETECTED_MEMORY_MB = detected_memory_mb()
SNAP_RESOURCE_PROFILE = str(config.get("snap_resource_profile", "unspecified") or "unspecified")
try:
    SNAP_PROFILE_MEM_MB = int(config.get("snap_profile_mem_mb", 0) or 0)
except (TypeError, ValueError):
    raise ValueError("snap_profile_mem_mb must be an integer number of MiB")

PLATFORM_RELEASE = platform.release()
PLATFORM_DESCRIPTION = platform.platform()
IS_WSL2 = "microsoft" in PLATFORM_RELEASE.lower() or "wsl" in PLATFORM_RELEASE.lower()

RESOURCE_WARNINGS = []
if DETECTED_MEMORY_MB and SNAP_PROFILE_MEM_MB and SNAP_PROFILE_MEM_MB > DETECTED_MEMORY_MB:
    RESOURCE_WARNINGS.append(
        "The selected resource profile requests more memory than the detected physical memory. "
        f"Detected memory: {DETECTED_MEMORY_MB / 1024:.1f} GiB; "
        f"Profile memory budget: {SNAP_PROFILE_MEM_MB / 1024:.1f} GiB. "
        "Use a smaller resource profile or increase available memory."
    )

############################################
# Optional notifications
############################################

ENABLE_TELEGRAM = bool(cfg("telegram", "enabled", False))
PREFLIGHT_ONLY = bool(config.get("snap_preflight_only", False))

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
# Shared parameters and per-rule resources
############################################

TRIMMOMATIC_PHRED   = cfg("trimmomatic", "phred", "-phred33")
TRIMMOMATIC_THREADS = int(cfg("resources", "trimmomatic_threads", cfg("trimmomatic", "threads", 2)))
LEADING             = int(cfg("trimmomatic", "leading", 3))
TRAILING            = int(cfg("trimmomatic", "trailing", 3))
SLIDINGWINDOW       = str(cfg("trimmomatic", "slidingwindow", "4:30"))
MINLEN              = int(cfg("trimmomatic", "minlen", 36))
ILLUMINACLIP_SE     = str(cfg("trimmomatic", "illuminaclip_se", "2:30:10"))
ILLUMINACLIP_PE     = str(cfg("trimmomatic", "illuminaclip_pe", "2:30:10"))

BWA_THREADS           = int(cfg("resources", "bwa_threads", 4))
SAMTOOLS_VIEW_THREADS = int(cfg("resources", "samtools_view_threads", 1))
SAMTOOLS_SORT_THREADS = int(cfg("resources", "samtools_sort_threads", 2))
HC_THREADS            = int(cfg("resources", "haplotypecaller_threads", 2))
MIN_MAPPING_QUALITY   = int(cfg("calling", "minimum_mapping_quality", 30))

# Per-job memory declarations. Workflow profiles can override these with set-resources.
TRIMMOMATIC_MEM_MB = int(cfg("resources", "trimmomatic_mem_mb", 2048))
BWA_MEM_MB         = int(cfg("resources", "bwa_mem_mb", 4096))
SORT_MEM_MB        = int(cfg("resources", "samtools_sort_mem_mb", 4096))
MARKDUP_MEM_MB     = int(cfg("resources", "markduplicates_mem_mb", 6144))
MERGE_MEM_MB       = int(cfg("resources", "merge_bam_mem_mb", 4096))
HC_MEM_MB          = int(cfg("resources", "haplotypecaller_mem_mb", 8192))
BCFTOOLS_MEM_MB    = int(cfg("resources", "bcftools_mem_mb", 2048))

MARKDUP_JAVA_EXTRA = str(cfg("java_options", "markduplicates_extra", "") or "")
HC_JAVA_EXTRA = str(cfg("java_options", "haplotypecaller_extra", "") or "")
GENOMICSDB_JAVA_EXTRA = str(cfg("java_options", "genomicsdbimport_extra", "-XX:ParallelGCThreads=2") or "")
GENOTYPE_JAVA_EXTRA = str(cfg("java_options", "genotypegvcfs_extra", "-XX:ParallelGCThreads=2") or "")

GENOMICSDB_BATCH_SIZE    = int(cfg("genomicsdb", "batch_size", 30))
GENOMICSDB_READER_THREADS = int(cfg("genomicsdb", "reader_threads", 1))
GENOMICSDB_MEM_MB        = int(cfg("genomicsdb", "import_mem_mb", 8192))
GENOTYPE_MEM_MB          = int(cfg("genomicsdb", "genotype_mem_mb", 12288))
MAX_ALTERNATE_ALLELES    = int(cfg("calling", "max_alternate_alleles", 6))
MAX_GENOTYPE_COUNT       = int(cfg("calling", "max_genotype_count", 1024))


def java_options_from_mem(mem_mb, extra=""):
    """Use ~75% of the declared job memory for Java heap, leaving native/off-heap headroom."""
    mem_mb = int(mem_mb)
    heap_mb = max(512, int(mem_mb * 0.75))
    return f"-Xmx{heap_mb}m" + (f" {extra.strip()}" if str(extra).strip() else "")


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

SUPPORTED_FASTQ_EXTENSIONS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
FASTQ_RE = re.compile(r"\.(?:fastq|fq)(?:\.gz)?$", flags=re.IGNORECASE)
R1_RE = re.compile(r"^(?P<unit>.+)_R1(?P<lane>_001)?(?P<ext>\.(?:fastq|fq)(?:\.gz)?)$", flags=re.IGNORECASE)
R2_RE = re.compile(r"^(?P<unit>.+)_R2(?P<lane>_001)?(?P<ext>\.(?:fastq|fq)(?:\.gz)?)$", flags=re.IGNORECASE)
INTERMEDIATE_FASTQ_RE = re.compile(
    r"(?:^|[._-])(trim(?:med)?|paired|unpaired)(?:[._-]|$)", flags=re.IGNORECASE
)

SAMPLE_GROUP_REGEX = config.get("sample_group_regex", r"^((?:SM)[0-9]{3})")
BATCH_TOKEN_RE = re.compile(r"^(?=.*[A-Za-z])(?=.*\d)[A-Za-z0-9]{9,16}$")

try:
    SAMPLE_GROUP_PATTERN = re.compile(str(SAMPLE_GROUP_REGEX)) if SAMPLE_GROUP_REGEX else None
except re.error as exc:
    raise ValueError(f"Invalid sample_group_regex: {SAMPLE_GROUP_REGEX}: {exc}") from exc
if SAMPLE_GROUP_PATTERN and SAMPLE_GROUP_PATTERN.groups < 1:
    raise ValueError("sample_group_regex must contain capture group 1 for the biological sample ID.")


def strip_fastq_suffix(filename):
    return FASTQ_RE.sub("", filename)


def get_sm_name(raw_unit):
    """Convert one FASTQ unit name to the biological sample (GATK SM)."""
    if SAMPLE_GROUP_PATTERN:
        m = SAMPLE_GROUP_PATTERN.search(raw_unit)
        if not m:
            raise ValueError(
                f"sample_group_regex did not match raw unit '{raw_unit}': {SAMPLE_GROUP_REGEX}. "
                "Capture group 1 must contain only the biological sample ID."
            )
        value = m.group(1)
        if not value:
            raise ValueError(f"sample_group_regex capture group 1 is empty for raw unit '{raw_unit}'.")
        return safe_id(value)

    if "__" in raw_unit:
        return safe_id(raw_unit.split("__", 1)[0])
    sm = re.sub(r"_S\d+_L\d{3}$", "", raw_unit, flags=re.IGNORECASE)
    sm = re.sub(r"_L\d{3}$", "", sm, flags=re.IGNORECASE)
    if "_" in sm:
        prefix, suffix = sm.rsplit("_", 1)
        if BATCH_TOKEN_RE.fullmatch(suffix):
            sm = prefix
    return safe_id(sm)


def _is_supported_fastq(name):
    lower = name.lower()
    return any(lower.endswith(ext) for ext in SUPPORTED_FASTQ_EXTENSIONS)


def _candidate_pair_key(name):
    m1 = R1_RE.match(name)
    if m1:
        return (m1.group("unit"), m1.group("lane") or "", "R1")
    m2 = R2_RE.match(name)
    if m2:
        return (m2.group("unit"), m2.group("lane") or "", "R2")
    return None


def discover_fastqs(fq_dir, layout):
    """Discover raw FASTQs with strict SE/PE validation.

    Supported: .fastq.gz, .fq.gz, .fastq, .fq.
    Files containing trimmed/paired/unpaired output tokens are ignored so workflow
    intermediates cannot be rediscovered as new raw units.
    """
    if not os.path.isdir(fq_dir):
        raise FileNotFoundError(f"FASTQ directory does not exist: {fq_dir}")

    entries = []
    ignored_intermediate = []
    for f in sorted(os.listdir(fq_dir)):
        path = os.path.join(fq_dir, f)
        if not os.path.isfile(path) or not _is_supported_fastq(f):
            continue
        if INTERMEDIATE_FASTQ_RE.search(strip_fastq_suffix(f)):
            ignored_intermediate.append(f)
            continue
        entries.append(f)

    if not entries:
        raise ValueError(
            f"No raw FASTQ files were found in {fq_dir}. Supported extensions: "
            + ", ".join(SUPPORTED_FASTQ_EXTENSIONS)
        )

    if layout == "SE":
        found = {}
        for f in entries:
            if R2_RE.match(f):
                continue
            m = R1_RE.match(f)
            raw_unit = m.group("unit") if m else strip_fastq_suffix(f)
            path = os.path.join(fq_dir, f)
            if raw_unit in found:
                raise ValueError(
                    f"Multiple SE FASTQs resolve to unit '{raw_unit}': {found[raw_unit]} and {path}"
                )
            found[raw_unit] = path
        if not found:
            raise ValueError(f"No SE FASTQ units were retained from {fq_dir}.")
        return found, ignored_intermediate

    # PE: index R1/R2 independently so the mates may use any supported extension.
    pairs = {}
    unrecognized = []
    for f in entries:
        key = _candidate_pair_key(f)
        if key is None:
            unrecognized.append(f)
            continue
        unit, lane, mate = key
        pair_key = (unit, lane)
        pairs.setdefault(pair_key, {})
        if mate in pairs[pair_key]:
            raise ValueError(f"Duplicate {mate} FASTQ for PE unit '{unit}{lane}' in {fq_dir}")
        pairs[pair_key][mate] = os.path.join(fq_dir, f)

    if unrecognized:
        preview = ", ".join(unrecognized[:10])
        raise ValueError(
            "PE mode requires filenames ending in _R1/_R2 (optionally _001) before the FASTQ extension. "
            f"Unrecognized raw FASTQ file(s): {preview}"
        )

    found = {}
    missing = []
    for (unit, lane), mates in sorted(pairs.items()):
        absent = [x for x in ("R1", "R2") if x not in mates]
        if absent:
            missing.append(f"{unit}{lane}: missing {'/'.join(absent)}")
            continue
        raw_unit = unit
        if raw_unit in found:
            raise ValueError(
                f"Multiple PE lane groups resolve to the same raw unit '{raw_unit}'. "
                "Include lane/run information before _R1/_R2 so each raw unit remains unique."
            )
        found[raw_unit] = (mates["R1"], mates["R2"])
    if missing:
        raise FileNotFoundError("Incomplete PE FASTQ pair(s):\n  - " + "\n  - ".join(missing))
    if not found:
        raise ValueError(f"No complete PE FASTQ pairs were found in {fq_dir}.")
    return found, ignored_intermediate


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
IGNORED_INTERMEDIATE_FASTQS = []
UNIT_FASTQ_MAP = {}
UNIT_META = {}
local_units = []
skipped_units = []

if WORKFLOW_STAGE in {"preprocess", "full"}:
    RAW_FASTQ_MAP, IGNORED_INTERMEDIATE_FASTQS = discover_fastqs(fastq_dir, READ_LAYOUT)

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

# Restrict the {unit} wildcard to retained raw-unit keys. Without this constraint,
# the generic SE output pattern ``{unit}.trim.fastq.gz`` can also match PE outputs
# such as ``<unit>_R1_paired.trim.fastq.gz`` by treating the suffix as part of
# the wildcard. That can route a PE target through trim_fastq_se and trigger a
# KeyError in the SE input function. Escaping keeps literal sample/unit IDs safe
# inside the regular expression.
if local_units:
    UNIT_REGEX = "(?:" + "|".join(
        re.escape(unit) for unit in sorted(local_units, key=lambda x: (-len(x), x))
    ) + ")"
else:
    # No preprocessing rules are instantiated for stages without local FASTQ units.
    # A never-match regex is still defined so the global constraint is valid.
    UNIT_REGEX = r"(?!)"

wildcard_constraints:
    unit=UNIT_REGEX

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

# Joint shard IDs are deterministic for a given reference selection + sharding config.
SHARD_RECORDS = list(JOINT_SHARDS)
SHARD_BY_ID = {x["shard"]: x for x in SHARD_RECORDS}
SHARD_BY_SEQUENCE = {
    sequence: rec["shard"]
    for rec in SHARD_RECORDS
    for sequence in rec["sequences"]
}


def read_contig_assignment(path):
    """Read sequence->partition assignments for distributed joint_shard execution."""
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
            if parts[0].lower() in {"contig", "chromosome", "lg", "sequence"}:
                continue
            if len(parts) < 2:
                raise ValueError(f"Malformed sequence assignment line {lineno}: {line}")
            sequence, partition = parts[0], parts[1]
            if sequence not in SELECTED_SEQUENCE_SET:
                # Excluded reference sequences do not need assignment and are ignored.
                if sequence in {name for name, _ in REFERENCE_CONTIGS}:
                    continue
                raise ValueError(f"Sequence in assignment is absent from reference .fai: {sequence}")
            if sequence in assignment:
                raise ValueError(f"Duplicate sequence in joint_contig_assignment_file: {sequence}")
            assignment[sequence] = partition
    missing = SELECTED_SEQUENCE_SET - set(assignment)
    if missing:
        raise ValueError(
            "Selected reference sequences missing from joint_contig_assignment_file: "
            + ", ".join(sorted(missing))
        )
    return assignment


JOINT_PARTITION = str(config.get("joint_partition", "") or "").strip()
CONTIG_ASSIGNMENT = {}
SHARD_PARTITION = {}
if WORKFLOW_STAGE == "joint_shard":
    assignment_file = config.get("joint_contig_assignment_file")
    if not assignment_file or not JOINT_PARTITION:
        raise ValueError("joint_shard requires joint_contig_assignment_file and joint_partition")
    CONTIG_ASSIGNMENT = read_contig_assignment(assignment_file)
    for rec in SHARD_RECORDS:
        parts = {CONTIG_ASSIGNMENT[seq] for seq in rec["sequences"]}
        if len(parts) != 1:
            raise ValueError(
                f"Balanced shard {rec['shard']} contains sequences assigned to multiple partitions: "
                + ", ".join(f"{seq}={CONTIG_ASSIGNMENT[seq]}" for seq in rec["sequences"])
                + ". Reassign sequences so every sequence within one balanced shard uses the same partition."
            )
        SHARD_PARTITION[rec["shard"]] = next(iter(parts))
    ACTIVE_SHARD_RECORDS = [
        rec for rec in SHARD_RECORDS if SHARD_PARTITION[rec["shard"]] == JOINT_PARTITION
    ]
    if not ACTIVE_SHARD_RECORDS:
        raise ValueError(f"No joint shards are assigned to joint_partition={JOINT_PARTITION}")
elif WORKFLOW_STAGE in {"joint", "full"}:
    ACTIVE_SHARD_RECORDS = list(SHARD_RECORDS)
else:
    ACTIVE_SHARD_RECORDS = []
ACTIVE_SHARDS = [x["shard"] for x in ACTIVE_SHARD_RECORDS]

# joint_finalize receives completed per-contig VCFs from one or more joint_shard servers.
external_per_contig = {}
if WORKFLOW_STAGE == "joint_finalize":
    for d in directories(config, "per_contig_dirs"):
        for p in sorted(Path(d).glob("*.vcf.gz")):
            if p.name.endswith(".vcf.gz.tbi"):
                continue
            shard = p.name[:-len(".vcf.gz")]
            if shard not in SHARD_BY_ID:
                continue
            ap = str(p.resolve())
            if shard in external_per_contig:
                raise ValueError(f"Duplicate per-contig VCF for shard {shard}: {external_per_contig[shard]} and {ap}")
            if not os.path.isfile(ap + ".tbi"):
                raise ValueError(f"Missing index for per-contig VCF: {ap}.tbi")
            external_per_contig[shard] = ap
    missing = [x["shard"] for x in SHARD_RECORDS if x["shard"] not in external_per_contig]
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
# Concise console summary (full details are written by rule preflight)
############################################

print(f"SNaP v1.0.4-rc1 | stage={WORKFLOW_STAGE} | layout={READ_LAYOUT} | ploidy={PLOIDY}")
print(
    f"Reference: {len(SELECTED_REFERENCE_CONTIGS)}/{len(REFERENCE_CONTIGS)} sequences selected "
    f"({REFERENCE_INTERVAL_MODE}); joint strategy={JOINT_SHARD_STRATEGY}, shards={len(SHARD_RECORDS)}"
)
if WORKFLOW_STAGE in {"preprocess", "full"}:
    print(
        f"FASTQ: {len(local_units)} retained units / "
        f"{len(set(UNIT_META[u]['sample'] for u in local_units))} biological samples; "
        f"ignored intermediates={len(IGNORED_INTERMEDIATE_FASTQS)}"
    )
elif WORKFLOW_STAGE == "gvcf":
    print(f"gVCF calling: {len(sample_names)} biological samples")
elif WORKFLOW_STAGE in {"joint", "joint_shard"}:
    extra = f", partition={JOINT_PARTITION}" if WORKFLOW_STAGE == "joint_shard" else ""
    print(f"Joint genotyping: {len(sample_names)} samples, {len(ACTIVE_SHARD_RECORDS)} active shards{extra}")
elif WORKFLOW_STAGE == "joint_finalize":
    print(f"Joint finalize: {len(external_per_contig)} shard VCFs")
print("Run target 'preflight' for the complete validation summary.")

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
    if PREFLIGHT_ONLY or not TELEGRAM_READY:
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
    label = "SNaP preflight" if PREFLIGHT_ONLY else f"SNaP ({WORKFLOW_STAGE})"
    msg = (
        f"[Snakemake] {label} started at "
        f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
    )
    print(msg)
    send_telegram(msg)


onsuccess:
    label = "SNaP preflight" if PREFLIGHT_ONLY else f"SNaP ({WORKFLOW_STAGE})"
    if PREFLIGHT_ONLY:
        result_note = preflight_summary_path
    elif WORKFLOW_STAGE == "preprocess":
        result_note = os.path.join(metadata_dir, f"{BATCH_ID}.units.tsv")
    elif WORKFLOW_STAGE == "gvcf":
        result_note = gvcf_dir
    elif WORKFLOW_STAGE == "joint_shard":
        result_note = PER_CONTIG_DIR
    else:
        result_note = os.path.join(results_dir, "final_combined_filtered_snps.vcf.gz")
    msg = (
        f"[Snakemake] {label} completed successfully at "
        f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        f"Output: {result_note}"
    )
    print(msg)
    send_telegram(msg)


onerror:
    label = "SNaP preflight" if PREFLIGHT_ONLY else f"SNaP ({WORKFLOW_STAGE})"
    msg = (
        f"[Snakemake] {label} failed at "
        f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        f"Please check logs under: {log_dir}"
    )
    print(msg)
    send_telegram(msg)

############################################
# Preflight metadata / validation report
############################################

rule preflight:
    output:
        summary = preflight_summary_path,
        reference = reference_metadata_path,
        intervals = selected_intervals_path,
        fastqs = fastq_metadata_path
    run:
        os.makedirs(metadata_dir, exist_ok=True)
        with open(output.intervals, "w") as out:
            for name, _ in SELECTED_REFERENCE_CONTIGS:
                out.write(f"{name}\n")

        with open(output.reference, "w") as out:
            out.write("sequence\tlength\tselected\texclusion_reason\tshard\n")
            for rec in REFERENCE_AUDIT:
                shard = SHARD_BY_SEQUENCE.get(rec["sequence"], ".")
                out.write(
                    f"{rec['sequence']}\t{rec['length']}\t"
                    f"{'TRUE' if rec['selected'] else 'FALSE'}\t"
                    f"{rec['exclusion_reason']}\t{shard}\n"
                )

        with open(output.fastqs, "w") as out:
            out.write("unit_key\tbiological_sample\traw_unit\tread_layout\tR1_or_SE\tR2\n")
            for unit in local_units:
                meta = UNIT_META[unit]
                fq = UNIT_FASTQ_MAP[unit]
                if READ_LAYOUT == "PE":
                    r1, r2 = fq
                else:
                    r1, r2 = fq, "."
                out.write(f"{unit}\t{meta['sample']}\t{meta['raw_unit']}\t{READ_LAYOUT}\t{r1}\t{r2}\n")

        samples_local = sorted({UNIT_META[u]["sample"] for u in local_units})
        sample_preview = []
        for raw_unit in sorted(RAW_FASTQ_MAP)[:10]:
            sample_preview.append(f"  {raw_unit} -> {get_sm_name(raw_unit)}")

        with open(output.summary, "w") as out:
            out.write("SNaP v1.0.4-rc1 preflight\n")
            out.write("======================\n\n")
            out.write(f"Stage: {WORKFLOW_STAGE}\n")
            out.write(f"Read layout: {READ_LAYOUT}\n")
            out.write(f"Ploidy: {PLOIDY}\n")
            out.write(f"Platform: {PLATFORM_DESCRIPTION}\n")
            out.write("\nResources\n---------\n")
            out.write(f"Detected CPU threads: {DETECTED_CPU_THREADS}\n")
            out.write(
                "Detected physical memory: "
                + (f"{DETECTED_MEMORY_MB / 1024:.1f} GiB\n" if DETECTED_MEMORY_MB else "unknown\n")
            )
            out.write(f"Selected profile: {SNAP_RESOURCE_PROFILE}\n")
            out.write(
                "Profile memory budget: "
                + (f"{SNAP_PROFILE_MEM_MB / 1024:.1f} GiB ({SNAP_PROFILE_MEM_MB} MiB)\n" if SNAP_PROFILE_MEM_MB else "not provided\n")
            )
            if IS_WSL2:
                out.write(
                    "WSL2 detected: memory visible to SNaP is the memory assigned to WSL2, "
                    "not necessarily the total physical RAM installed in Windows.\n"
                )
            if RESOURCE_WARNINGS:
                out.write("\nWARNING:\n")
                for warning in RESOURCE_WARNINGS:
                    out.write(f"- {warning}\n")
                if IS_WSL2:
                    out.write(
                        '- Consider using the "minimal" profile or increasing the WSL2 memory limit.\n'
                    )
                else:
                    out.write(
                        '- Consider using a smaller resource profile or increasing available RAM.\n'
                    )
            out.write("\nReference\n---------\n")
            out.write(f"FASTA: {reference_dir}\n")
            out.write(f"Sequences in FASTA: {len(REFERENCE_CONTIGS)}\n")
            out.write(f"Total reference length: {REFERENCE_TOTAL_LENGTH:,} bp\n")
            out.write(f"Interval mode: {REFERENCE_INTERVAL_MODE}\n")
            out.write(f"Selected sequences: {len(SELECTED_REFERENCE_CONTIGS)}\n")
            out.write(f"Excluded sequences: {len(REFERENCE_CONTIGS) - len(SELECTED_REFERENCE_CONTIGS)}\n")
            out.write(f"Selected length: {SELECTED_REFERENCE_LENGTH:,} bp\n")
            out.write(f"Joint sharding strategy: {JOINT_SHARD_STRATEGY}\n")
            out.write(f"Global joint shards: {len(SHARD_RECORDS)}\n")
            if WORKFLOW_STAGE in {"preprocess", "full"}:
                out.write("\nFASTQ\n-----\n")
                out.write("Supported extensions: .fastq.gz, .fq.gz, .fastq, .fq\n")
                out.write(f"Detected raw units: {len(RAW_FASTQ_MAP)}\n")
                out.write(f"Retained raw units: {len(local_units)}\n")
                out.write(f"Biological samples: {len(samples_local)}\n")
                out.write(f"Ignored intermediate FASTQs: {len(IGNORED_INTERMEDIATE_FASTQS)}\n")
                raw_paths = []
                for value in RAW_FASTQ_MAP.values():
                    if isinstance(value, (tuple, list)):
                        raw_paths.extend(value)
                    else:
                        raw_paths.append(value)
                extension_counts = {}
                for path in raw_paths:
                    lower = str(path).lower()
                    ext = next((x for x in SUPPORTED_FASTQ_EXTENSIONS if lower.endswith(x)), "unknown")
                    extension_counts[ext] = extension_counts.get(ext, 0) + 1
                if extension_counts:
                    out.write("FASTQ extension counts: " + ", ".join(
                        f"{ext}={extension_counts[ext]}" for ext in SUPPORTED_FASTQ_EXTENSIONS if ext in extension_counts
                    ) + "\n")
                out.write(f"sample_group_regex: {SAMPLE_GROUP_REGEX}\n")
                if sample_preview:
                    out.write("Sample parsing preview:\n" + "\n".join(sample_preview) + "\n")
            if WORKFLOW_STAGE == "gvcf":
                out.write(f"\nBiological samples entering HaplotypeCaller: {len(sample_names)}\n")
            if WORKFLOW_STAGE in {"joint", "joint_shard"}:
                out.write(f"\nSingle-sample gVCFs: {len(sample_names)}\n")
                out.write(f"Active joint shards: {len(ACTIVE_SHARD_RECORDS)}\n")
            if WORKFLOW_STAGE == "joint_finalize":
                out.write(f"\nPer-shard VCFs available: {len(external_per_contig)}\n")
            if RESOURCE_WARNINGS:
                out.write("\nValidation: PASS WITH WARNINGS\n")
            else:
                out.write("\nValidation: PASS\n")

        print(Path(output.summary).read_text())
        print(f"Preflight summary written to: {output.summary}")


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
    records = ACTIVE_SHARD_RECORDS if records is None else records
    shards = [x["shard"] for x in records]
    return (
        expand(os.path.join(PER_CONTIG_DIR, "{shard}.vcf.gz"), shard=shards)
        + expand(os.path.join(PER_CONTIG_DIR, "{shard}.vcf.gz.tbi"), shard=shards)
    )


def all_targets(wc):
    preflight_files = [preflight_summary_path, reference_metadata_path, selected_intervals_path, fastq_metadata_path]
    if WORKFLOW_STAGE == "full":
        return preflight_files + preprocess_targets(wc) + final_variant_targets(wc) + [os.path.join(metadata_dir, "full_inputs.tsv")]
    if WORKFLOW_STAGE == "preprocess":
        return preflight_files + preprocess_targets(wc)
    if WORKFLOW_STAGE == "gvcf":
        return (preflight_files + expand(os.path.join(gvcf_dir, "{sm}.g.vcf.gz"), sm=sample_names)
                + expand(os.path.join(gvcf_dir, "{sm}.g.vcf.gz.tbi"), sm=sample_names)
                + [os.path.join(metadata_dir, "gvcf_inputs.tsv")])
    if WORKFLOW_STAGE == "joint_shard":
        return preflight_files + per_contig_targets() + [
            os.path.join(metadata_dir, "joint_shard_inputs.tsv"),
            os.path.join(metadata_dir, "joint_contigs.tsv"),
            os.path.join(metadata_dir, "genomicsdb.sample_map.tsv"),
        ]
    if WORKFLOW_STAGE == "joint_finalize":
        return preflight_files + final_variant_targets(wc) + [os.path.join(metadata_dir, "joint_finalize_inputs.tsv")]
    return preflight_files + final_variant_targets(wc) + [os.path.join(metadata_dir, f"{WORKFLOW_STAGE}_inputs.tsv")]


rule all:
    input:
        all_targets
    default_target: True


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
            fq = se_fastq,
            preflight = preflight_summary_path
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
        resources:
            mem_mb = TRIMMOMATIC_MEM_MB
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
            R2 = pe_r2,
            preflight = preflight_summary_path
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
        resources:
            mem_mb = TRIMMOMATIC_MEM_MB
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
        resources:
            mem_mb = BWA_MEM_MB
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
        resources:
            mem_mb = SORT_MEM_MB
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
            java_option = lambda wc, resources: java_options_from_mem(resources.mem_mb, MARKDUP_JAVA_EXTRA)
        output:
            bam = os.path.join(marked_dir, "{unit}.marked.bam"),
            bai = os.path.join(marked_dir, "{unit}.marked.bam.bai"),
            metrics = os.path.join(marked_dir, "{unit}.marked_dup_metrics.txt")
        log:
            os.path.join(log_dir, "{unit}.markdup.log")
        benchmark:
            os.path.join(benchmark_dir, "{unit}.markdup.txt")
        resources:
            mem_mb = MARKDUP_MEM_MB
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
        resources:
            mem_mb = MERGE_MEM_MB
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
            reference_files = REFERENCE_INPUTS,
            preflight = preflight_summary_path,
            intervals = selected_intervals_path
        params:
            java_option = lambda wc, resources: java_options_from_mem(resources.mem_mb, HC_JAVA_EXTRA),
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
        resources:
            mem_mb = HC_MEM_MB
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
                -L {input.intervals:q} \
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
# Joint genotyping: GenomicsDB per deterministic shard
############################################


def sequences_for_shard(wc):
    return SHARD_BY_ID[wc.shard]["sequences"]


def interval_args_for_shard(wc):
    return " ".join(f"-L {shlex.quote(seq)}" for seq in sequences_for_shard(wc))


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
            reference_files = REFERENCE_INPUTS,
            preflight = preflight_summary_path
        output:
            sample_map = os.path.join(metadata_dir, "genomicsdb.sample_map.tsv"),
            contigs = os.path.join(metadata_dir, "joint_contigs.tsv")
        run:
            with open(output.sample_map, "w") as out:
                for sm in sorted(external_gvcfs):
                    out.write(f"{sm}\t{external_gvcfs[sm]}\n")
            with open(output.contigs, "w") as out:
                out.write("shard\tsequences\ttotal_length\tpartition\n")
                for rec in ACTIVE_SHARD_RECORDS:
                    partition = SHARD_PARTITION.get(rec["shard"], "ALL")
                    out.write(
                        f"{rec['shard']}\t{','.join(rec['sequences'])}\t{rec['length']}\t{partition}\n"
                    )


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
            java_option = lambda wc, resources: java_options_from_mem(resources.mem_mb, GENOMICSDB_JAVA_EXTRA),
            interval_args = interval_args_for_shard,
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
                {params.interval_args} \
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
            java_option = lambda wc, resources: java_options_from_mem(resources.mem_mb, GENOTYPE_JAVA_EXTRA),
            interval_args = interval_args_for_shard,
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
                {params.interval_args} \
                --max-alternate-alleles {MAX_ALTERNATE_ALLELES} \
                --max-genotype-count {MAX_GENOTYPE_COUNT} \
                --create-output-variant-index true \
                -O {output.vcf:q} \
                > {log:q} 2>&1
            """


if WORKFLOW_STAGE == "joint_finalize":
    rule write_finalize_manifest:
        input:
            preflight = preflight_summary_path,
            vcfs = lambda wc: [external_per_contig[x["shard"]] for x in SHARD_RECORDS],
            tbis = lambda wc: [external_per_contig[x["shard"]] + ".tbi" for x in SHARD_RECORDS]
        output:
            os.path.join(metadata_dir, "joint_finalize_inputs.tsv")
        run:
            with open(output[0], "w") as out:
                out.write("shard\tsequences\ttotal_length\tvcf\n")
                for rec in SHARD_RECORDS:
                    out.write(
                        f"{rec['shard']}\t{','.join(rec['sequences'])}\t{rec['length']}\t"
                        f"{external_per_contig[rec['shard']]}\n"
                    )


if WORKFLOW_STAGE in {"joint", "full", "joint_finalize"}:
    rule concat_joint_vcfs:
        input:
            preflight = preflight_summary_path,
            vcfs = lambda wc: (
                [external_per_contig[x["shard"]] for x in SHARD_RECORDS]
                if WORKFLOW_STAGE == "joint_finalize"
                else [os.path.join(PER_CONTIG_DIR, x["shard"] + ".vcf.gz") for x in SHARD_RECORDS]
            ),
            tbis = lambda wc: (
                [external_per_contig[x["shard"]] + ".tbi" for x in SHARD_RECORDS]
                if WORKFLOW_STAGE == "joint_finalize"
                else [os.path.join(PER_CONTIG_DIR, x["shard"] + ".vcf.gz.tbi") for x in SHARD_RECORDS]
            )
        output:
            vcf = os.path.join(results_dir, "allsample.vcf.gz"),
            tbi = os.path.join(results_dir, "allsample.vcf.gz.tbi")
        log:
            os.path.join(log_dir, "concat_joint_vcfs.log")
        resources:
            mem_mb = BCFTOOLS_MEM_MB
        shell:
            r"""
            # Balanced shards can contain non-adjacent reference sequences. Concatenate
            # all disjoint shard VCFs, then sort using the VCF header/reference order.
            {BCFTOOLS:q} concat -a -Ou {input.vcfs:q} 2> {log:q} \
            | {BCFTOOLS:q} sort -Oz -o {output.vcf:q} - >> {log:q} 2>&1
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
            {BCFTOOLS:q} annotate --set-id '%CHROM_%POS' {input:q} -Oz -o {output.vcf:q} > {log:q} 2>&1
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
            {BCFTOOLS:q} annotate --set-id '%CHROM_%POS' {input:q} -Oz -o {output.vcf:q} > {log:q} 2>&1
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
