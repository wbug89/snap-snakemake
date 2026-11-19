#!/usr/bin/env bash
set -uo pipefail

usage() {
  cat <<'EOF'
SNaP run wrapper (v1.0.3)

Usage:
  bash scripts/run_snap.sh --config CONFIG [--profile minimal|low-memory|workstation|server]
                           [--mode run|dry-run|preflight] [-- EXTRA_SNAKEMAKE_ARGS]

Examples:
  bash scripts/run_snap.sh --config config.full.yaml --profile workstation --mode preflight
  bash scripts/run_snap.sh --config config.full.yaml --profile workstation --mode dry-run
  bash scripts/run_snap.sh --config config.full.yaml --profile workstation --mode run
  bash scripts/run_snap.sh --config config.full.yaml --profile server --mode run -- --keep-going

Every invocation creates:
  <project_dir>/logs/runs/<timestamp>/
    run.log
    command.txt
    config.snapshot.yaml
    preflight_summary.txt
    environment.txt
    git_version.txt
    run_metadata.tsv
EOF
}

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

iso_now() {
  python - <<'PYTIME'
from datetime import datetime, timezone
print(datetime.now().astimezone().isoformat(timespec='seconds'))
PYTIME
}

CONFIG=""
PROFILE="low-memory"
MODE="run"
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG="$2"; shift 2 ;;
    --profile)
      PROFILE="$2"; shift 2 ;;
    --mode)
      MODE="$2"; shift 2 ;;
    --)
      shift; EXTRA_ARGS=("$@"); break ;;
    -h|--help)
      usage; exit 0 ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage >&2
      exit 2 ;;
  esac
done

if [[ -z "$CONFIG" ]]; then
  echo "ERROR: --config is required." >&2
  exit 2
fi
if [[ ! "$PROFILE" =~ ^(minimal|low-memory|workstation|server)$ ]]; then
  echo "ERROR: --profile must be minimal, low-memory, workstation, or server." >&2
  exit 2
fi
if [[ ! "$MODE" =~ ^(run|dry-run|preflight)$ ]]; then
  echo "ERROR: --mode must be run, dry-run, or preflight." >&2
  exit 2
fi

CONFIG_ABS="$(python - "$CONFIG" <<'PY'
import os, sys
print(os.path.abspath(os.path.expanduser(os.path.expandvars(sys.argv[1]))))
PY
)"
if [[ ! -f "$CONFIG_ABS" ]]; then
  echo "ERROR: config file not found: $CONFIG_ABS" >&2
  exit 2
fi

PROJECT_DIR="$(python - "$CONFIG_ABS" <<'PY'
import os, sys, yaml
with open(sys.argv[1]) as fh:
    cfg = yaml.safe_load(fh) or {}
value = cfg.get('project_dir')
if not value:
    raise SystemExit('config project_dir is missing')
print(os.path.abspath(os.path.expanduser(os.path.expandvars(str(value)))))
PY
)" || exit 2

RUN_ID="$(date '+%Y-%m-%d_%H%M%S')"
RUN_ROOT="${PROJECT_DIR}/logs/runs"
RUN_DIR="${RUN_ROOT}/${RUN_ID}"
mkdir -p "$RUN_DIR"
ln -sfn "$RUN_DIR" "${RUN_ROOT}/latest" 2>/dev/null || true

RUN_LOG="${RUN_DIR}/run.log"
METADATA="${RUN_DIR}/run_metadata.tsv"
PROFILE_DIR="${ROOT_DIR}/profiles/${PROFILE}"
PROFILE_FILE="${PROFILE_DIR}/profile.yaml"
SNAKEFILE="${ROOT_DIR}/Snakefile"

if [[ ! -f "$PROFILE_FILE" ]]; then
  echo "ERROR: resource profile not found: $PROFILE_FILE" >&2
  exit 2
fi

PROFILE_MEM_MB="$(python - "$PROFILE_FILE" <<'PYPROFILE'
import sys, yaml
with open(sys.argv[1]) as fh:
    cfg = yaml.safe_load(fh) or {}
resources = cfg.get('resources') or {}
value = resources.get('mem_mb')
if value is None:
    raise SystemExit('profile resources.mem_mb is missing')
print(int(value))
PYPROFILE
)" || exit 2

cp "$CONFIG_ABS" "${RUN_DIR}/config.snapshot.yaml"

{
  echo -e "key\tvalue"
  echo -e "run_id\t${RUN_ID}"
  echo -e "start_time\t$(iso_now)"
  echo -e "mode\t${MODE}"
  echo -e "profile\t${PROFILE}"
  echo -e "profile_mem_mb\t${PROFILE_MEM_MB}"
  echo -e "config\t${CONFIG_ABS}"
  echo -e "project_dir\t${PROJECT_DIR}"
  echo -e "status\tRUNNING"
} > "$METADATA"

{
  echo "SNaP VERSION: $(cat "${ROOT_DIR}/VERSION" 2>/dev/null || echo unknown)"
  echo "Git commit: $(git -C "$ROOT_DIR" rev-parse HEAD 2>/dev/null || echo not-a-git-checkout)"
  echo "Git describe: $(git -C "$ROOT_DIR" describe --tags --always --dirty 2>/dev/null || echo unavailable)"
  echo "Git branch: $(git -C "$ROOT_DIR" branch --show-current 2>/dev/null || echo unavailable)"
} > "${RUN_DIR}/git_version.txt"

{
  echo "date: $(iso_now)"
  echo "uname: $(uname -a 2>/dev/null || true)"
  echo "python: $(python --version 2>&1 || true)"
  echo "snakemake: $(snakemake --version 2>&1 || true)"
  echo "conda_env: ${CONDA_DEFAULT_ENV:-not-set}"
  echo
  echo "===== conda list ====="
  conda list 2>&1 || true
} > "${RUN_DIR}/environment.txt"

BASE_CMD=(snakemake -s "$SNAKEFILE" --configfile "$CONFIG_ABS" --workflow-profile "$PROFILE_DIR")
case "$MODE" in
  preflight)
    MAIN_TARGET="preflight"
    MAIN_CMD=("${BASE_CMD[@]}" \
      --config snap_resource_profile="$PROFILE" snap_profile_mem_mb="$PROFILE_MEM_MB" snap_preflight_only=true \
      preflight --forcerun preflight --quiet)
    ;;
  dry-run)
    MAIN_TARGET="all"
    MAIN_CMD=("${BASE_CMD[@]}" \
      --config snap_resource_profile="$PROFILE" snap_profile_mem_mb="$PROFILE_MEM_MB" \
      all --dry-run)
    ;;
  run)
    MAIN_TARGET="all"
    MAIN_CMD=("${BASE_CMD[@]}" \
      --config snap_resource_profile="$PROFILE" snap_profile_mem_mb="$PROFILE_MEM_MB" \
      all --rerun-incomplete)
    ;;
esac
MAIN_CMD+=("${EXTRA_ARGS[@]}")

{
  printf 'cd %q\n' "$ROOT_DIR"
  printf '%q ' "${MAIN_CMD[@]}"
  printf '\n'
} > "${RUN_DIR}/command.txt"

exec_logged() {
  (cd "$ROOT_DIR" && "$@") 2>&1 | tee -a "$RUN_LOG"
  local rc=${PIPESTATUS[0]}
  return "$rc"
}

# Dry-run output can be enormous for large cohorts. Keep the complete native
# Snakemake output in run.log, but show only the DAG summary / Job stats and
# important status or error lines on the console. This avoids relying on bare
# --quiet, which can suppress the Job stats summary in some Snakemake 9.x
# logging implementations.
exec_dryrun_logged() {
  (cd "$ROOT_DIR" && "$@") 2>&1 \
    | tee -a "$RUN_LOG" \
    | awk '
      /^Building DAG of jobs/ { print; next }
      /^Job stats:/ { in_stats=1; print; next }
      in_stats {
        print
        if ($0 ~ /^[[:space:]]*$/) in_stats=0
        next
      }
      /^Nothing to be done/ { print; next }
      /^This was a dry-run/ { print; next }
      /^(WorkflowError|RuleException|MissingInputException|MissingOutputException|SyntaxError|ValueError|ERROR:|Error:|Exception:)/ { print; next }
    '
  local rc=${PIPESTATUS[0]}
  if [[ $rc -ne 0 ]]; then
    echo "Dry-run failed. See the complete output in: ${RUN_LOG}"
    echo "Last 40 log lines:"
    tail -n 40 "$RUN_LOG"
  fi
  return "$rc"
}

{
  echo "SNaP v$(cat "${ROOT_DIR}/VERSION")"
  echo "Run ID: ${RUN_ID}"
  echo "Mode: ${MODE}"
  echo "Profile: ${PROFILE}"
  echo "Config: ${CONFIG_ABS}"
  echo "Run history: ${RUN_DIR}"
  echo
} | tee -a "$RUN_LOG"

# Record the explicit Snakemake target selected by the wrapper.
echo -e "target\t${MAIN_TARGET}" >> "$METADATA"

# Always refresh preflight before dry-run or execution so config/reference/sample summaries
# correspond exactly to this run.
if [[ "$MODE" != "preflight" ]]; then
  PREFLIGHT_CMD=("${BASE_CMD[@]}" \
    --config snap_resource_profile="$PROFILE" snap_profile_mem_mb="$PROFILE_MEM_MB" snap_preflight_only=true \
    preflight --forcerun preflight --quiet)
  echo "===== PREFLIGHT =====" | tee -a "$RUN_LOG"
  exec_logged "${PREFLIGHT_CMD[@]}"
  PREFLIGHT_RC=$?
  if [[ $PREFLIGHT_RC -ne 0 ]]; then
    echo -e "end_time\t$(iso_now)" >> "$METADATA"
    echo -e "status_final\tFAILED_PREFLIGHT" >> "$METADATA"
    echo -e "exit_code\t${PREFLIGHT_RC}" >> "$METADATA"
    exit "$PREFLIGHT_RC"
  fi
  if [[ -f "${PROJECT_DIR}/metadata/preflight_summary.txt" ]]; then
    cp "${PROJECT_DIR}/metadata/preflight_summary.txt" "${RUN_DIR}/preflight_summary.txt"
  fi
  echo >> "$RUN_LOG"
fi

MODE_LABEL="$(printf '%s' "$MODE" | tr '[:lower:]' '[:upper:]')"
echo "===== ${MODE_LABEL} =====" | tee -a "$RUN_LOG"
if [[ "$MODE" == "dry-run" ]]; then
  exec_dryrun_logged "${MAIN_CMD[@]}"
else
  exec_logged "${MAIN_CMD[@]}"
fi
RC=$?

if [[ -f "${PROJECT_DIR}/metadata/preflight_summary.txt" ]]; then
  cp "${PROJECT_DIR}/metadata/preflight_summary.txt" "${RUN_DIR}/preflight_summary.txt"
fi

if [[ $RC -eq 0 ]]; then
  case "$MODE" in
    run)
      if [[ "$MAIN_TARGET" == "all" ]]; then
        STATUS="SUCCESS"
        RESULT_MESSAGE="SNaP SUCCESS; target all satisfied; exit code 0"
      else
        STATUS="FAILED_TARGET_CHECK"
        RC=3
        RESULT_MESSAGE="SNaP FAILED; run mode did not target all; exit code ${RC}"
      fi
      ;;
    dry-run)
      STATUS="DRY_RUN_PASSED"
      RESULT_MESSAGE="SNaP DRY-RUN PASSED; target all; exit code 0"
      ;;
    preflight)
      STATUS="PREFLIGHT_PASSED"
      RESULT_MESSAGE="SNaP PREFLIGHT PASSED; exit code 0"
      ;;
  esac
else
  STATUS="FAILED"
  RESULT_MESSAGE="SNaP FAILED; exit code ${RC}"
fi
{
  echo -e "end_time\t$(iso_now)"
  echo -e "status_final\t${STATUS}"
  echo -e "exit_code\t${RC}"
} >> "$METADATA"

echo | tee -a "$RUN_LOG"
echo "$RESULT_MESSAGE" | tee -a "$RUN_LOG"
echo "Run log: ${RUN_LOG}" | tee -a "$RUN_LOG"
exit "$RC"
