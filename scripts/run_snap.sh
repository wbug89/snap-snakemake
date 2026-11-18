#!/usr/bin/env bash
set -uo pipefail

usage() {
  cat <<'EOF'
SNaP run wrapper (v1.0.2)

Usage:
  bash scripts/run_snap.sh --config CONFIG [--profile low-memory|workstation|server]
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
if [[ ! "$PROFILE" =~ ^(low-memory|workstation|server)$ ]]; then
  echo "ERROR: --profile must be low-memory, workstation, or server." >&2
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
SNAKEFILE="${ROOT_DIR}/Snakefile"

cp "$CONFIG_ABS" "${RUN_DIR}/config.snapshot.yaml"

{
  echo -e "key\tvalue"
  echo -e "run_id\t${RUN_ID}"
  echo -e "start_time\t$(iso_now)"
  echo -e "mode\t${MODE}"
  echo -e "profile\t${PROFILE}"
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
    MAIN_CMD=("${BASE_CMD[@]}" preflight --forcerun preflight --quiet --config snap_preflight_only=true)
    ;;
  dry-run)
    MAIN_CMD=("${BASE_CMD[@]}" --dry-run --quiet)
    ;;
  run)
    MAIN_CMD=("${BASE_CMD[@]}" --rerun-incomplete)
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

{
  echo "SNaP v$(cat "${ROOT_DIR}/VERSION")"
  echo "Run ID: ${RUN_ID}"
  echo "Mode: ${MODE}"
  echo "Profile: ${PROFILE}"
  echo "Config: ${CONFIG_ABS}"
  echo "Run history: ${RUN_DIR}"
  echo
} | tee -a "$RUN_LOG"

# Always refresh preflight before dry-run or execution so config/reference/sample summaries
# correspond exactly to this run.
if [[ "$MODE" != "preflight" ]]; then
  PREFLIGHT_CMD=("${BASE_CMD[@]}" preflight --forcerun preflight --quiet --config snap_preflight_only=true)
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
exec_logged "${MAIN_CMD[@]}"
RC=$?

if [[ -f "${PROJECT_DIR}/metadata/preflight_summary.txt" ]]; then
  cp "${PROJECT_DIR}/metadata/preflight_summary.txt" "${RUN_DIR}/preflight_summary.txt"
fi

if [[ $RC -eq 0 ]]; then
  STATUS="SUCCESS"
else
  STATUS="FAILED"
fi
{
  echo -e "end_time\t$(iso_now)"
  echo -e "status_final\t${STATUS}"
  echo -e "exit_code\t${RC}"
} >> "$METADATA"

echo | tee -a "$RUN_LOG"
echo "SNaP ${STATUS}; exit code ${RC}" | tee -a "$RUN_LOG"
echo "Run log: ${RUN_LOG}" | tee -a "$RUN_LOG"
exit "$RC"
