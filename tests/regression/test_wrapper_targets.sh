#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/snap-wrapper.XXXXXX")"
trap 'rm -rf "$TMP_ROOT"' EXIT
mkdir -p "$TMP_ROOT/bin" "$TMP_ROOT/project"

cat > "$TMP_ROOT/config.yaml" <<EOF_CFG
workflow_stage: full
read_layout: PE
batch_id: test
ploidy: 2
project_dir: $TMP_ROOT/project
reference: $TMP_ROOT/ref.fa
fastq_dir: $TMP_ROOT/fastq
sample_group_regex: '^((?:SM)[0-9]{3})'
EOF_CFG

cat > "$TMP_ROOT/bin/snakemake" <<'EOF_FAKE'
#!/usr/bin/env bash
if [[ "${1:-}" == "--version" ]]; then
  echo 9.26.1
  exit 0
fi
if printf '%s\n' "$@" | grep -qx -- '--dry-run'; then
  cat <<'OUT'
Building DAG of jobs...
Job stats:
job   count
----  -----
all       1
total     1

This was a dry-run (flag -n).
OUT
fi
exit 0
EOF_FAKE
chmod +x "$TMP_ROOT/bin/snakemake"

run_and_check() {
  local mode="$1" expected="$2"
  rm -rf "$TMP_ROOT/project/logs"
  PATH="$TMP_ROOT/bin:$PATH" bash "$ROOT_DIR/scripts/run_snap.sh" \
    --config "$TMP_ROOT/config.yaml" \
    --profile minimal \
    --mode "$mode" >/dev/null 2>&1
  local cmd="$TMP_ROOT/project/logs/runs/latest/command.txt"
  grep -q -- "$expected" "$cmd"
  grep -q 'snap_resource_profile=minimal' "$cmd"
  grep -q 'snap_profile_mem_mb=6000' "$cmd"
  echo "PASS: wrapper $mode -> $expected"
}

run_and_check preflight '-- preflight'
preflight_cmd="$TMP_ROOT/project/logs/runs/latest/command.txt"
# The preflight invocation is intentionally quiet.
grep -q -- '--quiet' "$preflight_cmd"
echo 'PASS: wrapper preflight retains --quiet'

run_and_check dry-run '-- all'
dryrun_cmd="$TMP_ROOT/project/logs/runs/latest/command.txt"
# Snakemake 9.x may hide Job stats under --quiet; dry-run must not use it.
if grep -q -- '--quiet' "$dryrun_cmd"; then
  echo 'FAIL: wrapper dry-run unexpectedly contains --quiet' >&2
  cat "$dryrun_cmd" >&2
  exit 1
fi
echo 'PASS: wrapper dry-run omits --quiet so Job stats remain available'

run_and_check run '-- all'

# Verify wrapper-supplied Snakemake options remain before the target delimiter.
rm -rf "$TMP_ROOT/project/logs"
PATH="$TMP_ROOT/bin:$PATH" bash "$ROOT_DIR/scripts/run_snap.sh" \
  --config "$TMP_ROOT/config.yaml" \
  --profile minimal \
  --mode run -- --keep-going >/dev/null 2>&1
cmd="$TMP_ROOT/project/logs/runs/latest/command.txt"
grep -Eq -- '--keep-going[[:space:]]+--[[:space:]]+all' "$cmd"
echo 'PASS: extra Snakemake options remain before -- all'
