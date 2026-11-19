#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PROFILE="${SNAP_TEST_PROFILE:-minimal}"

if ! command -v snakemake >/dev/null 2>&1; then
  echo "SKIP: snakemake is not available in PATH." >&2
  exit 77
fi

TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/snap-regression.XXXXXX")"
TEST_OK=0

cleanup() {
  local rc=$?
  trap - EXIT
  if [[ $rc -eq 0 && $TEST_OK -eq 1 ]]; then
    rm -rf "$TMP_ROOT"
  else
    echo >&2
    echo "Regression test failed." >&2
    echo "Temporary test files retained at:" >&2
    echo "  $TMP_ROOT" >&2
  fi
  exit "$rc"
}
trap cleanup EXIT

fail() {
  local message="$1"
  local output_file="${2:-}"
  echo "FAIL: $message" >&2
  if [[ -n "$output_file" && -f "$output_file" ]]; then
    echo >&2
    echo "--- captured output: $output_file ---" >&2
    cat "$output_file" >&2
    echo "--- end captured output ---" >&2
  fi
  return 1
}

assert_contains() {
  local pattern="$1" file="$2" label="$3"
  if ! grep -q -- "$pattern" "$file"; then
    fail "$label (missing: $pattern)" "$file"
  fi
}

assert_regex() {
  local pattern="$1" file="$2" label="$3"
  if ! grep -Eq -- "$pattern" "$file"; then
    fail "$label (pattern not found: $pattern)" "$file"
  fi
}

assert_not_contains() {
  local pattern="$1" file="$2" label="$3"
  if grep -Eq -- "$pattern" "$file"; then
    fail "$label (unexpected match: $pattern)" "$file"
  fi
}

make_reference() {
  local dir="$1"
  mkdir -p "$dir"
  cat > "$dir/ref.fa" <<'FA'
>chr1
AAAAAAAAAA
>chr2
CCCCCCCCCCCC
FA
  cat > "$dir/ref.fa.fai" <<'FAI'
chr1	10	6	10	11
chr2	12	23	12	13
FAI
  printf '@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10\n@SQ\tSN:chr2\tLN:12\n' > "$dir/ref.dict"
  for ext in amb ann bwt pac sa; do printf 'fixture\n' > "$dir/ref.fa.$ext"; done
}

write_config() {
  local cfg="$1" project="$2" fastq="$3" ref="$4" layout="$5"
  cat > "$cfg" <<EOF_CFG
workflow_stage: full
read_layout: ${layout}
batch_id: test
ploidy: 2
project_dir: ${project}
fastq_dir: ${fastq}
reference: ${ref}
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
EOF_CFG
}

run_case() {
  local name="$1" layout="$2" expected_trim="$3" forbidden_trim="$4" expected_count="$5"
  local case_dir="$TMP_ROOT/$name"
  local fq="$case_dir/fastq" project="$case_dir/project" refdir="$case_dir/ref"

  echo
  echo "========================================"
  echo "TEST: $name ($layout)"
  echo "========================================"

  mkdir -p "$fq" "$project"
  make_reference "$refdir"

  case "$name" in
    pe_one)
      printf '@r\nA\n+\nI\n' > "$fq/SM001_L001_R1.fastq.gz"
      printf '@r\nA\n+\nI\n' > "$fq/SM001_L001_R2.fastq.gz"
      ;;
    se_one)
      printf '@r\nA\n+\nI\n' > "$fq/SM001.fastq.gz"
      ;;
    pe_multi)
      for lane in L001 L002; do
        printf '@r\nA\n+\nI\n' > "$fq/SM001_${lane}_R1.fastq.gz"
        printf '@r\nA\n+\nI\n' > "$fq/SM001_${lane}_R2.fastq.gz"
      done
      ;;
    *) fail "Unknown case: $name" ;;
  esac

  local cfg="$case_dir/config.yaml"
  write_config "$cfg" "$project" "$fq" "$refdir/ref.fa" "$layout"

  if [[ "$name" == "pe_one" ]]; then
    local preflight_out="$case_dir/preflight.out"
    echo "[1/2] preflight"
    if ! bash "$ROOT_DIR/scripts/run_snap.sh" \
      --config "$cfg" \
      --profile "$PROFILE" \
      --mode preflight >"$preflight_out" 2>&1; then
      fail "$name preflight command failed" "$preflight_out"
    fi

    local preflight_latest="$project/logs/runs/latest"
    assert_contains '-- preflight' "$preflight_latest/command.txt" "$name wrapper target delimiter"
    assert_contains "Selected profile: $PROFILE" "$preflight_latest/preflight_summary.txt" "$name selected profile metadata"
    assert_contains 'Profile memory budget:' "$preflight_latest/preflight_summary.txt" "$name profile memory metadata"
    echo "PASS: preflight target/profile metadata"
  fi

  local out="$case_dir/dryrun.out"
  echo "[2/2] full dry-run"
  if ! bash "$ROOT_DIR/scripts/run_snap.sh" \
    --config "$cfg" \
    --profile "$PROFILE" \
    --mode dry-run >"$out" 2>&1; then
    fail "$name dry-run command failed" "$out"
  fi

  assert_contains 'Building DAG of jobs' "$out" "$name full DAG was not built"
  assert_contains 'Job stats:' "$out" "$name Job stats table"
  assert_regex "^[[:space:]]*${expected_trim}[[:space:]]+${expected_count}([[:space:]]|$)" "$out" "$name expected trimming job count"
  assert_regex '^[[:space:]]*all[[:space:]]+1([[:space:]]|$)' "$out" "$name all target count"

  if [[ "$name" == "pe_multi" ]]; then
    assert_regex '^[[:space:]]*bwa_mapping[[:space:]]+2([[:space:]]|$)' "$out" "$name BWA job count"
    assert_regex '^[[:space:]]*mark_duplicates[[:space:]]+2([[:space:]]|$)' "$out" "$name MarkDuplicates job count"
    assert_regex '^[[:space:]]*merge_marked_bam[[:space:]]+1([[:space:]]|$)' "$out" "$name merge job count"
    assert_regex '^[[:space:]]*gatk_haplotypecaller[[:space:]]+1([[:space:]]|$)' "$out" "$name HaplotypeCaller job count"
  fi

  assert_not_contains "^[[:space:]]*${forbidden_trim}[[:space:]]+[1-9]" "$out" "$name unexpectedly scheduled $forbidden_trim"
  assert_not_contains 'InputFunctionException in rule trim_fastq_se' "$out" "$name SE/PE wildcard collision"
  assert_not_contains '_R[12]_paired.*KeyError' "$out" "$name paired-output wildcard collision"

  local latest="$project/logs/runs/latest"
  assert_contains "Selected profile: $PROFILE" "$latest/preflight_summary.txt" "$name selected profile metadata after dry-run"
  assert_contains 'Profile memory budget:' "$latest/preflight_summary.txt" "$name profile memory metadata after dry-run"
  assert_contains '-- all' "$latest/command.txt" "$name dry-run all target delimiter"

  echo "PASS: $name"
}

echo "SNaP regression dry-run suite"
echo "Snakemake: $(snakemake --version 2>&1)"
echo "Profile: $PROFILE"
echo "Temporary root: $TMP_ROOT"

run_case pe_one PE trim_fastq_pe trim_fastq_se 1
run_case se_one SE trim_fastq_se trim_fastq_pe 1
run_case pe_multi PE trim_fastq_pe trim_fastq_se 2

TEST_OK=1
echo
echo "All SNaP regression dry-runs passed."
