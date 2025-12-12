#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

BIN="${ROOT_DIR}/target/release/HORSCAN"

SOURCE_MON="${SOURCE_MON:-${ROOT_DIR}/test/sim_chrX_mon.bed}"
TARGET_MON="${TARGET_MON:-${ROOT_DIR}/test/chm13_chrX_mon.bed}"
SOURCE_HOR="${SOURCE_HOR:-${ROOT_DIR}/test/sim_chrX_hor.bed}"
TARGET_HOR="${TARGET_HOR:-${ROOT_DIR}/test/chm13_chrX_hor.bed}"
OUT_PREFIX="${OUT_PREFIX:-${ROOT_DIR}/test/test_chrX}"

# Scoring (positive magnitudes; penalties are negated internally)
MON_SCORES=( ${MON_SCORES:-10 60 150 5} )
HOR_SCORES=( ${HOR_SCORES:-60 100 100 1} )

# Heuristics
HOR_PAIR_BAND="${HOR_PAIR_BAND:-0}"   # 0 disables
RELAXED="${RELAXED:-0}"               # 1 enables --relaxed
DEBUG="${DEBUG:-0}"                   # 1 enables --debug
FORCE_BUILD="${FORCE_BUILD:-0}"       # 1 rebuilds even if binary exists

usage() {
  cat <<EOF
Usage (from repo root):
  bash test/run_HORSCAN.sh

Optional environment overrides:
  OUT_PREFIX=... SOURCE_MON=... TARGET_MON=... SOURCE_HOR=... TARGET_HOR=...
  MON_SCORES="10 60 150 5"  HOR_SCORES="60 100 100 1"
  HOR_PAIR_BAND=0  RELAXED=0  DEBUG=0  FORCE_BUILD=0

Notes:
  - chrX test may need ~6.5 GiB RAM.
  - Outputs: <OUT_PREFIX>.alignment.tsv and <OUT_PREFIX>.hor.tsv
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

cd -- "${ROOT_DIR}"

if [[ ! -x "${BIN}" || "${FORCE_BUILD}" == "1" ]]; then
  echo "[HORSCAN] Building (cargo build --release) ..."
  cargo build --release
fi

if [[ ! -x "${BIN}" ]]; then
  echo "[HORSCAN] ERROR: Binary not found at: ${BIN}" >&2
  echo "[HORSCAN] Hint: ensure Rust toolchain is installed and build succeeds." >&2
  exit 1
fi

CMD=(
  "${BIN}"
  --source "${SOURCE_MON}"
  --target "${TARGET_MON}"
  --source-hor "${SOURCE_HOR}"
  --target-hor "${TARGET_HOR}"
  -o "${OUT_PREFIX}"
  -m "${MON_SCORES[@]}"
  -h "${HOR_SCORES[@]}"
  --hor-pair-band "${HOR_PAIR_BAND}"
)

if [[ "${RELAXED}" == "1" ]]; then
  CMD+=( --relaxed )
fi

if [[ "${DEBUG}" == "1" ]]; then
  CMD+=( --debug )
fi

echo "[HORSCAN] Running chrX test ..."
"${CMD[@]}"

echo "[HORSCAN] Done. Open visualize/horscan_viewer.html to inspect results."
