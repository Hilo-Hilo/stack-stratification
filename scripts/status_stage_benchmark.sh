#!/usr/bin/env zsh
set -euo pipefail

ROOT_DIR=${0:A:h:h}
DEFAULT_RUN_DIR="${ROOT_DIR}/runs/luad_stage_benchmark"
RUN_DIR="${DEFAULT_RUN_DIR}"

for ((i = 1; i <= $#; i++)); do
  if [[ "${@[$i]}" == "--workdir" && $((i + 1)) -le $# ]]; then
    CANDIDATE="${@[$((i + 1))]}"
    if [[ "${CANDIDATE}" == /* ]]; then
      RUN_DIR="${CANDIDATE}"
    else
      RUN_DIR="${ROOT_DIR}/${CANDIDATE}"
    fi
    break
  fi
done

PID_PATH="${RUN_DIR}/benchmark.pid"
LOG_PATH="${RUN_DIR}/benchmark.log"
METRICS_PATH="${RUN_DIR}/metrics.json"

if [[ -f "${PID_PATH}" ]]; then
  PID=$(cat "${PID_PATH}")
  if kill -0 "${PID}" 2>/dev/null; then
    echo "status=running pid=${PID}"
  else
    echo "status=not-running pid=${PID}"
  fi
else
  echo "status=no-pid"
fi

if [[ -f "${METRICS_PATH}" ]]; then
  echo "metrics=${METRICS_PATH}"
fi

if [[ -f "${LOG_PATH}" ]]; then
  echo "log=${LOG_PATH}"
  tail -n 20 "${LOG_PATH}"
fi
