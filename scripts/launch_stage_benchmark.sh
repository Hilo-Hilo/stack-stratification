#!/usr/bin/env zsh
set -euo pipefail

ROOT_DIR=${0:A:h:h}
DEFAULT_RUN_DIR="${ROOT_DIR}/runs/luad_stage_benchmark"
PYTHON_BIN="${STACK_STRAT_PYTHON:-${ROOT_DIR}/.venv311/bin/python}"

ARGS=("$@")
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

LOG_PATH="${RUN_DIR}/benchmark.log"
PID_PATH="${RUN_DIR}/benchmark.pid"

cd "${ROOT_DIR}"

mkdir -p "${RUN_DIR}" "${RUN_DIR}/tmp" "${RUN_DIR}/matplotlib"

export TMPDIR="${RUN_DIR}/tmp"
export MPLCONFIGDIR="${RUN_DIR}/matplotlib"
export LOKY_MAX_CPU_COUNT="${LOKY_MAX_CPU_COUNT:-1}"

nohup "${PYTHON_BIN}" "${ROOT_DIR}/scripts/run_luad_stage_benchmark.py" "${ARGS[@]}" > "${LOG_PATH}" 2>&1 &
echo $! > "${PID_PATH}"
echo "pid=$(cat "${PID_PATH}")"
echo "log=${LOG_PATH}"
