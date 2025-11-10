#!/usr/bin/env bash
# Wrapper that forwards TOPAS input files to the OpenTOPAS Docker image while
# mounting the current test directory and the TOPAS-nBio sources.
set -euo pipefail

if [ "$#" -lt 1 ]; then
	echo "Usage: $0 <topas-parameter-file> [additional arguments...]" >&2
	exit 1
fi

CONTAINER_ARGS=()
for arg in "$@"; do
	if [[ "$arg" = /* ]]; then
		CONTAINER_ARGS+=("$arg")
	else
		CONTAINER_ARGS+=("/simulations/$arg")
	fi
done

if [ -z "${TOPAS_NBIO_ROOT:-}" ]; then
	TOPAS_NBIO_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
fi

: "${TOPAS_DOCKER_LAUNCHER:=$TOPAS_NBIO_ROOT/docker-run}"
: "${G4DATA_DIR:=$HOME/Applications/GEANT4/G4DATA}"

if [ ! -x "$TOPAS_DOCKER_LAUNCHER" ]; then
	echo "Unable to execute docker launcher at $TOPAS_DOCKER_LAUNCHER" >&2
	exit 1
fi

CMD=( "$TOPAS_DOCKER_LAUNCHER" "-g4data=$G4DATA_DIR" "-extensions=$TOPAS_NBIO_ROOT" "${CONTAINER_ARGS[@]}" )

if [ -t 1 ]; then
	exec "${CMD[@]}"
fi

if command -v python3 >/dev/null 2>&1; then
	python3 - "$TOPAS_DOCKER_LAUNCHER" "-g4data=$G4DATA_DIR" "-extensions=$TOPAS_NBIO_ROOT" "${CONTAINER_ARGS[@]}" <<'PY'
import os
import pty
import sys

cmd = sys.argv[1:]
if not cmd:
	sys.exit("no command provided to pseudo-tty launcher")

raise SystemExit(pty.spawn(cmd))
PY
	exit $?
fi

LOG_FILE=$(mktemp "${TMPDIR:-/tmp}/topas-docker.XXXXXX")
script -q "$LOG_FILE" "${CMD[@]}"
status=$?
cat "$LOG_FILE"
rm -f "$LOG_FILE"
exit $status
