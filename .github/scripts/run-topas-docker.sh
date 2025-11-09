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

"$TOPAS_DOCKER_LAUNCHER" -g4data="$G4DATA_DIR" -extensions="$TOPAS_NBIO_ROOT" "${CONTAINER_ARGS[@]}"
