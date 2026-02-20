#!/bin/bash
set -e
# Postprocess runner: sync S3 prefix to local, accumulate run_*/outputs (exclude *.log),
# download user script, run it from the accumulated dir. User script writes results to
# the current directory (same dir as the script and accumulated data). New files created
# by the script are uploaded to S3 under SIM_DIR/postP_results/.
# Env: OUTPUT_BUCKET, SIM_DIR, SCRIPT_S3_URI, SCRIPT_FILENAME (optional, default: infer from S3 URI)

WORK_DIR="/work"
ACCUM_DIR="${WORK_DIR}/accumulated"
UPLOAD_DIR="${WORK_DIR}/upload_results"
mkdir -p "${ACCUM_DIR}"
cd "${WORK_DIR}"

echo "Postprocess runner: syncing s3://${OUTPUT_BUCKET}/${SIM_DIR}/ to ${WORK_DIR}/data/"
aws s3 sync "s3://${OUTPUT_BUCKET}/${SIM_DIR}/" "${WORK_DIR}/data/" --only-show-errors

echo "Accumulating run_*/outputs (excluding *.log) into ${ACCUM_DIR} with run_N_ prefix"
for run_dir in "${WORK_DIR}"/data/run_*/; do
  [ -d "$run_dir" ] || continue
  run_name="$(basename "$run_dir")"
  out_dir="${run_dir}outputs"
  [ -d "$out_dir" ] || continue
  for f in "$out_dir"/*; do
    [ -e "$f" ] || continue
    case "$f" in
      *.log) continue ;;
    esac
    b="$(basename "$f")"
    cp "$f" "${ACCUM_DIR}/${run_name}_${b}"
  done
done

echo "Downloading user script from ${SCRIPT_S3_URI}"
SCRIPT_FILENAME="${SCRIPT_FILENAME:-$(basename "${SCRIPT_S3_URI}")}"
aws s3 cp "${SCRIPT_S3_URI}" "${ACCUM_DIR}/${SCRIPT_FILENAME}" --only-show-errors

cd "${ACCUM_DIR}"
find . -type f -print | sort > "${WORK_DIR}/before.txt"

echo "Running user script: ${SCRIPT_FILENAME} (write results to current directory)"
if [[ "$SCRIPT_FILENAME" == *.py ]]; then
  python3 "${SCRIPT_FILENAME}" || { echo "User script failed"; exit 1; }
else
  chmod +x "${SCRIPT_FILENAME}" && ./"${SCRIPT_FILENAME}" || { echo "User script failed"; exit 1; }
fi

find . -type f -print | sort > "${WORK_DIR}/after.txt"
mkdir -p "${UPLOAD_DIR}"
comm -13 "${WORK_DIR}/before.txt" "${WORK_DIR}/after.txt" | while read -r path; do
  mkdir -p "${UPLOAD_DIR}/$(dirname "$path")"
  cp "./$path" "${UPLOAD_DIR}/$path"
done

echo "Syncing results to s3://${OUTPUT_BUCKET}/${SIM_DIR}/postP_results/"
aws s3 sync "${UPLOAD_DIR}/" "s3://${OUTPUT_BUCKET}/${SIM_DIR}/postP_results/" --only-show-errors
echo "Done."
