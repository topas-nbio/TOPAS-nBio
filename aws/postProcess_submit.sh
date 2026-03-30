#!/bin/bash
set -e

########################################################################################
# Postprocess submit: upload a local script to S3 and run it on AWS Batch.
# The container syncs OUTPUT_BUCKET/SIM_DIR (run_1/outputs, run_2/outputs, ...),
# accumulates files (excluding *.log) with run_N_ prefix, runs your script from
# that directory. Your script writes results to the current directory (same dir
# as the script); those new files are uploaded to SIM_DIR/postP_results/ on S3.
########################################################################################

########################################################################################
# Edit these variables for each run
########################################################################################

PROJECT_NAME="fricke_project"
RUN_DATE="2026-03-03"   # Date when the simulations were launched (must match topas_submit.sh run)
OUTPUT_BUCKET="topas-nbio-output"

LOCAL_SCRIPT="./postP.py"   # Local path to your postprocessing script
EXTRA_PIP_PACKAGES="numpy matplotlib"  # Optional: Python packages required by your script which will be pip installed in the postprocess container

########################################################################################
# AWS Commands
########################################################################################

SIM_DIR="projects/${PROJECT_NAME}/${RUN_DATE}"
SCRIPT_BASENAME="$(basename "$LOCAL_SCRIPT")"
SCRIPT_S3_URI="s3://${OUTPUT_BUCKET}/${SIM_DIR}/${SCRIPT_BASENAME}"

echo "Uploading script to ${SCRIPT_S3_URI}"
aws s3 cp "${LOCAL_SCRIPT}" "${SCRIPT_S3_URI}" --only-show-errors

echo "Submitting postprocess job (job-definition: topas-nbio-postprocess-job)"
aws batch submit-job \
  --job-name "postprocess-${PROJECT_NAME}-${RUN_DATE}" \
  --job-queue topas-nbio-queue \
  --job-definition topas-nbio-postprocess-job \
  --container-overrides "{
    \"environment\": [
      {\"name\": \"OUTPUT_BUCKET\",   \"value\": \"${OUTPUT_BUCKET}\"},
      {\"name\": \"SIM_DIR\",         \"value\": \"${SIM_DIR}\"},
      {\"name\": \"SCRIPT_S3_URI\",   \"value\": \"${SCRIPT_S3_URI}\"},
      {\"name\": \"SCRIPT_FILENAME\", \"value\": \"${SCRIPT_BASENAME}\"},
      {\"name\": \"EXTRA_PIP_PACKAGES\", \"value\": \"${EXTRA_PIP_PACKAGES}\"}
    ]
  }"
