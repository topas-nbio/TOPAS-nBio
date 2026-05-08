#!/bin/bash
set -e

########################################################################################
# Postprocess submit: upload a local script to S3 and run it on AWS Batch.
# The container syncs OUTPUT_BUCKET/SIM_DIR (run_1/outputs, run_2/outputs, ...),
# accumulates files (excluding *.log) and gives them a run_N_ prefix, then runs your
# script from that directory. 

# Ensure that your script is runnable on raw data files that are located in the same 
# directory as the script. See the following for an example: 
# https://topas-nbio.readthedocs.io/en/latest/Cloud/CloudPart2.html#submitting-simulations-with-topas-submit-sh

# The new files created by the script are uploaded to SIM_DIR/postP_results/ on S3.
########################################################################################

########################################################################################
# Edit these variables for each run
########################################################################################

PROJECT_NAME="fricke_project"
RUN_DATE="2026-03-03"   # Date when the simulations were launched (must match topas_submit.sh run)
OUTPUT_BUCKET="topas-nbio-output"

LOCAL_SCRIPT="./postP.py"   # Local path to your postprocessing script
EXTRA_PIP_PACKAGES="numpy matplotlib"  # Optional: Python packages required by your script which will be pip installed in the postprocess container

JOB_QUEUE="topas-nbio-queue"
POSTPROCESS_JOB_DEFINITION="topas-nbio-postprocess-job"

########################################################################################
# AWS Commands
########################################################################################

SIM_DIR="projects/${PROJECT_NAME}/${RUN_DATE}"
SCRIPT_BASENAME="$(basename "$LOCAL_SCRIPT")"
SCRIPT_S3_URI="s3://${OUTPUT_BUCKET}/${SIM_DIR}/${SCRIPT_BASENAME}"

echo "Uploading script to ${SCRIPT_S3_URI}"
aws s3 cp "${LOCAL_SCRIPT}" "${SCRIPT_S3_URI}" --only-show-errors

echo "Submitting postprocess job (job-definition: ${POSTPROCESS_JOB_DEFINITION})"
aws batch submit-job \
  --job-name "postprocess-${PROJECT_NAME}-${RUN_DATE}" \
  --job-queue "${JOB_QUEUE}" \
  --job-definition "${POSTPROCESS_JOB_DEFINITION}" \
  --container-overrides "{
    \"environment\": [
      {\"name\": \"OUTPUT_BUCKET\",   \"value\": \"${OUTPUT_BUCKET}\"},
      {\"name\": \"SIM_DIR\",         \"value\": \"${SIM_DIR}\"},
      {\"name\": \"SCRIPT_S3_URI\",   \"value\": \"${SCRIPT_S3_URI}\"},
      {\"name\": \"SCRIPT_FILENAME\", \"value\": \"${SCRIPT_BASENAME}\"},
      {\"name\": \"EXTRA_PIP_PACKAGES\", \"value\": \"${EXTRA_PIP_PACKAGES}\"}
    ]
  }"
