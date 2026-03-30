#!/bin/bash
set -e

########################################################################################
# Edit these variables for each run
########################################################################################

PROJECT_NAME="fricke_project"
RUN_DATE=$(date +%Y-%m-%d)
SIM_DIR_BASE="projects/${PROJECT_NAME}/${RUN_DATE}"

NUM_JOBS=2
JOB_NAME="topas-run-fricke"
INPUT_BUCKET="topas-nbio-input" # input and output bucket names should match batch-job-definition.json
OUTPUT_BUCKET="topas-nbio-output"

LOCAL_SIM_DIR="/Applications/TOPAS/public/TOPAS-nBio/examples/scorers/Fricke"
FILE_TO_RUN="FrickeIRT.txt"       # Main TOPAS parameter file

########################################################################################
#
########################################################################################

# For each job, create a job-specific copy of the input directory, append a random seed
# to the main parameter file, sync to a job-specific S3 prefix, and submit the job.
for i in $(seq 1 "$NUM_JOBS"); do
  # Job-specific temporary directory
  JOB_DIR="$(mktemp -d)"

  # Copy local simulation files into this job-specific directory
  cp -r "${LOCAL_SIM_DIR}/." "${JOB_DIR}/"

  # Generate a random seed for this job
  SEED=$RANDOM

  # Append the seed to the end of the main parameter file for this job
  {
    echo ""
    echo "i:Ts/Seed = ${SEED}"
  } >> "${JOB_DIR}/${FILE_TO_RUN}"

  # Job-specific SIM_DIR prefix in S3, e.g. projects/fricke_project/2025-02-18/run_1
  JOB_SIM_DIR="${SIM_DIR_BASE}/run_${i}"

  # Sync this job's inputs to its own S3 prefix
  aws s3 sync "${JOB_DIR}/" "s3://${INPUT_BUCKET}/${JOB_SIM_DIR}/" --only-show-errors

  # Clean up local temp directory
  rm -rf "${JOB_DIR}"

  # Submit the job. TOPAS_XVFB_DISPLAY prevents Xvfb conflicts when multiple jobs share an instance.
  aws batch submit-job \
    --job-name "${JOB_NAME}-${i}" \
    --job-queue topas-nbio-queue \
    --job-definition topas-nbio-job \
    --container-overrides "{
      \"environment\": [
        {\"name\": \"TOPAS_AWS_BATCH_MODE\", \"value\": \"1\"},
        {\"name\": \"TOPAS_XVFB_DISPLAY\",    \"value\": \":$((98 + i))\"},
        {\"name\": \"JOB_INDEX\",          \"value\": \"${i}\"},
        {\"name\": \"INPUT_BUCKET\",      \"value\": \"${INPUT_BUCKET}\"},
        {\"name\": \"OUTPUT_BUCKET\",     \"value\": \"${OUTPUT_BUCKET}\"},
        {\"name\": \"SIM_DIR\",           \"value\": \"${JOB_SIM_DIR}\"},
        {\"name\": \"FILE_TO_RUN\",       \"value\": \"${FILE_TO_RUN}\"}
      ]
    }"
done
