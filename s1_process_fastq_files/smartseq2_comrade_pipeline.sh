#!/bin/bash
#SBATCH -A ACCOUNT-SL2-CPU
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p partition
# set max wallclock time
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=user.name@cam.ac.uk
#SBATCH --job-name logfile
#SBATCH -e logfile.%N.%j.log


# ---------------------------------------------------------
#title           :smartseq2_pipeliner.sh
#author          :hpb29
#date            :20180424
#version         :1.1
#description     :Runs the full SmartSeq2 sequence data upstream
#                 analysis pipeline (FastQC > GSNAP > HTSeq) via
#                 SLURM.
#usage           :COPY this file to your RDS working area, edit
#                 variables within script & sbatch submit it.
# ----------------------------------------------------------


# ADJUST VARIABLES BELOW AS NEEDED
# =================================================================
CONFIG=ss2_config.json
SPECIES=human
WORKDIR=working_dir
OUTDIR=output_dir
# =================================================================


echo "Starting at `date`"


echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"

CONTAINER="singularity run containers/smartseq2_pipeline_comrades_versions.simg"


${CONTAINER} -c "${CONFIG}" -s "${SPECIES}" -w "${WORKDIR}" QC -o "${OUTDIR}"
${CONTAINER} -c "${CONFIG}" -s "${SPECIES}" -w "${WORKDIR}" ALIGN -o "${OUTDIR}"

${CONTAINER} -c "${CONFIG}" -s "${SPECIES}" -w "${OUTDIR}/SAM" QUANT -o "${OUTDIR}"
${CONTAINER} -c "${CONFIG}" -s "${SPECIES}" -w "${OUTDIR}/HTSeq" CONCAT


echo "Finished at `date`"
