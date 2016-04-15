#!/bin/bash
#SBATCH --job-name=serial-job
#SBATCH --time=48:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=shared
#SBATCH --mem=3000

. /work-zfs/lhc/cms/cmsset_default.sh
cd /scratch/groups/lhc/park/CMSSW_7_6_3/
eval $(scram ru -sh)
if [ ${SLURM_SUBMIT_DIR} ]; then cd ${SLURM_SUBMIT_DIR}; else cd - > /dev/null; fi || exit $?
echo "SLURM job running in: " `pwd`
./JHUGen Unweighted=0 VegasNc0=99999999 VegasNc1=99999999 VegasNc2=99999999 LHAPDF=NNPDF30_lo_as_0130/NNPDF30_lo_as_0130.info MReso=$1 $2

