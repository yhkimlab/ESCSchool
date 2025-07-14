#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32      # Cores per node
#SBATCH --mem=191000
#SBATCH --partition=g1           # Partition Name
##SBATCH --nodelist=n011
##
#SBATCH --job-name=PWSCFExample
##SBATCH --time=96:00:00           # Runtime: HH:MM:SS
#SBATCH -o tjob.%N.%j.out         # STDOUT
#SBATCH -e tjob.%N.%j.err         # STDERR
##

hostname
date

module purge
module add compiler/2022.1.0
module add mkl/2022.1.0
module add mpi/2021.6.0
module add QE/7.1


export OMP_NUM_THREADS=1
export RSH_COMMAND="/usr/bin/ssh -x"
export MPIARG="-genv I_MPI_DEBUG=5 "


# slurm: PWSCF FOR MULTI-NODES
# LOCATION OF INPUT FILE
## To use current dir, set INP_DIR="local", OUT_DIR="local"


#Remove scratch (YES/NO)
RMSCRATCH="YES"
##
################################################################
#-------------------------------------------
# setup. check backup files.
#-------------------------------------------

NPROCS="${SLURM_NTASKS}"
basedir="${SLURM_SUBMIT_DIR}"

runlog="cal.log"

cd ${basedir}



# run PWSCF
#--------------------------------------------------

RUNBIN="$HOME/Install/q-e-qe-7.4/bin/pw.x"
RUNBIN="$HOME/Install/q-e-qe-7.4/bin/pw.x"
RUNBIN2="$HOME/Install/q-e-qe-7.4/PP/src/summer_tuto.x"
RUNBIN2="$HOME/Install/q-e-qe-7.4/PP/src/summer_tuto.x"
PWSCF_INPUT="scf"
   mpirun -np $NPROCS ${MPIARG} ${RUNBIN} -nk 2  < ${PWSCF_INPUT}.in > ${PWSCF_INPUT}.out
PWSCF_INPUT="summ"
   mpirun -np $NPROCS ${MPIARG} ${RUNBIN2}  < ${PWSCF_INPUT}.in > ${PWSCF_INPUT}.out

#--------------------------------------------------
