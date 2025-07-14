# run PWSCF
#--------------------------------------------------

export LD_LIBRARY_PATH=/opt/conda/pkgs/liblapack-3.9.0-32_h7ac8fdf_openblas/lib:$LD_LIBRARY_PATH


RUNBIN="$HOME/install/q-e-qe-7.4/bin/pw.x"
RUNBIN="$HOME/install/q-e-qe-7.4/bin/pw.x"
PWSCF_INPUT="scf"
 ${RUNBIN}  < ${PWSCF_INPUT}.in > ${PWSCF_INPUT}.out
PWSCF_INPUT="nscf"
#   mpirun -np $NPROCS ${MPIARG} ${RUNBIN}  < ${PWSCF_INPUT}.in > ${PWSCF_INPUT}.out
PWSCF_INPUT="bands"
#   mpirun -np $NPROCS ${MPIARG} ${RUNBIN2}  < ${PWSCF_INPUT}.in > ${PWSCF_INPUT}.out

#--------------------------------------------------
