#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 24:00:00
module add CMake
module add GSL/2.5-GCC-7.3.0-2.30
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/tawabyx/lhapdf/lib/ ./penmp 0.05 1.00 0.01 pp 0 ${1} ${2} mom ${3} ${4} > ./data_JT/KCBK_asmom_parent_x0bk1_mu${4}/mult_asmom_KCBK_LHCb_pp_${1}_${2}_p=${3}.dat
