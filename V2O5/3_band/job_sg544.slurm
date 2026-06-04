#!/bin/bash
#SBATCH -J vasp
#SBATCH -N 2
#SBATCH --ntasks-per-node=32
#SBATCH -p kshcnormal

module purge
module load vasp-5.4.4-intelmpi2017_optcell

export MKL_DEBUG_CPU_TYPE=5 #加速代码
export MKL_CBWR=AVX2 #使cpu默认支持avx2
export I_MPI_PIN_DOMAIN=numa #内存位置与cpu位置绑定，加速内存读取。对于内存带宽要求高的计算提速明显

srun --mpi=pmi2 vasp_std
