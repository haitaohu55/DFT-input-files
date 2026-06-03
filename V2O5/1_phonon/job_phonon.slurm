#!/bin/bash
#SBATCH -J hht-GdI2
#SBATCH -p kshcnormal          # 换成只要CPU、不强制DCU的分区
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -o %j.out
#SBATCH -e %j.err

module purge
module load vasp-5.4.4-intelmpi2017_optcell

export MKL_DEBUG_CPU_TYPE=5 #加速代码
export MKL_CBWR=AVX2 #使cpu默认支持avx2
export I_MPI_PIN_DOMAIN=numa #内存位置与cpu位置绑定，加速内存读取。对于内存带宽要求高的计算提速明显

for i in {001..014}
do
    mkdir disp-$i
    cd disp-$i
    cp ../files/*  .
    cp ../POSCAR-$i POSCAR
    srun --mpi=pmi2 vasp_std
    cp vasprun.xml ../vasprun.xml-$i
    cd ..
done

