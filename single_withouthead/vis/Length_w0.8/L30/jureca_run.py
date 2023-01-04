import os
infiles = ['in1.run','in2.run','in3.run','in4.run','in5.run']

Num_of_Nodes = 1
Ntasks_per_node = 24

for infile in infiles:
    with open('job.slm','w') as fp:
        jobname = infile
        fp.write(r'#!/bin/bash -x'+'\n')
        fp.write(r'#SBATCH --job-name %s'%jobname+'\n')
        fp.write(r'#SBATCH --nodes=%d'%Num_of_Nodes+'\n')
        fp.write(r'#SBATCH --ntasks-per-node=%d'%Ntasks_per_node+'\n')
        fp.write(r'#SBATCH --output=mpi-out.%j'+'\n')
        fp.write(r'#SBATCH --error=mpi-err.%j'+'\n')
        fp.write(r'#SBATCH --time=24:00:00'+'\n')
        fp.write(r'#SBATCH --partition=batch'+'\n')
        fp.write(r'#SBATCH --mail-user mcj_buaa@qq.com'+'\n')
        fp.write(r'#SBATCH --mail-type ALL'+'\n')
        #fp.write(r'module load Architecture/KNL'+'\n')
        fp.write(r'module load Intel ParaStationMPI'+'\n')
        fp.write(r'srun -N %d  /p/project/cjics21/MO/lammps/src/lmp_iff -in %s'%(Num_of_Nodes,infile)+'\n')

    os.system(r'sbatch -A jiff44 job.slm')
    print 'job "%s" submitted!'%os.getcwd().split('/')[-1]


