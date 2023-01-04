import os
infiles = ['in-1.run','in-2.run','in-3.run','in-4.run']+['in-5.run','in-6.run','in-7.run','in-8.run']

Num_of_Nodes = 1
Ntasks_per_node = 24
partition = 'jics21'
for i,infile in enumerate(infiles):
    with open('job.slm','w') as fp:
        jobname = os.getcwd().split('/')[-1]+'-'+str(i+1)
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
        fp.write(r'module load Intel ParaStationMPI'+'\n')
        fp.write(r'srun -N %d  /p/project/cjics21/MO/lammps/src/lmp_iff -in %s'%(Num_of_Nodes,infile)+'\n')

    os.system(r'sbatch -A %s job.slm'%partition)
    print 'job "%s" submitted!'%jobname
