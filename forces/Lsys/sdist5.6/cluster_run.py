import os
infiles = ['in-0.run',]

Num_of_Nodes = 1
Ntasks_per_node = 12
partition = 'th2'
for infile in infiles:
    jobname = infile
    with open('job.cluster.slm','w') as fp:
        fp.write(r'#!/bin/bash -x'+'\n')
        fp.write(r'#SBATCH --job-name %s'%jobname+'\n')
        fp.write(r'#SBATCH --nodes=%d'%Num_of_Nodes+'\n')
        fp.write(r'#SBATCH --ntasks-per-node=%d'%Ntasks_per_node+'\n')
        fp.write(r'#SBATCH --output=mpi-out.%j'+'\n')
        fp.write(r'#SBATCH --error=mpi-err.%j'+'\n')
        fp.write(r'#SBATCH --partition=%s'%partition+'\n')
        fp.write(r'srun -N %d  ~/lammps/src/lmp_iff -in %s'%(Num_of_Nodes,infile)+'\n')

    os.system(r'sbatch job.cluster.slm')
    print 'job "%s" submitted!'%jobname


