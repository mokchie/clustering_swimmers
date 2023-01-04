import os
Num_of_Nodes = 1
Ntasks_per_node = 48
wall_time = '4:00:00'
partition = 'regular'
nr = 5
N = range(nr)
infiles = ['in%s.run'%i for i in N]
for infile in infiles:
    with open('job.cluster.slm','w') as fp:
        jobname = infile
        fp.write(r'#!/bin/bash'+'\n')
        fp.write(r'#SBATCH --job-name %s'%jobname+'\n')
        fp.write(r'#SBATCH --nodes=%d'%Num_of_Nodes+'\n')
        fp.write(r'#SBATCH --ntasks-per-node=%d'%Ntasks_per_node+'\n')
        fp.write(r'#SBATCH --output=mpi-out.%j'+'\n')
        fp.write(r'#SBATCH --error=mpi-err.%j'+'\n')
        fp.write(r'#SBATCH --partition=%s'%partition+'\n')
        if partition == 'bcam-exclusive':
            fp.write(r'#SBATCH --account=bcam-exclusive'+'\n')
        fp.write(r'#SBATCH --mem=%s'%48000+'\n')
        fp.write(r'#SBATCH --time=%s'%wall_time+'\n')

        fp.write(r'module load intel/2020a'+'\n')
        fp.write(r'module load FFTW/3.3.8-intel-2020a'+'\n')

        fp.write(r'srun /scratch/cmo/lmp_dipc -in %s'%(infile,)+'\n')
    os.system(r'sbatch job.cluster.slm')
    print( 'job "%s" submitted!'%jobname )


