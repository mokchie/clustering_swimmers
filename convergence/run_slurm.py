import os
Num_of_Nodes = 1
Ntasks_per_node = 16
wall_time = '00:40:00'
#partition = 'bcam-exclusive'
partition = 'regular'
#dirs_list = ['.',]
dirs_list = ['rho'+str(i) for i in (6.0,7.0,8.0,9.0,10.0,11.0,12.0)]
cwd = os.getcwd()
for dir in dirs_list:
    infile = 'in.run'
    if dir == '.':
        jobname = cwd.split('/')[-1]
    else:
        jobname = dir
    with open(dir+'/job.cluster.slm','w') as fp:
        fp.write(r'#!/bin/bash'+'\n')
        fp.write(r'#SBATCH --job-name %s'%jobname+'\n')
        fp.write(r'#SBATCH --nodes=%d'%Num_of_Nodes+'\n')
        fp.write(r'#SBATCH --ntasks-per-node=%d'%Ntasks_per_node+'\n')
        fp.write(r'#SBATCH --output=mpi-out.%j'+'\n')
        fp.write(r'#SBATCH --error=mpi-err.%j'+'\n')
        fp.write(r'#SBATCH --partition=%s'%partition+'\n')
        if partition == 'bcam-exclusive':
            fp.write(r'#SBATCH --account=bcam-exclusive'+'\n')
        fp.write(r'#SBATCH --mem=%s'%64000+'\n')
        fp.write(r'#SBATCH --time=%s'%wall_time+'\n')

        fp.write(r'module load intel/2020a'+'\n')
        fp.write(r'module load FFTW/3.3.8-intel-2020a'+'\n')

        fp.write(r'srun /scratch/cmo/lmp_dipc -in %s'%(infile,)+'\n')
    os.chdir(dir)
    os.system(r'sbatch job.cluster.slm')
    print( 'job "%s" submitted!'%jobname )
    os.chdir(cwd)

