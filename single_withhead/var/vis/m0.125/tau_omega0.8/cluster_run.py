import os
dirs_list = ['.',]
Num_of_Nodes = 1
Ntasks_per_node = 16
partition = 'jubio'
for dir in dirs_list:
    with open(dir+'/job.cluster.slm','w') as fp:
        if dir == '.':
            #jobname = '/'.join(os.getcwd().split('/')[-2:])
            jobname = os.getcwd().split('/')[-1]
        else:
            jobname = dir
        fp.write(r'#!/bin/bash -x'+'\n')
        fp.write(r'#SBATCH --job-name %s'%jobname+'\n')
        fp.write(r'#SBATCH --nodes=%d'%Num_of_Nodes+'\n')
        fp.write(r'#SBATCH --ntasks-per-node=%d'%Ntasks_per_node+'\n')
	#fp.write(r'#SBATCH --mincpus 8'+'\n')
	#fp.write(r'#SBATCH --ntasks 16'+'\n')
        fp.write(r'#SBATCH --output=mpi-out'+'\n')
        fp.write(r'#SBATCH --error=mpi-err'+'\n')
        fp.write(r'#SBATCH --time=24:00:00'+'\n')
        fp.write(r'#SBATCH --partition=%s'%partition+'\n')
#        fp.write(r'#SBATCH --mail-user mcj_buaa@qq.com'+'\n')
#        fp.write(r'#SBATCH --mail-type ALL'+'\n')
        fp.write(r'srun -N %d  ~/lammps/src/lmp_iff -in in.run'%Num_of_Nodes+'\n')

    if dir == '.':
        os.system(r'sbatch job.cluster.slm')
        print 'job "%s" submitted!'%os.getcwd().split('/')[-1]
    else:
        os.system(r'cd %s; sbatch job.cluster.slm; cd ..'%dir)
        print 'job "%s" submitted!'%dir


