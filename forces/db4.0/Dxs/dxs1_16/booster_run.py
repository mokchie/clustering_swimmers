import os
#dirs_list = ['b'+str(i) for i in (0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6)]
dirs_list = ['.']
Num_of_Nodes = 1
Ntasks_per_node = 64

for dir in dirs_list:
    with open(dir+'/job.slm','w') as fp:
        if dir == '.':
            #jobname = '/'.join(os.getcwd().split('/')[-2:])
            jobname = os.getcwd().split('/')[-1]
        else:
            jobname = dir
        fp.write(r'#!/bin/bash -x'+'\n')
        fp.write(r'#SBATCH --job-name %s'%jobname+'\n')
        fp.write(r'#SBATCH --nodes=%d'%Num_of_Nodes+'\n')
        fp.write(r'#SBATCH --ntasks-per-node=%d'%Ntasks_per_node+'\n')
        fp.write(r'#SBATCH --output=mpi-out.%j'+'\n')
        fp.write(r'#SBATCH --error=mpi-err.%j'+'\n')
        fp.write(r'#SBATCH --time=24:00:00'+'\n')
        fp.write(r'#SBATCH --partition=booster'+'\n')
        fp.write(r'#SBATCH --mail-user mcj_buaa@qq.com'+'\n')
        fp.write(r'#SBATCH --mail-type ALL'+'\n')
        fp.write(r'module load Architecture/KNL'+'\n')
        fp.write(r'module load Intel ParaStationMPI'+'\n')
        fp.write(r'srun -N %d  /p/project/cjics21/MO/lammps_booster/src/lmp_iff -in in.run'%Num_of_Nodes+'\n')

    if dir == '.':
        os.system(r'sbatch -A jics21 job.slm')
        print 'job "%s" submitted!'%os.getcwd().split('/')[-1]
    else:
        os.system(r'cd %s; sbatch -A jics21 job.slm; cd ..'%dir)
        print 'job "%s" submitted!'%dir


