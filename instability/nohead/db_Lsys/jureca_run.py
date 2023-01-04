import os

Num_of_Nodes = 1
Ntasks_per_node = 24
prj = 'jiff44'
#IND = [9.0,9.25,9.5,9.75,10.0,10.25,10.5,10.75,11.0,11.25,11.5,11.75,12.0]
IND = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]
IND = [6.25,6.5,6.75,7.25,7.5,7.75]
cwd = os.getcwd()
for ind in IND:
    dir = cwd+'/db'+str(ind)
    os.chdir(dir)
    infiles = ['in-1.run','in-2.run','in-3.run','in-4.run']
    for cn,infile in enumerate(infiles):
        with open('job.slm','w') as fp:
            jobname = dir.split('/')[-1]+'-'+str(cn)
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

        os.system(r'sbatch -A %s job.slm'%prj)
        print( 'job "%s" submitted!'%jobname )
    os.chdir(cwd)


