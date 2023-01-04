import os

dirs_list = ['.']
for dir in dirs_list:
    with open(dir+'/job.pbs','w') as fp:
        fp.write(r'#PBS -S /bin/bash'+'\n')
        fp.write(r'#PBS -e mpi.err'+'\n')
        fp.write(r'#PBS -o mpi.log'+'\n')
        if dir == '.':
            fp.write(r'#PBS -N %s'%os.getcwd().split('/')[-1]+'\n')
        else:
            fp.write(r'#PBS -N %s'%dir+'\n')
        fp.write(r'#PBS -l nodes=1:ppn=12'+'\n')
        fp.write(r'cd $PBS_O_WORKDIR'+'\n')
        fp.write(r'mpirun -np 12 /usr/users/iff_th2/cmo/lammps/src/lmp_iff -in ./in.run'+'\n')
    if dir == '.':
        os.system(r'qsub -q th2 job.pbs')
    else:
        os.system(r'cd %s; qsub -q th2 job.pbs; cd ..'%dir)
    print 'job "%s" submitted!'%dir
