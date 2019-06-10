import sys
import subprocess
import numpy as np
import os
import argparse

from mpi4py import MPI

from geometry import Geometry
from input import Input
from taskmanager import TaskManager

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    p_required = parser.add_argument_group('required arguments')
    p_required.add_argument("-g",type=str,required=True,
                        help='Path to geometry files.',nargs=2,
                        metavar=('XYZ_FILE','MOL_FILE'))

    p_required.add_argument("-rs",type=str,required=True,
                        help='Path to script template.',
                        metavar=('SCRIPT_FILE'))

    p_required.add_argument("-tdir",type=str,required=True,
                        help='Temporary root directory.',
                        metavar=('TMP_PATH'))

    p_optional = parser.add_argument_group('optional arguments')

    p_optional.add_argument('-EGH_out', nargs='?', type=argparse.FileType('w'),
                            help='EGH output file. (default: print to stdout)',
                            default=sys.stdout)

    p_optional.add_argument('-run_out', nargs='?', type=str,
                            help='Run script output file. (default: output.out)',
                            default="output.out")

    p_optional.add_argument('--calc_grad', action='store_true',
                        help='Calculate gradients.')
    p_optional.add_argument('--calc_hess', action='store_true',
                        help='Calculate hessians.')

    p_optional.add_argument('--z_sym', action='store_true',
                        help='Indicate that molecule has z symmetry.')

    p_optional.add_argument('-incr', type=float, default=1.e-3,
                        help='Incriment for geometry displacements. (default: 1.e-3)')

    #parse arguments
    args = parser.parse_args()

    geom_xyz_fname, geom_mol_fname = args.g   #geometry files

    template_script_fname = args.rs           #template script file
    tmp_dir = args.tdir                       #temporary directory

    EGH_out_fname = args.EGH_out              #EGH program output
    run_out_fname = args.run_out              #run script output

    calc_energy = True                        #calculate reference energy
    calc_force  = args.calc_grad              #calculate gradients
    calc_hess   = args.calc_hess              #calculate hessians

    z_sym = args.z_sym                        #z symmetry of a molecule
    incr  = float(args.incr)                  #displacements increment

    #initialize mpi variables
    mpi_comm = MPI.COMM_WORLD
    mpi_size = mpi_comm.Get_size()  #number of processes
    mpi_rank = mpi_comm.Get_rank()  #rank of current process
    mpi_master = 0                  #master process

    if mpi_rank == mpi_master:  #print only on master node
        print(
        "----------------------\n"
        "--- EGHrun program ---\n"
        "----------------------\n"
        "  Reference geometry will be loaded from the files:\n"
        "    %s and %s\n"
        "  Calculations will be performed in the directory:\n"
        "    %s\n"
        "  Template run script:\n"
        "    %s\n" %
                (geom_xyz_fname,geom_mol_fname,
                 tmp_dir,
                 template_script_fname))

    ref_geom = Geometry()
    ref_geom.read_xyz(geom_xyz_fname)
    ref_geom.read_mol(geom_mol_fname)

    #create task manager
    task_manager = TaskManager(tmp_dir,template_script_fname,run_out_fname)

    #prepare job list for calculations
    if mpi_rank == mpi_master:  #preparation is performed on master node
        #here we form full list of jobs
        job_list_full = []
        if calc_energy == True:  #create task for reference geometry
            job_list_full += task_manager.get_task_ref_geom(ref_geom)
        if calc_force == True:   #create list of tasks for gradients
            job_list_full += task_manager.get_task_list_grad(ref_geom,incr=incr,z_sym=z_sym)
        if calc_hess == True:    #create list of tasks for hessians
            job_list_full += task_manager.get_task_list_hess(ref_geom,incr=incr,z_sym=z_sym)

        #calculation of numbers of jobs per node
        n_jobs_all = len(job_list_full)                         #total number of jobs
        n_remain = n_jobs_all % mpi_size                        #number of unnisigned jobs
        n_main = int((n_jobs_all - n_remain) / mpi_size)        #number of jobs per process (without remaining)

        print("Total number of jobs: %s" % n_jobs_all)
        print("Number of processes requested: %s" % mpi_size)
        print("Job list size: %s, Remaining jobs: %s" % (n_main,n_remain))

        job_list = []
        n_local_jobs = np.zeros(mpi_size,int)
        i_job = 0

        #loop over processes
        for i_process in range(mpi_size):
            job_local = []

            #loop over block size
            for i_el in range(n_main):
                job_local.append(job_list_full[i_job])
                n_local_jobs[i_process] += 1
                i_job += 1

            #append remaining
            if i_process < n_remain:
                job_local.append(job_list_full[i_job])
                n_local_jobs[i_process] += 1
                i_job += 1

            print("Length of task list: %s for process number: %s" % (n_local_jobs[i_process],i_process+1))

            #append task to task list
            job_list.append(job_local)
    else:
        job_list = None

    #distribute tasks between processes
    local_job_list = mpi_comm.scatter(job_list,root=mpi_master)

    #perform calculations
    local_results = task_manager.calc(local_job_list)

    #collect results from all processes
    results_all = mpi_comm.gather(local_results)

    if mpi_rank == mpi_master:
        #analyze results
        out_str = task_manager.analyze(results_all,
            print_energy=calc_energy,print_grad=calc_force,print_hess=calc_hess)

        #print output
        EGH_out_fname.write(out_str)
        EGH_out_fname.close()