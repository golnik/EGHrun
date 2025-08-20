import sys
import subprocess
import numpy as np
import os
import argparse
import copy

from mpi4py import MPI

from geometry import Geometry
from input import Input
from taskmanager import TaskManager

bohr2A  = 0.529177
au2eV   = 27.211386
au2cm_1 = 219474.6313705
c = 137.035999

# Global error handler
def global_except_hook(exctype, value, traceback):
    import sys
    try:
        import mpi4py.MPI
        sys.stderr.write("\n*****************************************************\n")
        sys.stderr.write("Uncaught exception was detected on rank {}. \n".format(
            mpi4py.MPI.COMM_WORLD.Get_rank()))
        from traceback import print_exception
        print_exception(exctype, value, traceback)
        sys.stderr.write("*****************************************************\n\n\n")
        sys.stderr.write("\n")
        sys.stderr.write("Calling MPI_Abort() to shut down MPI processes...\n")
        sys.stderr.flush()
    finally:
        try:
            import mpi4py.MPI
            mpi4py.MPI.COMM_WORLD.Abort(1)
        except Exception as e:
            sys.stderr.write("*****************************************************\n")
            sys.stderr.write("Sorry, we failed to stop MPI, this process will hang.\n")
            sys.stderr.write("*****************************************************\n")
            sys.stderr.flush()
            raise e

sys.excepthook = global_except_hook

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

    p_optional.add_argument('-EGH_out', nargs='?', type=str,
                            help='EGH output file. (default: print to stdout)',
                            default="EGH.out")

    p_optional.add_argument('-run_out', nargs='?', type=str,
                            help='Run script output file. (default: output.out)',
                            default="output.out")

    p_optional.add_argument('-modes', nargs='?', type=str,
                            help='File with displacements along which gradients and hessians will be evaluted.')

    p_optional.add_argument('--calc_grad', action='store_true',
                        help='Calculate gradients.')
    p_optional.add_argument('--calc_hess', action='store_true',
                        help='Calculate hessians.')

    p_optional.add_argument('--z_sym', action='store_true',
                        help='Indicate that molecule has z symmetry.')

    p_optional.add_argument('-incr', type=float, default=1.e-2,
                        help='Incriment for geometry displacements. (default: 1.e-2)')

    p_optional.add_argument('--print_script_out', action='store_true',
                        help='Print output produced by script (if any).')

    #parse arguments
    args = parser.parse_args()

    geom_xyz_fname, geom_mol_fname = args.g   #geometry files

    template_script_fname = args.rs           #template script file
    tmp_dir = args.tdir                       #temporary directory

    EGH_out_fname = args.EGH_out              #EGH program output
    run_out_fname = args.run_out              #run script output

    modes_fname = args.modes                  #file with coord displacements

    calc_energy = True                        #calculate reference energy
    calc_force  = args.calc_grad              #calculate gradients
    calc_hess   = args.calc_hess              #calculate hessians

    z_sym = args.z_sym                        #z symmetry of a molecule
    incr  = float(args.incr)                  #displacements increment

    print_script_out=args.print_script_out          #print output from external script

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
        sys.stdout.flush()

    ref_geom = Geometry()
    ref_geom.read_xyz(geom_xyz_fname)
    ref_geom.read_mol(geom_mol_fname)

    n_coords = ref_geom.get_n_coords()

    #create list of "coordinates" that will be used for computing the derivatives
    coords = [] #array of geometries specifying the displacements
    dd     = [] #array of norms of the displacements' vectors

    #case of xyz displacements
    if modes_fname is None:
        for i_mode in range(n_coords):
            ngeom = copy.deepcopy(ref_geom)

            #set coordinates of all atoms to zero
            for i_coord in range(n_coords):
                ngeom.set_i_coord(i_coord,0.)

            val = 1.*incr
            ngeom.set_i_coord(i_mode,val)
            coords.append(ngeom)
    else:
        with open(modes_fname,'r') as file:
            for line in file:
                data = line.split()
                
                N = len(data)  #get number of elements in the line
                if N != (n_coords+1):
                    raise ValueError("Number of elements in each line of modes file must by 3*N + 1!")
                
                ngeom = copy.deepcopy(ref_geom)

                w = float(data[0])/au2cm_1  #convert normal modes frequency to au
                d = 1. #np.sqrt(1./(4. * np.pi**2 * c * w))

                #loop over coordinates in mode
                for i_coord in range(n_coords):
                    val = incr * d * float(data[i_coord+1])
                    ngeom.set_i_coord(i_coord,val)

                coords.append(ngeom)

    #generate norms of displacement vector
    masses = ngeom.atomic_masses
    Na = len(masses)    
    n_modes = len(coords)
    for i_mode in range(n_modes):
        res = 0.                
        indx = 0
        for ia in range(Na):
            for ixyz in range(3):
                res += (np.sqrt( masses[ia] ) * coords[i_mode].get_i_coord(indx))**2    # remove mass scaling from normal modes
                indx += 1
                
        norm = 2. * np.sqrt(res) / bohr2A
        dd.append(norm)

    #create task manager
    task_manager = TaskManager(tmp_dir,template_script_fname,run_out_fname,coords,dd)

    #prepare job list for calculations
    if mpi_rank == mpi_master:  #preparation is performed on master node
        #here we form full list of jobs
        job_list_full = []
        if calc_energy == True:  #create task for reference geometry
            job_list_full += task_manager.get_task_ref_geom(ref_geom)
        if calc_force == True:   #create list of tasks for gradients
            job_list_full += task_manager.get_task_list_grad(ref_geom)
        if calc_hess == True:    #create list of tasks for hessians
            job_list_full += task_manager.get_task_list_hess(ref_geom)

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

        print(
        "--------------------------\n"
        "--- Calculations start ---\n"
        "--------------------------\n")
        sys.stdout.flush()
    else:
        job_list = None

    #distribute tasks between processes
    local_job_list = mpi_comm.scatter(job_list,root=mpi_master)

    #perform calculations
    local_results, states = task_manager.calc(local_job_list,print_script_out)
    nstates = len(states)

    #collect results from all processes
    results_all = mpi_comm.gather(local_results)

    if mpi_rank == mpi_master:
        #analyze results
        energy, grad, hess = task_manager.analyze(results_all, nstates,
            print_energy=calc_energy,print_grad=calc_force,print_hess=calc_hess)

        for ist in range(nstates):
            state_str = states[ist]
            fname = "%s.%s" % (EGH_out_fname,state_str)

            with open(fname,'w') as outfile:
                if calc_energy == True:
                    outfile.write("$energy\n")
                    outfile.write("%15.8f\n" % energy[ist])

                if calc_force == True:
                    outfile.write("$gradient\n")
                    for i_mode in range(n_modes):
                        outfile.write("%15.8f" % grad[i_mode][ist])
                    outfile.write("\n")

                if calc_hess == True:
                    outfile.write("$hessian\n")
                    for i_mode in range(n_modes):
                        for j_mode in range(n_modes):
                            outfile.write("%15.8f" % hess[i_mode][j_mode][ist])
                        outfile.write("\n")