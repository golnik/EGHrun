import numpy as np
import copy
import subprocess
import os

from input import Input
from reader import Reader

class TaskManager(object):
    def __init__(self,tmp_dir,template_script_fname,run_out_fname,modes,dd):
        '''
        Args:
            modes:   array of geometry instances specifying displacements vector
                     along which gradient and hessian will be computed
            disps:   vector of displacements  
        '''
        self.tmp_dir = tmp_dir
        self.template_script_fname = template_script_fname
        self.run_out_fname = run_out_fname
        self.modes = modes
        self.n_modes = len(self.modes)
        self.dd = dd

    def get_task_ref_geom(self,geom):
        '''
        Return list containing one task with reference geometry.
        '''
        tmp_dir = os.path.join(self.tmp_dir,"ref_geom")
        task_type = 0
        res = [geom,tmp_dir,task_type]
        return [res]

    def get_task_list_grad(self,geom):
        '''
        Return list containing tasks to calculate gradients
        '''
        res_list = [] #results list

        counter = 0
        #loop over modes
        for i_mode in range(self.n_modes):
            #loop over plus and minus increment
            for sign in [-1,+1]:
                ngeom = copy.deepcopy(geom)

                #loop over coordinates within mode
                for i_coord in range(ngeom.get_n_coords()):
                    coord = ngeom.get_i_coord(i_coord)
                    coord += sign * self.modes[i_mode].get_i_coord(i_coord)
                    ngeom.set_i_coord(i_coord,coord)

                tmp_dir = os.path.join(self.tmp_dir,"grad_%s"%counter)

                task_type = 1
                res = [ngeom,tmp_dir,task_type,[i_mode,sign]]
                res_list.append(res)

                counter += 1

        return res_list

    def get_task_list_hess(self,geom):
        '''
        Return list containing tasks to calculate hessian
        '''
        res_list = [] #results list

        counter = 0
        #loop over modes
        for i_mode in range(self.n_modes):
            for j_mode in range(i_mode,self.n_modes):
                #loop over plus and minus increment
                for i_sign in [-1,+1]:
                    for j_sign in [-1,+1]:
                        ngeom = copy.deepcopy(geom)

                        #loop over coordinates within mode
                        for i_coord in range(ngeom.get_n_coords()):
                            coord = ngeom.get_i_coord(i_coord)
                            coord += i_sign * self.modes[i_mode].get_i_coord(i_coord)
                            ngeom.set_i_coord(i_coord,coord)

                        for j_coord in range(ngeom.get_n_coords()):
                            coord = ngeom.get_i_coord(j_coord)
                            coord += j_sign * self.modes[j_mode].get_i_coord(j_coord)
                            ngeom.set_i_coord(j_coord,coord)

                        tmp_dir = os.path.join(self.tmp_dir,"hess_%s"%counter)

                        task_type = 2
                        res = [ngeom,tmp_dir,task_type,[i_mode,i_sign,
                                                        j_mode,j_sign]]
                        res_list.append(res)

                        counter += 1

        return res_list

    def calc(self,task_list,print_out=False):
        '''
        This function runs calculations for specified task list.
        '''
        res_list = []   #list to store results
        nstates  = 0    #number of electronic states

        #loop over tasks in list
        for task in task_list:
            geom = task[0]      #geometry
            tmp_dir = task[1]   #directory to perform calculations
            task_type = task[2] #task type

            print("Calculation is started in tmp_dir: %s" % tmp_dir)

            #create tmp directory
            if not os.path.exists(tmp_dir):
                os.makedirs(tmp_dir)

            #specify outfile file with energy
            run_out_fname = os.path.join(tmp_dir,self.run_out_fname)

            #create input object
            inp_fname = os.path.join(tmp_dir,"input.inp")

            inp = Input(geom,tmp_dir,run_out_fname)
            inp.subst_template(self.template_script_fname,inp_fname)

            #initiate subprocess for calculations
            process = subprocess.Popen(['bash',inp_fname],
                                       stdin=subprocess.PIPE,
                                       stdout=subprocess.PIPE,
                                       bufsize=1,universal_newlines=True)
            #process.wait()
            while process.poll() is None:
                out = process.stdout.readline()
                if print_out==True:     #print output on-the-fly
                    print(out,end='')

            #here we parse run_output
            run_out_str = ''
            try:
                with open(run_out_fname,'r') as run_out_stream:
                    run_out_str = run_out_stream.read()
            except:
                raise Exception("File %s with energy output is not found." % run_out_fname)

            #read output file
            reader = Reader()
            energies, states = reader.get_energy(run_out_str)  #get calculated energies

            #result list
            res = [energies,task_type]

            #add aditional arguments to result list
            if task_type != 0:
                args = task[3]
                res.append(args)

            print("Calculation is performed. Results: %s" % res)

            res_list.append(res)

        return res_list, states

    def analyze(self,res_list,nstates,print_energy=True,print_grad=True,print_hess=True):
        energy_ref = np.zeros(nstates)

        if print_grad == True:
            Ep = np.zeros((self.n_modes,nstates))
            Em = np.zeros((self.n_modes,nstates))

        if print_hess == True:
            Upp = np.zeros((self.n_modes,self.n_modes,nstates))
            Ump = np.zeros((self.n_modes,self.n_modes,nstates))
            Upm = np.zeros((self.n_modes,self.n_modes,nstates))
            Umm = np.zeros((self.n_modes,self.n_modes,nstates))

        #loop over list of results
        for local_list in res_list:
            #loop over results in a particular list
            for result in local_list:
                energy = result[0]      #get energies
                task_type = result[1]   #get task type

                if task_type == 1:   #gradients
                    i_mode = result[2][0]
                    sign = result[2][1]

                    if(sign>0):
                        Ep[i_mode] = energy
                    else:
                        Em[i_mode] = energy
                elif task_type == 2: #hessians
                    i_mode = result[2][0]
                    i_sign = result[2][1]
                    j_mode = result[2][2]
                    j_sign = result[2][3]

                    if(i_sign>0 and j_sign>0):
                        Upp[i_mode][j_mode] = energy
                    elif(i_sign>0 and j_sign<0):
                        Upm[i_mode][j_mode] = energy
                    elif(i_sign<0 and j_sign>0):
                        Ump[i_mode][j_mode] = energy
                    elif(i_sign<0 and j_sign<0):
                        Umm[i_mode][j_mode] = energy
                else:
                    energy_ref = energy

        grad = np.zeros((self.n_modes,nstates))
        hess = np.zeros((self.n_modes,self.n_modes,nstates))
        
        #calculate gradients
        if print_grad == True:
            for i_mode in range(self.n_modes):
                grad[i_mode] = (Ep[i_mode] - Em[i_mode]) / self.dd[i_mode]

        #calculate hessians
        if print_hess == True:
            for i_mode in range(self.n_modes):
                for j_mode in range(i_mode,self.n_modes):
                    hess[i_mode][j_mode] = (Upp[i_mode][j_mode] - Ump[i_mode][j_mode]
                                           -Upm[i_mode][j_mode] + Umm[i_mode][j_mode]) / (self.dd[i_mode] * self.dd[j_mode])
                    #hessian matrix is symmetric
                    hess[j_mode][i_mode] = hess[i_mode][j_mode]

        return energy_ref, grad, hess