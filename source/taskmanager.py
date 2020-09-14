import numpy as np
import copy
import subprocess
import os

from input import Input
from reader import Reader

class TaskManager(object):
    def __init__(self,tmp_dir,template_script_fname,run_out_fname):
        self.tmp_dir = tmp_dir
        self.template_script_fname = template_script_fname
        self.run_out_fname = run_out_fname

    def get_task_ref_geom(self,geom):
        '''
        Return list containing one task with reference geometry.
        '''
        tmp_dir = os.path.join(self.tmp_dir,"ref_geom")
        task_type = 0
        res = [geom,tmp_dir,task_type]
        return [res]

    def get_task_list_grad(self,geom,incr=0.001,z_sym=False):
        '''
        Return list containing tasks to calculate gradients
        '''
        self.incr = incr
        self.n_coords = geom.get_n_coords()

        res_list = [] #results list

        counter = 0
        #loop over coordinates
        for i_coord in range(self.n_coords):
            i_rem = (i_coord+1) % 3
            #z symmetry conditions for gradients
            if z_sym == True:
                if i_rem == 0:
                    continue

            #loop over plus and minus increment
            for increment in [-incr,+incr]:
                ngeom = copy.deepcopy(geom)

                coord = ngeom.get_i_coord(i_coord)
                coord += increment
                ngeom.set_i_coord(i_coord,coord)

                tmp_dir = os.path.join(self.tmp_dir,"grad_%s"%counter)

                task_type = 1
                res = [ngeom,tmp_dir,task_type,[i_coord,int(np.sign(increment))]]
                res_list.append(res)

                counter += 1

        return res_list

    def get_task_list_hess(self,geom,incr=0.001,z_sym=False):
        '''
        Return list containing tasks to calculate hessian
        '''
        self.incr = incr
        self.n_coords = geom.get_n_coords()

        res_list = [] #results list

        counter = 0
        #loop over coordinates
        for i_coord in range(self.n_coords):
            i_rem = (i_coord+1) % 3
            for j_coord in range(i_coord,self.n_coords):
                j_rem = (j_coord+1) % 3

                #z symmetry conditions for hessians
                if z_sym == True:
                    if i_rem == 0 or j_rem == 0:
                        if i_rem == j_rem:  #if both coordinates are z-symmetric
                            if i_coord == j_coord:   #identical coordinates
                                tmp_dir = os.path.join(self.tmp_dir,"hess_%s"%counter)

                                #special case with no displacements
                                task_type = 2
                                res = [geom,tmp_dir,task_type,[i_coord,0,
                                                               j_coord,0]]
                                res_list.append(res)
                            else:                   #z coordinates of different atoms
                                ngeom = copy.deepcopy(geom)

                                #one of the coordinates should be with positive displacement
                                coord_i = ngeom.get_i_coord(i_coord)
                                coord_i += incr
                                ngeom.set_i_coord(i_coord,coord_i)

                                #another coordinate with negative displacement
                                coord_j = ngeom.get_i_coord(j_coord)
                                coord_j -= incr
                                ngeom.set_i_coord(j_coord,coord_j)

                                tmp_dir = os.path.join(self.tmp_dir,"hess_%s"%counter)

                                task_type = 2
                                res = [ngeom,tmp_dir,task_type,[i_coord,0,
                                                                j_coord,0]]
                                res_list.append(res)
                            counter += 1
                        else:           #if only one coordinate is z-symmetric, hessian is zero
                            continue

                #loop over plus and minus increment
                for i_incr in [-incr,+incr]:
                    for j_incr in [-incr,+incr]:
                        ngeom = copy.deepcopy(geom)

                        coord_i = ngeom.get_i_coord(i_coord)
                        coord_i += i_incr
                        ngeom.set_i_coord(i_coord,coord_i)

                        coord_j = ngeom.get_i_coord(j_coord)
                        coord_j += j_incr
                        ngeom.set_i_coord(j_coord,coord_j)

                        tmp_dir = os.path.join(self.tmp_dir,"hess_%s"%counter)

                        task_type = 2
                        res = [ngeom,tmp_dir,task_type,[i_coord,int(np.sign(i_incr)),
                                                        j_coord,int(np.sign(j_incr))]]
                        res_list.append(res)

                        counter += 1

        return res_list

    def calc(self,task_list,print_out=False):
        '''
        This function runs calculations for specified task list.
        '''
        res_list = []   #list to store results

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
            energy = reader.get_energy(run_out_str)  #get calculated energy

            #result list
            res = [energy,task_type]

            #add aditional arguments to result list
            if task_type != 0:
                args = task[3]
                res.append(args)

            print("Calculation is performed. Results: %s" % res)

            res_list.append(res)

        return res_list

    def analyze(self,res_list,print_energy=True,print_grad=True,print_hess=True):
        energy_ref = 0.

        if print_grad == True:
            Ep = np.zeros(self.n_coords)
            Em = np.zeros(self.n_coords)

        if print_hess == True:
            Upp = np.zeros((self.n_coords,self.n_coords))
            Ump = np.zeros((self.n_coords,self.n_coords))
            Upm = np.zeros((self.n_coords,self.n_coords))
            Umm = np.zeros((self.n_coords,self.n_coords))

        #loop over list of results
        for local_list in res_list:
            #loop over results in a particular list
            for result in local_list:
                energy = result[0]      #get energy
                task_type = result[1]   #get task type

                if task_type == 1:   #gradients
                    i_coord = result[2][0]
                    sign = result[2][1]

                    if(sign>0):
                        Ep[i_coord] = energy
                    else:
                        Em[i_coord] = energy
                elif task_type == 2: #hessians
                    i_coord = result[2][0]
                    i_sign = result[2][1]
                    j_coord = result[2][2]
                    j_sign = result[2][3]

                    if(i_sign>0 and j_sign>0):
                        Upp[i_coord][j_coord] = energy
                    elif(i_sign>0 and j_sign<0):
                        Upm[i_coord][j_coord] = energy
                    elif(i_sign<0 and j_sign>0):
                        Ump[i_coord][j_coord] = energy
                    elif(i_sign<0 and j_sign<0):
                        Umm[i_coord][j_coord] = energy
                    elif(i_sign==0 and j_sign==0):  #special case, z symmetry
                        Upm[i_coord][j_coord] = energy
                        Ump[i_coord][j_coord] = energy
                else:
                    energy_ref = energy

        #geometries are in angstrom. so we sould convert them to au.
        if print_grad == True or print_hess == True:
            bohr2A = 0.529177
            A2bohr = 1.88973
            dx = 2. * self.incr * A2bohr

        #calculate gradients
        if print_grad == True:
            grad = np.zeros(self.n_coords)
            for i_coord in range(self.n_coords):
                grad[i_coord] = (Ep[i_coord] - Em[i_coord]) / dx

        #calculate hessians
        if print_hess == True:
            hess = np.zeros((self.n_coords,self.n_coords))
            for i_coord in range(self.n_coords):
                for j_coord in range(i_coord,self.n_coords):
                    hess[i_coord][j_coord] = (Upp[i_coord][j_coord] - Ump[i_coord][j_coord]
                                             -Upm[i_coord][j_coord] + Umm[i_coord][j_coord]) / dx**2
                    #hessian matrix is symmetric
                    hess[j_coord][i_coord] = hess[i_coord][j_coord]

        #output
        out_str = ''

        if print_energy == True:
            out_str += "$energy\n  "
            out_str += "%s\n" % energy_ref

        if print_grad == True:
            out_str += "$gradient\n  "
            for i_coord in range(self.n_coords):
                out_str += "%s " % grad[i_coord]
            out_str += "\n"

        if print_hess == True:
            out_str += "$hessian\n"
            for i_coord in range(self.n_coords):
                out_str += "  "
                for j_coord in range(self.n_coords):
                    out_str += "%s " % hess[i_coord][j_coord]
                out_str += "\n"

        return out_str