#!/usr/bin/ipython

# Different tests of environments
#--------------------------------
# Global system environment
# #!/usr/bin/python
# #!/usr/bin/ipython

# # Virtualenv environment
# #!/localhome/david/scripts/vasp/AFM/AFM_plot/bin/ipython 

# pipenv evironment
# #!/localhome/david/.local/share/virtualenvs/AFM-GusOYeYZ/bin/ipython

import os
import sys
import re
import subprocess
from distutils.util import strtobool
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from scipy.integrate import  simps
from matplotlib.widgets import Button, Slider
from matplotlib.colors import Colormap

import warnings
warnings.filterwarnings("ignore")

tip_name = "AFM tip"
prefix = "wip"

#  # KKKK--------------------------------------------------------------------------- 
#  def print_forces():
#      global dist_interval
#      ff = open('temp_forces.dat', 'w')
#      fd = open('temp_dist.dat', 'w')
#      for i in range(len(dist_interval)):
#          ff.write(str(-force_spl(dist_interval[i]))+"\n")
#          fd.write(str(dist_interval[i])+"\n")
#      ff.close()
#      fd.close()
#  #-------------------------------------------------------------------------------- 

# Experimental parameters
#------------------------
A = 1        # Amplitud of the oscillation of the experimental cantiliever (in Ang)
# w_0 = 23000   # Natural resonance frequency of the cantilever (in Hz) (f in reality)
w_0 = 23000/(2*np.pi)   # Natural resonance frequency of the cantilever  Times 2pi (in Hz)
#w_0 = 150   # Natural resonance frequency of the cantilever (in Hz)
#k = 40      # Spring constant (N/m)
k = 40*0.0624      # Spring constant (eV/ang**2)
#k = 4        # Spring constant (nN/Ang)
#k = 4*0.0624        # Spring constant (eV/Ang**2)
omega_exp = -5

show_graphs = True
save_graphs = True
show_afmimage = True
save_afmimage = True
autosave = True

#      GP  Description
#      1   On As atom 1sst layer
#      2   -.
#      3    .
#      4    .
#      5    .
#      6    .- 2nd layer
#      7    .
#      8    .
#      9   -.
#      10  On Ga atom 1sst layer

# X,Y coordinates of the grid surface
grid = {1: (0.00000, 0.00000),
        2: (2.02290, 0.00000),
        3: (0.00000, 1.12780),
        4: (2.02290, 1.12780),
        5: (0.00000, 2.25560),
        6: (2.02290, 2.25560),
        7: (0.00000, 3.38340),
        8: (2.02290, 3.38340),
        9: (0.00000, 4.51120),
        10: (2.02290, 4.51120)}

afmimage = {1: [], # Object containing (x,y) grid points, [Z_i] and the [ delta_w ] points 
         2: [],
         3: [],
         4: [],
         5: [],
         6: [],
         7: [],
         8: [],
         9: [],
         10:[]}
afmforces = {1: [], # Object containing (x,y) grid points, [Z_i] and the force curve at each grid point 
          2: [],
          3: [],
          4: [],
          5: [],
          6: [],
          7: [],
          8: [],
          9: [],
          10:[]}

imageAtW = {1: [], # Object containing the  (x,y) grid points and the Z value at a given delta W
            2: [],
            3: [],
            4: [],
            5: [],
            6: [],
            7: [],
            8: [],
            9: [],
            10:[]}
afm_image = []

# Used to repeat the AFM image
x_interval = 2.02290 # Distance betwwen x grid points
y_interval = 1.12780 # Distance betwwen x grid points
nxgp = 2             # Number of original X Grid points
nygp = 5             # Number of original Y Grid points
xrep = 6             # Number of repetition of X Grid points
yrep = 3             # Number of repetition of Y Grid points


def string_to_bool (string):
    return bool(strtobool(str(string)))


def get_atom_index (kind):
    """
    Returns the atom index with higher Z if kind == tip or atom index of the
    one with lowest Z if kind == surface
    TODO: So far only returns the biggest and the lowest Z value regardless of the atom kind.
          - Can be improved by ignoring H in the index selection
    Reads the POSCAR (for VASP 5.x)
    """
    poscar = "POSCAR.1"

    with open (poscar,'r') as f:
        # Skipping the comment
        f.readline()

        acell = float(f.readline())
        cell = []
        for i in  range (0,3):
            cell.append(f.readline().split())
        cell = np.array(cell)         # Cell will be an array of strings
        cell = cell.astype(np.float)  # Transform cell into an array of floats

        # read kind of atoms
        at_kinds = f.readline().split()

        # read number of each kind of atoms
        no_at_each_kind = f.readline().split()
        total_atoms = 0
        for i in  no_at_each_kind:
            total_atoms = total_atoms + int(i)
        f.readline()
        coord_type = f.readline()

        # If coordinates are in direct coordinates, transform to cartesians
        transform = False
        if coord_type[0].lower() == "d":
            transform = True

        idx = 1
        coordinates = []
        constraints = []
        while idx <= total_atoms:
            line = f.readline()
            if line != "\n":
                coordinates_i = line.split()[0:3]; constraints_i = line.split()[3:7]
                if transform:
                    coordinates_i = transform_coordiantes ( cell, coordinates_i, acell)
                coordinates.append([idx, coordinates_i])
                constraints.append([idx, constraints_i])
                idx +=1

        atom_index = 1
        z = coordinates[0][1][2]
        for atom in coordinates:
            if kind == 'tip':
                if atom[1][2] > float(z):
                    atom_index = atom[0]
                    z = atom[1][2]
            else:
                if atom[1][2] < float(z):
                    atom_index = atom[0]
                    z = atom[1][2]
        # print (kind + " index: ", atom_index)

        return  atom_index


def get_grid_point ():
    grid_point = os.getcwd().split("_")[-1]
    return grid_point


def error_critical (msg = "Unknown critical error"):
    sys.exit("ERROR: " + msg )


def calculate_distance (coordinates, at_1, at_2 ):
    for i in range (0, len(coordinates)):
        idx = coordinates[i][0]
        if idx == at_1:
            coord_at_1 = coordinates[i][1]
        if idx == at_2:
            coord_at_2 = coordinates[i][1]
    distance = abs(float(coord_at_1[2]) - float(coord_at_2[2]))
    return distance


def transform_coordiantes (cell, coordinates, acell = 1.0):
    """ Transform coordinates from Direct to Cartesian """ 

    # Convert the values to numpy arrays
    cell = np.array(cell)
    coordinates = np.array(coordinates)
    coordinates = coordinates.astype(np.float)
    new_coordinates = []

    if coordinates.ndim == 1:
        x = coordinates[0] * cell[0,:] * acell
        y = coordinates[1] * cell[1,:] * acell
        z = coordinates[2] * cell[2,:] * acell
        new_coordinates = x + y + z
        return new_coordinates
    else:
        # TODO: This block does not work.
        # If coordinates is an array of coordinates of several atoms, the function
        # does not work
        for i in range (0,len(coordinates)):
            x = coordinates[i,0] * cell[0,:] * acell
            y = coordinates[i,1] * cell[1,:] * acell
            z = coordinates[i,2] * cell[2,:] * acell
            coord_i = x + y + z
            new_coordinates.append(coord_i)
        return new_coordinates


def read_poscar (poscar):
    """ Reads the POSCAR (for VASP 5.x)"""

    with open (poscar,'r') as f:
        # Skipping the comment
        f.readline()

        acell = float(f.readline())
        cell = []
        for i in  range (0,3):
            cell.append(f.readline().split())
        cell = np.array(cell)         # Cell will be an array of strings
        cell = cell.astype(np.float)  # Transform cell into an array of floats

        # read kind of atoms
        at_kinds = f.readline().split()

        # read number of each kind of atoms
        no_at_each_kind = f.readline().split()
        total_atoms = 0
        for i in  no_at_each_kind:
            total_atoms = total_atoms + int(i)

        f.readline()
        coord_type = f.readline()
        # If coordinates are in direct coordinates, transform to cartesians
        transform = False
        if coord_type[0].lower() == "d":
            transform = True
            # print "coordinates in file", poscar, "are Direct. Will be transformed to cartesian"

        idx=1
        coordinates = []
        constraints = []
        while idx <= total_atoms:
            line = f.readline()
            if line != "\n":
                coordinates_i = line.split()[0:3]; constraints_i = line.split()[3:7]
                if transform:
                    coordinates_i = transform_coordiantes ( cell, coordinates_i, acell)
                coordinates.append([idx, coordinates_i])
                constraints.append([idx, constraints_i])
                idx +=1

        return  coordinates, constraints


def calc_force (distances, energies, method = "direct", pol_rank = 12):
    global f_interp
    if method == "poly":
        # Fit Distance-Energy curve to a polynomy of rank <pol_rank>
        # Overwriting pol_rank value (shall I keep it in this way?)
        pol_rank = 12
        p = np.poly1d(np.polyfit(distances, energies,pol_rank))
        # Compute the -derivative of the polynomy (the force)
        pp = -np.polyder(p)
        force = pp
        f_interp = force
        #return p, force

    if method=="pol_interp4":
        force = []
        for i in range (4, len(distances)-4):
            section_x = distances[i-4:i+5]
            section_y = energies[i-4:i+5]
            # Overwriting pol_rank value (shall I keep it in this way?)
            pol_rank = 5
            # Fit section of Distance-Energy curve to a polynomy of rank <pol_rank>
            p = np.poly1d(np.polyfit(section_x, section_y,pol_rank))
            # Compute the -derivative of the polynomy (the force)
            pp = -np.polyder(p)
            f_i = np.polyval(pp, distances[i])
            force.append(f_i)
        f_interp = force
        #return force

    if method == "pol_interp3":
        force = []
        for i in range (3, len(distances)-3):
            section_x = distances[i-3:i+4]
            section_y = energies[i-3:i+4]
            # Overwriting pol_rank value (shall I keep it in this way?)
            pol_rank = 5
            # Fit section of Distance-Energy curve to a polynomy of rank <pol_rank>
            p = np.poly1d(np.polyfit(section_x, section_y,pol_rank))
            # Compute the -derivative of the polynomy (the force)
            pp = -np.polyder(p)
            f_i = np.polyval(pp, distances[i])
            force.append(f_i)
        f_interp = force
        #return force

    if method == "pol_interp2":
        force = []
        for i in range (2, len(distances)-2):
            section_x = distances[i-2:i+3]
            section_y = energies[i-2:i+3]
            # Overwriting pol_rank value (shall I keep it in this way?)
            pol_rank = 5
            # Fit section of Distance-Energy curve to a polynomy of rank <pol_rank>
            p = np.poly1d(np.polyfit(section_x, section_y,pol_rank))
            # Compute the -derivative of the polynomy (the force)
            pp = -np.polyder(p)
            f_i = np.polyval(pp, distances[i])
            force.append(f_i)
        f_interp = force
        #return force

    if method == "pol_interp":
        force = []
        for i in range (1, len(distances)-1):
            section_x = distances[i-1:i+2]
            section_y = energies[i-1:i+2]
            # Overwriting pol_rank value (shall I keep it in this way?)
            pol_rank = 5
            # Fit section of Distance-Energy curve to a polynomy of rank <pol_rank>
            p = np.poly1d(np.polyfit(section_x, section_y,pol_rank))
            # Compute the -derivative of the polynomy (the force)
            pp = -np.polyder(p)
            f_i = np.polyval(pp, distances[i])
            force.append(f_i)
        f_interp = force
        #return force

    if method == "num1":
        # Compute the forces via numerical derivative
        dE = []
        for i in range (0, len(distances)-1):
            dE_i =  -(energies[i] - energies[i+1]) / (distances[i] - distances[i+1])
            dE.append(dE_i)
        dE_fit =  np.poly1d (np.polyfit(distances[:-1],dE,pol_rank))
        force = dE
        f_interp = force
        #return force

    if method == "num2":
        # Using (f(x-h) - f(x+h)) / 2h
        dE_2 = []
        for i in range (1, len(distances)-1):
            dE_i =  -(energies[i-1] - energies[i+1]) / (distances[i-1] - distances[i+1])
            dE_2.append(dE_i)
        dE_2_fit =  np.poly1d (np.polyfit(distances[1:-1],dE_2,pol_rank))
        force = dE_2
        f_interp = force
        #return force

    if method == "spline":
        global force_spl
        global e_fit_spl
        e_fit_spl = None
        force_spl = None
        # Fit the curve E(z) to a spline
        e_fit_spl = UnivariateSpline(distances, energies, s=spl_smoothing)
        # Compute the first derivative of the spline
        force_spl = e_fit_spl.derivative(1) # Remember to use -force_spl later!
        #return force, e_fit


def read_outcar(filename = 'OUTCAR'):
    """Read OUTCAR type file.
    Reads unitcell, atom positions, energies, and forces from the OUTCAR file.
    CAREFUL: does not explicitly read constraints (yet?)
    Based on "No recuerdo el nombre del Author"'s script
    """

    if isinstance(filename, str):
        f = open(filename)
    else: # Assume it's a file-like object
        f = filename
    data    = f.readlines()

    natoms  = 0
    images  = []
    # atoms   = Atoms(pbc = True)
    energy  = 0
    atoms   = []
    forces  = []
    species = []
    symbols = []
    species_num = []

    for n,line in enumerate(data):
        if 'VRHFIN' in line:
            temp = line.split('=')
            species.append(temp[1][0:2].strip(':'))
        if 'ions per type' in line:
            temp = line.split()
            for ispecies in range(len(species)):
                species_num += [int(temp[ispecies+4])]
                natoms += species_num[-1]
                for iatom in range(species_num[-1]): symbols += [species[ispecies]]
        if 'direct lattice vectors' in line:
            cell = []
            for i in range(3):
                temp = data[n+1+i].split()
                cell += [[float(temp[0]), float(temp[1]), float(temp[2])]]
        if 'FREE ENERGIE OF THE ION-ELECTRON SYSTEM' in line:
            energy = float(data[n+2].split()[4])
        if 'POSITION          ' in line:
            # Restart here the forces and atoms list in order to take just the last geometry
            forces = []
            atoms  = []
            for iatom in range(natoms):
                temp    = data[n+2+iatom].split()
                atoms  += [ (symbols[iatom], [float(temp[0]),float(temp[1]),float(temp[2])])]
                forces += [[float(temp[3]),float(temp[4]),float(temp[5])]]
            images = [(atoms, forces)]
            # Uncomment to store all the iterations
            # images += [(atoms, forces)]
    # TODO: Check that the OUTCAR has converged

    write_geom_and_forces (natoms, atoms, forces)

    return atoms, forces


def write_geom_and_forces (natoms,atoms,forces):

    jmol_script = 'jmolscript: vectors on; vectors 2; set  vectorscale 2; set percentVdwAtom 30;  set bondradiusmilliangstroms 120;'
    fout = open('temp.xyz', 'a')
    fout.write(str(natoms) + "\n")
    fout.write(jmol_script + "\n")
    for i in range (0,len(atoms)):
        # Sorry for the long line
        fout.write('%s  %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f \n' % ( atoms[i][0], atoms[i][1][0],  atoms[i][1][1], atoms[i][1][2], forces[i][0], forces[i][1], forces[i][2]))
    fout.close()

## TODO: Understand this: ##
##      Not sure what I wanted to do with this, but rememeber it took me ages to understand this structure that I've already forgotten ##
    # # If we store all the images, instead of the last one, <n> will be the no. of each iteration
    # n=1
    # print "images: ", len(images)
    # print "images[0][1] shape: ", images[0][0][0+natoms*(n-1):natoms+natoms*(n-1)]

    # images[0][0] = [( atom, [x , y , z])]  # atom and coordinates
    # images[0][1] = [[Fx , Fy , Fz]]        # Forces on atom
    # return images


def sum_forces (outcar, constraints):

    atoms, forces = read_outcar (outcar)
    vertical_force=0
    for i in range (0, len(forces)):
        if constraints[i][1][2].lower() == "f":
            # print atoms[1] , forces[i], constraints[i][1]
            # print 'forces: ', forces[i][2],  constraints[i][1][2]
            vertical_force += forces[i][2]
    return vertical_force


def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]


def set_interval ():
    '''Set a E fitting interval between z_min and z_max'''

        # Select points between height (z) interval
    global dist_interval
    global ener_interval
    global outcar_interval

    # Delete previously set interval
    del dist_interval[:]
    del ener_interval[:]
    del outcar_interval[:]

    dist = []
    ener = []
    outc = []
    for i,z in enumerate(distances):
        if z >= z_min and z<= z_max:
            dist.append(distances[i])
            ener.append(energies[i])
            outc.append(outcar_files[i])
    dist_interval = dist
    ener_interval = ener
    outcar_interval = outc


def plot_outcar_labels():
    '''Plot outcar labels'''

    # Add labels to energy points
    idx_label = 0
    for i in range (len(ener_interval)):
        if (idx_label % label_interval == 0 ):
            ax1.text(dist_interval[i]+0.2,ener_interval[i],outcar_interval[i].replace("OUTCAR.",""),fontsize=8)
            ax2.text(dist_interval[i]+0.2,-force_spl(dist_interval[i]),outcar_interval[i].replace("OUTCAR.",""),fontsize=8)
            ax3.text(dist_interval[i]+0.2,ener_interval[i],outcar_interval[i].replace("OUTCAR.",""),fontsize=8)
        idx_label+=1


def replot_graphs():
    ''' Replot graphs'''

    # Clear previous data but not axes
    ax1.cla()
    ax2.cla()
    ax3.cla()

    # Padding 
    plt.tight_layout(pad=1.1)
    ax1.margins(0.05, 0.2)
    ax2.margins(0.05, 0.2)
    ax3.margins(0.05, 0.2)

    #plot_name = prefix+"_E_"+str(grid_point)+".pdf"
    ax1.plot(dist_interval, ener_interval, 'o', label="Energies", lw=1)
    ax1.plot(dist_interval, e_fit_spl(dist_interval), '--', label="E fit", lw=1)

    #plot_name = prefix+"_F_"+str(grid_point)+".pdf"
    #ax2.plot(dist_interval, -force(dist_interval), 'k-',  label="Force", lw=1)
    # Computed via with spline
    ax2.plot(dist_interval, -force_spl(dist_interval), 'k-',  label="Force", lw=1)

    #plot_name = prefix+"delta_w_"+str(grid_point)+".pdf"
    ax3.plot (h,delta_w, 'm-', label = 'omega', lw=1)
    #ax3.plot (h,delta_w, 'm-', lw=1.5)

    # Energy
    ax1.set_title('Energy', fontsize=10)
    ax1.set_xlabel (r'z ($\AA$)')
    ax1.set_ylabel ('Energy (ev)')

    # Force
    ax2.set_title('Force vs vertical distance', fontsize=10)
    ax2.set_xlabel (r'z ($\AA$)')
    ax2.set_ylabel (r'Force(ev/$\AA$)')

    # Omega
    ax3.set_title('Frequency shift', fontsize=10)
    ax3.set_xlabel (r'h ($\AA$)')
    ax3.set_ylabel (r'$\Delta \omega$ (Hz)')

    #plot_name = prefix+"_allgraphs_"+str(grid_point)+".pdf"
    if add_labels: plot_outcar_labels()
    if show_graphs: plt.draw()


def integrate_forces():
    '''Calculate the integral of the forces'''

    global w
    global h
    global delta_w

    # Empty the values before computing the integrals
    w = []  # w => (w_exp/w_0)**2
    h = []  # h : z + and - the apmplitud
    delta_w = [] # (w_exp - w_0)

    w_vs_h = []
    min_z = dist_interval[0]
    max_z = dist_interval[-1]
    # Integrate at every z value
    for hh in dist_interval:
        if hh-A < min_z: continue  # Check that we are not integrating out of the fiting boundaries
        if hh+A > max_z: continue

        # NOTE: With this xx value, got good results
        xx = np.linspace (hh-A, hh+A, int_steps)

        # NOTE: With this xx value all integrals are 0
        #  xx = hh + A*np.sin(phi)

        #f_sum  = -force_spl(hh + np.sin(phi))* np.sin(phi)*(2*np.pi/int_steps)
        # f_int_sum=sum(f_sum)
        f  = -force_spl(hh + A*np.sin(phi))* np.sin(phi)

        # KKKK Temp file to write out the forces and read in into fortran code
        fla = open('fla.dat_'+str(hh), 'w')
        fla.write('# xx    hh+Asin(phi)     f(h+Asin(phi))sin(phi)    phi  \n')
        for i, fi in enumerate(f):
            fla.write(str(xx[i]) + "\t" + str(hh + A*np.sin(phi[i])) + "\t" + str(fi) + "\t" + str(phi[i]) + '\n')

        f_int = simps(f,phi)
        fla.write("# Integral: "  + str(f_int))
        fla.close()
        function =1.0-1.0/(np.pi*k*A)*f_int

        w_vs_h += [(hh, function)]

    for i in range(len(w_vs_h)):
        h.append(w_vs_h[i][0]);  w.append(w_vs_h[i][1])

    delta_w = []
    for i, w_i in enumerate(w):
        temp = np.sqrt(w_i)*w_0
        delta_w.append(temp-w_0)

def read_autosave(ifile, grid_point):
    ''' Read autosaved file. Note ifile is an already read file as f.readlines()'''

    for n,line in enumerate(ifile):
        if 'grid_point' in line and not 'end' in line:
            if int(grid_point) == int(line.split()[1]):
                # range 4: Max number of options expected in this block
                for i in range(4):
                    option = ifile[n+i]
                    print option
                    if 'end_grid_point' in option:
                        break
                    elif 'smoothing' in option:
                        global spl_smoothing
                        spl_smoothing = float(option.split()[1])
                    elif 'z_min' in option:
                        global z_min
                        z_min = float(option.split()[1])
                    elif 'z_max' in option:
                        global z_max
                        z_max = float(option.split()[1])

def write_autosave(ofile, grid_point):
    global spl_smoothing
    global z_min
    global z_max

    print grid_point
    ofile.write("grid_point  " + str(grid_point) +'\n')
    ofile.write("  z_min  " + str(z_min) + '\n')
    ofile.write("  z_max  " + str(z_max) + '\n')
    ofile.write("  smoothing  " + str(spl_smoothing) + '\n')
    ofile.write("end_grid_point " + '\n')

def read_input(filename='inp.afm'):
    ''' Read input file. Default inp.afm'''

    if isinstance(filename, str):
        f = open(filename)
    else: # Assume it's a file-like object
        f = filename
    data  = f.readlines()

    for n,line in enumerate(data):
        if 'smoothing' in line:
            global spl_smoothing
            spl_smoothing = float(line.split()[1])
        elif 'z_min' in line:
            global z_min
            z_min = float(line.split()[1])
        elif 'z_max' in line:
            global z_max
            z_max = float(line.split()[1])
        elif 'omega_exp' in line:
            global omega_exp
            omega_exp = float(line.split()[1])
        elif 'datfile' in line:
            global outdatfile
            outdatfile = str(line.split()[1])
        elif 'tip_name' in line:
            global tip_name
            tip_name = line[9::]
        elif 'add_labels' in line:
            global add_labels
            add_labels = string_to_bool(line.split()[1])
        elif 'label_interval' in line:
            global label_interval
            label_interval = int(line.split()[1])
        elif 'prefix' in line:
            global prefix
            prefix = str(line.split()[1])
        elif 'grid_points' in line:
            global plot_grid_points
            plot_grid_points = []
            for i in line.split()[1::]:
                plot_grid_points.append(int(i))
        elif 'show_graphs' in line:
            global show_graphs
            show_graphs = string_to_bool(line.split()[1])
        elif 'show_afmimage' in line:
            global show_afmimage
            show_afmimage = string_to_bool(line.split()[1])
        elif 'autosave' in line:
            global autosave
            autosave = string_to_bool(line.split()[1])
        else:
            print ("Warning! Keyword "+ line + " not recognized")


def addCurveToGridPoint(gp, h, curve):
    global afmimage
    global grid

    #                       GP(x,y) coord,   h vs delta_w 
    #                                |       |     |
    afmimage[int(gp)].append([grid[int(gp)], h, curve])

def addForcesToGridPoint(gp, h, curve):
    global afmforces
    global grid

    #                       GP(x,y) coord,    h vs forces 
    #                                |        |     |
    afmforces[int(gp)].append([grid[int(gp)], h, -curve])


def get_Z_at_delta_w (omega,h, curve):
    f = np.array(curve)-omega
    ffit = UnivariateSpline(h, f, s=0)

    roots = ffit.roots()

    if roots.size > 1:
        #root = max(roots) # Take the Delta_w from the "atractive" part => Gives no contrast
        root = min(roots)  # Take the Delta_w from the "repulsive" part => Gives contrast 
    elif roots.size == 1:
        root = roots
    else:
        print ("Error! No Z found at this delta Omega. Please increase your Delta Omega and try again")
        return
    return root


# Get atom numbers to use in order to measure Z
atom_tip  = get_atom_index ("tip")
atom_surf = get_atom_index ("surface")
# atom_tip  = 204
# atom_surf = 108

if len(sys.argv) > 1:
    infile = sys.argv[1]
else:
    infile = 'inp.afm'

# Get the list of directories (grid points)
dir_list= [d for d in os.listdir('.') if re.match(r'grid_point_[0-9]+', d)]

# Set default parameters
# ----------------------
# General fitting and ploting parameters
z_min = 10.0
z_max = 100.0
spl_smoothing = 0.0005
outdatfile = 'Evsz.dat'
using_dat_file = False
add_labels = False
label_interval = 5

# Overwrite default parameters with an input file
read_input (infile)

# Save initial values for reuse at the begining of each grid point
z_min_init = z_min
z_max_init = z_max
spl_smooth_init = spl_smoothing

if_save=None
if autosave:  # Autosave contains Grid Point specific fitting values
    # Load previously existing autosave file
    try:
        if_save = open(prefix+"_autosave.afm").readlines()
    except:
        pass
    of_save = open(prefix+"_autosave.wip", 'w')

for directory in dir_list:

    energies = []
    distances = []
    vertical_forces = []

    # Start from initial parameters
    z_min = z_min_init
    z_max = z_max_init
    spl_smoothing = spl_smooth_init

    os.chdir(directory)
    grid_point = get_grid_point()

    # Print only the desired grid points
    if not int(grid_point) in plot_grid_points:
        os.chdir("..")
        continue
    else:
        print ("Plotting grid point: "+str(grid_point))
        print ("-----------------------")

    # Get the data in the dir and skip if nothing found
    if not outdatfile in os.listdir('.'):
        outcar_files = [f for f in os.listdir('.') if re.match(r'OUTCAR.[0-9]+', f)]
        if not outcar_files:
            # print ("No OUTCAR files found in grid point " + grid_point)
            print ("WARNING: No data available in grid point: " + grid_point)
            os.chdir("..")
            continue

    # Look for E vs Z dat file in the current dir
    if outdatfile in os.listdir('.'):
        print (outdatfile + " found in " + directory )
        data=np.genfromtxt(outdatfile, dtype=None, names=True)
        distances = data['Distance']
        energies  = data['Energy']
        outcar_files = data['Filename']
        using_dat_file=True
    # If not E vs Z data file found,then parse the OUTCARs
    else:
        # Sort the list in a human readable style
        outcar_files.sort(key=natural_keys)
        # Get the list of files
        poscar_files = [f for f in os.listdir('.') if re.match(r'POSCAR.[0-9]+', f)]
        # Sort the list in a human readable style
        poscar_files.sort(key=natural_keys)

        idx_en=0
        for outcar in outcar_files:
            idx_en=idx_en+1
            # Energy when sigma -> 0
            cmd = "grep 'energy  w' " + outcar +  " | tail -n 1 | awk '{print $7}'"
            energy=subprocess.check_output(cmd,shell=True)
            #energies.append([idx, float(energy.split()[0])])
            energies.append(float(energy.split()[0]))

        idx_d=0
        for poscar in poscar_files:
            idx_d=idx_d+1
            # coordinates: : [ idx, ( x y z) ]
            # constraints: : [ idx, ( T T T) ]
            coordinates, constraints = read_poscar(poscar)
            distance = calculate_distance(coordinates, atom_tip, atom_surf)
            #distances.append([idx, distance])
            vertical_force = sum_forces (outcar_files[idx_d-1], constraints)
            vertical_forces.append(vertical_force)
            distances.append(distance)

        # Checking consistency on number of files
        if idx_en != idx_d:
            print "Number of OUTCARS: ", idx_en
            print "Number of POSCARs/CONTCARs: ", idx_d
            error_critical ("Number of OUTCARs and CONTCARs/POSCARs does not match")

        os.rename('temp.xyz', 'Tip-surface_forces.xyz')

        # Save E vs z data in a file (to avoid parsing over and over again OUTCAR files)
        e_file=open(outdatfile,'w')
        e_file.write ("Distance      Energy    Filename\n")
        for i in range(len(distances)):
            e_file.write("%g\t%g\t%s\n" % (distances[i], energies[i], outcar_files[i]))

    # Override initial parameters for grid point specific parameters if requested
    if if_save: read_autosave(if_save, grid_point)
    title=tip_name + "Grid point: "+ str(grid_point)

    # Displacing energies up to 0
    max_ener = max(energies)
    for i in range(0,len(energies)):
        energies[i] = energies[i] - max_ener

    # Sorting the points from shorter to longer distances
    if distances[0] > distances[-1]:
        distances = distances[::-1]
        energies = energies[::-1]
        outcar_files=outcar_files[::-1]
        if not using_dat_file:
            poscar_files=poscar_files[::-1]

    dist_interval = []
    ener_interval = []
    outcar_interval = []
    force_spl = None
    e_fit_spl = None
    f_interp = None

    # Select points between height (z) interval
    set_interval()

    # Fit the E(z) to a spline and compute its first derivative
    # Note that force will be a spline
    force_method="spline"
    calc_force (dist_interval, ener_interval, force_method)


# # ---------------------------------------------------------
# # KKK TEMP: WRITING FORCES
#     print_forces()
# # ---------------------------------------------------------


    # Transform Force(z) to Delta Frequency(h)
    #-----------------------------------------
    # NOTE! Check frequency if it is 2pi/T or just 1/T 
    int_steps = 1001       # Integration steps
    #int_steps = 101         # Integration steps
    phi = np.linspace(0,2*np.pi,int_steps) # Angle

    w = []  # w ==> (w_exp/w_0)**2
    h = []  # h ==>  z +- A : range [min(z)+A , max(z)-A]
    delta_w = []
    integrate_forces()

    # GENERATING  GRAPHS
    #-------------------
    plt.style.use('fivethirtyeight')
    #plt.title(title, {'fontsize':'18'})

    fig = plt.figure()
    fig.suptitle(tip_name+'Grid Point: '+ grid_point,  fontweight='bold')
    # Energy graph
    ax1  = fig.add_subplot(221)
    # Force graph
    ax2  = fig.add_subplot(223)
    # Frequency graph
    ax3  = fig.add_subplot(224)

    out_format='svg'
    plot_name = prefix+"_allgraphs_"+str(grid_point)+"." + out_format

    # PLOTTING STUFFS
    #----------------
    plt.grid(True)
    # Plot all graphs
    replot_graphs()

    # Print the legend, and locate it bottom right
    ax1.legend(loc='best')
    #ax2.legend(loc='best')
    #ax3.legend(loc='best')

    #plt.savefig(plot_name, bbox_inches="tight")
    if save_graphs: plt.savefig(plot_name, dpi=1200,  bbox_inches="tight")

    # INTERACTIVE PLOTTING
    #---------------------

    def on_change_z_min (val):
        global z_min

        z_min = val
        set_interval()
        calc_force(dist_interval, ener_interval, 'spline')
        integrate_forces()
        replot_graphs()
        # print z_min
    def on_change_z_max (val):
        global z_max

        z_max = val
        set_interval()
        calc_force(dist_interval, ener_interval, 'spline')
        integrate_forces()
        replot_graphs()
        # print z_min

    def on_change_splsmooth (val):
        global spl_smoothing

        spl_smoothing = val
        set_interval()
        calc_force(dist_interval, ener_interval, 'spline')
        integrate_forces()
        replot_graphs()
        #print spl_smoothing

    #fig.subplots_adjust(bottom=0.2, left=0.1)

    #Placing the sliders
    slider_z_ax1 = plt.axes([0.50,0.90,0.4, 0.02])
    slider_z_ax2 = plt.axes([0.50,0.85,0.4, 0.02])
    slider_z_ax3 = plt.axes([0.50,0.80,0.4, 0.02])

    # Slider: Z_min
    slider_zmin = Slider (slider_z_ax1, "Z min", distances[0],distances[-4], valinit=z_min, color='#AAAAAA')
    slider_zmin.on_changed(on_change_z_min)

    # Slider: Z_max
    slider_zmax = Slider (slider_z_ax2, "Z max", distances[4],distances[-1]+1, valinit=z_max, color='#AAAAAA')
    slider_zmax.on_changed(on_change_z_max)

    # Slider: Spline smoothing
    #slider_smooth = Slider (slider_z_ax3, "Smooth", 0.000001,0.0001, valinit=spl_smoothing, color='#AAAAAA')
    slider_smooth = Slider (slider_z_ax3, "Smooth", 0.0,0.0001, valinit=spl_smoothing, color='#AAAAAA')
    slider_smooth.on_changed(on_change_splsmooth)

    fig1 = plt.gcf()
    if show_graphs: plt.show()
    if save_graphs: fig1.savefig(plot_name, bbox_inches="tight")

    plt.close()

    # Once fitting parameters have been adjusted, store the delta_w curve at each grid point 
    addCurveToGridPoint(grid_point, h, delta_w)
    addForcesToGridPoint(grid_point, dist_interval, force_spl(dist_interval))

    os.chdir("..")

    if autosave: write_autosave(of_save, grid_point)

# Close Aotosaving file
if autosave:
    of_save.close()
    os.rename(prefix+"_autosave.wip", prefix+"_autosave.afm")

def create_AFM_image():
    '''Create the NC-AFM image from the previous spectroscopies'''
    global imageAtW
    global afm_image
    # global x_interval
    # global y_interval
    global nxgp
    global nygp
    global xrep
    global yrep


    # Empty the previous image
    imageAtW = {1: [], # Object containing the  (x,y) grid points and the Z value at a given delta W
                2: [],
                3: [],
                4: [],
                5: [],
                6: [],
                7: [],
                8: [],
                9: [],
                10:[]}

    print ("--------------------\n")
    for gp in range (1,11):
        z_i = get_Z_at_delta_w (omega_exp, afmimage[gp][0][1], afmimage[gp][0][2])
        imageAtW[gp].append((grid[gp], z_i))

    # Fill the matrix (X,Y,Z)
    image_xyz = []
    for i in range(xrep):
        for j in range(yrep):
            for gp in [1, 3, 5, 7, 9]:
                image_xyz.append((imageAtW[gp][0][0][0] + x_interval*nxgp*i, imageAtW[gp][0][0][1] + y_interval*nygp*j, imageAtW[gp][0][1]))
        for j in range(yrep):
            for gp in [2, 4, 6, 8, 10]:
                image_xyz.append((imageAtW[gp][0][0][0] + x_interval*nxgp*i, imageAtW[gp][0][0][1] + y_interval*nygp*j, imageAtW[gp][0][1]))
    afm_image = np.array(image_xyz)



def plot_AFM_image ():

    # TODO: Generate finner grid by interpolating! wIdea: griddata
    x = afm_image[:,0]
    y = afm_image[:,1]
    z = afm_image[:,2]
    Z = z.reshape(nxgp*xrep, nygp*yrep)
    ZT = Z.T
    #plt.set_cmap('Greys')
    plt.set_cmap('afmhot')
    #ax2.imshow(ZT, aspect='auto', origin='lower', interpolation='gaussian', extent = ( x.min(), x.max(), y.min(), y.max()))
    ax2.grid(False)
    #ax2.imshow(ZT, aspect='equal', origin='lower', interpolation='gaussian')
    #ax2.imshow(ZT, aspect='equal', origin='lower', interpolation='lanczos')
    ax2.imshow(ZT, aspect='equal', origin='lower', interpolation='bicubic')
    #ax2.imshow(ZT, aspect='auto', origin='lower',  extent = ( x.min(), x.max(), y.min(), y.max()))
    # fig.colorbar(im)
    # if show_afmimage: plt.show()

def on_change_w (val):
    global omega_exp

    ax2.cla()
    omega_exp = val
    create_AFM_image()
    plot_AFM_image()
    plt.draw()

# Plot all the Delta omega curves togheter
fig = plt.figure()
# Plot of all the curves of delta w vs Z
ax1  = fig.add_subplot(221)
# Plot of full image at a omega_exp
ax2  = fig.add_subplot(122)

# Plot all force curves for each Grid Point
ax3  = fig.add_subplot(223)

color = [ 
    '#e34545',
    '#de27e5',
    '#5d4fe5',
    '#4fb9e5',
    '#10be3d',
    '#3f5d36',
    '#7d4804',
    '#ff783c',
    '#c70707',
    '#560606'
]

for gp in range (1,11):
    x = afmimage[gp][0][1]
    y = afmimage[gp][0][2]
    label = "GP " + str(gp)
    ax1.plot(x, y, '-', color=color[gp-1], label=label, lw=1.3)
    x = afmforces[gp][0][1]
    y = afmforces[gp][0][2]
    ax3.plot(x, y, '-', color=color[gp-1], label=label, lw=1.3)

# Slider: W
slider_z_ax1 = plt.axes([0.60,0.90,0.3, 0.02])
slider_w = Slider (slider_z_ax1, r"$\Delta \omega$", -100, -1, valinit=omega_exp, color='#AAAAAA')
slider_w.on_changed(on_change_w)

create_AFM_image ()
plot_AFM_image ()

ax3.legend(loc='lower right', ncol=2, numpoints=5)
# Once plt.show() is called, a new blank figure is created, so plt.savefig () will be a blank figure after plt.show()
# Get Current Figure in fig2 to save it whenever we want as plt.gcf()
fig2 = plt.gcf()
plt.show()
#plot_name = prefix+"_AFMImage"+".pdf"
plot_name = prefix+"_AFMImage"+".svg"
if save_afmimage: fig2.savefig(plot_name, bbox_inches="tight")

plt.close()
