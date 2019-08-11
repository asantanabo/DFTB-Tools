#!/usr/bin/env python

import numpy as np
import os, sys, math, re, subprocess
import optparse
import argparse
import itertools 
from itertools import izip_longest

# Obtaining parameters from the user

parser = optparse.OptionParser()

parser.add_option("-f", "--file", dest="mol.cif",
                  help="mol.cif unit cell ")

parser.add_option('--nx',
    action="store", dest="nx",
    help="supercell x direction", default="1")

parser.add_option('--ny',
    action="store", dest="ny",
    help="supercell y direction", default="1")

parser.add_option('--nz',
    action="store", dest="nz",
    help="supercell z direction", default="1")

options, args = parser.parse_args()

def get_fractional_to_cartesian_matrix(a, b, c, alpha, beta, gamma,
                                       angle_in_degrees=True):
    """
    Return the transformation matrix that converts fractional coordinates to
    cartesian coordinates.
    Parameters
    ----------
    a, b, c : float
        The lengths of the edges.
    alpha, gamma, beta : float
        The angles between the sides.
    angle_in_degrees : bool
        True if alpha, beta and gamma are expressed in degrees.
    Returns
    -------
    r : array_like
        The 3x3 rotation matrix. ``V_cart = np.dot(r, V_frac)``.
    """
    if angle_in_degrees:
        alpha = np.deg2rad(alpha)
        beta = np.deg2rad(beta)
        gamma = np.deg2rad(gamma)
    cosa = np.cos(alpha)
    sina = np.sin(alpha)
    cosb = np.cos(beta)
    sinb = np.sin(beta)
    cosg = np.cos(gamma)
    sing = np.sin(gamma)
    volume = 1.0 - cosa**2.0 - cosb**2.0 - cosg**2.0 + 2.0 * cosa * cosb * cosg
    volume = np.sqrt(volume)
    r = np.zeros((3, 3))
    r[0, 0] = a
    r[0, 1] = b * cosg
    r[0, 2] = c * cosb
    r[1, 1] = b * sing
    r[1, 2] = c * (cosa - cosb * cosg) / sing
    r[2, 2] = c * volume / sing
    return r.transpose()

def chunker(iterable, chunksize):
    return map(None,*[iter(iterable)]*chunksize)

# stored parameters from the user in the repetition of the unit cell

NX = int(options.nx)
NY = int(options.ny)
NZ = int(options.nz)

# Getting the atomic positions from the mol.cif file and saved in the
# file moltmp

ofile = open("moltmp", "w")
print_flag = False
with open('mol.cif', 'r') as f:
     for line in f:    
         if "_atom_site_fract_z" in line:
            print_flag = True
         if print_flag:
             for lines in f: 
                 ofile.write(lines)
ofile.close();


# obtaining the number of atoms and the coordinates from the moltmp file

element = []; x=[]; y=[]; z=[];

with open("moltmp") as g:
    for line in g:
        elem, xval, yval, zval = line.split()
        element.append(str(elem)); x.append(float(xval)); y.append(float(yval)); z.append(float(zval))

nratoms = len(element)
new_nratoms = int(nratoms*NX*NY*NZ)
#print new_nratoms

# preparing new atoms for the supercell

element_rescaled= []; x_rescaled = []; y_rescaled = []; z_rescaled = []

for i in xrange(0,len(x),1):
    element_rescaled.append(str(element[i])), x_rescaled.append(float(x[i]/NX)); y_rescaled.append((y[i]/NY)); z_rescaled.append((z[i]/NZ))

# creating the translations based on th enumber of supercells

x_super = []; y_super = []; z_super=[];

# making a supercell of (NX,0,0); (0,NY,0); (0,0,NZ)

for i in xrange(0,NX):
    x_super.append(i/(NX*1.0))
for j in xrange(0,NY):
    y_super.append(j/(NY*1.0))  
for k in xrange(0,NZ):
    z_super.append(k/(NZ*1.0))

#print x_super,y_super,z_super

re_unit_cell = zip(x_rescaled, y_rescaled, z_rescaled)
translations = [];

# Creating all possible combination of translations with the given NX,NY,NZ
for x in itertools.product(x_super,y_super,z_super):
      translations.append(x)

a = np.array(re_unit_cell)
b = np.array(translations)

loh = [];

# Creating the coordinates of the supercell

for i in xrange(0,len(re_unit_cell)):
    for j in xrange(0,len(translations)):
        loh.append(a[i]+b[j])

# Flattening the previous list of arrays 
list_super = np.concatenate(loh).ravel()

#all_coord = np.array(list_super).reshape(new_nratoms,3)

endlist = int(new_nratoms/(1.0*len(element_rescaled)))

#print element_rescaled
#print endlist
#print len(element_rescaled)
a = chunker(list_super, 3)
list_atoms = chunker(a,endlist)

final = [];

for i in xrange(0,len(element_rescaled)):
   for j in xrange(0,endlist):
       ele = [element_rescaled[i], list_atoms[i][j][0], list_atoms[i][j][1], list_atoms[i][j][2]]
       final.append(ele)

with open('supercelltmp','w') as f :       
  for i in xrange(0,len(final)):
        f.write("%s %12.5e %12.5e %12.5e\n" % (final[i][0],final[i][1],final[i][2],final[i][3]))
    
f.close();

# header of the xyz file to convert to gen
nrofatoms = open("nratoms","w")
nrofatoms.write('%g \n' %(new_nratoms))
nrofatoms.write('\n') 
nrofatoms.close();

with open('try.xyz', 'w+') as f:
    subprocess.Popen(["cat", "nratoms", "supercelltmp"], stdout=f)

subprocess.call('xyz2gen try.xyz',shell=True)


# Creating the supercell.gen file

with open('mol.cif', 'r') as search:
    for line in search:
        line = line.rstrip()
        if "_cell_length_a" in line:
            atmp=str(line)
        elif "_cell_length_b" in line:
            btmp=str(line)
        elif "_cell_length_c" in line:
            ctmp=str(line)
        elif "_cell_angle_alpha" in line:
            alpatmp=str(line)
        elif "_cell_angle_beta" in line:
            betatmp=str(line)
        elif "_cell_angle_gamma" in line:
            gammatmp=str(line)

# obtining the values of the lattice from the .cif file            
cell_a = re.findall(r"[-+]?\d*\.*\d+", atmp)
cell_b = re.findall(r"[-+]?\d*\.*\d+", btmp)
cell_c = re.findall(r"[-+]?\d*\.*\d+", ctmp)
cell_alpha = re.findall(r"[-+]?\d*\.*\d+", alpatmp)
cell_beta = re.findall(r"[-+]?\d*\.*\d+", betatmp)
cell_gamma= re.findall(r"[-+]?\d*\.*\d+", gammatmp)

# Converting them into numbers for using them in the fractcor code
a=NX*float(cell_a[0]); b=NY*float(cell_b[0]); c=NZ*float(cell_c[0]); alpha=float(cell_alpha[0]); beta=float(cell_beta[0]); gamma=float(cell_gamma[0])

lattice_cell = get_fractional_to_cartesian_matrix(a,b,c,alpha,beta,gamma)
translational = np.zeros(3)  # ading the three zeros of the translational
mat_temp = np.matrix(lattice_cell) # matrix of the cell
mat = np.vstack([translational,mat_temp]) # mergin zeros and lattice cell for mol.gen
    
with open('unit_cell.tmp','w') as f:
    for line in mat:
        np.savetxt(f, line, fmt='%12.5f')
      
with open('local.tmp', 'w') as f:
    subprocess.Popen(["cat", "try.gen", "unit_cell.tmp"], stdout=f)
    
# Changing the header of the local

newline = [len(final),"F"]
head = open('header.tmp', 'w')
head.write("%d %s \n" % (newline[0],newline[1]))
head.close();

subprocess.call(["sed -i '1d' local.tmp"], shell=True)

with open('supercell.gen', 'w+') as f:
    subprocess.Popen(["cat", "header.tmp", "local.tmp"], stdout=f)

subprocess.call(["rm -r *.tmp"], shell=True)
subprocess.call(["rm -r *tmp"], shell=True)
subprocess.call(["rm -r try.*"], shell=True)
subprocess.call(["rm -r nratoms"], shell=True)

