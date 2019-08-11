#/usr/bin/env python 

import sys, math, subprocess
import numpy as np
import scipy as sp
from itertools import * 
from itertools import izip_longest


###################################################
# This program reorganizes the hessian.out and    #
# print it in a square version. Then, change it   #
# to normal mode coordinates and get the          #
# eigenvectors.                                   #
###################################################

try:
    infile = sys.argv[1]; outfile = sys.argv[2]
except:
    print "Usage:",sys.argv[0], "hessian.out mol.xyz"; sys.exit(1)

ifile = open(infile, 'r')   # open the hessian to read it
ofile = open(outfile,  'r')   # open file with the eigenvectors

def chunker(iterable, chunksize):
    return map(None,*[iter(iterable)]*chunksize)

test = [];

for line in ifile:
     line = line.strip()
     if not line:
         continue
     test.append(line.split())

a=[];
for i in xrange(0,len(test)-1,2):
         a.append(list(chain(test[i],test[i+1])))

c=[];
for i in xrange(0,len(a)):
    c.append(np.array(a[i]))

d = np.array(c)
k = [];

for i in xrange(0,len(c)):
    k.append(np.vstack(c[i].reshape(len(c[i]),1)))

final = np.concatenate(k,axis=0).flatten()

linear = open('linear.tmp', 'w')
for i in xrange(0,len(final)):
      s = str(final[i])
      linear.write(s+"\n")

linear.close();


### creating hessianmatrix.dat
num_atoms=70

temp = np.loadtxt('linear.tmp')
hessian = temp.reshape(-1,num_atoms*3)

print hessian[0],hessian[1]

########################################################################################################################################################
# We will read the hessianmat.dat and transform it into mass weighted cordinates skiprows make the trick not to read the first line in hessianmat.dat  #
########################################################################################################################################################

# Dictionary with all masses of the elements to construct the mass weigthed coordinates

masses = {'H': 1.00079,'He': 4.0026,'Li': 6.941 ,'Be': 9.0122, 'B': 10.811, 
'C': 12.017, 'N': 14.0067, 'O': 15.9994, 'F': 18.9984, 'Ne': 20.1797, 
'Na': 22.9897,'Mg': 24.305, 'Al': 26.9815, 'Si': 28.0855, 'P': 30.9738,
'S': 32.065, 'Cl': 35.453, 'K': 39.0983, 'Ar': 39.948, 'Ca': 40.078, 
'Sc': 44.9559, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938,
'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 64.38,
'Ga': 69.723, 'Ge': 72.630, 'As': 74.922, 'Se':78.970, 'Br': 79.904, 
'Kr': 83.798, 'Rb':85.468, 'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224, 
'Nb': 92.906,'Mo': 95.95, 'Tc': 98.00, 'Ru': 101.07, 'Rh': 102.91, 
'Pd':106.42, 'Ag': 107.87, 'Cd':112.41, 'In': 114.82, 'Sn': 118.71, 
'Sb': 121.76, 'Te': 127.60, 'I': 126.90, 'Xe': 131.29, 'Cs': 132.91,
'Ba': 137.33, 'Hf': 178.49, 'Ta': 180.95, 'W': 183.84, 'Re': 186.21,
'Os': 190.23, 'Ir': 192.22,'Pt': 195.08,'Au': 196.97,'Hg': 200.59,
'Ti': 204.38,'Pb': 207.2,'Bi': 208.98,'Po': 209, 'At': 210, 'Rn': 222,
'Fr': 223, 'Ra': 226, 'Rf': 265, 'Db': 268, 'Sg': 271, 'Bf':207, 'Hs': 277,
'Mt':276,'Ds':281,'Rg':280,'Cn':285,'Nh':286,'Fl': 289, 'Mc': 289, 
'Lv':293,'Ts':294,'Og':294}

with open(outfile) as f:
    f.next()                             # Skipping the first line of the mol.xyz
    f.next()                             # Skipping the second line of the mol.xyz
    elements = [] ; atoms = []  
    for line in f:                       # Converting into Lists to get the masses of the atoms in atomic mass units        
     pair = line.split()                 # splitting the lists 
     if len(pair) > 1:
           atoms.append(str(pair[0]))
           elements.append(float(masses[str(pair[0])]))

elements = np.repeat(np.array(elements), 3)
el = 1 / np.sqrt(elements) 
red_masses_mat = np.outer(el, el)


#a = np.loadtxt('hessianmatrix.dat', dtype=np.float64, skiprows=1)
#a = a.reshape((3*atom,3*atom), order="C")    # We have to change the number in reshape for atom*atom. 

                                              
#res = np.multiply(red_masses_mat,a)
#final = res * 26424605.2076   # changing units to obtain frequencies in cm-1

#v, w = np.linalg.eigh(final)

#print sp.sqrt(v)
