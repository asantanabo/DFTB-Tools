#!/usr/bin/env python3
import sys 
import glob
import shutil
import random
import numpy as np
from itertools import islice
import os

# Function created to move many files

def moveFilesByType(source, name_folder, extension):
   for basename in os.listdir(source):
     if basename.endswith('.' + extension):
        pathname = os.path.join(source, name_folder)
        shutil.copy(basename,pathname)

def removeFilesByType(source, extension):
    for basename in os.listdir(source):
        if basename.endswith('.' + extension):
           pathname = os.path.join(source, basename)
           os.remove(pathname)


         
# opening the MD file from DFTB
ifile=open(sys.argv[1],"r")

# Splitting into arrays of numatoms
natoms = int(ifile.readline())
all_coord= []

for line in ifile:
     if line.startswith('MD iter:'):
        all_coord.append(''.join(islice(ifile,natoms)))
     
#creating a list of random numbers to take configurations
listofnum = []
numconf = 1000
for x in range (0,numconf):
    listofnum.append(random.randint(0, len(all_coord)))

selconform = []
for i in range(0,len(listofnum)):
     selconform.append(all_coord[listofnum[i]])

for index, line in enumerate(selconform):
     with open('conformer_{}.txt'.format(index), 'w') as output:
          output.write(str(natoms))
          output.write('\n')
          output.write('\n')
          output.write(line)
    
# Creating directories for each conformer
current_dir = os.getcwd()
dir_new = current_dir + "/analysis_pes"

if not os.path.exists(dir_new):
    os.mkdir(dir_new)
    print("Directory ",dir_new,"Created")
else:    
    print("Directory",dir_new,"already exists")

# Moving all files to conformer analysis folder
moveFilesByType(current_dir, 'analysis_pes', 'txt')

# Removing all files from the current directory 
removeFilesByType(current_dir, 'txt')

# Creating many different folders to do SP calculations

basenames = [os.path.splitext(f)[0] for f in os.listdir(dir_new)]

for name in basenames:
    os.mkdir(os.path.join(dir_new,name))

# Counting the directories inside the dir_new folder and moving them to new folders
# and remnoving the .txt from the previous step

directories = [next(os.walk(dir_new))]
initial = []
final = []

for i in range(0,len(directories[0][1])):
   initial.append(os.path.join(directories[0][0], directories[0][1][i] + ".txt"))
   final.append(os.path.join(directories[0][0],directories[0][1][i]))
   
for i in range(0,len(initial)):
   shutil.move(initial[i],final[i])

# changing to analysis_pes folder looping all directories and moving to the mol.xyz and mol.gen 

os.chdir(final[0])
