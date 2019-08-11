#!/usr/bin/env python3
import sys 
import numpy as np


# Function created to read a file
ifile=open(sys.argv[1],"r")
timestep= [] ; angle_1 = []; angle_2 = [];
for line in ifile:                              
     pair = line.split()                 
     if len(pair) > 1:
           timestep.append(float(pair[0]))
           angle_1.append(float(pair[1]))
           angle_2.append(float(pair[2]))


# Applyitng statistics for list_1
angle_1 = np.array(angle_1)
average_1 = np.mean(angle_1)
stand_1 = np.std(angle_1, dtype=np.float64)
print (average_1)
print(stand_1)


# Applying staistics for list_2
angle_2 = np.array(angle_2)
average_2 = np.mean(angle_2)
stand_2 = np.std(angle_2, dtype=np.float64)
print(average_2)
print(stand_2)

#covarance

total = np.stack((angle_1,angle_2),axis=0)
print(np.cov(total))
