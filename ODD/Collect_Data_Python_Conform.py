#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from matplotlib import pyplot as plt
import numpy as np
import collections
from collections import defaultdict
import sys
import os

import matplotlib

#makeLaTex available
# =============================================================================
params = {
     #'text.latex.preamble': ['\\usepackage{gensymb}'],
     'image.origin': 'lower',
     'image.interpolation': 'nearest',
     'image.cmap': 'gray',
     'axes.grid': False,
     'savefig.dpi': 150,  # to adjust notebook inline plot size
     'axes.labelsize': 8, # fontsize for x and y labels (was 10)
     'axes.titlesize': 8,
     'font.size': 8, # was 10
     'legend.fontsize': 6, # was 10
     'xtick.labelsize': 8,
     'ytick.labelsize': 8,
     'text.usetex': True,
     'figure.figsize': [3.39, 2.10],
     'font.family': 'serif',
 } 
matplotlib.rcParams.update(params)

##############################################################################
#Collate Data in two dictionaries
##############################################################################

Type = 'MutInf' #Type of data to be plotted (options: 'Purity', 'MutInf')

os.chdir(sys.path[0] + '/Data') #Change working directory to Data folder

#names of txt files in the folder containing the data
Hamiltonians = 'DD_Hamiltonians.txt'
States = 'DD_Init_States.txt'
Data_distance_unitary = 'DD_Data_dist_unit_opt.txt'
Data_max_eig = 'DD_Data_eig_opt.txt'
Operations_dist_unitary = 'DD_Operations_Optimized_unitary.txt'
Operations_max_eig = 'DD_Operations_Optimized_eig.txt'

                                   

#load data from the files

#Load Hamiltonians
HamDict = {}
Hamil = []
count = 0
f = open(Hamiltonians , "r")
for line in f:
    if (line[0] != '#') and (len(line.rstrip("\n")) >0):
        if line[0] == '[':
            line = line.rstrip("\n").split(' ')
            identifier = int(line[0][1:-1])
        else:
            count +=1
            line = line.rstrip("\n").split(' ')
            for item in line:
                Hamil.append(complex(item.replace('i', 'j')))
            if count ==4:
                HamDict[identifier] = np.array(Hamil).reshape(4,4)
                count = 0
                Hamil = []
                
#Load Environment States
StateDict = {}
State = []
count = 0
f = open(States, "r")
for line in f:
    if (line[0] != '#') and (len(line.rstrip("\n")) >0):
        if line[0] == '[':
            line = line.rstrip("\n").split(' ')
            identifier = int(line[0][1:-1])
        else:
            count +=1
            line = line.rstrip("\n").split(' ')
            for item in line:
                State.append(complex(item.replace('i', 'j')))
            if count ==2:
                StateDict[identifier] = np.array(State).reshape(2,2)
                count = 0
                State = []



#Load best operations from distance to unitaries algorithm
OpDictDistUnit = defaultdict(list)
Op = []
count = 0
id_count = 0
f = open(Operations_dist_unitary, "r")
for line in f:
    if (line[0] != '#') and (len(line.rstrip("\n")) >0):
        if line[0] == '[':
            line = line.rstrip("\n").split(' ')
            identifier = int(line[0][1:-1])
            
        else:
            count +=1
            line = line.rstrip("\n").split(' ')
            for item in line:
                Op.append(complex(item.replace('i', 'j')))
            if count >= 12:
                OpDictDistUnit[identifier].append(np.array(Op[0:16]).reshape(4,4))
                OpDictDistUnit[identifier].append(np.array(Op[16:32]).reshape(4,4))
                OpDictDistUnit[identifier].append(np.array(Op[32:]).reshape(4,4))
                count = 0
                Op = []
                
                
#Load best operations from max Eig algorithm
OpDictEig = defaultdict(list)
Op = []
count = 0
id_count = 0
f = open(Operations_max_eig, "r")
for line in f:
    if (line[0] != '#') and (len(line.rstrip("\n")) >0):
        if line[0] == '[':
            line = line.rstrip("\n").split(' ')
            identifier = int(line[0][1:-1])
            
        else:
            count +=1
            line = line.rstrip("\n").split(' ')
            for item in line:
                Op.append(complex(item.replace('i', 'j')))
            if count >= 12:
                OpDictEig[identifier].append(np.array(Op[0:16]).reshape(4,4))
                OpDictEig[identifier].append(np.array(Op[16:32]).reshape(4,4))
                OpDictEig[identifier].append(np.array(Op[32:]).reshape(4,4))
                count = 0
                Op = []


#load numerical data from the files

#from max eig algorithm
DataMaxEig = collections.defaultdict(dict)

f = open(Data_max_eig, "r")
for line in f:
    if line[0] != '#':
        line = line.rstrip("\n").split(' ')
        DataMaxEig[int(line[0])]['DelT'] = float(line[1])
        DataMaxEig[int(line[0])]['PurOpt'] = float(line[2])
        DataMaxEig[int(line[0])]['PurStand'] = float(line[3])
        DataMaxEig[int(line[0])]['MutInfOpt'] = float(line[4])
        DataMaxEig[int(line[0])]['MutInfStand'] = float(line[5])
        DataMaxEig[int(line[0])]['Hamiltonian'] = HamDict[int(line[0])]
        DataMaxEig[int(line[0])]['EnvState'] = StateDict[int(line[0])]
        DataMaxEig[int(line[0])]['Op1'] = OpDictEig[int(line[0])][0]
        DataMaxEig[int(line[0])]['Op2'] = OpDictEig[int(line[0])][1]
        DataMaxEig[int(line[0])]['Op3'] = OpDictEig[int(line[0])][2]

DataDistUnit = collections.defaultdict(dict)
            
f = open(Data_distance_unitary, "r")
for line in f:
    if line[0] != '#':
        line = line.rstrip("\n").split(' ')
        DataDistUnit[int(line[0])]['DelT'] = float(line[1])
        DataDistUnit[int(line[0])]['PurOpt'] = float(line[2])
        DataDistUnit[int(line[0])]['PurStand'] = float(line[3])
        DataDistUnit[int(line[0])]['MutInfOpt'] = float(line[4])
        DataDistUnit[int(line[0])]['MutInfStand'] = float(line[5])
        DataDistUnit[int(line[0])]['Hamiltonian'] = HamDict[int(line[0])]
        DataDistUnit[int(line[0])]['EnvState'] = StateDict[int(line[0])]
        DataDistUnit[int(line[0])]['Op1'] = OpDictDistUnit[int(line[0])][0]
        DataDistUnit[int(line[0])]['Op2'] = OpDictDistUnit[int(line[0])][1]
        DataDistUnit[int(line[0])]['Op3'] = OpDictDistUnit[int(line[0])][2]

##############################################################################
#Sort and bin data so that it can be plotted
##############################################################################

#For Data obtained from MaxEig algorithm
PlotDataMaxEigPur = defaultdict(list)
PlotDataStandPur = defaultdict(list)
PlotDataMaxEigMutInf = defaultdict(list)
PlotDataStandMutInf = defaultdict(list)
Del = DataMaxEig[1]['DelT']


for k in np.sort(np.array(list(DataMaxEig.keys()))):
    if abs(Del - DataMaxEig[k]['DelT']) < 0.00001:
        PlotDataMaxEigPur[Del].append(DataMaxEig[k]['PurOpt'])
        PlotDataStandPur[Del].append(DataMaxEig[k]['PurStand'])
        PlotDataMaxEigMutInf[Del].append(DataMaxEig[k]['MutInfOpt'])
        PlotDataStandMutInf[Del].append(DataMaxEig[k]['MutInfStand'])
    else:
        Del = DataMaxEig[k]['DelT']
        PlotDataMaxEigPur[Del].append(DataMaxEig[k]['PurOpt'])
        PlotDataStandPur[Del].append(DataMaxEig[k]['PurStand'])
        PlotDataMaxEigMutInf[Del].append(DataMaxEig[k]['MutInfOpt'])
        PlotDataStandMutInf[Del].append(DataMaxEig[k]['MutInfStand'])
        
#Create Dictionary with averaged values
for key in PlotDataMaxEigPur.keys():
    PlotDataMaxEigPur[key] = np.average(np.array(PlotDataMaxEigPur[key]))
    PlotDataStandPur[key] = np.average(np.array(PlotDataStandPur[key]))
    PlotDataMaxEigMutInf[key] = np.average(np.array(PlotDataMaxEigMutInf[key]))
    PlotDataStandMutInf[key] = np.average(np.array(PlotDataStandMutInf[key]))
    

#For Data obtained from Distance to unitaries algorithm
PlotDataDistUnitPur = defaultdict(list)
PlotDataDistUnitMutInf = defaultdict(list)
Del = DataDistUnit[1]['DelT']


for k in np.sort(np.array(list(DataDistUnit.keys()))):
    if abs(Del - DataDistUnit[k]['DelT']) < 0.00001:
        PlotDataDistUnitPur[Del].append(DataDistUnit[k]['PurOpt'])
        PlotDataDistUnitMutInf[Del].append(DataDistUnit[k]['MutInfOpt'])
    else:
        Del = DataDistUnit[k]['DelT']
        PlotDataDistUnitPur[Del].append(DataDistUnit[k]['PurOpt'])
        PlotDataDistUnitMutInf[Del].append(DataDistUnit[k]['MutInfOpt'])
        
#Create Dictionary with averaged values
for key in PlotDataDistUnitPur.keys():
    PlotDataDistUnitPur[key] = np.average(np.array(PlotDataDistUnitPur[key]))
    PlotDataDistUnitMutInf[key] = np.average(np.array(PlotDataDistUnitMutInf[key]))


#Store x- and y- Data in ascending order with respect to time

#Sort the keys of the two Dictionaries      
keys_sorted1 = np.sort(np.array(list(PlotDataMaxEigPur.keys())))
keys_sorted2 = np.sort(np.array(list(PlotDataDistUnitPur.keys())))

#Fill arrays with purity and MutInf in ascending time order
PurOptarr1 = []
PurStandarr1 =[]
MutInfOptarr1 = []
MutInfStandarr1 = []

for key1 in keys_sorted1:
    PurOptarr1.append(PlotDataMaxEigPur[key1])
    PurStandarr1.append(PlotDataStandPur[key1])
    MutInfOptarr1.append(PlotDataMaxEigMutInf[key1])
    MutInfStandarr1.append(PlotDataStandMutInf[key1])

#Convert to array for easier manipulation
PurOptarr1 = np.array(PurOptarr1)  
PurStandarr1 = np.array(PurStandarr1)
MutInfOptarr1 = np.array(MutInfOptarr1)
MutInfStandarr1 = np.array(MutInfStandarr1)
DelT = np.array(keys_sorted1)


PurOptarr2 = []
PurStandarr2 =[]
MutInfOptarr2 = []
MutInfStandarr2 = []

for key2 in keys_sorted2:
    PurOptarr2.append(PlotDataDistUnitPur[key2])
    MutInfOptarr2.append(PlotDataDistUnitMutInf[key2])

#Convert to array for easier manipulation
PurOptarr2 = np.array(PurOptarr2)  
PurStandarr2 = np.array(PurStandarr2)
MutInfOptarr2 = np.array(MutInfOptarr2)
MutInfStandarr2 = np.array(MutInfStandarr2)


##############################################################################
#Plot Data
##############################################################################


#Plot main figure
fig, ax1 = plt.subplots()

if Type == 'MutInf':
    ax1.plot(np.log10(DelT),  MutInfStandarr1, 'r', label = 'Standard DD Sequence')
    ax1.plot(np.log10(DelT), MutInfOptarr1,'g', label = 'Optimized largest Eigenvalue')
    ax1.plot(np.log10(DelT), MutInfOptarr2,'b', label = 'Optimized distance to unitaries')
    plt.title('Mutual Information')
else:
    if Type == 'Purity':
        ax1.plot(np.log10(DelT),  PurStandarr1, 'r', label = 'Standard DD Sequence')
        ax1.plot(np.log10(DelT), PurOptarr1,'g', label = 'Optimized largest Eigenvalue')
        ax1.plot(np.log10(DelT), PurOptarr2,'b', label = 'Optimized distance to unitaries')
        plt.title('Purity')

ax1.set_xlabel(r'$\log(\Delta \textsf{T})$', labelpad=3.2)
ax1.set_ylabel(r'$\textsf{Mutual Information}$',labelpad=10.2)
ax1.legend(loc=0)
plt.grid()



plt.show()













































