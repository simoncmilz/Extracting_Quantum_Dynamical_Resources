#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from matplotlib import pyplot as plt
import numpy as np
import collections
from collections import defaultdict
import sys
import os


import matplotlib


Type = 'MutInf' #Type of data to be plotted. Options: 'Purity', 'MutInf'


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

os.chdir(sys.path[0] + '/Data')	#set directory to the Data folder

##############################################################################
#Collate Data in two dictionaries
##############################################################################



#names of txt files in the folder containing the data
Hamiltonians = 'DD_Hamiltonians_Conc.txt'
States = 'DD_Init_States_Conc.txt'
Data_distance_unitary = 'DD_Data_dist_unit_opt_Conc.txt'
Data_max_eig = 'DD_Data_eig_opt_Conc.txt'
Data_Pur_opt = 'DD_Data_Purity_opt_Conc.txt'
Operations_dist_unitary = 'DD_Operations_Optimized_unitary_Conc.txt'
Operations_max_eig = 'DD_Operations_Optimized_eig_Conc.txt'
Operations_max_Pur = 'DD_Operations_Optimized_Purity_Conc.txt'
Lindbladians = 'DD_Lindbladians_Conc.txt'
                                   

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

#Load Lindbladians
LindDict = {}
Lind = []
count = 0
f = open(Lindbladians , "r")
for line in f:
    if (line[0] != '#') and (len(line.rstrip("\n")) >0):
        if line[0] == '[':
            line = line.rstrip("\n").split(' ')
            identifier = int(line[0][1:-1])
        else:
            count +=1
            line = line.rstrip("\n").split(' ')
            for item in line:
                Lind.append(complex(item.replace('i', 'j')))
            if count ==16:
                LindDict[identifier] = np.array(Lind).reshape(16,16)
                count = 0
                Lind = []


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
            if count >= 60:
                OpDictDistUnit[identifier].append(np.array(Op[0:16]).reshape(4,4))
                OpDictDistUnit[identifier].append(np.array(Op[16:32]).reshape(4,4))
                OpDictDistUnit[identifier].append(np.array(Op[32:48]).reshape(4,4))
                OpDictDistUnit[identifier].append(np.array(Op[48:64]).reshape(4,4))
                OpDictDistUnit[identifier].append(np.array(Op[64:80]).reshape(4,4))
                OpDictDistUnit[identifier].append(np.array(Op[80:96]).reshape(4,4))
                OpDictDistUnit[identifier].append(np.array(Op[96:112]).reshape(4,4))
                OpDictDistUnit[identifier].append(np.array(Op[112:128]).reshape(4,4))
                OpDictDistUnit[identifier].append(np.array(Op[128:144]).reshape(4,4))
                OpDictDistUnit[identifier].append(np.array(Op[144:160]).reshape(4,4))
                OpDictDistUnit[identifier].append(np.array(Op[160:176]).reshape(4,4))
                OpDictDistUnit[identifier].append(np.array(Op[176:192]).reshape(4,4))
                OpDictDistUnit[identifier].append(np.array(Op[192:208]).reshape(4,4))
                OpDictDistUnit[identifier].append(np.array(Op[208:224]).reshape(4,4))
                OpDictDistUnit[identifier].append(np.array(Op[224:]).reshape(4,4))
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
            if count >= 60:
                OpDictEig[identifier].append(np.array(Op[0:16]).reshape(4,4))
                OpDictEig[identifier].append(np.array(Op[16:32]).reshape(4,4))
                OpDictEig[identifier].append(np.array(Op[32:48]).reshape(4,4))
                OpDictEig[identifier].append(np.array(Op[48:64]).reshape(4,4))
                OpDictEig[identifier].append(np.array(Op[64:80]).reshape(4,4))
                OpDictEig[identifier].append(np.array(Op[80:96]).reshape(4,4))
                OpDictEig[identifier].append(np.array(Op[96:112]).reshape(4,4))
                OpDictEig[identifier].append(np.array(Op[112:128]).reshape(4,4))
                OpDictEig[identifier].append(np.array(Op[128:144]).reshape(4,4))
                OpDictEig[identifier].append(np.array(Op[144:160]).reshape(4,4))
                OpDictEig[identifier].append(np.array(Op[160:176]).reshape(4,4))
                OpDictEig[identifier].append(np.array(Op[176:192]).reshape(4,4))
                OpDictEig[identifier].append(np.array(Op[192:208]).reshape(4,4))
                OpDictEig[identifier].append(np.array(Op[208:224]).reshape(4,4))
                OpDictEig[identifier].append(np.array(Op[224:]).reshape(4,4))
                count = 0
                Op = []
                
#Load best operations from max Purity algorithm
OpDictPur = defaultdict(list)
Op = []
count = 0
id_count = 0
f = open(Operations_max_Pur, "r")
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
            if count >= 60:
                OpDictPur[identifier].append(np.array(Op[0:16]).reshape(4,4))
                OpDictPur[identifier].append(np.array(Op[16:32]).reshape(4,4))
                OpDictPur[identifier].append(np.array(Op[32:48]).reshape(4,4))
                OpDictPur[identifier].append(np.array(Op[48:64]).reshape(4,4))
                OpDictPur[identifier].append(np.array(Op[64:80]).reshape(4,4))
                OpDictPur[identifier].append(np.array(Op[80:96]).reshape(4,4))
                OpDictPur[identifier].append(np.array(Op[96:112]).reshape(4,4))
                OpDictPur[identifier].append(np.array(Op[112:128]).reshape(4,4))
                OpDictPur[identifier].append(np.array(Op[128:144]).reshape(4,4))
                OpDictPur[identifier].append(np.array(Op[144:160]).reshape(4,4))
                OpDictPur[identifier].append(np.array(Op[160:176]).reshape(4,4))
                OpDictPur[identifier].append(np.array(Op[176:192]).reshape(4,4))
                OpDictPur[identifier].append(np.array(Op[192:208]).reshape(4,4))
                OpDictPur[identifier].append(np.array(Op[208:224]).reshape(4,4))
                OpDictPur[identifier].append(np.array(Op[224:]).reshape(4,4))
                count = 0
                Op = []


#load numerical data from the files

#from max eig algorithm
DataMaxEig = collections.defaultdict(dict)

f = open(Data_max_eig, "r")
for line in f:
    if line[0] != '#' and (len(line.rstrip("\n")) >0):
        line = line.rstrip("\n").split(' ')
        DataMaxEig[int(line[0])]['DelT'] = float(line[1])
        DataMaxEig[int(line[0])]['PurOpt'] = float(line[2])
        DataMaxEig[int(line[0])]['PurStand'] = float(line[3])
        DataMaxEig[int(line[0])]['PurBase'] = float(line[4])
        DataMaxEig[int(line[0])]['MutInfOpt'] = float(line[5])
        DataMaxEig[int(line[0])]['MutInfStand'] = float(line[6])
        DataMaxEig[int(line[0])]['MutInfBase'] = float(line[7])
        DataMaxEig[int(line[0])]['Hamiltonian'] = HamDict[int(line[0])]
        DataMaxEig[int(line[0])]['EnvState'] = StateDict[int(line[0])]
        DataMaxEig[int(line[0])]['Lindbladian'] = LindDict[int(line[0])]        
        DataMaxEig[int(line[0])]['Op1'] = OpDictEig[int(line[0])][0]
        DataMaxEig[int(line[0])]['Op2'] = OpDictEig[int(line[0])][1]
        DataMaxEig[int(line[0])]['Op3'] = OpDictEig[int(line[0])][2]
        DataMaxEig[int(line[0])]['Op4'] = OpDictEig[int(line[0])][3]
        DataMaxEig[int(line[0])]['Op5'] = OpDictEig[int(line[0])][4]
        DataMaxEig[int(line[0])]['Op6'] = OpDictEig[int(line[0])][5]
        DataMaxEig[int(line[0])]['Op7'] = OpDictEig[int(line[0])][6]
        DataMaxEig[int(line[0])]['Op8'] = OpDictEig[int(line[0])][7]
        DataMaxEig[int(line[0])]['Op9'] = OpDictEig[int(line[0])][8]
        DataMaxEig[int(line[0])]['Op10'] = OpDictEig[int(line[0])][9]
        DataMaxEig[int(line[0])]['Op11'] = OpDictEig[int(line[0])][10]
        DataMaxEig[int(line[0])]['Op12'] = OpDictEig[int(line[0])][11]
        DataMaxEig[int(line[0])]['Op13'] = OpDictEig[int(line[0])][12]
        DataMaxEig[int(line[0])]['Op14'] = OpDictEig[int(line[0])][13]
        DataMaxEig[int(line[0])]['Op15'] = OpDictEig[int(line[0])][14]

#from optimization of distance to the unitaries
DataDistUnit = collections.defaultdict(dict)
            
f = open(Data_distance_unitary, "r")
for line in f:
    if line[0] != '#' and (len(line.rstrip("\n")) >0):
        line = line.rstrip("\n").split(' ')
        DataDistUnit[int(line[0])]['DelT'] = float(line[1])
        DataDistUnit[int(line[0])]['PurOpt'] = float(line[2])
        DataDistUnit[int(line[0])]['PurStand'] = float(line[3])
        DataDistUnit[int(line[0])]['MutInfOpt'] = float(line[5])
        DataDistUnit[int(line[0])]['MutInfStand'] = float(line[6])
        DataDistUnit[int(line[0])]['Hamiltonian'] = HamDict[int(line[0])]
        DataDistUnit[int(line[0])]['EnvState'] = StateDict[int(line[0])]
        DataDistUnit[int(line[0])]['Lindbladian'] = LindDict[int(line[0])] 
        DataDistUnit[int(line[0])]['Op1'] = OpDictDistUnit[int(line[0])][0]
        DataDistUnit[int(line[0])]['Op2'] = OpDictDistUnit[int(line[0])][1]
        DataDistUnit[int(line[0])]['Op3'] = OpDictDistUnit[int(line[0])][2]
        DataDistUnit[int(line[0])]['Op4'] = OpDictDistUnit[int(line[0])][3]
        DataDistUnit[int(line[0])]['Op5'] = OpDictDistUnit[int(line[0])][4]
        DataDistUnit[int(line[0])]['Op6'] = OpDictDistUnit[int(line[0])][5]
        DataDistUnit[int(line[0])]['Op7'] = OpDictDistUnit[int(line[0])][6]
        DataDistUnit[int(line[0])]['Op8'] = OpDictDistUnit[int(line[0])][7]
        DataDistUnit[int(line[0])]['Op9'] = OpDictDistUnit[int(line[0])][8]
        DataDistUnit[int(line[0])]['Op10'] = OpDictDistUnit[int(line[0])][9]
        DataDistUnit[int(line[0])]['Op11'] = OpDictDistUnit[int(line[0])][10]
        DataDistUnit[int(line[0])]['Op12'] = OpDictDistUnit[int(line[0])][11]
        DataDistUnit[int(line[0])]['Op13'] = OpDictDistUnit[int(line[0])][12]
        DataDistUnit[int(line[0])]['Op14'] = OpDictDistUnit[int(line[0])][13]
        DataDistUnit[int(line[0])]['Op15'] = OpDictDistUnit[int(line[0])][14]
        
    
#from optimization of distance to the putiry of the resulting channel
DataPurOpt = collections.defaultdict(dict)
            
f = open(Data_Pur_opt, "r")
for line in f:
    if line[0] != '#' and (len(line.rstrip("\n")) >0):
        line = line.rstrip("\n").split(' ')
        DataPurOpt[int(line[0])]['DelT'] = float(line[1])
        DataPurOpt[int(line[0])]['PurOpt'] = float(line[2])
        DataPurOpt[int(line[0])]['PurStand'] = float(line[3])
        DataPurOpt[int(line[0])]['MutInfOpt'] = float(line[5])
        DataPurOpt[int(line[0])]['MutInfStand'] = float(line[6])
        DataPurOpt[int(line[0])]['Hamiltonian'] = HamDict[int(line[0])]
        DataPurOpt[int(line[0])]['EnvState'] = StateDict[int(line[0])]
        DataPurOpt[int(line[0])]['Lindbladian'] = LindDict[int(line[0])] 
        DataPurOpt[int(line[0])]['Op1'] = OpDictPur[int(line[0])][0]
        DataPurOpt[int(line[0])]['Op2'] = OpDictPur[int(line[0])][1]
        DataPurOpt[int(line[0])]['Op3'] = OpDictPur[int(line[0])][2]
        DataPurOpt[int(line[0])]['Op4'] = OpDictPur[int(line[0])][3]
        DataPurOpt[int(line[0])]['Op5'] = OpDictPur[int(line[0])][4]
        DataPurOpt[int(line[0])]['Op6'] = OpDictPur[int(line[0])][5]
        DataPurOpt[int(line[0])]['Op7'] = OpDictPur[int(line[0])][6]
        DataPurOpt[int(line[0])]['Op8'] = OpDictPur[int(line[0])][7]
        DataPurOpt[int(line[0])]['Op9'] = OpDictPur[int(line[0])][8]
        DataPurOpt[int(line[0])]['Op10'] = OpDictPur[int(line[0])][9]
        DataPurOpt[int(line[0])]['Op11'] = OpDictPur[int(line[0])][10]
        DataPurOpt[int(line[0])]['Op12'] = OpDictPur[int(line[0])][11]
        DataPurOpt[int(line[0])]['Op13'] = OpDictPur[int(line[0])][12]
        DataPurOpt[int(line[0])]['Op14'] = OpDictPur[int(line[0])][13]
        DataPurOpt[int(line[0])]['Op15'] = OpDictPur[int(line[0])][14]


##############################################################################
#Sort and bin data so that it can be plotted
##############################################################################

#For Data obtained from MaxEig algorithm
PlotDataMaxEigPur = defaultdict(list)
PlotDataStandPur = defaultdict(list)
PlotDataBasePur = defaultdict(list)
PlotDataMaxEigMutInf = defaultdict(list)
PlotDataStandMutInf = defaultdict(list)
PlotDataBaseMutInf = defaultdict(list)

#Purities of the best operations
PurityOpEig = defaultdict(list)

Del = DataMaxEig[1]['DelT']

print(len(np.array(list(DataMaxEig.keys()))))
for k in np.sort(np.array(list(DataMaxEig.keys()))):
    if abs(Del - DataMaxEig[k]['DelT']) < 0.00001:
        PlotDataMaxEigPur[Del].append(DataMaxEig[k]['PurOpt'])
        PlotDataStandPur[Del].append(DataMaxEig[k]['PurStand'])
        PlotDataBasePur[Del].append(DataMaxEig[k]['PurBase'])
        PlotDataMaxEigMutInf[Del].append(DataMaxEig[k]['MutInfOpt'])
        PlotDataStandMutInf[Del].append(DataMaxEig[k]['MutInfStand'])
        PlotDataBaseMutInf[Del].append(DataMaxEig[k]['MutInfBase'])
        for Op in OpDictEig[k]:
            PurityOpEig[Del].append(np.trace(np.dot(Op,Op)/(np.trace(Op)**2)))
            
    else:
        Del = DataMaxEig[k]['DelT']
        PlotDataMaxEigPur[Del].append(DataMaxEig[k]['PurOpt'])
        PlotDataStandPur[Del].append(DataMaxEig[k]['PurStand'])
        PlotDataBasePur[Del].append(DataMaxEig[k]['PurBase'])
        PlotDataMaxEigMutInf[Del].append(DataMaxEig[k]['MutInfOpt'])
        PlotDataStandMutInf[Del].append(DataMaxEig[k]['MutInfStand'])
        PlotDataBaseMutInf[Del].append(DataMaxEig[k]['MutInfBase'])
        for Op in OpDictEig[k]:
            PurityOpEig[Del].append(np.trace(np.dot(Op,Op)/(np.trace(Op)**2)))
        
#Create Dictionary with averaged values
for key in PlotDataMaxEigPur.keys():
    PlotDataMaxEigPur[key] = np.average(np.array(PlotDataMaxEigPur[key]))
    PlotDataStandPur[key] = np.average(np.array(PlotDataStandPur[key]))
    PlotDataBasePur[key] = np.average(np.array(PlotDataBasePur[key]))
    
    PurityOpEig[key] = np.average(np.array(PurityOpEig[key]))
    
    PlotDataMaxEigMutInf[key] = np.average(np.array(PlotDataMaxEigMutInf[key]))
    PlotDataStandMutInf[key] = np.average(np.array(PlotDataStandMutInf[key]))
    PlotDataBaseMutInf[key] = np.average(np.array(PlotDataBaseMutInf[key]))
    

#For Data obtained from Distance to unitaries algorithm
PlotDataDistUnitPur = defaultdict(list)
PlotDataDistUnitMutInf = defaultdict(list)
PurityOpDistUnit = defaultdict(list)
Del = DataDistUnit[1]['DelT']


for k in np.sort(np.array(list(DataDistUnit.keys()))):
    if abs(Del - DataDistUnit[k]['DelT']) < 0.00001:
        PlotDataDistUnitPur[Del].append(DataDistUnit[k]['PurOpt'])
        PlotDataDistUnitMutInf[Del].append(DataDistUnit[k]['MutInfOpt'])
        for Op in OpDictDistUnit[k]:
            PurityOpDistUnit[Del].append(np.trace(np.dot(Op,Op)/(np.trace(Op)**2)))
    else:
        Del = DataDistUnit[k]['DelT']
        PlotDataDistUnitPur[Del].append(DataDistUnit[k]['PurOpt'])
        PlotDataDistUnitMutInf[Del].append(DataDistUnit[k]['MutInfOpt'])
        for Op in OpDictDistUnit[k]:
            PurityOpDistUnit[Del].append(np.trace(np.dot(Op,Op)/(np.trace(Op)**2)))
        
#Create Dictionary with averaged values
for key in PlotDataDistUnitPur.keys():
    PlotDataDistUnitPur[key] = np.average(np.array(PlotDataDistUnitPur[key]))
    PlotDataDistUnitMutInf[key] = np.average(np.array(PlotDataDistUnitMutInf[key]))
    
    PurityOpDistUnit[key] = np.average(np.array(PurityOpDistUnit[key]))
    

#For Data obtained from Purity Optimization
PlotDataPurOptPur = defaultdict(list)
PlotDataPurOptMutInf = defaultdict(list)
PurityOpsPur = defaultdict(list)
Del = DataPurOpt[1]['DelT']


for k in np.sort(np.array(list(DataPurOpt.keys()))):
    if abs(Del - DataDistUnit[k]['DelT']) < 0.00001:
        PlotDataPurOptPur[Del].append(DataPurOpt[k]['PurOpt'])
        PlotDataPurOptMutInf[Del].append(DataPurOpt[k]['MutInfOpt'])
        for Op in OpDictPur[k]:
            PurityOpsPur[Del].append(np.trace(np.dot(Op,Op)/(np.trace(Op)**2)))
    else:
        Del = DataDistUnit[k]['DelT']
        PlotDataPurOptPur[Del].append(DataPurOpt[k]['PurOpt'])
        PlotDataPurOptMutInf[Del].append(DataPurOpt[k]['MutInfOpt'])
        for Op in OpDictPur[k]:
            PurityOpsPur[Del].append(np.trace(np.dot(Op,Op)/(np.trace(Op)**2)))
        
#Create Dictionary with averaged values
for key in PlotDataDistUnitPur.keys():
    PlotDataPurOptPur[key] = np.average(np.array(PlotDataPurOptPur[key]))
    PlotDataPurOptMutInf[key] = np.average(np.array(PlotDataPurOptMutInf[key]))
    
    PurityOpsPur[key] = np.average(np.array(PurityOpsPur[key]))


#Store x- and y- Data in ascending order with respect to time

#Sort the keys of the two Dictionaries      
keys_sorted1 = np.sort(np.array(list(PlotDataMaxEigPur.keys())))
keys_sorted2 = np.sort(np.array(list(PlotDataDistUnitPur.keys())))
keys_sorted3 = np.sort(np.array(list(PlotDataPurOptPur.keys())))

#Fill arrays with purity and MutInf in ascending time order
PurOptarr1 = []
PurStandarr1 =[]
PurBasearr1 = []

MutInfOptarr1 = []
MutInfStandarr1 = []
MutInfBasearr1 = []

PurOps1 = []
PurOps2 = []
PurOps3 = []

for key1 in keys_sorted1:
    PurOptarr1.append(PlotDataMaxEigPur[key1])
    PurStandarr1.append(PlotDataStandPur[key1])
    PurBasearr1.append(PlotDataBasePur[key1])
    
    MutInfOptarr1.append(PlotDataMaxEigMutInf[key1])
    MutInfStandarr1.append(PlotDataStandMutInf[key1])
    MutInfBasearr1.append(PlotDataBaseMutInf[key1])
    
    PurOps1.append(PurityOpEig[key1])

#Convert to array for easier manipulation
PurOptarr1 = np.array(PurOptarr1)  
PurStandarr1 = np.array(PurStandarr1)
PurBasearr1 = np.array(PurBasearr1)

MutInfOptarr1 = np.array(MutInfOptarr1)
MutInfStandarr1 = np.array(MutInfStandarr1)
MutInfBasearr1 = np.array(MutInfBasearr1)

DelT = np.array(keys_sorted1)


PurOptarr2 = []
PurStandarr2 =[]
MutInfOptarr2 = []
MutInfStandarr2 = []

for key2 in keys_sorted2:
    PurOptarr2.append(PlotDataDistUnitPur[key2])
    MutInfOptarr2.append(PlotDataDistUnitMutInf[key2])
    PurOps2.append(PurityOpDistUnit[key2])

#Convert to array for easier manipulation
PurOptarr2 = np.array(PurOptarr2)  
PurStandarr2 = np.array(PurStandarr2)
MutInfOptarr2 = np.array(MutInfOptarr2)
MutInfStandarr2 = np.array(MutInfStandarr2)

PurOptarr3 = []
PurStandarr3 =[]
MutInfOptarr3 = []
MutInfStandarr3 = []

for key3 in keys_sorted3:
    PurOptarr3.append(PlotDataPurOptPur[key3])
    MutInfOptarr3.append(PlotDataPurOptMutInf[key3])
    PurOps3.append(PurityOpsPur[key3])

#Convert to array for easier manipulation
PurOptarr3 = np.array(PurOptarr3)  
PurStandarr3 = np.array(PurStandarr3)
MutInfOptarr3 = np.array(MutInfOptarr3)
MutInfStandarr3 = np.array(MutInfStandarr3)


##############################################################################
#Plot Data
##############################################################################

#Plot main figure
fig, ax1 = plt.subplots()

if Type == 'MutInf':
    ax1.plot(np.log10(DelT),  MutInfStandarr1, 'r', label = 'Standard DD Sequence')
    ax1.plot(np.log10(DelT), MutInfOptarr1,'g', label = 'Optimized largest Eigenvalue')
    ax1.plot(np.log10(DelT), MutInfOptarr2,'b', label = 'Optimized distance to unitaries')
    ax1.plot(np.log10(DelT), MutInfOptarr3,'black', label = 'Optimized Purity')
    ax1.plot(np.log10(DelT), MutInfBasearr1,'y', label = 'Do-nothing protocol')
    #ax1.plot(np.log10(DelT), PurOps1, label = 'purity of operations Eig')
    #ax1.plot(np.log10(DelT), PurOps2, label = 'purity of operations DistUnit')
    #ax1.plot(np.log10(DelT), PurOps3, label = 'purity of operations Purity')
    ax1.set_ylabel(r'$\textsf{Mutual Information}$',labelpad=10.2)
    plt.title('Mutual Information')
else:
    if Type == 'Purity':
        ax1.plot(np.log10(DelT),  PurStandarr1, 'r', label = 'Standard DD Sequence')
        ax1.plot(np.log10(DelT), PurOptarr1,'g', label = 'Optimized largest Eigenvalue')
        ax1.plot(np.log10(DelT), PurOptarr2,'b', label = 'Optimized distance to unitaries')
        ax1.plot(np.log10(DelT), PurOptarr3,'black', label = 'Optimized Purity')
        ax1.plot(np.log10(DelT), PurBasearr1,'y', label = 'do-nothing protocol')
        #ax1.set_ylabel(r'$\textsf{Purity}$',labelpad=10.2)
        plt.title('Purity')
ax1.set_xlabel(r'$\log(\Delta \textsf{T})$', labelpad=3.2)
ax1.legend(loc=0)
plt.grid()


plt.show()













































