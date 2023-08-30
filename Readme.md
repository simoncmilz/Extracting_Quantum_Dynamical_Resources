## Code to accompany: *[Extracting Quantum Dynamical Resources: Consumption of Non-Markovianity for Noise Reduction](https://arxiv.org/abs/2110.02613)* 
#### Graeme D. Berk, Simon Milz, Felix A. Pollock, Kavan Modi

This is a repository for the code used in the article "*Extracting Quantum Dynamical Resources: Consumption of Non-Markovianity for Noise Reduction*, Graeme D. Berk, Simon Milz, Felix A. Pollock, Kavan Modi, [arXiv:2110.02613 [quant-ph]]([https://arxiv.org/abs/2110.03233](https://arxiv.org/abs/2110.02613))".

All code is written in matlab (for the creation of the data) and python (for visualisation) and requires:
- [cvx](http://cvxr.com/cvx/) - a free Matlab Software for Disciplined Convex Programming 
- [MOSEK](https://www.mosek.com) - a software package for solving mathematical optimization problems (under the free personal academic license)

- [numpy](https://numpy.org/) - the fundamental package for scientific computing with Python
- [matplotlib](https://matplotlib.org/) - a comprehensive library for creating static, animated, and interactive visualizations in Python.

**TL;DR**: To plot the data for **ODD**, run ODD/Collect_Data_Python_Conform.py. In said file, the type of data that is to be plotted (purity or mutual information) can be chosen by changing the variable T=**Type**.

To plot the data for **MODD**, run MODD/Collect_Data_Python_Conform_Concatinaed.py. In said file, the type of data that is to be plotted (purity or mutual information) can be chosen by changing the variable T=**Type**.


###### This repository consists of the following:

#### CODE

##### Main programs (see [arXiv:2110.02613](https://arxiv.org/abs/2110.02613) for details)

###### ODD (single round optimization)

- [ODD/DD_Sequence_unitary.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/DD_Sequence_unitary.m):
given a **Comb**, the number of slots **N**, the dimension of the considered system **dim**, and an **exitTolerance** that fixes when the see-saw algorithm stops, this program computes (locally) optimal decoupling operations by **minimizing the distance of the resulting channel to the set of unitaries**. This is achieved via a see-saw SDP that optimizes the operation in each slot individually, and stops once the relative change of the resulting distance is lower than exitTolerance.

    - As an output, the program yields

      **PurityOpt, PurityStand**: The purity of the (Choi matrix of the) resulting channel, both for the starting CPTP maps and the optimized ones

      **MutInfOpt, MutInfStand**: The input-output mutual information of the (Choi matrix of the) resulting channel, both for the starting CPTP maps and the optimized 
ones
      **Purities**: Array containing the purities of the 

      **OpUnit1, OpUnit2, OpUnit3**: The purities of the employed CPTP maps

    Currently, the program is only fully functional for **N=3**, to change this, only the output of OpUnit1, OpUnit2, OpUnit3 would have to be expanded

- [ODD/DD_Sequence_max_eig.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/DD_Sequence_max_eig.m):
given a **Comb**, the number of slots **N**, the dimension of the considered system **dim**, and an **exitTolerance** that fixes when the see-saw algorithm stops, this program computes (locally) optimal decoupling operations by **maximizing the eigenvalue of the resulting channe**. This is achieved via a see-saw SDP that iteratively optimizes the operation in each slot and the projector on the largest eigenvalue of the resulting channel, and stops once the relative change of the eigenvalue is lower than exitTolerance.

    - As an output, the program yields

      **PurityOpt, PurityStand**: The purity of the (Choi matrix of the) resulting channel, both for the starting CPTP maps and the optimized ones

      **MutInfOpt, MutInfStand**: The input-output mutual information of the (Choi matrix of the) resulting channel, both for the starting CPTP maps and the optimized 
ones
      **Purities**: Array containing the purities of the 

      **OpEig1, OpEig2, OpEig3**: The purities of the employed CPTP maps

  Currently, the program is only fully functional for **N=3**, to change this, only the output of OpUnit1, OpUnit2, OpUnit3 would have to be expanded

- [ODD/DD_Sequence_Print.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/DD_Sequence_Print.m):
given a number of slots **N**, system dimension **dim**, number of times **n**, number of samples **NPerTime**, and exitTolerance, this program computes ideal uncoupling sequences for NPerTime x n Combs and writes all data into the respective .txt files

  Currently, the program is only fully functional for **N=3**, to change this, the size of the respective dictionaries has to be made variable, and the programs [ODD/DD_Sequence_unitary.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/DD_Sequence_unitary.m) and [DD_Sequence_max_eig.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/DD_Sequence_max_eig.m) would have to be adjusted accordingly

- [ODD/Collect_Data_Python_Conform.py](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/Collect_Data_Python_Conform.py):
Python script that plots the data contained in the ODD/Data folder. Type of data to be plotted (options: 'Purity', 'MutInf') can be chosen by changing the variable **Type** in the script. 

###### MODD (multi-round optimization)

- [MODD/Concatenated_DD.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Concatenated_DD.m):
Given a number of rounds **m**, a number of slots **N** per round, a Hamiltonian **Hamil**, an initial environment state **Init_Env**, a step size **delT**, an **exitTolerance**, and the type of **optimization**, this function computes optimal MODD decoupling operations for the given parameters and returns the mutual information and purity of the resulting channel, as well as the optimized operations and the used Linbladian. As optimization methods, it accepts **eigenvalue** (for maximization of the eigenvalue of the resulting channel), **dist_unit** (for minimization of the resulting channel to the set of unitaries) as well as **purity** (for maximization of the purity of the resulting channel).

    - As an output, the program yields

      **MutInfOpt**: The mutual information of the (Choi matrix of the) resulting channel after decoupling

      **PurityOpt**: Purity of the (Choi matrix of the) resulting channel after decoupling
ones
      **Op**: Array containing the optimal decoupling operations 

      **Lindblad**: Lindbladian used for the dynamics

  Currently, the program is only fully functional for **N=3** and system dimension **d=2**

- [MODD/DD_Sequence_Print_Concatenated.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/DD_Sequence_Print_Concatenated.m):
  Given a number of rounds **m**, a number of slots **N**, a system dimension **dim**, a number of times **n**, a number **NPerTime** of iterations per time, and an **exitTolerance**, this function computes optimized DD sequences for three different optimization techniques (maximization of the eigenvalue of the resulting channel, minimization of the resulting channel to the set of unitaries, and maximization of the purity of the resulting channel) and writes all produced data to the respective .txt files ([MODD/Data/DD_Data_eig_opt_Conc.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Data/DD_Data_eig_opt_Conc.txt), [MODD/Data/DD_Data_dist_unit_opt_Conc.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Data/DD_Data_dist_unit_opt_Conc.txt), [MODD/Data/DD_Hamiltonians_Conc.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Data/DD_Hamiltonians_Conc.txt), [MODD/Data/DD_Operations_Optimized_eig_Conc.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Data/DD_Operations_Optimized_eig_Conc.txt), [MODD/Data/DD_Operations_Optimized_unitary_Conc.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Data/DD_Operations_Optimized_unitary_Conc.txt), [MODD/Data/DD_Init_States_Conc.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Data/DD_Init_States_Conc.txt), [MODD/Data/DD_Data_Purity_opt_Conc.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Data/DD_Data_Purity_opt_Conc.txt), [MODD/Data/DD_Lindbladians_Conc.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Data/DD_Lindbladians_Conc.txt), [MODD/Data/DD_Operations_Optimized_Purity_Conc.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Data/DD_Operations_Optimized_Purity_Conc.txt))

- [MODD/data_printer_conc.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/data_printer_conc.m)
  Helper function that handles the data printing for [MODD/DD_Sequence_Print_Concatenated.m](). Necessary due to parallelization of the for loops in [MODD/DD_Sequence_Print_Concatenated.m]().

- [MODD/Dist_Unit_Opt_conc.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Dist_Unit_Opt_conc.m):
Given a **Comb**, an initial environment state **Init_Env**, a number of slots **N**, a number of rounds **m** for MODD, a dimension **d**, and an **exitTolerance**, this function returns

    - **OpUnit**: The sequence of ideal operations from MODD
    - **OpFinal**: The sequence of operations optimizing the second round in the MODD protocol (i.e., optimization over the slots that are still open after the first round)
    - **Res_Comb**: The resulting Comb after the first round of MODD
  
  Optimization is via the minimization of the distance to the set of unitaries. Currently, it has mostly been tested for **N=3**, **m=4**, and **dim=2**. Generalization to other parameters is straight forward but requires changes in all functions this function depends on.

- [MODD/Dist_Unit_One_It.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Dist_Unit_One_It.m):
Given a **Comb**, a number of slots **N**, a system dimension **dim**, an initial environment state **Init_Env**, and an **exitTolerance**, this function finds an ideal DD sequence (via minimization to the set of unitaries) and computes the resulting environment state.

So far, it has mostly been tested for **N=3** and **dim=2**. Generalization to other parameters is straight forward but requires changes in all functions this function depends on.

- [MODD/Eigenvalue_Opt_conc.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Eigenvalue_Opt_conc.m):
Same as [MODD/Dist_Unit_Opt_conc.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Dist_Unit_Opt_conc.m) but optimization via maximization of the maximal eigenvalue of the resulting channel.

- [MODD/Eigenvalue_Opt_One_It.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Eigenvalue_Opt_One_It.m):
Same as [MODD/Dist_Unit_One_It.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Dist_Unit_One_It.m) but optimization via maximization of the maximal eigenvalue of the resulting channel.

- [MODD/BaselineSequence_Conc.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/BaselineSequence_Conc.m):
Given a two-qubit Lindbladian (in vectorized form) **Lindblad**, a qubit state **Init_Env** and a step size **delT**, this function returns the the resulting purity and mutial information of the channel resulting from doing-nothing in a 15-slot comb (used as a baseline to compare the decoupling data to).

  - As an output, the program yields
  
    **PurityBase**: Purity of the (Choi matrix of) the resulting channel

    **MutInfBase**: Input-Output mutual information of the (Choi matrix of) the resulting channel

  Currently only set up to work for **d=2** (can be easily changed by making the dimension a parameter)

- [MODD/StandardSequence_Conc.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/StandardSequence_Conc.m):
Given a Lindbladian **Lindblad**, an initial environment state **Init_Env** and a step size **delT**, this function returns the Purity and mutual information for the resulting channel of a 15 step comb if a standard DD sequence is employed (used as a baseline to compare DD protocols to).

Currently only set up to work for **d=2**.

- [MODD/Collect_Data_Python_Conform_Concatinaed.py](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Collect_Data_Python_Conform_Concatinaed.py):
Python script that plots the data contained in the MODD/Data folder. Type of data to be plotted (options: 'Purity', 'MutInf') can be chosen by changing the variable **Type** in the script. 

##### Auxiliary Functions (both in ODD/ and MODD/
  
- [MaxEnt.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/MaxEnt.m):
given a dimension **dim**, an unnormalized maximally entangled state of  dimension dim is returned

- [PartTr.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/PartTr.m):
given a state **rho**, an array **sys** and and array **dim**, this function performs the partial transpose of rho (defined on the systems with dimensions given in dim) over the systems given in sys. Code taken from Tony Cubitt's [homepage](https://www.dr-qubit.org/matlab.html).

- [quantum_entr.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/quantum_entr.m)):
  Given a matrix **X** and optimization parameters **(m,k)** this function returns the von-Neumann entropy of X. m and k do not have to be specified explicitly. Based on https: *Semidefinite approximations of the matrix logarithm* by Hamza Fawzi, James Saunderson and Pablo A. Parrilo, [arXiv:1705.00812](https://arxiv.org/abs/1705.00812)

- [quantum_rel_entr.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/quantum_rel_entr.m)):
given two matrices **A** and **B** optimization parameters **(m,k)**, this function returns the quantum relative entropy between A and B. m and k do not have to be specified explicitly. Based on https: *Semidefinite approximations of the matrix logarithm* by Hamza Fawzi, James Saunderson and Pablo A. Parrilo, [arXiv:1705.00812](https://arxiv.org/abs/1705.00812)

- [Rand_comb.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/Rand_comb.m)):
given a number of slots **N**, a two-qubit Hamiltonian **Ham** and a step size delT, this function samples an initial two-qubit state and returns the resulting comb for the given Hamiltonian and additional Markovian noise (detailed in [arXiv:2110.02613](https://arxiv.org/abs/2110.02613)) and the step size. Currently, this function **only works for two-qubit Hamiltonians**.

- [Rand_state.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/Rand_state.m)):
given system size **n**, this function returns a random (according to Haar measure) pure quantum state.

- [randU.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/randU.m)):
given **n**, this function returns a random (according to Haar measure) n x n unitary. Code taken from Tony Cubitt's [homepage](https://www.dr-qubit.org/matlab.html).

- [syspermute.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/syspermute.m)):
given a matrix **p**, a permutation **perm** and dimensions of subsystems **dim**, this function returns a permutation of p (defined on susbsystems of dimensions given by dim) according to the permutation perm. Code taken from Tony Cubitt's [homepage](https://www.dr-qubit.org/matlab.html).

- [TnProduct.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/TnProduct.m)):
given a **list of input objects**, this function returns their tensor product. Code taken from Tony Cubitt's [homepage](https://www.dr-qubit.org/matlab.html).

- [TrX.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/TrX.m)):
given a matrix **p**, a list of systems **sys** and a list **dim** of subsystem dimensions, this function returns the partial trace of p (defined on subsystems wih dimensions given by dim) over the systems given by sys. Code taken from Tony Cubitt's [homepage](https://www.dr-qubit.org/matlab.html).

- [chanconv.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/chanconv.m)):
Given a matrix **E**, this function can change its representation by choosing **from** and **to** for the the represenation it is in and the one it is supposed to be transformed to. Code taken from Tony Cubitt's [homepage](https://www.dr-qubit.org/matlab.html).

- [Rand_comb_Var_In.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Rand_comb_Var_In.m)):
  Given a number of slots **N**, a Hamiltonian **N**, and a step size **delT**, this function returns an N-slot comb, built from the Hamiltonian and the noise model of arXiv:2110.02613 that has an open initial and final environment line.

  Currently only works for **d=2**, but extension to higher system dimensions is trivially possible.

#### DATA FILES

###### ODD (single round optimization)

- [ODD/Data/DD_Data_dist_unit_opt.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/Data/DD_Data_dist_unit_opt.txt):
  Text file containing data created by [ODD/DD_Sequence_Print.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/DD_Sequence_Print.m) for single round (i.e., three slots ODD optimization) via minimization of the distance to the unitaries. Contains the following columns:

  **Identifier**: Unique identifier for every line for easier further handling of the data
  
  **DeltaT**: Employed temporal step size in the creation of the data
  
  **PurityOpt**: Purity of the final channel for the optimized DD sequence
  
  **PurityStand**: Purity of the final channel under standard DD
  
  **MutInfOpt**:  Mutual Information of the final channel for the optimized DD sequence

  **MutInfStand**: Mutual Information of the final channel under standard DD
  
  **Purity1, Purity2, Purity3**: Purity of the three employed decoupling sequences

- [ODD/Data/DD_Data_eig_opt.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/Data/DD_Data_eig_opt.txt):
    Same as [ODD/Data/DD_Data_dist_unit_opt.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/Data/DD_Data_dist_unit_opt.txt), but for optimization of the maximum eigenvalue of the resulting channel. Created by [ODD/DD_Sequence_Print.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/DD_Sequence_Print.m).

- [ODD/Data/DD_Hamiltonians.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/Data/DD_Hamiltonians.txt):
  Text File containing the Hamiltonians used for the optimization of the DD sequences. Created by [ODD/DD_Sequence_Print.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/DD_Sequence_Print.m). Contains:

  **Identifier**: Unique identifier to match Hamiltonian to run of simulation
  
  **DelT**: Employed temporal step size

  **Hamiltonian**: 4x4 Hamiltonian matrix that was used in the run of the simulation

- [ODD/Data/DD_Init_States.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/Data/DD_Init_States.txt):
  Text File containing the initial states used for the optimization of the DD sequences. Created by [ODD/DD_Sequence_Print.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/DD_Sequence_Print.m). Contains:
  
  **Identifier**: Unique identifier to match Hamiltonian to run of simulation
  
  **DelT**: Employed temporal step size

  **State**: 2x2 initial environment state that was used in the run of the simulation

- [ODD/Data/DD_Operations_Optimized_unitary.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/Data/DD_Operations_Optimized_unitary.txt):
 Text File containing the optimized operations for DD Sequences, from minimization to the set of unitaries of the resulting channel. Created by [ODD/DD_Sequence_Print.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/DD_Sequence_Print.m). Contains:

**Identifier**: Unique identifier to match sequence to run of simulation

**Op1, Op2, Op3**: Three optimal operations for the given run of the simulation

- [DataFilesODD/DD_Operations_Optimized_eig.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/Data/DD_Operations_Optimized_eig.txt):
 Same as [ODD/Data/DD_Operations_Optimized_unitary.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/Data/DD_Operations_Optimized_unitary.txt), but for optimization of the maximum eigenvalue of the resulting channel. Created by [ODD/DD_Sequence_Print.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/ODD/DD_Sequence_Print.m).


###### MODD (multi-round optimization)

- [MODD/Data/DD_Data_dist_unit_opt_Conc.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Data/DD_Data_dist_unit_opt_Conc.txt):
  
Text file containing data created by [MODD/DD_Sequence_Print_Concatenated.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/DD_Sequence_Print_Concatenated.m) for multi-round round (i.e., 4 times three slots MODD optimization) via minimization of the distance to the unitaries. Contains the following columns:

  **Identifier**: Unique identifier for every line for easier further handling of the data 
  
  **DeltaT**: Employed temporal step size in the creation of the data
  
  **PurityOpt**: Purity of the final channel for the optimized DD sequence
  
  **PurityStand**: Purity of the final channel under standard DD
  
  **MutInfOpt**:  Mutual Information of the final channel for the optimized DD sequence

  **MutInfStand**: Mutual Information of the final channel under standard DD
 
  **MutInfBase**: Mutual Information of the final channel under do-nothing operations

- [MODD/Data/DD_Data_eig_opt_Conc.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Data/DD_Data_eig_opt_Conc.txt):
Same as [MODD/Data/DD_Data_dist_unit_opt_Conc.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Data/DD_Data_dist_unit_opt_Conc.txt), but for optimization of the maximum eigenvalue of the resulting channel. Created by [MODD/DD_Sequence_Print_Concatenated.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/DD_Sequence_Print_Concatenated.m).
 
- [MODD/Data/DD_Data_Purity_opt_Conc.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Data/DD_Data_Purity_opt_Conc.txt):
Same as [MODD/Data/DD_Data_dist_unit_opt_Conc.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Data/DD_Data_dist_unit_opt_Conc.txt), but for optimization of the purity of the resulting channel (generated using an algorithm for purity maximization: J. Morris, *in preparation* (2023)).

- [MODD/Data/DD_Operations_Optimized_unitary_Conc.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Data/DD_Operations_Optimized_unitary_Conc.txt):
  Text File containing the optimized operations for DD Sequences, from minimization to the set of unitaries of the resulting channel via MODD. Created by [MODD/DD_Sequence_Print_Concatenated.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/DD_Sequence_Print_Concatenated.m). Contains:

**Identifier**: Unique identifier to match sequence to run of simulation

**Op1, ..., Op15**: 15 optimal operations for the given run of the simulation

- [MODD/Data/DD_Operations_Optimized_eig_Conc.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Data/DD_Operations_Optimized_eig_Conc.txt):
  Same as [DataFilesMODD/DD_Operations_Optimized_unitary_Conc.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Data/DD_Operations_Optimized_unitary_Conc.txt), but for optimization of the maximum eigenvalue of the resulting channel. Created by [MODD/DD_Sequence_Print_Concatenated.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/DD_Sequence_Print_Concatenated.m).

- [MODD/Data/DD_Operations_Optimized_Purity_Conc.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Data/DD_Operations_Optimized_Purity_Conc.txt):
Same as [DataFilesMODD/DD_Operations_Optimized_unitary_Conc.txt](), but for optimization of the purity of the resulting channel (generated using an algorithm for purity maximization: J. Morris, *in preparation* (2023)).

- [MODD/Data/DD_Hamiltonians_Conc.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Data/DD_Hamiltonians_Conc.txt):
  Text File containing the Hamiltonians used for the optimization of the DD sequences in MODD scenario. Created by [MODD/DD_Sequence_Print_Concatenated.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/DD_Sequence_Print_Concatenated.m). Contains:

  **Identifier**: Unique identifier to match Hamiltonian to run of simulation

  **Hamiltonian**: 4x4 Hamiltonian matrix that was used in the run of the simulation

- [MODD/Data/DD_Lindbladians_Conc.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Data/DD_Lindbladians_Conc.txt):
Text File containing the Lindbladians used for the optimization of the DD sequences in MODD scenario. Created by [MODD/DD_Sequence_Print_Concatenated.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/DD_Sequence_Print_Concatenated.m). Contains:
**Identifier**: Unique identifier to match Hamiltonian to run of simulation

**Hamiltonian**: 16x16 Lindbladian in vectorized form.

- [MODD/Data/DD_Init_States_Conc.txt](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/Data/DD_Init_States_Conc.txt):
Text File containing the initial states used for the optimization of the DD sequences. Created by [MODD/DD_Sequence_Print_Concatenated.m](https://github.com/simoncmilz/Extracting_Quantum_Dynamical_Resources/blob/main/MODD/DD_Sequence_Print_Concatenated.m). Contains:
  
  **Identifier**: Unique identifier to match Hamiltonian to run of simulation

  **State**: 2x2 initial environment state that was used in the run of the simulation

