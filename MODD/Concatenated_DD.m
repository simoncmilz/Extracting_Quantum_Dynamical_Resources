function[MutInfOpt, PurityOpt, Op, Lindblad] = Concatenated_DD(m, N,Hamil,Init_Env, delT, exitTolerance, optimization)
%Function that finds the ideal decoupling sequence for m rounds of
%decoupling (each comb of the round corresponds to N=3 slots). The dimension 
%of the system is set to 2. Hamil is the employed system-environment hamiltonian
%and delT the step size between slots. The 'optimization' argument
%determines what type of optimization is used.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Preliminaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim = 2; %dimension of the system

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute ideal sequences for a given comb according to the optimization
%method that is given
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Produce Comb with open start and end line
[Comb, Lindblad] = Rand_comb_Var_In(N, Hamil, delT);

%Possible optimizations: 
%
%Maximization of the largest eigenvalue of the resulting channel
%
%Minimization of the distance to the unitaries
%
%Maximization of the purity of the resulting channel

if  strcmp(optimization,'eigenvalue')
    [Op, OpFinal, Res_Comb] = Eigenvalue_Opt_conc(Comb, Init_Env, N, m, dim, exitTolerance);
end
if strcmp(optimization,'dist_unit')
    [Op, OpFinal, Res_Comb] = Dist_Unit_Opt_conc(Comb, Init_Env, N, m, dim, exitTolerance);
end
if strcmp(optimization, 'purity')
    [Op, OpFinal, Res_Comb] = Purity_Opt_conc(Comb, Init_Env, N, m, dim, exitTolerance);
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute relevant parameters, like mutual info of resulting channel and so
%on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
%Preliminaries
%%%%%%%%%%%%%%%%%%

%Compute channel resulting from optimized operations

%Create Sequence corresponding to collection of optimal operations
SeqOpt = 1;  %Optimized Sequence
Purities = []; %Array for the Purities of the optimal maps in the sequence

%Make Sequence of the best operations and store the purities of the
%operations
for j = 1:N
    SeqOpt = kron(SeqOpt,OpFinal{j});
    Purities(j) = trace(OpFinal{j}*OpFinal{j})/(trace(OpFinal{j})^2);
end

%Compute resulting Channel and normalize
ChannOpt = TrX(Res_Comb*transpose(TnProduct(eye(dim),SeqOpt,eye(dim))), [2], [dim dim^(2*N) dim]);
ChannOpt = ChannOpt/trace(ChannOpt);

%Compute resulting purity and mutual information

%Purity
PurityOpt = trace(ChannOpt*ChannOpt);

%Mutual information
MutInfOpt = -quantum_entr(ChannOpt) + quantum_entr(TrX(ChannOpt,[2],[dim,dim]))  + quantum_entr(TrX(ChannOpt,[1],[dim,dim]));

end

