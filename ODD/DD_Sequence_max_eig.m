function [PurityOpt, PurityStand, MutInfOpt, MutInfStand, Purities, OpEig1, OpEig2, OpEig3]= DD_Sequence_max_eig(Comb, N, dim, exitTolerance)

%Function that optimizes the sequence of decoupling operations for a given
%comb with N slots, where the system size is d for all slots. Works via a
%see saw algorithm that maximizes the maximal eigenvalue of the resulting channel. 
%
%Arguments:
%
%Comb: Comb for which the best sequence is to be found
%N: Number of slots
%dim: system dimension
%
%Requires:
%MaxEnt
%randChann
%randV, chanconv
%TrX

%%%%%%%%%%%%%%%%%%
%Preliminaries
%%%%%%%%%%%%%%%%%%

MaxEntSt = MaxEnt(dim);   %Choi of identity for optimization


%if dimension is 2 and N is 3, start optimization with DD sequence
if N==3 & dim==2
    sigmaX = [[0 1];[1 0]];
    sigmaZ = [[1 0];[0 -1]];

    %Choi matrix of the Pauli matrices
    SX = kron(sigmaX,eye(2))*MaxEnt(2)*kron(sigmaX,eye(2));
    SZ = kron(sigmaZ,eye(2))*MaxEnt(2)*kron(sigmaZ,eye(2));

    DD_Seq = {SX,SZ,SX};
    DD_Seq2 = DD_Seq; %store for later use

%otherwise start with random sequence
    else   
    %Sample initial channels for DD sequence
        DD_Seq = cell(N, 1);
        for k = 1 : N
            DD_Seq{k} = randChan(dim,'choi');
        end
    DD_Seq2= DD_Seq; %Store for later use
end

%Array with the respective spatial dimensions
dims = dim^2*ones([1,N]);
dims = [dim,dims,dim];

%%%%%%%%%%%%%%%%%%
%Main
%%%%%%%%%%%%%%%%%%

%Looping parameters
exit = 0;
Vals1 = 100; %Initial dummy value for exit condition check
counter = 0; %counter for additional exit condition
counterExit = 20;
depth = 10; %depth for largest eigenvalue optimization
EigenVal1 = 100; %initial dummy parameter for eigenvalue optimization
relEigExit = 0.0001; %Exit parameter for optimization of eigenvalues


%Compute channel that stems from standard DD_sequence to have good initial
%guess for projector E onto largest eigenvalue

%tensor the channels in DD sequence
Seq = 1;
for j = 1:N
    Seq = kron(Seq,DD_Seq{j});
end

%add identities at beginning and end to match dimension of the comb
Seq = TnProduct(eye(dim),Seq,eye(dim));
OverallSys = 2:N+1; %systems that get contracted over

%Contract comb with sequence, compute eigenvector of largest eigenvalue
Chan = TrX(mtimes(Comb,transpose(Seq)),[OverallSys],dims);
[E,lam] = eigs(Chan,1);
E = kron(E',E);


%Loop over optimization as long as the relative change still exceeds
%exitTolerance or is below counterExit
while ~exit
    counter = counter +1;
    %one round of optimization over each operation of the sequence
    %individually:
    for i = 1:N
        %Create tensor product of operations in DD_Seq, replace the one to be
        %optimzed by identity matrix
        Seq=1;
        for j = 1:N
            if j ~= i
            Seq = kron(Seq,DD_Seq{j});
            else
                Seq = kron(Seq,eye(dim^2));
            end
        end
        %add identities at beginning and end to match dimension of the comb
        Seq = TnProduct(eye(dim),Seq,eye(dim));

        %array of systems over which we contract
        sys = 1:N+1;
        sys(i+1) = [];
        sys(1)=[];

        %Contract with Comb to obtain one-slot comb
        One_slot = TrX(mtimes(Comb,transpose(Seq)),sys,dims);

        %Optimize over the remaining open slot by maximizing the largest
        %eigenvalue via see-saw
        count2 = 0;
        relEigChange = relEigExit+1; %reset value
        while relEigChange > relEigExit
            count2 = count2 + 1; 
            
            cvx_begin quiet
            %variables for the CPTP map in the slot
            variable Lambda(dim^2,dim^2) hermitian semidefinite;

            %Criterion: Maximize projection in direction of E
            maximize real(trace(E*TrX(mtimes(One_slot,TnProduct(eye(dim),Lambda,eye(dim))),[2],[dim,dim^2,dim])))

            subject to
                %partial trace condition
                TrX(Lambda,[2],[dim, dim]) == eye(dim);
            cvx_end
            
            %Optimization over projectors
            cvx_begin quiet
            %variable for the projector for largest eigenvalue
            variable E(dim^2,dim^2) hermitian semidefinite;

            %Criterion: maximize Projection
            maximize real(trace(E*TrX(mtimes(One_slot,TnProduct(eye(dim),Lambda,eye(dim))),[2],[dim,dim^2,dim])))

            subject to
                %'Projector' property
                trace(E) == 1;
            cvx_end
        
        %compute relative change of eigenvalues
        EigenVal2 = real(trace(E*TrX(mtimes(One_slot,TnProduct(eye(dim),Lambda,eye(dim))),[2],[dim,dim^2,dim])));
        relEigChange = norm((EigenVal2 - EigenVal1)/EigenVal1);
        EigenVal1 = EigenVal2;
        end
        %update DD_Sequence
        DD_Seq{i} = transpose(Lambda);
    end
    %Check if exit condition is fulfilled:
    Vals2 = EigenVal2;
    relValChange = norm((Vals2-Vals1)/Vals1);
    if relValChange < exitTolerance
        exit = 1;
    else
        if counter > 20
            exit = 1;
        end
    end
    Vals1 = Vals2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Check if resulting DD_Sequence outperforms standard one/starting one in
%terms of the mutual information between input and output

SeqOpt = 1;  %Optimized Sequence
SeqStand = 1;  %Standard Sequence
Purities = []; %Array for the Purities of the optimal maps in the sequence

for j = 1:N
    SeqOpt = kron(SeqOpt,DD_Seq{j});
    Purities(j) = trace(DD_Seq{j}*DD_Seq{j})/(trace(DD_Seq{j})^2);
    SeqStand = kron(SeqStand,DD_Seq2{j});
end

%Pad out with identities
SeqOpt = TnProduct(eye(dim),SeqOpt,eye(dim));
SeqStand = TnProduct(eye(dim),SeqStand,eye(dim));

%Compute resulting Channel and normalize
ChannOpt = TrX(mtimes(Comb,transpose(SeqOpt)),[OverallSys],dims);
ChannOpt = ChannOpt/trace(ChannOpt);

ChannStand = TrX(mtimes(Comb,transpose(SeqStand)),[OverallSys],dims);
ChannStand = ChannStand/trace(ChannStand);

%Compute mutual information between input and output
MutInfOpt = -quantum_entr(ChannOpt) + quantum_entr(TrX(ChannOpt,[2],[dim,dim]))  + quantum_entr(TrX(ChannOpt,[1],[dim,dim]));
MutInfStand = -quantum_entr(ChannStand) + quantum_entr(TrX(ChannStand,[2],[dim,dim]))  + quantum_entr(TrX(ChannStand,[1],[dim,dim]));

%Compute purities:
PurityOpt = trace(ChannOpt*ChannOpt);
PurityStand = trace(ChannStand*ChannStand);

%Return best operations
OpEig1 = DD_Seq{1};
OpEig2 = DD_Seq{2};
OpEig3 = DD_Seq{3};


end

