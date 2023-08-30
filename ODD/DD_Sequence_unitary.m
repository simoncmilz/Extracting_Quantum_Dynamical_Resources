function [PurityOpt, PurityStand, MutInfOpt, MutInfStand, Purities, OpUnit1, OpUnit2, OpUnit3]= DD_Sequence_unitary(Comb, N, dim, exitTolerance)

%Function that optimizes the sequence of decoupling operations for a given
%comb with N slots, where the system size is d for all slots. Works via a
%see saw algorithm that minimizes the distance to the identity channel in 
%every iteration. Important: Algorithm also optimizes over initial and
%final time, i.e., altogether over sequence of N+2 operations, thus aiming
%to minimize the distance to the set of unitary operations.
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


%if dimension is 2 and N is 3, start optimization with standard DD sequence
if N==3 & dim==2
    sigmaX = [[0 1];[1 0]];
    sigmaZ = [[1 0];[0 -1]];
    
    %Choi of Pauli matrices
    SX = kron(sigmaX,eye(2))*MaxEnt(2)*kron(sigmaX,eye(2));
    SZ = kron(sigmaZ,eye(2))*MaxEnt(2)*kron(sigmaZ,eye(2));

    %Add initial and final channel
    DD_Seq = {MaxEnt(2),SX,SZ,SX,SZ};
    DD_Seq2 = DD_Seq;   %Copy for later use

    %otherwise start with random sequence
    else   
        %Sample initial channels for DD sequence
        DD_Seq = cell(N+2, 1) ;
        for k = 1 : N+2
           DD_Seq{k} = randChan(dim,'choi');
        end

        DD_Seq2= DD_Seq;    %Copy for later use
end

%Array with the respective spatial dimensions of the slots and initial
%input/final output
dims = dim^2*ones([1,N+2]);
dims = [dim,dims,dim];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exit = 0;
Vals = [10, 0]; %Initial dummy values for exit condition check
counter = 0; %loop counter for additional exit condition

%add identity matrices at beginning and end to match dimension of the comb
%and the sequence of channels
Comb2 = Comb; %Copy for later use
Comb = TnProduct(eye(dim),Comb,eye(dim));

%Loop over optimization as long as the relative change still exceeds
%exitTolerance
while ~exit
    counter = counter +1;
    %one round of optimization over each operation of the sequence
    %individually:
    for i = 1:N+2
        %Create tensor product of operations in DD_Seq, replace the one to be
        %optimzed by identity matrix
        Seq=1;
        for j = 1:N+2
            if j ~= i
            Seq = kron(Seq,DD_Seq{j});
            else
                Seq = kron(Seq,eye(dim^2));
            end
        end
        
        %If optimization over first CPTP map
        if i == 1
            %array of systems over which we contract
            sys = 2:N+2;
            %dimensions of subsystems
            dims = dim^2*ones([1,N+2]);
            dims(N+2) = dim;
            dims(N+3) = dim;
            %Contract with Comb to obtain one-slot comb
            One_slot = TrX(mtimes(Comb,PartTr(Seq,sys,dims)),sys,dims);
            
            %Optimize over the remaining open slot          
            cvx_begin quiet
            %variables for the CPTP map in the slot
            variable Lambda(dim^2,dim^2) hermitian semidefinite;

            %Criterion: Minimize distance of the resulting comb to the identity
            %channel
            minimize norm(TrX(mtimes(TnProduct(One_slot),TnProduct(Lambda,eye(dim))),[2],[dim,dim,dim]) - MaxEntSt)
            
            subject to
                %partial trace condition
                TrX(Lambda,[2],[dim, dim]) == eye(dim);
            cvx_end
            
        end
        
        %If optimization over CPTP map in the 'bulk'
        if i ~= N+2 & i~=1
            %array of systems over which we contract
            sys = 2:N+3;
            sys(i) = [];
            %dimensions of subsystems
            dims = [dim,dim,dim^2*ones([1,N]),dim,dim];
            %Contract with Comb to obtain one-slot comb
            One_slot = TrX(mtimes(Comb,transpose(Seq)),sys,dims);
            
            %Optimize over the remaining open slot          
            cvx_begin quiet
            %variables for the CPTP map in the slot
            variable Lambda(dim^2,dim^2) hermitian semidefinite;

            %Criterion: Minimize distance of the resulting comb to the identity
            %channel
            minimize norm(TrX(mtimes(One_slot,TnProduct(eye(dim),Lambda,eye(dim))),[2],[dim,dim^2,dim]) - MaxEntSt)

            subject to
                %partial trace condition
                TrX(Lambda,[2],[dim, dim]) == eye(dim);
            cvx_end
        end
        
        %If optimization over final CPTP map    
        if i == N+2
            %array of systems over which we contract
            sys = 2:N+2;
            %dimensions of subsystems
            dims = dim^2*ones([1,N+1]);
            dims = [dim,dim,dims];
            %Contract with Comb to obtain one-slot comb
            One_slot = TrX(mtimes(Comb,transpose(Seq)),sys,dims);
            
            %Optimize over the remaining open slot          
            cvx_begin quiet
            %variables for the CPTP map in the slot
            variable Lambda(dim^2,dim^2) hermitian semidefinite;

            %Criterion: Minimize distance of the resulting comb to the identity
            %channel
            minimize norm(TrX(mtimes(TnProduct(One_slot),TnProduct(eye(dim),Lambda)),[2],[dim,dim,dim]) - MaxEntSt)
            
            subject to
                %partial trace condition
                TrX(Lambda,[2],[dim, dim]) == eye(dim);
            cvx_end
        end
        
        %update DD_Sequence
        DD_Seq{i} = transpose(Lambda);
    
    end
        
        

%Check if exit condition is fulfilled:
Vals(2) = norm(TrX(mtimes(TnProduct(One_slot),TnProduct(eye(dim),Lambda)),[2],[dim,dim,dim]) - MaxEntSt);
relValChange = norm((Vals(2)-Vals(1))/Vals(1));
if relValChange < exitTolerance
    exit = 1;
else
    if counter > 20
        exit = 1;
    end
end
Vals(1) = Vals(2);
end


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
OverallSys = 2:N+1;
dims = dim^2*ones([1,N]);
dims = [dim,dims,dim];


ChannOpt = TrX(mtimes(Comb2,transpose(SeqOpt)),[OverallSys],dims);
ChannStand = TrX(mtimes(Comb2,transpose(SeqStand)),[OverallSys],dims);

%Normalize for meaningful entropic statements
ChannOpt = ChannOpt/trace(ChannOpt);
TrX(ChannOpt/trace(ChannOpt),[2],[dim,dim]);
ChannStand = ChannStand/trace(ChannStand);

%Compute mutual information between input and output
MutInfOpt = -quantum_entr(ChannOpt) + quantum_entr(TrX(ChannOpt,[2],[dim,dim]))  + quantum_entr(TrX(ChannOpt,[1],[dim,dim]));
MutInfStand = -quantum_entr(ChannStand) + quantum_entr(TrX(ChannStand,[2],[dim,dim]))  + quantum_entr(TrX(ChannStand,[1],[dim,dim]));

%Compute purities:
PurityOpt = trace(ChannOpt*ChannOpt);
PurityStand = trace(ChannStand*ChannStand);

%Return best operations
OpUnit1 = DD_Seq{1};
OpUnit2 = DD_Seq{2};
OpUnit3 = DD_Seq{3};

end 
