function [Env_new, DD_Seq] = Dist_Unit_One_It(Comb, N, dim, Init_Env, exitTolerance)

%Function that optimizes the sequence of decoupling operations for a given
%comb with N slots, where the system size is d for all slots.
%Works via optimizing three operations at a time (N will basically always
%be equal to 3), and then doing this for m times, with a final optimization
%over the remaining open slots.
%
%Returns:
%
%Env_new: New environment state after sequence
%DD_Seq: Sequence of decoupling operations
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

    SX = kron(sigmaX,eye(2))*MaxEnt(2)*kron(sigmaX,eye(2));
    SZ = kron(sigmaZ,eye(2))*MaxEnt(2)*kron(sigmaZ,eye(2));

    %Add initial and final channel
    DD_Seq = {MaxEntSt,SX,SZ,SX,SZ};

    %otherwise start with random sequence
    else   
        %Sample initial channels for DD sequence
        DD_Seq = cell(N+2, 1) ;
        for k = 1 : N+2
           DD_Seq{k} = randChan(dim,'choi');
        end

end

%Array with the respective spatial dimensions
dims = dim^2*ones([1,N+2]);

%%%%%%%%%%%%%%%%%%
%Main
%%%%%%%%%%%%%%%%%%

%Compute initial comb by contracting it with the initial environment state
%and tracing out the environment
Res_Comb = TrX(Comb*TnProduct(eye(dim),transpose(Init_Env), eye(dim^(2*N+2))), [2 4], [dim dim dim^(2*N+1) dim]);
%pad out to match dimensions of the sequence of operations
Res_Comb = TnProduct(eye(dim),Res_Comb,eye(dim));

%Looping Parameters
exit = 0;
Vals = [10, 0]; %Initial values for exit condition check
counter = 0; %counter for additional exit condition


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
            Seq = kron(Seq,full(DD_Seq{j}));
            else
                Seq = kron(Seq,eye(dim^2));
            end
        end
         
        
        if i == 1
            %array of systems over which we contract
            sys = 2:N+2;
            %dimensions of subsystems
            dims = dim^2*ones([1,N+1]);
            dims(N+2) = dim;
            dims(N+3) = dim;
            %Contract with Comb to obtain channel comb
            One_slot = TrX(Res_Comb*PartTr(Seq,sys,dims),sys,dims);
            
            %Optimize over the remaining open slot          
            cvx_begin quiet
            %variables for the CPTP map in the slot
            variable Lambda(dim^2,dim^2) hermitian semidefinite;

            %Criterion: Minimize distance of the resulting comb to the identity
            %channel
            minimize norm(TrX(One_slot*TnProduct(PartTr(Lambda,[2], [dim dim]),eye(dim)),[2],[dim,dim,dim]) - MaxEntSt)
            
            subject to
                %partial trace condition
                TrX(Lambda,[2],[dim, dim]) == eye(dim);
            cvx_end
        end 
        
        if i ~= N+2 & i~=1
            %array of systems over which we contract
            sys = 2:N+3;
            sys(i) = [];
            %dimensions of subsystems
            dims = [dim,dim,dim^2*ones([1,N]),dim,dim];
            %Contract with Comb to obtain one-slot comb
            One_slot = TrX(Res_Comb*PartTr(Seq, sys, dims),sys,dims);
            
            %Optimize over the remaining open slot          
            cvx_begin quiet
            %variables for the CPTP map in the slot
            variable Lambda(dim^2,dim^2) hermitian semidefinite;

            %Criterion: Minimize distance of the resulting comb to the identity
            %channel
            minimize norm(TrX(mtimes(One_slot,TnProduct(eye(dim),transpose(Lambda),eye(dim))),[2],[dim,dim^2,dim]) - MaxEntSt)

            subject to
                %partial trace condition
                TrX(Lambda,[2],[dim, dim]) == eye(dim);
            cvx_end
        end
        
            
        if i == N+2
            %array of systems over which we contract
            sys = 2:N+2;
            %dimensions of subsystems
            dims = dim^2*ones([1,N+1]);
            dims = [dim,dim,dims];
            %Contract with Comb to obtain one-slot comb
            One_slot = TrX(Res_Comb*PartTr(Seq, sys, dims),sys,dims);
            
            %Optimize over the remaining open slot          
            cvx_begin quiet
            %variables for the CPTP map in the slot
            variable Lambda(dim^2,dim^2) hermitian semidefinite;

            %Criterion: Minimize distance of the resulting comb to the identity
            %channel
            minimize norm(TrX(One_slot*TnProduct(eye(dim),PartTr(Lambda,[1],[dim dim])),[2],[dim,dim,dim]) - MaxEntSt)
            
            subject to
                %partial trace condition
                TrX(Lambda,[2],[dim, dim]) == eye(dim);
            cvx_end
        end
        
        %update DD_Sequence
        DD_Seq{i} = Lambda;
    
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


%Compute resulting environment state
%Make sequence of operations
%tensor the channels in DD sequence without last operation
Seq = 1;
for j = 1:N+1
    Seq = kron(Seq,full(DD_Seq{j}));
end
Seq2 = Seq;

%Tensor sequence with final environment state and identity to match comb
Seq = TnProduct(Init_Env,Seq,eye(dim^2));
%Reshuffle to get spaces in right order
Seq = syspermute(Seq,[2 3 1 4], [dim dim dim dim^(2*N + 2)]);

%Systems over which to contract
SysContr = [1 2 3 4 5];

%Dimensions of all systems
Dims = [dim dim dim dim^(2*N) dim dim];

%Contract original comb with sequence
Env_new = 1/dim*TrX(TnProduct(eye(2),Comb)*PartTr(Seq, SysContr, Dims),SysContr, Dims);

DD_Seq(1) = [];
DD_Seq(N+1) = [];


end 
