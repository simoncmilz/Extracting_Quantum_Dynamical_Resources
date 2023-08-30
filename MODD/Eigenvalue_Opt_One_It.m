function [Env_new, DD_Seq]= Eigenvalue_Opt_One_It(Comb, N, dim, Init_Env, exitTolerance)

%Function that optimizes the sequence of decoupling operations for a given
%comb with N slots, where the system size is d for all slots.
%Works via optimizing three operations at a time (N will basically always
%be equal to 3), and then doing this for m times, with a final optimization
%over the remaining open slots.

%Returns:
%
%Env_new: New environment state after sequence
%DD_Seq: Sequence of decoupling operations
%
%Arguments:
%
%Comb: Comb for which the best sequence is to be found. Has to open lines
%in beginning and end
%N: Number of slots
%dim: system dimension
%Init_Env: Initial environment state
%exitTolerance: Tolerance at which to end the computation.
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


%if dimension is 2 and N is 4, start optimization with DD sequence
if N==3 & dim==2
sigmaX = [[0 1];[1 0]];
sigmaZ = [[1 0];[0 -1]];

SX = kron(sigmaX,eye(2))*MaxEnt(2)*kron(sigmaX,eye(2));
SZ = kron(sigmaZ,eye(2))*MaxEnt(2)*kron(sigmaZ,eye(2));

DD_Seq = {SX,SZ,SX,SZ};

%otherwise start with random sequence
else   
%Sample initial channels for DD sequence
DD_Seq = cell(N+1, 1);
for k = 1 : N+1
   DD_Seq{k} = randChan(dim,'choi');
end
end



%%%%%%%%%%%%%%%%%%
%Main
%%%%%%%%%%%%%%%%%%

%Compute initial comb by contracting it with the initial environment state
%and tracing out the environment
Res_Comb = TrX(Comb*TnProduct(eye(dim),transpose(Init_Env), eye(dim^(2*N+2))), [2 4], [dim dim dim^(2*N+1) dim]);


%Looping parameters
exit = 0;
Vals1 = 100; %Initial values for exit condition check
counter = 0; %counter for additional exit condition
counterExit = 20;
depth = 10; %depth for largest eigenvalue optimization
EigenVal1 = 100; %initial parameters for eigenvalue optimization
relEigExit = 0.0001; %Exit parameter for optimization of eigenvalues


%Compute channel that stems from standard DD_sequence to have good initial
%guess for projector E

%tensor the channels in DD sequence
Seq = 1;
for j = 1:N+1
    Seq = kron(Seq,DD_Seq{j});
end

%add identity at beginning to sequence
Seq = TnProduct(eye(dim),Seq);
%Add identity to end of the comb to match dimension
Res_Comb = TnProduct(Res_Comb,eye(dim));

%Contract Comb with sequence
OverallSys = 2:N+2; %systems that get contracted over

%Array with the respective spatial dimensions
dims = dim^2*ones([1,N]);
dims = [dim,dims, dim, dim];
%Contract comb with sequence, compute eigenvector of largest eigenvalue
Chan = TrX(Res_Comb*PartTr(Seq, [OverallSys],dims),[OverallSys],dims);
[E,lam] = eigs(Chan,1);
E = kron(E',E);


%Loop over optimization as long as the relative change still exceeds
%exitTolerance or is below counterExit
while ~exit
    counter = counter +1;
    %one round of optimization over each operation of the sequence
    %individually:
    for i = 1:N+1
        %Create tensor product of operations in DD_Seq, replace the one to be
        %optimzed by identity matrix
        Seq=1;
        for j = 1:N+1
            if j ~= i
            Seq = kron(Seq,full(DD_Seq{j}));
            else
                Seq = kron(Seq,eye(dim^2));
            end
        end
        %add identity at beginning to match dimension of the comb
        Seq = TnProduct(eye(dim),Seq);

        
        %Contract with Comb to obtain one-slot comb
        %array of systems over which we contract
        sys = 2:N+2;
        sys(i) = [];
        %Contraction
        One_slot = TrX(Res_Comb*PartTr(Seq, sys, dims),sys,dims);
        
        %Optimize over the remaining open slot by maximizing the largest
        %eigenvalue via see-saw
        count2 = 0;
        relEigChange = relEigExit+1; %reset value
        while relEigChange > relEigExit & count2 < 40
            count2 = count2 + 1;
            
            %Distinguish between whether operation to be optimized over is 
            %last or not
            
            %last
            if i == N+1
                cvx_begin quiet
                %variables for the CPTP map after the resultant channel
                variable Lambda(dim^2,dim^2) hermitian semidefinite;

                %Criterion: Maximize projection in direction of E
                maximize real(trace(E*TrX(One_slot*TnProduct(eye(dim),PartTr(Lambda,[1], [dim dim])),[2],[dim,dim,dim])))

                subject to
                    %partial trace condition
                    TrX(Lambda,[2],[dim, dim]) == eye(dim);
                cvx_end
            
            %not last
            else
                cvx_begin quiet
                %variables for the CPTP map in the slot
                variable Lambda(dim^2,dim^2) hermitian semidefinite;

                %Criterion: Maximize projection in direction of E
                maximize real(trace(E*TrX(mtimes(One_slot,TnProduct(eye(dim),transpose(Lambda),eye(dim))),[2],[dim,dim^2,dim])))
            
                subject to
                    %partial trace condition
                    TrX(Lambda,[2],[dim, dim]) == eye(dim);
                cvx_end
            end
            
            
            %Optimization over projectors
            
            %Distinction between whether last operation or not
            
            %last
            if i == N+1
                cvx_begin quiet
                %variable for the projector for largest eigenvalue
                variable E(dim^2,dim^2) hermitian semidefinite;

                %Criterion: maximize Projection
                maximize real(trace(E*TrX(One_slot*TnProduct(eye(dim),PartTr(Lambda,[1],[dim dim])),[2],[dim,dim,dim])))

                subject to
                    %'Projector' property
                    trace(E) == 1;
                cvx_end
            
            %not last
            else
                cvx_begin quiet
                %variable for the projector for largest eigenvalue
                variable E(dim^2,dim^2) hermitian semidefinite;

                %Criterion: maximize Projection
                maximize real(trace(E*TrX(One_slot*TnProduct(eye(dim),transpose(Lambda),eye(dim)),[2],[dim,dim^2,dim])))

                subject to
                    %'Projector' property
                    trace(E) == 1;
                cvx_end
            end
        
        %compute relative change of eigenvalues
        if i==N+1
            EigenVal2 = real(trace(E*TrX(One_slot*TnProduct(eye(dim),PartTr(Lambda,[1],[dim dim])),[2],[dim,dim,dim])));
        else
            EigenVal2 =  real(trace(E*TrX(One_slot*TnProduct(eye(dim),transpose(Lambda),eye(dim)),[2],[dim,dim^2,dim])));
        end
        
        relEigChange = norm((EigenVal2 - EigenVal1)/EigenVal1);
        EigenVal1 = EigenVal2;
        end
        %update DD_Sequence
        DD_Seq{i} = Lambda;
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

%Compute resulting environment state
%Make sequence of operations
%tensor the channels in DD sequence without last operation
Seq = 1;
for j = 1:N
    Seq = kron(Seq,full(DD_Seq{j}));
end 

%Tensor sequence with initial environment state and identities to match comb
Seq = TnProduct(eye(dim),Init_Env,Seq,eye(dim^2));

%Systems over which to contract
SysContr = 1:N+3;

%Dimensions of all systems
Dims = dim^2*ones([1,N]);
Dims = [dim, dim, Dims, dim, dim];

%Contract original comb with sequence
Env_new = 1/dim*TrX(Comb*transpose(Seq),SysContr, Dims);
DD_Seq(N+1) = [];

end

