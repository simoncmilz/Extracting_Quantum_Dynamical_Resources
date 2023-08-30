function [PurityStand, MutInfStand] = StandardSequence_Conc(Lindblad, Init_Env, delT)

%Function that computes the purity and mutial information for the case of
%standard decoupling operations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Preliminaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim = 2;

%standard decoupling operations
sigmaX = [[0 1];[1 0]];
sigmaZ = [[1 0];[0 -1]];

%Vectorized versions thereof 
SX = TnProduct(sigmaX,eye(dim),transpose(sigmaX),eye(dim));
SZ = TnProduct(sigmaZ,eye(dim),transpose(sigmaZ),eye(dim));

%Unnormalized maxEnt state vector on system/environment (both have same dim)
MaxEntSt = zeros([1,dim]);
for i=1:dim
    MaxEntSt((i-1)*dim + i) = 1;
end

%Compute overall resulting system-environment map for standard sequence
L = expm(delT*Lindblad);
SysEnvVec = (L*SX*L*SZ*L*SX*L)*(SZ*L*SX*L*SZ*L*SX*L)*(SZ*L*SX*L*SZ*L*SX*L)*(SZ*L*SX*L*SZ*L*SX*L);

%Vectorize initial environment state
EnvStateVec = TnProduct(Init_Env,eye(dim))*transpose(MaxEntSt);

%Reshuffle matrices so that system corresponds to first two subsystems
SysEnvVec = syspermute(SysEnvVec,[1 3 2 4], [dim dim dim dim]);

%Multiply with initial environment state and trace out final environment
SysVec = TnProduct(eye(dim^2),MaxEntSt)*SysEnvVec*TnProduct(eye(dim^2),EnvStateVec);

%Compute Channels in Choi and normalize
ChannStand = TnProduct(eye(dim^2),MaxEntSt)*TnProduct(eye(dim),SysVec,eye(dim))*TnProduct(transpose(MaxEntSt),eye(dim^2));
ChannStand = ChannStand/trace(ChannStand);

%Compute purity and mutual informaion

%Purity
PurityStand = trace(ChannStand*ChannStand);

%Mutual information
MutInfStand = -quantum_entr(ChannStand) + quantum_entr(TrX(ChannStand,[2],[dim,dim]))  + quantum_entr(TrX(ChannStand,[1],[dim,dim]));

end
