function [PurityBase, MutInfBase] = BaselineSequence_Conc(Lindblad, Init_Env, delT)

%Function that computes the purity and mutial information for the case of
%do-nothing operations
dim = 2;
%Unnormalized maxEnt state vector on system/environment (both have same dim)
MaxEntSt = zeros([1,dim^2]);
for i=1:dim
    MaxEntSt((i-1)*dim + i) = 1;
end

%Compute overall resulting system-environment map for base sequence
L = expm(delT*Lindblad);
SysEnvVecBase = L^16;

%Vectorize initial environment state
EnvStateVec = TnProduct(Init_Env,eye(dim))*transpose(MaxEntSt);

%Reshuffle matrices so that system corresponds to first two subsystems
SysEnvVecBase = syspermute(SysEnvVecBase,[1 3 2 4], [dim dim dim dim]);

%Multiply with initial environment state and trace out final environment
SysVecBase = TnProduct(eye(dim^2),MaxEntSt)*SysEnvVecBase*TnProduct(eye(dim^2),EnvStateVec);

%Compute Channel in Choi and normalize
ChannBase = TnProduct(eye(dim^2),MaxEntSt)*TnProduct(eye(dim),SysVecBase,eye(dim))*TnProduct(transpose(MaxEntSt),eye(dim^2));
ChannBase = ChannBase/trace(ChannBase);

%Compute Mutual info and purity of the resulting channel

%Purity
PurityBase = trace(ChannBase*ChannBase);

%Compute mutual informations
MutInfBase = -quantum_entr(ChannBase) + quantum_entr(TrX(ChannBase,[2],[dim,dim]))  + quantum_entr(TrX(ChannBase,[1],[dim,dim]));

end
