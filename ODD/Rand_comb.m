function [ResComb, UsedInit]= Rand_comb(N, Ham, delT)

%Function that returns a comb on a system of dimension dim with N slots.
%Comb is created by using the input Hamiltonian and, together
%with the dissipator Graeme uses, propagate it for time delT and slot it all
%together via a link product. 
%fprintf('Total number of inputs = %d\n',nargin);
%system dimension (could become additional variable in the future
dim = 2;
n = dim^2;

%Create unnormalized maxEnt state vector on system and environment
%(required to switch between representations
MaxEntSt = zeros([1,dim^4]);
for i=1:dim^2
    MaxEntSt((i-1)*dim^2 + i) = 1;
end


%matrices for the dissipator
Sx = [[0 1];[1 0]];
Sy = [[0 -1j];[1j 0]];

SigmaP = Sx + 1j*Sy;
SigmaM = Sx - 1j*Sy;


%Create Lindbladian in Liouville Form
L = -1j*TnProduct(Ham,eye(dim^2)) + 1j*TnProduct(eye(dim^2), transpose(Ham)) + (0.3)^2*(TnProduct(SigmaM, eye(dim), transpose(SigmaP),eye(dim)) ...
    - 1/2*(TnProduct(SigmaP*SigmaM,eye(dim^3)) + TnProduct(eye(dim^2),transpose(SigmaP*SigmaM), eye(dim))));

%Compute system-environment map in Liouville picture
Lam = expm(L*delT);


%Transform to Choi Rep by correctly multiplying with MaxEnt from left and
%right
Choi = TnProduct(eye(dim^4),MaxEntSt)*TnProduct(eye(dim^2),Lam,eye(dim^2))*TnProduct(transpose(MaxEntSt),eye(dim^4));


%Link the Chois of the resulting maps to build a comb

%Contract first Choi with random environment state
UseSt = Rand_state(dim);
rho = TnProduct(eye(dim),UseSt,eye(dim^2));
Comb = TrX(Choi*rho,[2],[dim,dim,dim,dim]);

for i = 1:N
    %Link by padding out and permuting the new map correctly
    Comb = TnProduct(Comb, eye(dim^3));
    NewChoi = syspermute(Choi,[2,1,3,4], [dim, dim, dim, dim]);
    NewChoi = TnProduct(eye(dim^(2*i)),NewChoi);
    
    dims = dim*ones([1,2*i+4]); %dimensions of subsystems
    
    %Multiply, partially transpose and trace out correct space (i.e.,
    %perform link product:
    Comb = TrX(PartTr(Comb, [2*i+1], dims)*NewChoi, [2*i+1],dims);
end

%trace out final environment
dims = dim*ones([1,2*N+3]);
Comb = TrX(Comb,[2*N+3],dims);

ResComb = Comb;
UsedInit = UseSt;
end
