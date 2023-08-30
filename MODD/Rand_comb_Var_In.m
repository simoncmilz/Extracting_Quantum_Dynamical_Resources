function [ResComb, Lindblad]= Rand_comb_Var_In(N, Hamil, delT)

%Function that returns a comb on a system of dimension dim with N slots.
%Comb is created by taking the given Hamiltonian Hamil on dim^2 and, together
%with the dissipator of arXiv:2110.02613, propagate it for time delT and slot it all
%together via a link product. The initial and final environment line are
%deliberately left open

%system dimension (could become additional variable in the future
dim = 2;

%Create unnormalized maxEnt state vector on system/environment
%(required to switch between representations)
MaxEntSt = zeros([1,dim^4]);
for i=1:dim^2
    MaxEntSt((i-1)*dim^2 + i) = 1;
end

%normalize Hamiltonian
Ham = Hamil/norm(Hamil);

%matrices for the dissipator
Sx = [[0 1];[1 0]];
Sy = [[0 -1j];[1j 0]];

SigmaP = Sx + 1j*Sy;
SigmaM = Sx - 1j*Sy;


%Create Lindbladian in Liouville Form
L = -1j*TnProduct(Ham,eye(dim^2)) + 1j*TnProduct(eye(dim^2), transpose(Ham)) + (0.1^2)*(TnProduct(SigmaP, eye(dim), transpose(SigmaM),eye(dim)) ...
    - 1/2*(TnProduct(SigmaM*SigmaP,eye(dim^3)) + TnProduct(eye(dim^2),transpose(SigmaM*SigmaP), eye(dim))));

%Compute system-environment map in Liouville picture
Lam = expm(L*delT);


%Transform to Choi Rep by correctly multiplying with MaxEnt from left and
%right
Choi = TnProduct(eye(dim^4), MaxEntSt)*TnProduct(eye(dim^2),Lam,eye(dim^2))*TnProduct(transpose(MaxEntSt), eye(dim^4));


%Link the Chois of the resulting maps to build a comb

%First chunk of the comb
Comb = Choi;


for i = 1:N
    %Link by padding out and permuting the new map correctly
    Comb = TnProduct(Comb, eye(dim^3));
    NewChoi = syspermute(Choi,[2,1,3,4], [dim, dim, dim, dim]);
    NewChoi = TnProduct(eye(dim^(2*i+1)),NewChoi);
    
    dims = dim*ones([1,2*i+5]); %dimensions of subsystems
    
    %Multiply, partially transpose and trace out correct space (i.e.,
    %perform link product:
    Comb = TrX(PartTr(Comb, [2*i+2], dims)*NewChoi, [2*i+2],dims);
end

%Return Comb with open initial and final environment line
ResComb = Comb; 
Lindblad = L;

end
