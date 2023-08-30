function x = MaxEnt(dim)

%Program that returns an unnormalized maximally entangled state on the copy of a Hilbert
%space of dimension dim

vec = zeros(1,dim^2);

for i = 1:dim
    vec(dim*(i-1)+i) = 1;
end
vec;
x = TnProduct(vec,transpose(vec));
end