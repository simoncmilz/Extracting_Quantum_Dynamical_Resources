function x = Rand_state(n)
%Returns a random (according to Haar measure) pure state of size n.

FiducState = zeros(n);
FiducState(1,1) = 1;
U = randU(n);	%Sample random unitary
initial = mtimes(mtimes(U,FiducState),U');

x=initial;
end
