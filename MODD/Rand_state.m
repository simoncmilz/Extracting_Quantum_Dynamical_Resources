function x = Rand_state(n)
% Computes a random superchannel by means of an SDP. 
% To create different Superchannels for different calls, the SDP minimizes
% the distance of the superchannel to a randomly chosen pure (8x8) state,
% drawn according to a Haar measure.

%Random state (used for randomization of the superchannel
FiducState = zeros(n);
FiducState(1,1) = 1;
U = randU(n);
initial = mtimes(mtimes(U,FiducState),U');

x=initial;
end