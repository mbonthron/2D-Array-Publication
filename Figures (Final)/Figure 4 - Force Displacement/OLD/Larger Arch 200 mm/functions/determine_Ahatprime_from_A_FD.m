function [Ahatprime] = determine_Ahatprime_from_A_FD(A,data)
%DETERMINE_A0_FROM_A0 Summary of this function goes here
%   Calculates the total first order from of a
%% Load from data
N = data.N;
N_modes = data.N_modes;
modes_to_skip = data.modes_to_skip;


preserve_rows = setdiff(1:2*N*N_modes,[modes_to_skip (N*N_modes+modes_to_skip)]);

Ahat = A(preserve_rows,:);

Ahatprime = Ahat;
midIdx = length(Ahat)/2;

Ahatprime([1 midIdx]) = [];

end