function [Ahat] = determine_Ahat_from_A(A,data)
%DETERMINE_A0_FROM_A0 Summary of this function goes here
%   Calculates the total first order from of a
%% Load from data
N = data.N;
N_modes = data.N_modes;
modes_to_skip = data.modes_to_skip;


preseve_rows = setdiff(1:2*N*N_modes,[modes_to_skip (N*N_modes+modes_to_skip)]);

Ahat = A(preseve_rows,:);

end