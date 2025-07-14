function [Ahat] = determine_Ahat_from_A(A,data)
%DETERMINE_A0_FROM_A0 Summary of this function goes here
%   Calculates the total first order from of a
%% Load from data
N_modes = data.N_modes;

preserve_rows = [2:N_modes N_modes+2:2*N_modes]';

Ahat = A(preserve_rows,:);

end