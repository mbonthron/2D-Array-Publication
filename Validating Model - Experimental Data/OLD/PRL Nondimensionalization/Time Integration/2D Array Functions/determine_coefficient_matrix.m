function [data] = determine_coefficient_matrix(data)
%DETERMINE_STARTING_VALS Summary of this function goes here
%   Detailed explanation goes here
%% Load values from data
N_modes = data.N_modes;
eta     = data.eta;

coeff = zeros(N_modes+1,N_modes+1);

for i = 1:N_modes
    coeff(i,i) = 1;
    coeff(i,end) = sin(i*pi*eta);
end

coeff(end,1:N_modes) = sin(eta*pi*(1:N_modes));

data.coeff_matrix = coeff;
end