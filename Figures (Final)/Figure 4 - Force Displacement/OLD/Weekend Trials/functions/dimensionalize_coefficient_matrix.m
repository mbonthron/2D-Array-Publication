function [data] = dimensionalize_coefficient_matrix(data)
%DIMENSIONALIZE_COEFFICIENT_MATRIX Summary of this function goes here
%   Detailed explanation goes here
coefficient_matrix = data.coeff_matrix;

rho = data.rho;
AA  = data.AA;
eta = data.eta;

L   = data.L;

N       = data.N;
N_modes = data.N_modes;

% Add in the modal mass:
coefficient_matrix(1:N*N_modes,1:N*N_modes) = rho*AA*eye(N*N_modes,N*N_modes);

% Add in the terms for the forces
force_vector = 2*sin(eta*pi*(1:N_modes))/L';

coefficient_matrix(1:N_modes,end) = force_vector;

data.coeff_matrix = coefficient_matrix;
end

