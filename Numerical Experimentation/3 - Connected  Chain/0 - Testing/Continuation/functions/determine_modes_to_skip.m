function [data] = determine_modes_to_skip(data)
%DETERMINE_MODES_TO_SKIP From the coefficient matrix, determines a set of
%   modes to ignore which leaves the problem linearly dependent
%   INPUTS
%   ===================================================
%   coeff_matrix: coefficient matrix of the system
%   C: Number of constraint equations (number of modes to remove)
%   N: Number of arches
%   N_modes: Number of modes per arch
%
%   OUTPUTS
%   ===================================================
%   modes_to_skip: a vector which contains the index of modes to skip to be
%       left with a linearly independent system
%%
N = data.N;
N_modes = data.N_modes;
C = 2*N-data.V; %double check line
coeff_matrix = data.coeff_matrix;
C_rows = coeff_matrix(N_modes*N+data.V+1:N_modes*N+data.V+C,1:N*N_modes); 

[~,modes_to_skip] = rref(C_rows(:,3:end));
data.modes_to_skip = modes_to_skip + 2;

end