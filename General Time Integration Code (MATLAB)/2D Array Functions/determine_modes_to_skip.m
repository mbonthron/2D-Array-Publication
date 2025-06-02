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
coeff_matrix_modes = data.coeff_matrix_modes;
C = data.constraint_count;
N = data.N;
N_modes = data.N_modes;

%%
last_C_rows = coeff_matrix_modes(end-(C-1):end,1:N*N_modes); 

[~,modes_to_skip] = rref(last_C_rows(:,3:end));
data.modes_to_skip = modes_to_skip + 2;

end

