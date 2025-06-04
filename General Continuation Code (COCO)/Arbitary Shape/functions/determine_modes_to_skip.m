function [modes_to_skip] = determine_modes_to_skip(coeff_matrix,C,N,N_modes)
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
last_C_rows = coeff_matrix(end-(C-1):end,1:N*N_modes); 

[~,modes_to_skip] = rref(last_C_rows(:,3:end));
modes_to_skip = modes_to_skip + 2;

end

