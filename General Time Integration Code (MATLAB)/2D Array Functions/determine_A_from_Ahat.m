function [A] = determine_A_from_Ahat(Ahat,data)
%DETERMINE_A0_FROM_A0 Summary of this function goes here
%   Calculates the total first order from of a
%% Load from data
coeff_matrix = data.coeff_matrix_modes;
N = data.N;
N_modes = data.N_modes;
modes_to_skip = data.modes_to_skip;
C  = data.constraint_count;

% Recover the lost variables
last_C_rows = coeff_matrix(end-(C-1):end,1:N*N_modes); 
LHS = last_C_rows(:,modes_to_skip);
RHS = -1*last_C_rows(:,setdiff(1:N*N_modes,modes_to_skip));

Ahat = Ahat';

missingvals = (LHS\RHS)*Ahat(1:end/2,:);
Dmissingvals = (LHS\RHS)*Ahat(end/2+1:end,:);

% Produce the 'full' Ahat matrix which can be used in dVdAN
A = Ahat;
clear Ahat;
shift_modes = N*N_modes;
[M_size, N_size] = size(A);
C = length(modes_to_skip);
A_new = zeros(M_size + 2*C, N_size);

insert_locs = false(M_size + 2*C, 1);
insert_locs(modes_to_skip) = true;
insert_locsD = false(M_size + 2*C, 1);
insert_locsD(modes_to_skip+shift_modes) = true;

% Create final matrix
A_new(insert_locs, :) = missingvals;
clear missingvals;
A_new(insert_locsD, :) = Dmissingvals;
clear Dmissingvals;
A_new(~(insert_locs| insert_locsD), :) = A;

A = A_new';
end