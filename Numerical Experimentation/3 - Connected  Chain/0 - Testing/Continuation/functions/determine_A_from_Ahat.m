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
for i = 1:C
    mode = modes_to_skip(i);
    A = [A(1:mode-1,:); missingvals(i,:) ; A(mode:end,:)];
end
shift_modes = N*N_modes;    % Do the same for the derivative terms
for i = 1:C
    mode = modes_to_skip(i);
    A = [A(1:shift_modes+mode-1,:); Dmissingvals(i,:) ; A(shift_modes+mode:end,:)];
end

A = A';
end