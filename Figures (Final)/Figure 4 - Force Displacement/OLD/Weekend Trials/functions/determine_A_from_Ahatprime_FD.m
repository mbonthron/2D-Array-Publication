function [A] = determine_A_from_Ahatprime_FD(Ahatprime,data,t)
%DETERMINE_A0_FROM_A0 Summary of this function goes here
%   Calculates the total first order from of a
%% Load from data
coeff_matrix    = data.coeff_matrix_modes;
N               = data.N;
N_modes         = data.N_modes;
modes_to_skip   = data.modes_to_skip;
C               = data.constraint_count;

eta   = data.eta;
alpha = data.alpha;

%% Solve for a1 and a1dot from imposed displacement 
sin_terms = sin((2:N_modes) * eta * pi)*Ahatprime(1:N_modes-1,:);
a1 = (1/sin(eta*pi))*(data.initial_height - alpha*t - sin_terms);

sin_terms = sin((2:N_modes) * eta * pi)*Ahatprime((end/2+1):(end/2+1+N_modes-2),:);
a1dot = (1/sin(eta*pi))*(-t - sin_terms);

Ahat = [a1;Ahatprime(1:end/2,:);a1dot;Ahatprime((end/2+1):end,:)];

%% Recover the lost variables
last_C_rows = coeff_matrix(end-(C-1):end,1:N*N_modes); 
LHS = last_C_rows(:,modes_to_skip);
RHS = -1*last_C_rows(:,data.setdiffmodes);

% Ahat = Ahat';

missingvals = (LHS\RHS)*Ahat(1:end/2,:);
Dmissingvals = (LHS\RHS)*Ahat(end/2+1:end,:);

% Produce the 'full' Ahat matrix which can be used in dVdAN
A = Ahat;
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
A_new(insert_locsD, :) = Dmissingvals;
A_new(~(insert_locs| insert_locsD), :) = A;

A = A_new;
end