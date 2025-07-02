function [Q] = determine_force(A,data)
%%
N_modes         = data.N_modes;
coeff_matrix    = data.coeff_matrix;
beta            = data.beta;

[m,n] = size(A);

%% Use A to solve the full system
% Construct the vector composing of dv/daN
RHS = zeros(N_modes+1,m);
dVdaNvec = single_arch_dVdaN(A,data);

%
RHS(1:N_modes,:) = dVdaNvec;

%Add in damping to RHS
RHS = RHS + [beta.*A(:,N_modes+1:end)' ; zeros(1,m)];

%% Final Inversion
solution = -(inv(coeff_matrix))*RHS;

Q = solution(end,:)';

end

