function [solution] = determine_M_Q(t,A,data)
%DETERMINE_M_Q Summary of this function goes here
%   Detailed explanation goes here

beta = data.beta;
N = data.N;
N_modes = data.N_modes;
%% Construct the vector composing of dv/dan
RHS = zeros(height(data.coeff_matrix),length(t));
dVdaNvec = arbitrary_grid_dVdaN(A, data);

RHS(1:N_modes*N,:) = dVdaNvec';

% Add in damping to RHS
RHS = RHS + [beta.*A(:,N_modes*N+1:end)' ; zeros(height(RHS) - N_modes*N,length(t)) ] ;


% Include term relating first and second derivatives for imposed disp
displacement_count = 0;
for i = 1:length(data.impose_displacement_at)
    RHS(end - length(data.impose_displacement_at) + 1 + displacement_count,:) = 0;
    displacement_count = displacement_count + 1;
end

%%
solution = -(inv(data.coeff_matrix))*RHS;
end
