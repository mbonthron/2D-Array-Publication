function [solution] = determine_M_Q(t,A,data)
%DETERMINE_M_Q Summary of this function goes here
%   Detailed explanation goes here

beta = data.beta;
N = data.N;
Q = data.force_magnitude;

%% Construct the vector composing of dv/dan
RHS = zeros(height(data.coeff_matrix),length(t));
dVdaNvec = arbitrary_grid_dVdaN(A, data);

RHS(1:3*N,:) = dVdaNvec';

% Add in damping to RHS
RHS = RHS + [beta.*A(:,3*N+1:end)' ; zeros(height(RHS) - 3*N,length(t)) ] ;


% Include term relating first and second derivatives for imposed disp
displacement_count = 0;
for i = 1:length(data.impose_displacement_at)
    RHS(end - length(data.impose_displacement_at) + 1 + displacement_count,:) = data.b_vector(i)*data.displacement_omega^2*cos(data.displacement_omega*t);
    displacement_count = displacement_count + 1;
end

%%
solution = -(inv(data.coeff_matrix))*RHS;
end
