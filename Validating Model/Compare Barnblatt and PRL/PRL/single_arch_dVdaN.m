function [dVdt] = single_arch_dVdaN(A,data)
%   ARBITRARY_GRID_DVDT Determines the nonlinear terms associated with zeroth derivative
%   terms
%   INPUTS
%   ===================================================
%   A: Vector (or matrix) describing the state space of the system
%   N: Number of arches in the system
%   b_vector: Vector describing the b value of each arch from elastic
%       buckling
%   e_vector: Vector describing the e value associated with plastic
%       deformation
%   t_vector: Vector describing the thickness of each arch
%   N_modes: Number of modes describing each arch
%
%   OUTPUTS
%   ===================================================
%   dVdt: Since the system is in first order form, partial derivative with
%       respect to time of each variable

%%
b = data.b;
e = data.e;
N_modes = data.N_modes;

[m,n] = size(A);

indices = 1:N_modes;


dVdt = zeros(N_modes,1);

% Elastic Arch
deltaL = 1 + (b/2)^2;
sum1 = sum(((1:1:N_modes).^2.*A(indices).'.^2).');

% Write the equation for the first mode
j = 1;
dVdt(1) = (j^4).*A(indices(j)) - (j^2).*(deltaL - 1/2*e*A(indices(1))  - 0.25.*(sum1)).*(A(indices(j))+e);

% Write the equation for all subsequent modes
for j = 2:N_modes
    dVdt(j) = (j^4).*A(indices(j)) - (j^2).*(deltaL - 1/2*e*A(indices(1))  - 0.25.*(sum1)).*A(indices(j));
end

end


