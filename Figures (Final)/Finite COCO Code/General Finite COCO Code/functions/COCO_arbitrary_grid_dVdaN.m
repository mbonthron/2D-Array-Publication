function [dVdt] = COCO_arbitrary_grid_dVdaN(A,data)
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
N = data.N;
b_vector = data.b_vector;
e_vector = data.e_vector;
t_vector = data.t_vector;
N_modes = data.N_modes;

[m,n] = size(A);
dVdt = zeros(N_modes*N,n);


if n == 1
    % If only given a single input
    for i = 1:N
        indices = (N_modes*i-(N_modes-1)):(N_modes*i);
        % Elastic Arch
        deltaL = (b_vector/2).^2;
        sum1 = sum(((1:1:N_modes).^2.*A(indices).'.^2).');
        
        % Write the equation for the first mode
        j = 1;
        dVdt(N_modes*i-(N_modes-j)) = (1/12)*(t_vector.^2)*(j^4).*A(indices(j)) - (j^2).*(deltaL   - 0.25.*(sum1)).*(A(indices(j)));

        % Write the equation for all subsequent modes
        for j = 2:N_modes
            dVdt(N_modes*i-(N_modes-j)) = (1/12)*(t_vector.^2)*(j^4).*A(indices(j)) - (j^2).*(deltaL   - 0.25.*(sum1)).*A(indices(j));

        end

    end
else
    % Parameterized input
    dVdt = zeros(N_modes*N,n);
    for i = 1:N
        % Elastic Arch
        indices = (N_modes*i-(N_modes-1)):(N_modes*i);
        % Elastic Arch
        deltaL = (b_vector/2).^2;
        sum1 = sum(((1:1:N_modes).^2.*A(indices,:).'.^2).');

        for j = 1:N_modes
            dVdt(N_modes*i-(N_modes-j),:) = (1/12)*(t_vector.^2)*(j^4).*A(indices(j),:) - (j^2).*(deltaL  - 0.25.*(sum1)).*A(indices(j),:);
        end
    end
end