function dAdt = single_arch_ODE(t,Ahat,data)
%ARBITRARY_GRID_ODE Performs time integration on an arbitrary system of
%arches
%   INPUTS
%   ===================================================
%   t:  time
%
%   A:  Vector describing the state space of the system
%
%   data
%
%   OUTPUTS
%   ===================================================
%   dAdt: Since the system is in first order form, partial derivative with
%       respect to time of each variable

%% Load in values
N_modes         = data.N_modes;
coeff_matrix    = data.coeff_matrix;
beta            = data.beta;

% Imposed Displacement
% Set up empty vector used to describe the system
dAdt = zeros(2*(N_modes-1),1);

% Relate first and second derivative
dAdt(1:N_modes-1,:) = Ahat(N_modes:end,:);

% Recover A from Ahat
A = determine_A_from_Ahat(t,Ahat,data);


%% Use A to solve the full system
% Construct the vector composing of dv/daN
RHS = zeros(N_modes+1,1);
dVdaNvec = single_arch_dVdaN(A,data);

RHS(1:N_modes) = dVdaNvec;

%Add in damping to RHS
RHS = RHS + [beta.*A(N_modes+1:end) ; 0];

solution = -(inv(coeff_matrix))*RHS;

% Move the needed rows into the solution matrix
dAdt(N_modes:end,:) = solution(2:end-1);

end