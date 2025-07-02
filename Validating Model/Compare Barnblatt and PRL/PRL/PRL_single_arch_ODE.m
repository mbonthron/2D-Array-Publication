function dAdt = PRL_single_arch_ODE(t,Ahat,data)
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
N_modes         = 3;
coeff_matrix    = [0 1 ; 1 0];
beta            = 1.5;

%% Final Inversion
if data.imposed
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
elseif ~ data.imposed
    % Free - NOTE AHAT NOW REFERS TO A
    A = Ahat;

    % Set up empty vector used to describe the system
    dAdt = zeros(2*N_modes,1);
    
    % Relate first and second derivative
    dAdt(1:N_modes,:) = A(N_modes+1:end,:);
    
    %% Use A to solve the full system
    % Construct the vector composing of dv/daN
    dVdaNvec = single_arch_dVdaN(A,data);
    
    RHS = dVdaNvec;
    
    %Add in damping to RHS
    RHS = RHS + [beta.*A(N_modes+1:end)];

    % No inversion needed for RHS
    dAdt(N_modes+1:end,:) = -RHS;
end

end