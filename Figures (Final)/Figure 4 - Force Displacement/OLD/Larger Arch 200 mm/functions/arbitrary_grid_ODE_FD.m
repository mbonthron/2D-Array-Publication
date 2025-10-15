function dAdt = arbitrary_grid_ODE_FD(t,Ahatprime,data)
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
%% Load Data
coeff_matrix        = data.coeff_matrix;
coeff_matrix_modes  = data.coeff_matrix_modes;

%b_vector = data.b_vector;
%e_vector = data.e_vector;
%t_vector = data.t_vector;

impose_displacement_at  = data.impose_displacement_at;
displacement_omega      = data.displacement_omega;

force_vector            = data.force_magnitude;
force_eta               = data.force_eta;
force_omega             = data.force_omega;

impose_rotation_at      = data.impose_rotation_at;
rotation_omega          = data.rotation_omega;
rotation_mag            = data.rotation_mag;

beta = data.beta;
N = data.N;
N_modes = data.N_modes;

C  = data.constraint_count;
modes_to_skip = data.modes_to_skip;
%%
% Set up empty vector used to describe the system
dAdt = zeros(size(Ahatprime));

% Relate first and second derivative
dAdt(1:end/2,:) = Ahatprime(end/2+1:end,:);

%%
% Recover the lost variable
A = determine_A_from_Ahatprime_FD(Ahatprime,data,t);


%% Use A to solve the full system
% Construct the vector composing of dv/daN
RHS = zeros(length(coeff_matrix),1);
dVdaNvec = arbitrary_grid_dVdaN(A,data);

RHS(1:N_modes*N) = dVdaNvec;

%Add in damping to RHS
RHS = RHS + [beta.*A(N_modes*N+1:end) ; zeros(size(RHS,1) - N_modes *N,1)];

%% Add force to RHS if there are any forces
imposed_force_at = find(force_vector(:,1) ~=0, 1);
modes = (1:N_modes)';
if ~isempty(imposed_force_at)
    eta_vals = force_eta';
    omega_vals = force_omega';
    Q_vals = force_vector;
    
    % Precompute sine(j * eta * pi) for all eta and j
    sin_modes_eta_pi = sin(modes  .* pi .* (eta_vals));  % size: [N_modes x num_forces]
    
    % Compute time factor depending on omega
    time_factor = ones(1, length(eta_vals));
    nonzero_omega = omega_vals ~= 0;
    time_factor(nonzero_omega) = sin(omega_vals(nonzero_omega) * t);
    
    % Compute contributions
    contributions = Q_vals' .* time_factor .* sin_modes_eta_pi;  % [N_modes x num_forces]
    
    RHS(1:N_modes*N) = RHS(1:N_modes*N) + contributions(:);
end
% clear sin_modes_eta_pi;
%% Include term relating first and second derivatives for imposed disp
% imposed_displacement_at = find(impose_displacement_at ~= 0);
% displacement_count = length(imposed_displacement_at);
% if ~isempty(imposed_displacement_at)
%     % omega_vals = displacement_omega(imposed_displacement_at);
%     % RHS(N*N_modes+2*N+1:N*N_modes+2*N+displacement_count) = data.initial_height.*omega_vals.^2.*cos(omega_vals*t);
%     RHS(N*N_modes+2*N+1:N*N_modes+2*N+displacement_count) = 0;
% end


%% Final Inversion
solution = -(inv(coeff_matrix))*RHS;

% Move the needed rows into the solution matrix
dAdt(end/2+1:end,:) = solution(data.keep_row,:);
% clear keep_rows coeff_matrix;
end