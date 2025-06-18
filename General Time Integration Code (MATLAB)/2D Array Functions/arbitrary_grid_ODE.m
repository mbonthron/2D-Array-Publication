function dAdt = arbitrary_grid_ODE(t,Ahat,data)
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
coeff_matrix = data.coeff_matrix;
coeff_matrix_modes = data.coeff_matrix_modes;

%b_vector = data.b_vector;
%e_vector = data.e_vector;
%t_vector = data.t_vector;

impose_displacement_at = data.impose_displacement_at;
displacement_omega = data.displacement_omega;

force_vector = data.force_magnitude;
force_eta = data.force_eta;
force_omega = data.force_omega;

impose_rotation_at = data.impose_rotation_at;
rotation_omega = data.rotation_omega;
rotation_mag = data.rotation_mag;

beta = data.beta;
N = data.N;
%V = data.V;
N_modes = data.N_modes;

C  = data.constraint_count;
modes_to_skip = data.modes_to_skip;
%%
% Set up empty vector used to describe the system
dAdt = zeros(2*(N*N_modes-C),1);

% Relate first and second derivative
dAdt(1:N*N_modes-C,:) = Ahat(N*(N_modes)-C+1:end,:);

%%
% Recover the lost variable
last_C_rows = coeff_matrix_modes(end-(C-1):end,1:N*N_modes); 
clear coeff_matrix_modes;
LHS = last_C_rows(:,modes_to_skip);
RHS = -1*last_C_rows(:,setdiff(1:N*N_modes,modes_to_skip));
clear last_C_rows;

missingvals = (LHS\RHS)*Ahat(1:end/2,:);
Dmissingvals = (LHS\RHS)*Ahat(end/2+1:end,:);
clear LHS;

% Produce the 'full' A matrix which can be used in dVdAN
A = Ahat;
clear Ahat;
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
clear missingvals;
A_new(insert_locsD, :) = Dmissingvals;
clear Dmissingvals;
A_new(~(insert_locs| insert_locsD), :) = A;
A = A_new;
clear A_new insert_locs insert_locsD;



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
clear sin_modes_eta_pi;
%% Include term relating first and second derivatives for imposed disp
imposed_displacement_at = find(impose_displacement_at ~= 0);
displacement_count = length(imposed_displacement_at);
if ~isempty(imposed_displacement_at)
    omega_vals = displacement_omega(imposed_displacement_at);
    RHS(N*N_modes+2*N+1:N*N_modes+2*N+displacement_count) = data.initial_height.*omega_vals.^2.*cos(omega_vals*t);
   
end

%% Include term relating first and second derivates for imposed rotation
imposed_rotation_at = find(impose_rotation_at ~= 0);
if ~isempty(imposed_rotation_at)
    rotation_count = length(imposed_rotation_at);
    omega_vals = rotation_omega(imposed_rotation_at);
    mag_vals = rotation_mag(imposed_rotation_at);

    RHS(N*N_modes+2*N+displacement_count+1:N*N_modes+2*N+displacement_count+rotation_count) = mag_vals.*data.initial_angle.*omega_vals.^2.*cos(omega_vals*t);
end

%% Final Inversion
solution = -(inv(coeff_matrix))*RHS;

% Move the needed rows into the solution matrix
keep_rows = setdiff(1:N*N_modes,modes_to_skip);
dAdt(N*N_modes-C+1:end,:) = solution(keep_rows);
clear keep_rows coeff_matrix;
end