function dAdt = arbitrary_grid_ODE(t,A,data)
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
b_vector = data.b_vector;
e_vector = data.e_vector;
t_vector = data.t_vector;

impose_displacement_at = data.impose_displacement_at;
displacement_omega = data.displacement_omega;

force_vector = data.force_magnitude;
force_eta = data.force_eta;
force_omega = data.force_omega;

impose_rotation_at = data.impose_rotation_at;
rotation_omega = data.rotation_omega;

beta = data.beta;
N = data.N;
V = data.V;
N_modes = data.N_modes;

%%
% Set up empty vector used to describe the system
dAdt = zeros(2*N*N_modes,1);

% Construct the vector composing of dv/daN
RHS = zeros(length(coeff_matrix),1);

%% Add force to RHS if there are any forces
imposed_force_at = find(force_vector(:,1) ~=0);
modes = (1:N_modes)';
if ~isempty(imposed_force_at)
    eta_vals = force_eta;
    omega_vals = force_omega;
    Q_vals = force_vector;
    
    % Precompute sine(j * eta * pi) for all eta and j
    sin_modes_eta_pi = sin(modes  .* pi .* (eta_vals));  % size: [N_modes x num_forces]
    
    % Compute time factor depending on omega
    time_factor = ones(1, length(imposed_force_at));
    nonzero_omega = omega_vals ~= 0;
    time_factor(nonzero_omega) = sin(omega_vals(nonzero_omega) * t);
    
    % Compute contributions
    contributions = Q_vals' .* time_factor .* sin_modes_eta_pi;  % [N_modes x num_forces]
    
    RHS(1:N_modes*N) = contributions(:);
end

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

    if data.rotation_type == "const"
        % If constant rotation i.e. dw/dx = angle*(1-omega*t)
        for i = 1:length(rotation_count)
            % Use the coefficient matrix to find one mode to solve for in
            % terms of the other modes
            coeff_row = coeff_matrix(end-rotation_count+i,1:N*N_modes);
            idx_to_solve_for = find(coeff_row~=0,1,"first");
           
            % f(t)
            ft = data.initial_angle(i)*(1 - omega_vals(i)*t);

            % df(t)
            dft = data.initial_angle(i)*(-omega_vals(i));

            % Solve for the index - zero derivative
            solved0 = (1/coeff_row(idx_to_solve_for))*(ft - coeff_row*A(1:N*N_modes) + coeff_row(idx_to_solve_for)*A(idx_to_solve_for));
            A(idx_to_solve_for) = solved0;

            % Solve for the index - first derivative
            solved1 = (1/coeff_row(idx_to_solve_for))*(dft - coeff_row*A(N*N_modes+1:end) + coeff_row(idx_to_solve_for)*A(N*N_modes+idx_to_solve_for));
            A(N*N_modes + idx_to_solve_for) = solved1;  

        end
    else
        % Otherwise assume the rotation 
        RHS(N*N_modes+2*N+displacement_count+1:N*N_modes+2*N+displacement_count+rotation_count) = data.initial_angle.*omega_vals.^2.*cos(omega_vals*t);
    end
end

% Relate First and Second Derivatives
dAdt(1:N*N_modes) = A(N*N_modes+1:end);

%%
dVdaNvec = arbitrary_grid_dVdaN(A,data);

RHS(1:N_modes*N) = RHS(1:N_modes*N) + dVdaNvec;

%Add in damping to RHS
RHS = RHS + [beta.*A(N_modes*N+1:end) ; zeros(size(RHS,1) - N_modes *N,1)];

%%
solution = -(inv(coeff_matrix))*RHS;
dAdt(N_modes*N+1:end) = solution(1:N_modes*N);
end