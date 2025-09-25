t= 1;
A= data.A0;

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

% Relate first and second derivative
dAdt(1:N*N_modes) = A(N*N_modes+1:end);

% Construct the vector composing of dv/daN
RHS = zeros(length(coeff_matrix),1);
dVdaNvec = arbitrary_grid_dVdaN(A,data);

RHS(1:N_modes*N) = dVdaNvec;

%Add in damping to RHS
RHS = RHS + [beta.*A(N_modes*N+1:end) ; zeros(size(RHS,1) - N_modes *N,1)];

%% Add force to RHS if there are any forces
%  05/24/2025 - Can we vectorize this code to add all forces in
%  simultaneously without a 'for' loop?
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
    
    
    
%     for i = 1:length(imposed_force_at)
%         % Add the forces to the corresponding part of RHS at the
%         % corresponding eta values
%         Q = force_vector(i);
%         omega = force_omega(i);
%         eta   = force_eta(i);
%         for j = 1:N_modes
% 
%             
%             if omega == 0
%                 % Constant force in the direction
%                 RHS(N_modes*i-(N_modes-j)) = RHS(N_modes*i-(N_modes-j)) + Q*sin(j*eta*pi);
%             else
%                 % Oscillatory force in the direction
%                 RHS(N_modes*i-(N_modes-j)) = RHS(N_modes*i-(N_modes-j)) + Q*sin(omega*t)*sin(j*eta*pi);
%             end
% 
%         end
%     end
end