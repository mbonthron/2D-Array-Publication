function dAdt = arbitrary_grid_ODE(t,A,coeff_matrix,b_vector,e_vector,t_vector,displacement_vector,force_vector,eta_vector,omega_vector,beta,N,N_modes)
%ARBITRARY_GRID_ODE Performs time integration on an arbitrary system of
%arches
%   INPUTS
%   ===================================================
%   t:  time
%
%   A:  Vector describing the state space of the system
%
%   coeff_matrix: Coefficient matrix describing how the arches are
%       connected to one another
%   b_vector: Vector describing the b value of each arch from elastic
%       buckling
%   e_vector: Vector describing the e value associated with plastic
%       deformation
%   t_vector: Vector describing the thickness of each arch
%
%   displacement_vector: Vector which describes where displacement is
%       imposed on each arch. 0 if no displacement is imposed
%       Positive values correspond to a downward motion (local coordinate)
%       Negative values correspond to an upward motion (local coordinate)
%   force_vector: Vector describing the magnitude of force on each arch
%
%   eta_vector: Vector describing the location the force is applied on each
%       arch
%   omega_vector: Frequency of excitation of forces and displacements
%       corresponding to each arch
%       If omega_vector = 0, but the corresponding force is nonzero, assume
%       a constant force in the direction
%   beta: Modal damping coefficient
%
%   N: Number of arches in the system
%
%   N_modes: Number of modes used to describe each arch
%
%   OUTPUTS
%   ===================================================
%   dAdt: Since the system is in first order form, partial derivative with
%       respect to time of each variable

% Set up empty vector used to describe the system
dAdt = zeros(2*N*N_modes,1);

% Relate first and second derivative
dAdt(1:N*N_modes) = A(N*N_modes+1:end);

% Construct the vector composing of dv/daN
RHS = zeros(length(coeff_matrix),1);
dVdaNvec = arbitrary_grid_dVdaN_integration(A,N,b_vector,e_vector,t_vector,N_modes);

RHS(1:N_modes*N) = dVdaNvec;

%Add in damping to RHS
RHS = RHS + [beta.*A(N_modes*N+1:end) ; zeros(height(RHS) - N_modes *N,1)];

% Add force to RHS if there are any forces
imposed_force_at = find(force_vector(:,1) ~=0);
if ~isempty(imposed_force_at)
    for i = imposed_force_at
        % Add the forces to the corresponding part of RHS at the
        % corresponding eta values
        for j = 1:N_modes
            Q = force_vector(i);
            omega = omega_vector(i);
            eta   = eta_vector(i);
            
            if omega == 0
                % Constant force in the direction
                RHS(N_modes*i-(N_modes-j)) = RHS(N_modes*i-(N_modes-j)) + Q*sin(j*eta*pi);
            else
                % Oscillatory force in the direction
                RHS(N_modes*i-(N_modes-j)) = RHS(N_modes*i-(N_modes-j)) + Q*sin(omega*t)*sin(j*eta*pi);
            end

        end
    end
end

% Include term relating first and second derivatives for imposed disp
displacement_count = 0;
imposed_displacement_at = find(displacement_vector ~= 0);

if ~isempty(imposed_displacement_at)
    for i = imposed_displacement_at
        %Since sin is odd, will be negative when user inputs negative
        %number
        eta = displacement_vector(i);
        initial_height = e_vector(i)*sin(pi*eta) + b_vector(i)*sin(pi*eta);
        omega = omega_vector(i);
    
        RHS(end - length(imposed_displacement_at) + 1 + displacement_count) = ...
                initial_height*omega^2*cos(omega*t);

        % if eta > 0
        %     % Downward motion in local coordinates
        %     RHS(end - length(imposed_displacement_at) + 1 + displacement_count) = ...
        %         initial_height*omega^2*cos(omega*t);
        % else
        %     % Upward motion in local coordiante
        %     RHS(end - length(imposed_displacement_at) + 1 + displacement_count) = ...
        %         initial_height*omega^2*cos(omega*t);
        % end


        
        displacement_count = displacement_count + 1;
    end
end


solution = -(inv(coeff_matrix))*RHS;
dAdt(N_modes*N+1:end) = solution(1:N_modes*N);
end
