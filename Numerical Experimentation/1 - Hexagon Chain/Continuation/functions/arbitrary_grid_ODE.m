function dAdt = arbitrary_grid_ODE(A,p,coeff_matrix,N,N_modes,modes_to_skip)
%ARBITRARY_GRID_ODE Performs time integration on an arbitrary system of
%arches
%   INPUTS
%   ===================================================
%   A:  2*(3*N_modes - 3) state space describing the system
%
%   p:  p(1) = b, rise of the arch 
%       p(2) = t, thickness of the arches   
%
%   coeff_matrix: Coefficient matrix describing the system
%
%   N: Number of arches in the system
%
%   N_modes: Number of modes used to describe each arch
%   
%   modes_to_skip: Modes removed from Ahat to ensure linear independence
%       of the state variables
%
%   OUTPUTS
%   ===================================================
%   dAdt: Since the system is in first order form, partial derivative with
%       respect to time of each variable
%%
% Load the parameters
b = p(1,:);
t = p(2,:);

% Constraint Count
C = length(modes_to_skip);

% Set up empty vector used to describe the system
dAdt = zeros(size(A));

% Relate first and second derivative
dAdt(1:N*N_modes-C,:) = A(N*(N_modes)-C+1:end,:);

% Recover the lost variables
last_C_rows = coeff_matrix(end-(C-1):end,1:N*N_modes); 
LHS = last_C_rows(:,modes_to_skip);
RHS = -1*last_C_rows(:,setdiff(1:N*N_modes,modes_to_skip));

missingvals = (LHS\RHS)*A(1:end/2,:);
Dmissingvals = (LHS\RHS)*A(end/2+1:end,:);

% Produce the 'full' Ahat matrix which can be used in dVdAN
Ahat = A;
for i = 1:C
    mode = modes_to_skip(i);
    Ahat = [Ahat(1:mode-1,:); missingvals(i,:) ; Ahat(mode:end,:)];
end
shift_modes = N*N_modes;    % Do the same for the derivative terms
for i = 1:C
    mode = modes_to_skip(i);
    Ahat = [Ahat(1:shift_modes+mode-1,:); Dmissingvals(i,:) ; Ahat(shift_modes+mode:end,:)];
end

% Use dVdANvec function to solve for RHS of 
dVdaNvec = arbitrary_grid_dVdaN(Ahat,N,b,t,N_modes);

% Construct the vector composing of dv/daN
RHS = zeros(length(coeff_matrix),width(A));

RHS(1:N_modes*N,:) = dVdaNvec;

%Add in damping to RHS
beta = 0.00001;
RHS = RHS + [beta.*Ahat(N_modes*(N)+1:end,:) ; zeros(height(RHS) - N*(N_modes),width(A))];

% Solve the system
solution = -(inv(coeff_matrix))*RHS;

% Move the needed rows into the solution matrix
keep_rows = setdiff(1:N*N_modes,modes_to_skip);
dAdt(N*N_modes-C+1:end,:) = solution(keep_rows,:);

end