function dAdt = COCO_arbitrary_grid_ODE(Ahat,p,data)
%ARBITRARY_GRID_ODE Performs time integration on an arbitrary system of
%arches

%%
modes_to_skip   = data.modes_to_skip;
C               = data.constraint_count;
coeff_matrix    = data.coeff_matrix;
N               = data.N;
N_modes         = data.N_modes;

% Load the parameters
b = p(1,:);
t = p(2,:);

data.b_vector = b;
data.t_vector = t;




% Set up empty vector used to describe the system
dAdt = zeros(size(Ahat));

% Relate first and second derivative
dAdt(1:N*N_modes-C,:) = Ahat(N*(N_modes)-C+1:end,:);

% Recover the lost variables
A = determine_A_from_Ahat(Ahat',data)';

% Use dVdANvec function to solve for RHS of 
dVdaNvec = COCO_arbitrary_grid_dVdaN(A,data);

% Construct the vector composing of dv/daN
RHS = zeros(length(coeff_matrix),width(A));

RHS(1:N_modes*N,:) = dVdaNvec;

%Add in damping to RHS
beta = 0.00001;
RHS = RHS + [beta.*A(N_modes*(N)+1:end,:) ; zeros(height(RHS) - N*(N_modes),width(A))];

% Solve the system
solution = -(inv(coeff_matrix))*RHS;

% Move the needed rows into the solution matrix
keep_rows = setdiff(1:N*N_modes,modes_to_skip);
dAdt(N*N_modes-C+1:end,:) = solution(keep_rows,:);

end