yy = sqrt(3)/2;

points = [0 0;
         0.5 yy;
         1 0];

%%
data.points = points;
data.V = length(points);

data.shape_name = 'Super Simple Triangle';

data = determine_adjacency_matrix(data);

%% I don't know if I need this
data.e_vector = 0*ones(data.N,1);

data.impose_displacement_at = zeros(data.N,1);
data.impose_rotation_at = zeros(data.V,1);

data = determine_coefficient_matrix(data);
data = determine_modes_to_skip(data);

plot_system_once(zeros(data.N*data.N_modes),data);