points = [0 0;
          1 -1;
          1 0;
          0 -1;
          2 -1];

data.N_modes = 3;

data.points = points;

data = determine_adjacency_matrix(data);

data.adjacency_matrix(1,4) = 0;
data.adjacency_matrix(4,1) = 0;
data.N = data.N - 1;

data = determine_coefficient_matrix(data);
data = determine_modes_to_skip(data);

data.plot_grids = true;
plot_grid(data,true)