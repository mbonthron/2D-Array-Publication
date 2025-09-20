xx = cos(pi/6);
yy = sin(pi/6);
points = [-xx yy;
          0 0;
          xx yy;
          0 -1];

data.N_modes = 3;

data.points = points;

data = determine_adjacency_matrix(data);

data = determine_coefficient_matrix(data);
data = determine_modes_to_skip(data);

data.plot_grids = true;
plot_grid(data,true)