xx = cos(pi/6);
yy = sin(pi/6);
points = [0 0;
          0 -1;
          -xx yy;
          xx yy;
          -xx -1-yy
          xx -1-yy];

theta = -pi/2;
R = [cos(theta) -sin(theta);
     sin(theta) cos(theta)];

points = points*R;

data.N_modes = 5;
data.V = size(points,1);

data.points = points;

data = determine_adjacency_matrix(data);
data = determine_coefficient_matrix(data);
data = determine_modes_to_skip(data);

data.plot_grids = true;
plot_grid(data,true)