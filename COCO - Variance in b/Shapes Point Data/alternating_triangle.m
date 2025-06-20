%% Load in the Points
yy = sqrt(3)/2;

points = [0 0;
         0.5  yy;
         0.5 -yy;
          1    0;
         1.5  yy;
         1.5 -yy;
         2     0
         2.5  yy;
         2.5 -yy;
         3     0;
         3.5  yy;
         3.5 -yy;];

data.points = points;
data.V = length(points);

%% Things for clean up
nodes_to_remove = [3, 5, 9, 11];
nodes_to_remove2 = [data.N_cells*8+2];

nodes_to_hold = [];
arches_to_displace = [];
nodes_to_rotate = [];
arches_to_force_positive = [];
arches_to_force_negative = [];

connections_to_remove = [];
data.shape_name = 'Alternating Triangle';