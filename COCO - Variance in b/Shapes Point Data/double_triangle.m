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
nodes_to_remove = [];
nodes_to_remove2 = (data.N_cells*12)+[2 3];

connections_to_remove = [2 5
                        5 8
                        11 8
                        11 2
                        3 6
                        6 9
                        9 12
                        12 3];
data.shape_name = 'Double Triangle';

nodes_to_hold = [2 3 data.N_cells*12-1 data.N_cells*12];
arches_to_displace = [90];
nodes_to_rotate = [];
arches_to_force_positive = [];
arches_to_force_negative = [];