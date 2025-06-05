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
         3.5 -yy;
         4     0;];

data.points_finite = points;
data2.points = points;

data2 = determine_adjacency_matrix(data2);
data2 = remove_connection(data2);

data.adjacency_matrix_finite = data2.adjacency_matrix;
