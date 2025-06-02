yy = sqrt(3)/2;
number_of_elements = 10;

points0 = [0 0;
         0.5 yy;
         1 0;
         1.5 yy;
         2 0;
         2.5 yy;
          3 0;
          3.5 yy
          4 0
          4.5 yy];

points = points0;

for i = 2:number_of_elements
    points_to_add = points0(3:end,:) + [4*(i-1)*ones(8,1) zeros(8,1)];
    points = [points ; points_to_add];
end

data.points = points;
data.V = length(points);
