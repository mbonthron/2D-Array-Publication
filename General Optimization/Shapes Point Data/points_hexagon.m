xx = sqrt(3)/2;
number_of_elements = 2;

points = [0 2;
          0 1;
          0 0;
         xx 1.5;
         xx 0.5;
        2*xx 2;
        2*xx 1;
        2*xx 0
        3*xx 1.5;
        3*xx 0.5];

%% Add in the points for the number of cells
new_cell = points;
for i = 2:number_of_elements
    shift = [4*(i-1)*xx*ones(10,1) zeros(10,1)];
    points = [points ; new_cell + shift];
end

%% Add in the final elements
endx = points(end,1);
too_add = [endx+xx 2;
           endx+xx 1;
           endx+xx 0];

points = [points;too_add];

data.points = points;
data.V = length(points);