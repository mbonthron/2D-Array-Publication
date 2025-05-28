square_count_m = 3; % Number of Horizontal Squares
square_count_n = 3; % Number of Vertical Squares

points = [0 0;
          1 0;
          0 1;
          1 1];

% Make the Horizontal Boxes
for i = 2:square_count_m
    point1 = [i 0];
    point2 = [i 1];

    points = [points; point1 ; point2];
end

% Add the Vertical Boxes
for i = 2:square_count_n
    yval = i;
    xvals = 0:square_count_m;
    
    points_to_add = [xvals' yval*ones(square_count_m+1,1)];

    points = [points; points_to_add];
end

data.points = points;
data.V = length(points);