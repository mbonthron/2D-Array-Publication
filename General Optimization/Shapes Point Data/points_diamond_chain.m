yy = sqrt(3)/2;

% 2 Super Cells
points = [0    0;
         yy  0.5;
         yy -0.5;
       2*yy    0;
       3*yy 0.5;
       3*yy -0.5;];

% 4 Super Cells
% points = [0    0;
%          yy  0.5;
%          yy -0.5;
%        2*yy    0;
%        3*yy 0.5;
%        3*yy -0.5;
%        4*yy    0;
%        5*yy 0.5;
%        5*yy -0.5;
%        6*yy    0;
%        7*yy  0.5;
%        7*yy -0.5];

data.points = points;
data.V = length(points);