square_count = 5;

points = [0 0;
          1 0;
          0 1;
          1 1];

for i = 2:square_count
    point1 = [i 0];
    point2 = [i 1];

    points = [points; point1 ; point2];
end