rhombus_count = 2;

yy = sqrt(3)/2;
points = [0 0 ; 1 0 ; 0.5 yy ; 1.5 yy];

for i = 2:rhombus_count
    point1 = [i 0];
    point2 = [i+.5 yy];

    points = [points; point1 ; point2];
end