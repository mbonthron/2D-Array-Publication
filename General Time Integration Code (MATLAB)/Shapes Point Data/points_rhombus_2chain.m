rhombus_count = 2;

xx = sqrt(3)/2;
points = [0 0 ; 0 1 ; xx 0.5 ; xx 1.5];

for i = 2:rhombus_count
    if mod(i,2) == 0
        % "down" rhombus
        point1 = [i*xx 0];
        point2 = [i*xx 1];
    else
        % "up" rhombus
        point1 = [i*xx 0.5];
        point2 = [i*xx 1.5];
    end
    points = [points; point1 ; point2];
end

data.points = points;
data.V = length(points);