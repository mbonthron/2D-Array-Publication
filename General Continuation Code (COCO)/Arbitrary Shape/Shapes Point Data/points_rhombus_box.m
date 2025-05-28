horizontal_count = 3;
vertical_count = 3;

yy = sqrt(3)/2;

row1x = 0:1:horizontal_count;
row1y = zeros(size(row1x));

row2x = 0.5:1:horizontal_count+0.5;
row2y = yy*ones(size(row2x));

xs = [row1x' ; row2x'];
ys = [row1y' ; row2y'];

points = [xs ys];

for i = 2:vertical_count
    rowix_start = 0.5*i;
    rowix_end   = horizontal_count + 0.5*i;

    rowix = rowix_start:1:rowix_end;
    rowiy = yy*i*ones(size(rowix));
    
    add_in = [rowix' rowiy'];

    points = [points; add_in];
end