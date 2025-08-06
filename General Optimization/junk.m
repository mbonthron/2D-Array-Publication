% Original vector
v = [1;1;2;2;3;3;3;4;5;4;5;6;6;8;8;8;9;9];

% Mapping matrix: [newValue, oldValue]
map = [4 3;
       5 4;
       6 5;
       7 6;
       8 7;
       9 8;
      10 9;
       1 10;
       2 11;
       3 12];

% Apply replacements
for i = 1:size(map,1)
    v(v == map(i,2)) = map(i,1);
end

disp(v)