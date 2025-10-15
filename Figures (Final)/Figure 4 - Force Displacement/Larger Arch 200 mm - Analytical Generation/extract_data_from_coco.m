clear; close all
[midptstab1,midptunst1,b1] = get_data_from_coco('dogbone1');
[midptstab1_1,midptunst1_1,b1_1] = get_data_from_coco('dogbone1-1');

[midptstab1_8,midptunst1_8,b1_8] = get_data_from_coco('dogbone1-8');
[midptstab1_9,midptunst1_9,b1_9] = get_data_from_coco('dogbone1-9');

save("Dogbone Bifurcation.mat")
