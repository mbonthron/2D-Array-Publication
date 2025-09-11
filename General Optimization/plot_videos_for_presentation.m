count_array = [201 1552];

time_array = 2000*ones(size(count_array));

force_sign = [1 1];

for i = 1:length(count_array)
    run_and_plot(count_array(i),0.1,time_array(i),force_sign(i))
end