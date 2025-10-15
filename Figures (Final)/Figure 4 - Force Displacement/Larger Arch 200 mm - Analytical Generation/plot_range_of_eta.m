eta_vector1 = 0.1:0.05:0.9;

figure(1); hold on


eta_vector2 = eta_vector1(1:17);
for i = eta_vector2
    idx = find(i==eta_vector1,1,"first");
    run_name = 'S eta ' + string(idx);
    load(run_name)

    plot(displacement*1000,Q,"k--","linewidth",2)
end

ylim([0 5])