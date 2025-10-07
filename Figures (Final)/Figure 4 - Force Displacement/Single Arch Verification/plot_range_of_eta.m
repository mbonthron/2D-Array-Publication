eta_vector = 0.1:0.05:0.9;

% figure(1); clf; hold on

for i = 1:length(eta_vector)
    run_name = 'C eta ' + string(i);
    load(run_name)

    plot(displacement*1000,Q,"linewidth",2)
end

ylim([0 10])