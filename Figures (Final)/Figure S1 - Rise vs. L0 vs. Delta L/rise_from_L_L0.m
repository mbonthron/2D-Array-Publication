% clear
L = 100;            % Distance between hinges[mm]


L0_prime_vector = linspace(100.01,102,25);
b_vector  = [];
for i = 1:length(L0_prime_vector)
    L0_prime = L0_prime_vector(i);
    b_vector(i) = fzero(@(b) integral(@(x) sqrt(1+(pi/L*b*cos(pi*x/L)).^2),0,L)-L0_prime,5);

    fprintf("%.2f mm -> b = %.3f mm\n",[L0_prime_vector(i) b_vector(i)])
end


deltaL = L0_prime_vector - L;
%%
f = figure(1); clf; hold on
f.Units = "inches";
f.Position(3:4) = 1.05*[1.35 0.75];

plot(deltaL,b_vector)
plot(deltaL,-1*b_vector)

plot(deltaL,0*deltaL)

% ylim([0 20])

xlabel("$\Delta L$ - [mm]")
ylabel("$b$ - [mm]")
set(gca,'FontSize',10)
