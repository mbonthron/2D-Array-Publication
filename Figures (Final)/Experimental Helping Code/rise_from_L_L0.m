clear
L = 100;            % Distance between hinges[mm]
L0_prime_nominal = 102.42;  % Unbuckled length [mm]

variation = [707 702 702];

for i = 1:length(variation)
    L0_prime = L0_prime_nominal * variation(i) / variation(end);
    b = fzero(@(b) integral(@(x) sqrt(1+(pi/L*b*cos(pi*x/L)).^2),0,L)-L0_prime,5);

    fprintf("%.0f Pixels -> b = %.3f mm\n",[variation(i) b])
end