% clear
L = 100;            % Distance between hinges[mm]


L0_prime_nominal = [104.84 104.48 103.80];
b_vector  = [];
for i = 1:length(L0_prime_nominal)
    L0_prime = L0_prime_nominal(i);
    b_vector(i) = fzero(@(b) integral(@(x) sqrt(1+(pi/L*b*cos(pi*x/L)).^2),0,L)-L0_prime,5);

    fprintf("%.2f mm -> b = %.3f mm\n",[L0_prime_nominal(i) b_vector(i)])
end

%
idx_nominal = 3;
b_vector ./ b_vector(idx_nominal)