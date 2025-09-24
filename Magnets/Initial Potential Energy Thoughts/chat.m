clear; clc; close all;

% -------- Parameters you set --------
Br     = 1.2;       % Tesla, remanence (N52 NdFeB ~1.2 T)
D      = 0.01;      % m, diameter of magnets
L      = 0.01;      % m, length of track magnets
Lp     = 0.01;      % m, length of probe magnet
a      = 0.02;      % m, spacing between track magnets
h      = 0.005;     % m, vertical stand-off distance of probe magnet
N      = 10;        % number of magnets each side of probe to sum

% -------- Derived quantities --------
mu0 = 4*pi*1e-7;
V   = pi*(D/2)^2*L;     % volume of track magnet
m   = (Br/mu0)*V;       % dipole moment of each magnet
qm  = m/L;              % equivalent magnetic pole strength
qp  = qm;               % assume probe magnet has same strength

% -------- Probe position sweep --------
xspan = linspace(-3*a, 3*a, 1000);   % probe position over ~two periods
U = zeros(size(xspan));            % potential energy

% -------- Compute energy --------
for ix = 1:length(xspan)
    x = xspan(ix);
    Utot = 0;
    % loop over track magnets
    for n = -N:N
        x_n = n*a;  % position of nth magnet center
        % charges of track magnet
        for alpha = [+1, -1]  % +1 = north, -1 = south face
            q_n = ((-1)^n)*alpha*qm;   % alternating polarity
            x_n_alpha = x_n + alpha*L/2;
            % probe charges
            for beta = [+1, -1]
                q_p = beta*qp;
                x_p_beta = x + beta*Lp/2;
                r = sqrt((x_p_beta - x_n_alpha)^2 + h^2);
                Utot = Utot + (mu0/(4*pi))*q_n*q_p / r;
            end
        end
    end
    U(ix) = Utot;
end

% -------- Plot results --------
U = U - min(U); % shift baseline
figure;
plot(xspan*1e3, U*1e6, 'LineWidth', 1.5);
xlabel('Probe position x (mm)');
ylabel('Potential energy (µJ)');
title('Bistable Magnetic Potential Energy Landscape');
grid on;

% Mark wells and barrier
% [~, locs] = findpeaks(-U); % minima
% [~, barriers] = findpeaks(U); % maxima
% hold on;
% plot(xspan(locs)*1e3, U(locs)*1e6, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
% plot(xspan(barriers)*1e3, U(barriers)*1e6, 'kx', 'MarkerSize', 10, 'LineWidth', 1.5);
% legend('Energy landscape', 'Wells', 'Barrier');