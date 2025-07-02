%% PETG
% rho = 1270;     % [kg/m^3]
% EE   = 4.3*1e9;  % [N/m^2]

%% Steel
rho = 7930;     % [kg/m^3]
EE   = 210*1e9;  % [N/m^2]

% %% Elite 32
rho = 1;        % [kg/m^3]
EE   = 0.234e6;  % [N/m^2]  


%% ===================================================================
%  Dimensional Constants
L = 100 / 1000; % [m] Length
t = 1 / 1000;   % [m] Thicknesss
w = 11 / 1000;  % [m] Width

b = 10 / 1000;  % [m] Rise
alpha = 1/1000; % [m/s] Indentor Speed

beta = 0.5;     % [Ns/m^2] Damping


%% ===================================================================
%  Derived Constants
II = 1/12*w*t^3;
AA = w*t;
rr = sqrt(II/AA);

%% ===================================================================
%  PRL Non-Dimensionalization
PRL = struct();
PRL.b = b / rr;
PRL.beta = beta*L^2/(pi^2*sqrt(rho*AA*EE*II));
PRL.alpha = alpha * L^2/pi^2*sqrt(AA*rho/EE/II);

PRL.Vbar_to_V = 2*L^2/(pi^3*EE*II*rr);

%% ===================================================================
%  Barenblatt Non-Dimensionalization
BB = struct();
BB.b = b * pi / L;
BB.t = t * pi / L;
BB.beta = beta * pi / (L*sqrt(rho*EE));
BB.alpha = alpha*sqrt(rho/EE);

BB.Vbar_to_V = 1 / (EE*L^3);