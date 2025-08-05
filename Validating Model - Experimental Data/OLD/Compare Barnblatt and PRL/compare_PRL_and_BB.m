%% ========= Physical Dimensional Quantities =========
% PETG
% material = "petg";
% rho = 1270;     % [kg/m^3]
% EE   = 4.3*1e9;  % [N/m^2]

% Steel
material = "steel";
rho = 7930;     % [kg/m^3]
EE   = 210*1e9;  % [N/m^2]

% Elite 32
% material = "elite32";
% rho = 1070;        % [kg/m^3]
% EE   = 0.234e6;  % [N/m^2]  


%  Dimensional Constants
L = 120 / 1000; % [m] Length
t = 0.01*25.4 / 1000;   % [m] Thicknesss
w = 0.25*25.4 / 1000;  % [m] Width

b = 0.5*9.101 / 1000;  % [m] Rise

alpha = 0.999/1000;     % [m/s] Indentor Speed
eta = 0.5;              % [UL] Location of Indentor

beta = 2.5;     % [Ns/m^2] Damping

%  Derived Constants
II = 1/12*w*t^3;
AA = w*t;
rr = sqrt(II/AA);


N_modes = 5;
file_string = material + " b = "+b*1000 + "mm t = "+t*1000+"mm.mat";



%% ========= PRL Non-dimensionalization =========
PRL = struct();
PRL.file_string = file_string;
PRL.b = b / rr;
PRL.e = 0;

PRL.rr = rr;

PRL.beta = beta*L^2/(pi^2*sqrt(rho*AA*EE*II));
PRL.alpha = alpha/rr * L^2/pi^2*sqrt(AA*rho/EE/II);
PRL.eta = eta;
PRL.N_modes = N_modes;

PRL.Vbar_to_V = (pi^3*EE*II*rr)/(2*L^2);
PRL.tbar_to_t = L^2/pi^2*sqrt(AA*rho/EE/II);
PRL.Qbar_to_Q = pi^4*EE*II*rr / (2*L^3);

%% ========= Barenblatt Non-dimensionalization =========
BB = struct();
BB.file_string = file_string;
BB.b = b * pi / L;
BB.e = 0;
BB.t = t * pi / L;
BB.beta = beta * pi / (L*sqrt(rho*EE));
BB.alpha = alpha*pi/L*sqrt(rho/EE);
BB.eta = eta;
BB.N_modes = N_modes;

BB.Vbar_to_V = 1 / (EE*L^3);


%% Run PRL Time Integration Code
PRL = run_PRL(PRL);

%% Run BB Time Integration Code
% BB = run_BB(BB);