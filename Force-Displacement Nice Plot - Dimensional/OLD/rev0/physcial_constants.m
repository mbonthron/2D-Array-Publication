%% Dimensional Parameters
L         = 100 / 1000;      % Length [m]
t         = 1 / 1000;     % Thickness [m]
w         = 14.18 / 1000;    % Width [m]
rise      = 10 / 1000;       % Initial Rise [m]
EE        = 2.3e9;           % Young's Modulus [N/m^2]

rho = 1270;                  % Volumetric Density [kg/m^3]

indentor_speed_mms = 0.2;    % Indentor Speed mm/s

% Derived parameters
II = 1/12 * w * t^3;
AA = t*w;

%% Save into data structure
data.rho = rho;
data.AA  = AA;
data.II  = II;
data.EE  = EE;
data.L   = L;

data.bending = EE*II*pi^4/L^4;
data.stretch = AA*EE*pi^2/L^2;
