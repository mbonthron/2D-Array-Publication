%% Dimensional Parameters
L         = 100/1000;      % Length [m]
t         = 0.75/1000;     % Thickness [m]
w         = 14.18 / 1000;    % Width [m]
rise      = 7.5 / 1000;       % Initial Rise [m]
EE        = 1.9e9;           % Young's Modulus [N/m^2]

rho = 1270;                  % Volumetric Density [kg/m^3]

% Derived parameters
II = 1/12 * w * t^3;
AA = t*w;

%% Save into data structure
data.rho = rho;
data.AA  = AA;
data.II  = II;
data.EE  = EE;
data.L   = L;

data.t    = t;
data.t_ND = t*pi/L;

data.rise = rise;
data.rise_ND = rise*pi/L;

data.bending = EE*II*pi^4/L^4;
data.stretch = AA*EE*pi^2/L^2;
