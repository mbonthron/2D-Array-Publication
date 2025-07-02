%% PETG
% rho = 1270;     % [kg/m^3]
% EE   = 4.3*1e9;  % [N/m^2]

%% Steel
rho = 7930;     % [kg/m^3]
EE   = 210*1e9;  % [N/m^2]

% %% Elite 32
% rho = 1;        % [kg/m^3]
% EE   = 0.234e6;  % [N/m^2]  


%% ===================================================================
%  Dimensional Constants
L = 100 / 1000; % [m] Length
w = 11 / 1000;  % [m] Width
alpha = 1/1000; % [m/s] Indentor Speed
beta = 0.5;     % [Ns/m^2] Damping

b = linspace(1,20,25)/1000;     % [m] Rise
t = linspace(0.254,1,25)/1000;   % [m] Thickness

[B T] = meshgrid(b,t);

%% ===================================================================
%  Derived Constants
II = 1/12*w*T.^3;
AA = w*T;
rr = sqrt(II./AA);

%% PRL
PRLB = B ./ rr;

surf(B*1000,T*1000,PRLB)
xlabel('Rise - [mm]')
ylabel('Thickness - [mm]')
zlabel('b - [ND]')

%% PRL Rubber Arch
b = 2 *sqrt(3);