%% Clear Everything so there are no stragglers
clear; close all

%% Add the Paths to the Required Functions
restoredefaultpath
startup
addpath('Single Arch Snap Through - Barenblatt\')
addpath('Single Arch Snap Through - Dimensional\')
addpath('Single Arch Snap Through - PRL\')

%% Load in data


%% Initialize data
data = struct();
data.L = 100 / 1000;    % Length [m]
data.t = 0.75 / 1000;   % Thickness [m]
data.w = 14.18 / 1000;  % Width [m]
data.rise = 10 / 1000;  % Initial Rise [m]
data.EE = 2.3e9;        % Young's Modulus [N/m^2]
data.rho = 1270;        % Volumetric Density [kg/m^3]
data.indentor_speed_mms = 0.2;      % Indentor Speed mm/s
data.beta_PRL = .25;


%% Run numerical results

schemes_to_run = ['B', 'D', 'P'];
% Run for 3 different locations
eta_vals = [.500000000000];
for eta = eta_vals
    data.eta = eta;

    for i = 1:length(schemes_to_run)
        switch(schemes_to_run(i))
            case 'B'
                % Make directory and cd into it
                if ~exist("data Barenblatt", 'dir')
                    mkdir("data Barenblatt");
                end
                cd("data Barenblatt");

                % run data
                snap_through_analytical_B(data);

                % cd out
                cd("..\")
            case 'D'
                if ~exist("data Dimensional", 'dir')
                    mkdir("data Dimensional");
                end
                cd("data Dimensional");

                snap_through_analytical_D(data);

                cd("..\");
            case 'P'

                if ~exist("data PRL", 'dir')
                    mkdir("data PRL");
                end
                cd("data PRL");

                snap_through_analytical_P(data);

                cd("..\");
            otherwise
                error("Wrong Scheme Identification");
        end
    end
end

%% Plot each relevant combination
