%% Clear Everything so there are no stragglers
clear; close all

%% Add the Paths to the Required Functions
%restoredefaultpath
startup
addpath('Single Arch Snap Through - Barenblatt\')
addpath('Single Arch Snap Through - Dimensional\')
addpath('Single Arch Snap Through - PRL\')
addpath('HausdorffDist\')


%% Load in data
bvals = [10]; %mm
tvals = [.75]; %mm
schemes_to_run = ['B', 'D', 'P'];
% Run for 3 different locations
eta_vals = [.500000000000];

total_cell = cell(length(bvals), length(tvals), length(eta_vals), 3, 2);
hausdorff_cell = cell(length(bvals), length(tvals), length(eta_vals));
num_samples = 1000;

%% Initialize data
data = struct();
data.L = 100 / 1000;    % Length [m]
data.EE = 2.3e9;        % Young's Modulus [N/m^2]
data.rho = 1270;        % Volumetric Density [kg/m^3]
data.indentor_speed_mms = 0.2;      % Indentor Speed mm/s
data.beta_PRL = .25;
data.w = 14.18 / 1000;  % Width [m]

calculate_hausdorff = 0;
plot_experimental = 1;
data.save_each = 1; %Save each inner for loop
data.clear_each = 0; %Clear each inner for loop

nowTime = datetime('now');
data.timeStr = string(datestr(nowTime, 'yyyy-mm-dd_HH-MM-SS'));

for b_idx = 1:length(bvals)
    b = bvals(b_idx);
    for t_idx = 1:length(tvals)
        t = tvals(t_idx);

        data.t = t / 1000;   % Thickness [m]
        data.rise = b / 1000;  % Initial Rise [m]


        %% Run numerical results
        temp_cell = validation(data, schemes_to_run, eta_vals, num_samples);

        % Save to bigger cell
        total_cell(b_idx, t_idx, :, :, :) = reshape(temp_cell, [1 1 length(eta_vals) 3 2]);

    end
end

save()
%% Plot each relevant combination
if calculate_hausdorff
    for b_idx = 1:length(bvals)
        b = bvals(b_idx);
        for t_idx = 1:length(tvals)
            t = tvals(t_idx);
            for eta_idx = 1:length(eta_vals)
                point_set1 = [total_cell{b_idx, t_idx, eta_idx, 1, 1} total_cell{b_idx, t_idx, eta_idx, 1, 2}];
                point_set2 = [total_cell{b_idx, t_idx, eta_idx, 3, 1} total_cell{b_idx, t_idx, eta_idx, 3, 2}];

                hausdorff_cell{b_idx, t_idx, eta_idx} = HausdorffDist(point_set1, point_set2);

            end
        end
    end


    % Suppose data_cell is A x B x C
    A = size(hausdorff_cell,1);
    B = size(hausdorff_cell,2);
    C = size(hausdorff_cell,3);

    % Convert to full 3D numeric array
    data_mat = cell2mat(hausdorff_cell);           % yields A×B×C matrix
    data_mat = reshape(data_mat, A, B, C);    % ensure shape is correct

    % Choose slice indices
    c_idx = 1;
    b_idx = 1;
    a_idx = 1;

    % Extract slices
    xy_slice = data_mat(:,:,c_idx);
    xz_slice = squeeze(data_mat(:,b_idx,:));
    yz_slice = squeeze(data_mat(a_idx,:,:));

    % Plot
    figure;
    tiledlayout(1,3);

    nexttile;
    imagesc(bvals, tvals, xy_slice); colorbar; title('b-t Slice');

    nexttile;
    imagesc(bvals, eta_vals, xz_slice); colorbar; title('b-eta Slice');

    nexttile;
    imagesc(tvals, eta_vals, yz_slice); colorbar; title('t-eta Slice');
end
%% Plot experimental
if plot_experimental
    color_matrix = [    0.1190    0.5472    0.6160
        0.4984    0.1386    0.4733
        0.9597    0.1493    0.3517
        0.3404    0.2575    0.8308
        0.5853    0.8407    0.5853
        0.2238    0.2543    0.5497
        0.7513    0.8143    0.9172
        0.2551    0.2435    0.2858
        0.5060    0.9293    0.7572
        0.6991    0.3500    0.7537
        0.8909    0.1966    0.3804
        0.9593    0.2511    0.5678];

    counter = 1;

    %b_array = ["3.75" "5.625" "7.5" "9.375" "11.25" "15"];
    b_array = bvals;
    L_array = ["L = 100mm"];

    figure(4);
    cd("data Experimental\")

    for b = 1:length(b_array)
        for L = 1:length(L_array)
            file_names = dir("b = "+b_array(b) + "mm " + L_array(L) + "*.xlsx");

            if ~isempty(file_names)
                for i = 1:length(file_names)
                    A = readtable(file_names(i).name,ReadVariableNames=false);
                    position = A.Var2;
                    force = A.Var3;

                    startidx = find(force>0.05);

                    position = position(startidx:end);
                    force = force(startidx:end);

                    if i == 1
                        plot(position-position(1),force,'color',color_matrix(counter,:),"DisplayName","b = "+b_array(b)+"mm "+L_array(L),'LineWidth',2)
                    else
                        plot(position-position(1),force,'color',color_matrix(counter,:),"HandleVisibility","off",'LineWidth',2)
                    end

                end
            end
            counter = counter + 1;
        end
    end

    legend()
    cd("..\")

end