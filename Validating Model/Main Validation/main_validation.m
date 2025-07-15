%% Clear Everything so there are no stragglers
clear; close all

%% Add the Paths to the Required Functions
%restoredefaultpath
%startup
addpath('Single Arch Snap Through - Barenblatt\')
addpath('Single Arch Snap Through - Dimensional\')
addpath('Single Arch Snap Through - PRL\')
addpath('HausdorffDist\')


%% Load in data
bvals = [9.3]; %mm
tvals = [1.27]; %mm
schemes_to_run = ['B'];
% Run for 3 different locations
eta_vals = [.500000000000 .49 .475];

total_cell = cell(length(bvals), length(tvals), length(eta_vals), 3, 2);
hausdorff_cell = cell(length(bvals), length(tvals), length(eta_vals));
num_samples = 1000;

%% Initialize data
data = struct();
data.L = 100 / 1000;    % Length [m]
data.EE = 2.0e9;        % Young's Modulus [N/m^2]
data.rho = 1270;        % Volumetric Density [kg/m^3]
data.indentor_speed_mms = .4;      % Indentor Speed mm/s
data.beta_PRL = .25;
data.w = 14.18 / 1000;  % Width [m]

calculate_hausdorff = 0;
plot_experimental = 1;
data.save_each = 0; %Save each inner for loop
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

    color_matrix = [color_matrix .5*ones(12,1)];

    counter = 1;

    %b_array = ["3.75" "5.625" "7.5" "9.375" "11.25" "15"];
    b_array = bvals;
    L_array = ["L = 100mm"];

    figure(4);
    cd("data Experimental\")

    for b = 1:length(b_array)
        for L = 1:length(L_array)
            file_names = dir("b = "+b_array(b) + "mm " + L_array(L) + "*.xlsx");
            upsideDown = 0;
            backwards = 0;
            first_flag_up = 0;
            first_flag_back = 0;
            first_else = 0;
            if ~isempty(file_names)
                for i = 1:length(file_names)
                    color_mod = 0;
                    legend_mod = "";
                    A = readtable(file_names(i).name,ReadVariableNames=false);
                    position = A.Var2;
                    force = A.Var3;

                    if contains(lower(file_names(i).name), "upside")
                        upsideDown = 1;
                        color_mod = 1;
                        if first_flag_up == 0
                            first_flag_up = 1;
                        end
                        legend_mod = " upside down";
                    elseif contains(lower(file_names(i).name), "back")
                        backwards = 1;
                        color_mod = 2;
                        if first_flag_back == 0
                            first_flag_back = 1;
                        end
                        legend_mod = " backwards";
                    else
                        if first_else == 0
                            first_else = 1;
                        end
                    end
                        
                    startidx = find(force>0.05);

                    position = position(startidx:end);
                    force = force(startidx:end);

                    if (first_else == 1) || (first_flag_up == 1) || (first_flag_back == 1)
                        plot(position-position(1),force,'color',color_matrix(counter+color_mod,:),"DisplayName","b = "+b_array(b)+"mm "+L_array(L)+ legend_mod,'LineWidth',2)
                    else
                        plot(position-position(1),force,'color',color_matrix(counter+color_mod,:),"HandleVisibility","off",'LineWidth',2)
                    end

                    %activate only once
                    if first_flag_up == 1
                       first_flag_up = -1;
                    end
                    if first_flag_back == 1
                       first_flag_back = -1;
                    end
                    if first_else == 1
                       first_else = -1;
                    end

                end
            end
            counter = counter + 1 + upsideDown+backwards;
        end
    end

    legend()
    cd("..\")
    if ~exist("Images/"+data.timeStr, 'dir')
        mkdir("Images/"+data.timeStr);
    end
    print(gcf, "Images/"+data.timeStr+"/"+"b = "+data.rise*1000+"mm L =" +data.L*1000 +"mm EE = "+data.EE/(1e9)+"e9 beta = "+data.beta_PRL+" - Force Displacement.png", '-dpng', '-r600');
end