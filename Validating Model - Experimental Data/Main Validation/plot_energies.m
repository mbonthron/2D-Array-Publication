%% Plot experimental
color_matrix = [...                   % Color matrix for plotting
    0.9649    0.0344    0.2575
    0.1576    0.4387    0.8407
    0.9706    0.3816    0.2543
    0.9572    0.7655    0.8143
    0.4854    0.7952    0.2435
    0.8003    0.1869    0.9293
    0.1419    0.4898    0.3500
    0.4218    0.4456    0.1966
    0.9157    0.6463    0.2511
    0.7922    0.7094    0.6160
    0.9595    0.7547    0.4733
    0.6557    0.2760    0.3517
    0.0357    0.6797    0.8308
    0.8491    0.6551    0.5853
    0.9340    0.1626    0.5497
    0.6787    0.1190    0.9172
    0.7577    0.4984    0.2858
    0.7431    0.9597    0.7572
    0.3922    0.3404    0.7537
    0.6555    0.5853    0.3804
    0.1712    0.2238    0.5678
    0.7060    0.7513    0.0759
    0.0318    0.2551    0.0540
    0.2769    0.5060    0.5308
    0.0462    0.6991    0.7792
    0.0971    0.8909    0.9340
    0.8235    0.9593    0.1299
    0.6948    0.5472    0.5688
    0.3171    0.1386    0.4694
    0.9502    0.1493    0.0119];

color_matrix = [color_matrix .5*ones(length(color_matrix),1)];   % Add transparency column

counter = 1;

%b_array = ["3.75" "5.625" "7.5" "9.375" "11.25" "15"];
b_array = bvals;                   % Use bvals for labeling
L_array = ["L = " + data.L*1000+"mm"];  % Span label

figure(4); clf; hold on; legend()
cd("data Experimental\")           % Change directory to experimental data

figure(401); clf; hold on; legend();

if ~check_file_names
    b_array = b_array(1);         % Use only first values if not checking names
    L_array = L_array(1);
    tvals = tvals(1);
end
for b = 1:length(b_array)
    for L = 1:length(L_array)
        for t=1:length(tvals)
            if check_file_names
                file_names = dir("b = "+num2str(b_array(b)) + "mm " + L_array(L) + " t = "+num2str(tvals(t)) + "mm *.xlsx");
            else
                file_names = dir("*.xlsx");     % Get all .xlsx files

            end
            upsideDown = 0; backwards = 0; flipped = 0; unscrewed = 0;
            first_flag_up = 0; first_flag_back = 0;
            first_flag_flip = 0; first_flag_unscrewed = 0;
            first_else = 0;

            file_name_string_cell = {          % Define expected keywords
                "default", 0;
                "upside", 0;
                "back", 0;
                "flipped", 0;
                "unscrew", 0;
                "1hole", 0;
                "2hole", 0;
                "3hole", 0;
                "4hole", 0;
                "5hole", 0;
                "1rein", 0;
                "2rein", 0;
                "3rein", 0;
                "4rein", 0;
                "5rein", 0;
                };

            if ~isempty(file_names)
                for i = 1:length(file_names)
                    matches = arrayfun(@(p) contains(string(lower(file_names(i).name)), lower(p)), string(file_name_string_cell(:,1)));  % Check file name keywords
                    matched_rows = file_name_string_cell(matches, :);    % Get matching labels
                    color_mod = find(matches ,1,'first')-1;              % Find matching index

                    if isempty(color_mod)
                        color_mod = 0;    % Default to 0 if no match
                    end

                    legend_mod = "";     % Initialize legend modifier
                    A = readtable(file_names(i).name,ReadVariableNames=false);  % Read file
                    position = A.Var2;   % Displacement
                    force = A.Var3;      % Force


                    if file_name_string_cell{color_mod+1,2} == 0
                        file_name_string_cell{color_mod+1,2} = 1;   % Mark as first appearance
                        switch(color_mod)  % Set legend label
                            case 0, legend_mod = "";
                            case 1, legend_mod = " upside down";
                            case 2, legend_mod = " backwards";
                            case 3, legend_mod = " top plate flipped";
                            case 4, legend_mod = " undid and retightened";
                            case 5, legend_mod = " 1 holed arch";
                            case 6, legend_mod = " 2 holed arch";
                            case 7, legend_mod = " 3 holed arch";
                            case 8, legend_mod = " 4 holed arch";
                            case 9, legend_mod = " 5 holed arch";
                            case 10, legend_mod = " 1 reinforcement arch";
                            case 11, legend_mod = " 2 reinforcement arch";
                            case 12, legend_mod = " 3 reinforcement arch";
                            case 13, legend_mod = " 4 reinforcement arch";
                            case 14, legend_mod = " 5 reinforcement arch";
                            otherwise, legend_mod = "";
                        end
                    end
                    startidx = find(force>0.06);  % Trim noisy beginning

                    position = position(startidx:end);  % Trim data
                    force = force(startidx:end);

                    energy = cumtrapz(position-position(1),force);


                    if file_name_string_cell{color_mod+1,2} == 1
                        figure(4)
                        plot(position-position(1),force,'color',color_matrix(counter+color_mod,:),"DisplayName","b = "+b_array(b)+"mm "+L_array(L)+ " t = "+tvals(t) + "mm"+ legend_mod,'LineWidth',2)
                        figure(401); hold on
                        plot(position-position(1),energy,'color',color_matrix(counter+color_mod,:),"DisplayName","b = "+b_array(b)+"mm "+L_array(L)+ " t = "+tvals(t) + "mm"+ legend_mod,'LineWidth',2)

                    else
                        figure(4)
                        plot(position-position(1),force,'color',color_matrix(counter+color_mod,:),"HandleVisibility","off",'LineWidth',2)
                        figure(401); hold on
                        plot(position-position(1),energy,'color',color_matrix(counter+color_mod,:),"HandleVisibility","off",'LineWidth',2)

                    end

    
                    


                    if file_name_string_cell{color_mod+1,2} == 1
                        file_name_string_cell{color_mod+1,2} = -1;  % Mark as plotted
                    end
                end
            end
            counter = counter + sum(cell2mat(file_name_string_cell(:,2))==-1);  % Advance color index
        end
    end
end
title(sprintf(title_mod));        % Set figure title
legend()                          % Show legend
cd("..\")                         % Return to previous directory
if ~exist("Images/"+data.timeStr, 'dir')
    mkdir("Images/"+data.timeStr);   % Make directory if it doesn't exist
end
print(gcf, "Images/"+data.timeStr+"/"+"b = "+data.rise*1000+"mm L =" +data.L*1000 +"mm EE = "+data.EE/(1e9)+"e9 beta = "+data.beta_PRL+" - Force Displacement.png", '-dpng', '-r600');  % Save plot

%%
f = figure(4);
f.Position(3:4) = [900 420];
leg = legend();
leg.FontSize = 10;
leg.Location = "southoutside";
leg.NumColumns = 2;
xlabel('Displacement - [mm]')
ylabel('Force - [N]')
grid()

f = figure(401);
f.Position(3:4) = [900 420];
leg = legend();
leg.FontSize = 10;
leg.Location = "southoutside";
leg.NumColumns = 2;
xlabel('Displacement - [mm]')
ylabel('Energy - [Nmm]')
grid()