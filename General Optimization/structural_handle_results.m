index_idx       = find(strcmp(results(1,:),'Index'));
binary_idx      = find(strcmp(results(1,:),'Binary'));
connected_idx   = find(strcmp(results(1,:),'Connected'));
stable_idx      = find(strcmp(results(1,:),'Stable Found'));
asym_idx        = find(strcmp(results(1,:),'Asymmetric Wells'));
vhigh_idx       = find(strcmp(results(1,:),'VHigh'));
vlow_idx        = find(strcmp(results(1,:),'VLow'));
thigh_idx       = find(strcmp(results(1,:),'thetaHigh'));
tlow_idx        = find(strcmp(results(1,:),'thetaLow'));
archd_idx       = find(strcmp(results(1,:),'Arch Displacement'));
prohib_idx      = find(strcmp(results(1,:),'Prohibiting Walls'));


[m,n] = size(results);



%%
connected_count = 0;
asymmetric_count = 0;
prohibiting_walls_count =  0;

unconnected_systems = [];
prohib_walls_systems = [];

for count = 0:(m-2)
    % Count the number of connected systems
    if results{count+2,connected_idx} == 1
        connected_count = connected_count + 1;
    else
        unconnected_systems = [unconnected_systems count];
    end

    % Count the number of system with asymmetric wells
    if results{count+2,asym_idx} == 1 
        asymmetric_count = asymmetric_count + 1;
    end

    % Count the number of samples without prohibiting walls
    if results{count+2,prohib_idx} == 0
        prohibiting_walls_count = prohibiting_walls_count + 1;
    end
end

%% Raw Count
total_cases = 2^11;


per_connected = connected_count / total_cases;
per_asym      = asymmetric_count / total_cases;
per_no_walls  = prohibiting_walls_count / total_cases;


clc
fprintf("%0.f Total Cases\n====================================\n",total_cases)
fprintf("Connected:         %.0f \t %.2f percent\n",[connected_count per_connected*100])
fprintf("Asymmetric:        %.0f   \t %.2f percent\n",[asymmetric_count per_asym*100])
fprintf("No prohibit walls: %.0f   \t %.2f percent\n",[prohibiting_walls_count per_no_walls*100])


%% Running Count
connected_count2 = 0;
asymmetric_count2 = 0;
prohibiting_walls_count2 =  0;

check_out = [];

for count = 0:(m-2)
    % Check if the system is connected
    if results{count+2,connected_idx} == 1
        connected_count2 = connected_count2 + 1;

        % Count the number of system with asymmetric wells
        if results{count+2,asym_idx} == 1 
            asymmetric_count2 = asymmetric_count2 + 1;

            % Count the number of samples without prohibiting walls
            if results{count+2,prohib_idx} == 0
                prohibiting_walls_count2 = prohibiting_walls_count2 + 1;
                
                check_out = [check_out ; count];
            else
                prohib_walls_systems = [prohib_walls_systems count];
            end
        end
    end
end

per_connected2 = connected_count2 / total_cases;
per_asym2      = asymmetric_count2 / total_cases;
per_no_walls2  = prohibiting_walls_count2 / total_cases;


fprintf("\n\n%0.f Running Count\n====================================\n",total_cases)
fprintf("Connected:         %.0f \t %.2f percent\n",[connected_count2               per_connected2*100])
fprintf("Asymmetric:        %.0f   \t %.2f percent\n",[asymmetric_count2            per_asym2*100])
fprintf("No prohibit walls: %.0f   \t %.2f percent\n",[prohibiting_walls_count2     per_no_walls2*100])
