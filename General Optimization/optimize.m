function trans_percent_matrix = optimize(data, bpoints, betavals, sigmavals)
% Save data being iterated over
data_Orig1 = deepCopyStruct(data);
trans_percent_matrix = zeros(length(bpoints),length(betavals));
%% Inside for loop for each b
for b = bpoints
    for sigma_idx = 1:length(sigmavals)
        sigma = sigmavals(sigma_idx);
        data = data_Orig1;
        % Run time integration for each b
        % Pattern periodic into long chain %Michael
        data.b = b;
        data = initialize_time_integration(data);
        rng('default');
        data.b_vector = data.b_vector.*normrnd(1,sigma,[1,data.N])';
        data.sigma = sigma;
        if length(sigmavals) > 1
            figure(767); clf; hold on;
            histogram(data.b_vector, 20);
            ylabel('# of Arches');
            xlabel('b');
            data.file_name_sigma =    data.timeStr + "\"+data.shape_name + ", b = " + num2str(data.b) + ", NumCells = "+ num2str(data.N_cells) + ", t = "+num2str(data.t)+", sigma = " +num2str(sigma);
            sgtitle(data.shape_name + ", b = " + num2str(data.b) + ", NumCells = "+ num2str(data.N_cells) + ", t = "+num2str(data.t)+", sigma = " +num2str(sigma));
            print(gcf, "TransitionExcel\"+data.file_name_sigma + " - bHist.png", '-dpng', '-r600');

            figure(757); clf; hold on;
            scatter(1:data.N,data.b_vector,25,'filled', 'o');
            ylabel('b val');
            xlabel('Arch #');
            sgtitle(data.shape_name + ", b = " + num2str(data.b) + ", NumCells = "+ num2str(data.N_cells) + ", t = "+num2str(data.t)+", sigma = " +num2str(sigma));
            print(gcf, "TransitionExcel\"+data.file_name_sigma + " - bvArch.png", '-dpng', '-r600');

        end

        % Time integration and mitigate edge effects
        % Take a look at the initial condition
        %[~] = plot_system_once(data.A0,data);

        % Prepare for Time Integration
        data.beta = .1;
        T_end = 200;
        data = determine_coefficient_matrix(data);
        data = determine_starting_vals(data);
        data = determine_modes_to_skip(data);
        data.A0hat = determine_Ahat_from_A(data.A0,data);
        [t,Ahat] = ode45(@(t,A) arbitrary_grid_ODE(t,A,data),[0 T_end],data.A0hat);

        %% Look at end of time integration
        A = determine_A_from_Ahat(Ahat, data);
        data.plot_labels = 1;
        %[~] = plot_system_once(A(end,:)',data);
        data.plot_labels = 0;

        data.A0hat = Ahat(end,:)';
        clear Ahat;
        data.A0hat(data.N*data.N_modes+1:end) = 0;

        %% Determine which arch to force/displace %MICHAEL QUESTION

        % For now will try hardcoding into init_shape?
        % hold some nodes near edge stationary

        data.impose_rotation_at(data.nodes_to_hold) = 1;
        data.rotation_omega(data.nodes_to_hold) = 0.0013;
        data.rotation_mag(data.nodes_to_hold) = 0;

        data.impose_displacement_at(data.arches_to_displace) = 0.5;
        data.displacement_omega(data.arches_to_displace) = 2*pi/T_end;

        data.impose_displacement_at(data.arches_to_displace) = 0.5;
        data.displacement_omega(data.arches_to_displace) = 2*pi/T_end;

        data.impose_rotation_at(data.nodes_to_rotate) = 1;
        data.rotation_omega(data.nodes_to_rotate) = 2*pi/T_end;
        data.rotation_mag(data.nodes_to_rotate) = 2.5;

        data.force_eta(data.arches_to_force_positive) = .5;
        data.force_omega(data.arches_to_force_positive) = 0;
        data.force_magnitude(data.arches_to_force_positive) = .08;

        data.force_eta(data.arches_to_force_negative) = .5;
        data.force_omega(data.arches_to_force_negative) = 0;
        data.force_magnitude(data.arches_to_force_negative) = -.08;

        % data.impose_rotation_at(2) = 1;
        % data.rotation_omega(2) = 0.0013;
        % data.rotation_mag(2) = 0;
        %

        %
        % data.impose_rotation_at(46) = 1;
        % data.rotation_omega(46) = 0.0013;
        % data.rotation_mag(46) = 0;
        %
        % data.impose_rotation_at(47) = 1;
        % data.rotation_omega(47) = 0.0013;
        % data.rotation_mag(47) = 0;
        %


        % Remake data with new actuation
        data = determine_coefficient_matrix(data);
        data = determine_starting_vals(data);
        data = determine_modes_to_skip(data);
        %data.A0hat = determine_Ahat_from_A(data.A0,data);

        % We have tiled state with edge effects, now need to compare transition
        % Save data being iterated over
        data_Orig2 = deepCopyStruct(data);
        % Recover A
        %A_Orig2 = determine_A_from_Ahat(Ahat,data);

        % for loop for different beta values
        %% Run Time Integration for each Beta
        for beta = betavals
            %A = A_Orig2;
            data = data_Orig2;
            %% Prepare for Time Integration
            data.beta = beta;
            T_end = 2*min([2000*sqrt(beta/.1)-400,10000]);

            %% Run Time Integration
            [t,Ahat] = ode45(@(t,A) arbitrary_grid_ODE(t,A,data),[0 T_end],data.A0hat);

            %% Recover A
            A = determine_A_from_Ahat(Ahat,data);
            clear Ahat;


            %% Visualize the Results
            data.frames = 100;
            tinterp = linspace(t(1),t(end),data.frames);
            Ainterp = interp1(t,A,tinterp);
            clear A t;

            if ~exist("Videos\"+ data.timeStr + "\", 'dir')
                mkdir("Videos\"+data.timeStr + "\");
            end
            data.file_name = data.timeStr + "\"+ data.shape_name + " b = " + num2str(b) + " beta = " + num2str(data.beta) + " NumCells = "+ num2str(data.N_cells) + " t = "+num2str(data.t);
            data.file_name_trans = data.timeStr + "\"+ data.shape_name + " beta = " + num2str(data.beta) + " NumCells = "+ num2str(data.N_cells);

            if isfield(data, 'plot_videos') && data.plot_videos
                plot_system_over_time(tinterp,Ainterp,data)
            end

            % Determine if a transition occurred and save info (boolean? or distance of wave?)
            % Within each unit cell, calc potential energy, see which ones went from
            % high to low (within say 10% and within 3 unit/super cells)
            b_idx = round(bpoints,4) == round(b,4);
            beta_idx = round(betavals,4) == round(beta,4);

            trans_percent_matrix(b_idx,beta_idx) = snap_through_question(tinterp,Ainterp,data);
        end

        if length(sigmavals) > 1
            sigma_plot_vals(sigma_idx) = trans_percent_matrix(b_idx,1);
        end

    end
    if length(sigmavals) > 1
        figure(7676); clf; hold on;
        plot(sigmavals*100, sigma_plot_vals);
        ylabel('Percentage of Wave Propagation');
        xlabel('Sigma [% of b]');
        data.file_name_sigma =   data.timeStr + "\"+ data.shape_name + ", b = " + num2str(data.b) + ", NumCells = "+ num2str(data.N_cells) + ", t = "+num2str(data.t);
        sgtitle(data.shape_name + ", b = " + num2str(data.b) + ", NumCells = "+ num2str(data.N_cells) + ", t = "+num2str(data.t));
        print(gcf, "TransitionExcel\"+data.file_name_sigma + " - Sigma.png", '-dpng', '-r600');
    end

end
end

