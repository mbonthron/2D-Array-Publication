function plot_vals_small = validation(data, schemes_to_run, eta_vals, num_points)
save_each = data.save_each;
clear_each = data.clear_each;
plot_vals = cell(length(eta_vals),3, 2);
plot_vals_small = cell(length(eta_vals),3, 2);
for eta_idx = length(eta_vals)
    data.eta = eta_vals(eta_idx);
    for i = 1:length(schemes_to_run)
        switch(schemes_to_run(i))
            case 'B'
                % Make directory and cd into it
                if ~exist("data Barenblatt", 'dir')
                    mkdir("data Barenblatt");
                end
                cd("data Barenblatt");

                % run data
                [disp_dimensional_B, Rxn_dimensional_B] = snap_through_analytical_B(data);
                plot_vals{eta_idx, 1, 1} = disp_dimensional_B;
                plot_vals{eta_idx, 1, 2} = Rxn_dimensional_B;

                disp_dimensional_B_q = linspace(min(disp_dimensional_B), max(disp_dimensional_B), num_points);
                Rxn_dimensional_B_q = interp1(disp_dimensional_B, Rxn_dimensional_B, disp_dimensional_B_q, 'linear');
                plot_vals_small{eta_idx, 1, 1} = disp_dimensional_B_q;
                plot_vals_small{eta_idx, 1, 2} = Rxn_dimensional_B_q;

                % cd out
                cd("..\")
            case 'D'
                if ~exist("data Dimensional", 'dir')
                    mkdir("data Dimensional");
                end
                cd("data Dimensional");

                [disp_dimensional_D, Rxn_dimensional_D] = snap_through_analytical_D(data);
                plot_vals{eta_idx, 2, 1} = disp_dimensional_D;
                plot_vals{eta_idx, 2, 2} = Rxn_dimensional_D;

                disp_dimensional_D_q = linspace(min(disp_dimensional_D), max(disp_dimensional_D), num_points);
                Rxn_dimensional_D_q = interp1(disp_dimensional_D, Rxn_dimensional_D, disp_dimensional_D_q, 'linear');
                plot_vals_small{eta_idx, 2, 1} = disp_dimensional_D_q;
                plot_vals_small{eta_idx, 2, 2} = Rxn_dimensional_D_q;

                cd("..\");
            case 'P'

                if ~exist("data PRL", 'dir')
                    mkdir("data PRL");
                end
                cd("data PRL");

                [disp_dimensional_P, Rxn_dimensional_P] = snap_through_analytical_P(data);
                plot_vals{eta_idx, 3, 1} = disp_dimensional_P;
                plot_vals{eta_idx, 3, 2} = Rxn_dimensional_P;

                disp_dimensional_P_q = linspace(min(disp_dimensional_P), max(disp_dimensional_P), num_points);
                Rxn_dimensional_P_q = interp1(disp_dimensional_P, Rxn_dimensional_P, disp_dimensional_P_q, 'linear');
                plot_vals_small{eta_idx, 3, 1} = disp_dimensional_P_q;
                plot_vals_small{eta_idx, 3, 2} = Rxn_dimensional_P_q;

                cd("..\");
            otherwise
                error("Wrong Scheme Identification");
        end
    end
    if save_each
        figure(4);
        % title("Transition Percent " + data.shape_name + " beta = " + num2str(beta) + " NumCells = "+ num2str(data.N_cells));
        % xlabel("t");
        % ylabel("b");
        % clim([0 100]);
        % a=colorbar;
        % a.Label.String = 'Transition Percent %';
        % imagesc(tvals/pi, bpoints/pi, squeeze(trans_percent_tensor(:,beta_idx,:)))
        if ~exist("Images/"+data.timeStr, 'dir')
            mkdir("Images/"+data.timeStr);
        end
        print(gcf, "Images/"+data.timeStr+"/"+"b = "+data.rise*1000+"mm L =" +data.L*1000 +"mm EE = "+data.EE/(1e9)+"e9 beta = "+data.beta_PRL+" eta =" + data.eta + " - Force Displacement.png", '-dpng', '-r600');
    end
    if clear_each
        figure(4); clf;
    end
end
end

