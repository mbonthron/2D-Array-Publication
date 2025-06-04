function [data] = get_mode_shape_from_coco(data)

cd('data');
rhombus_runs = dir([data.shape_name '*']);
cd ..

for runNum = 1:length(rhombus_runs)
    % for i = [11]
    run_name = rhombus_runs(i).name;

    %% Things Needed

    N_modes = data.N_modes;

    % Determine the adjacency matrix and number of arches
    adjacency_matrix = data.adjacency_matrix;
    N = data.N;

    % Visualize the point and connection between the nodes
    plot_grid(data,1);

    % Determin the coefficient matrix and number of constraints of the system
    data = determine_coefficient_matrix(data);
    data = determine_modes_to_skip(data);

    %% Load the data from coco for the mode shapes
    bd = coco_bd_read(run_name);
    UZ = coco_bd_labs(run_name, 'UZ');

    % Constraint length
    C = length(modes_to_skip);

    bcrits = zeros(1,length(UZ));
    A = zeros(2*(N*N_modes-C),length(UZ));

    %%
    for k = 1:length(UZ)
        bcrits(k) = coco_bd_val(bd,UZ(k),'b');
        A(:,k) = coco_bd_val(bd,UZ(k),'x');
    end

    %% Recover the missing modes from the system
    % Recover the lost variables
    coeff_matrix = data.coeff_matrix;
    C_rows = coeff_matrix(N_modes*N+data.V+1:N_modes*N+data.V+C,1:N*N_modes);

    LHS = C_rows(:,modes_to_skip);
    RHS = -1*C_rows(:,setdiff(1:N*N_modes,modes_to_skip));

    missingvals = (LHS\RHS)*A(1:end/2,:);
    Dmissingvals = (LHS\RHS)*A(end/2+1:end,:);

    % Produce the 'full' Ahat matrix which can be used in dVdAN
    Ahat = A;
    for i = 1:C
        mode = modes_to_skip(i);
        Ahat = [Ahat(1:mode-1,:); missingvals(i,:) ; Ahat(mode:end,:)];
    end
    shift_modes = N*N_modes;    % Do the same for the derivative terms
    for i = 1:C
        mode = modes_to_skip(i);
        Ahat = [Ahat(1:shift_modes+mode-1,:); Dmissingvals(i,:) ; Ahat(shift_modes+mode:end,:)];
    end

    %% Plot the system at each UZ point
    for i = 1:4
        f = plot_system_once(Ahat(:,i),data);

        figure
        copyobj(allchild(f),gcf);

        lbl = gcf().Number;
        text(min(xlim),max(ylim),num2str(lbl),'HorizontalAlignment','center','FontSize',14,'FontWeight','bold')
        axis off
        V_vector = calculate_energy(Ahat(:,i)',N,bcrits(i)*ones(N,1)',0.1*pi*ones(N,1)',N_modes);
        title(sprintf("%.6f",sum(V_vector)))


        figure(9899); hold on
        scatter3(Ahat(1,i),Ahat(2,i),bcrits(i),200,"MarkerFaceColor",'c',"MarkerEdgeColor","k")
        text(Ahat(1,i),Ahat(2,i),bcrits(i),num2str(lbl),'HorizontalAlignment','center','FontSize',14,'FontWeight','bold')
    end

    %% Save the Solutions as AO
    for k = 1:length(UZ)
        A0p = Ahat(:,k);
        b_val = round(bcrits(k)/pi,2);

        % Check if energy larger
        save("b = "+ num2str(b_val) +" pi.mat","A0p")
    end
end