run_name = 'rhombus_direct_run3';

%% Things Needed
run('points_rhombus_direct')

N_modes = 3;

% Determine the adjacency matrix and number of arches
adjacency_matrix = determine_adjacency_matrix(points);
adjacency_matrix = add_periodicity(adjacency_matrix);
N = sum(triu(adjacency_matrix,1) ==1,'all');

% Visualize the point and connection between the nodes
plot_grid(adjacency_matrix,points,1);

% Determin the coefficient matrix and number of constraints of the system
[coeff_matrix,constraint_count] = determine_coefficient_matrix(adjacency_matrix,0,N,N_modes);
modes_to_skip = determine_modes_to_skip(coeff_matrix,constraint_count,N,N_modes);

%% Load the data from coco for the mode shapes
bd = coco_bd_read(run_name);
UZ = coco_bd_labs(run_name, 'UZ'); 

% Constraint length
C = length(modes_to_skip);

bcrits = zeros(1,length(UZ));
A = zeros(2*(N*N_modes-C),4);

%%
for k = 1:4
    bcrits(k) = coco_bd_val(bd,UZ(k),'b');
    A(:,k) = coco_bd_val(bd,UZ(k),'x');
end

%% Recover the missing modes from the system 
% Recover the lost variables
last_C_rows = coeff_matrix(end-(C-1):end,1:N*N_modes); 
LHS = last_C_rows(:,modes_to_skip);
RHS = -1*last_C_rows(:,setdiff(1:N*N_modes,modes_to_skip));

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
    f = plot_system_once(Ahat(:,i),N,N_modes,zeros(N,1),adjacency_matrix,points);

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
A0p = Ahat(:,1);
save("b = 0.20 pi.mat","A0p")

A0p = Ahat(:,2);
save("b = 0.15 pi.mat","A0p")

A0p = Ahat(:,3);
save("b = 0.10 pi.mat","A0p")

A0p = Ahat(:,4);
save("b = 0.05 pi.mat","A0p")