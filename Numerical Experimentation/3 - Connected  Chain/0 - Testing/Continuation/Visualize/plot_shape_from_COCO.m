function [] = plot_shape_from_COCO(run_name,data)
%PLOT_SHAPE_FROM_COCO Given a run from coco, plots the points defined by UZ
%   points in the coco run
%%
N = data.N;
N_modes = data.N_modes;
adjacency_matrix = data.adjacency_matrix;
points = data.points;
coeff_matrix = data.coeff_matrix;
modes_to_skip = data.modes_to_skip;


%% Load the data from coco for the mode shapes
bd = coco_bd_read(run_name);
UZ = coco_bd_labs(run_name, 'UZ'); 

% Constraint length
C = data.constraint_count;

bcrits = zeros(1,length(UZ));
Ahat = zeros(2*(N*N_modes-C),length(UZ));

%%
for k = 1:length(UZ)
    bcrits(k) = coco_bd_val(bd,UZ(k),'b');
    Ahat(:,k) = coco_bd_val(bd,UZ(k),'x');
end

%% Recover the missing modes from the system 
A = determine_A_from_Ahat(Ahat',data)';


%% Plot the system at each UZ point
for i = 1:length(UZ)
    f = plot_system_once(A(:,i),data);

    figure
    copyobj(allchild(f),gcf);

    lbl = gcf().Number;
    text(min(xlim),max(ylim),num2str(lbl),'HorizontalAlignment','center','FontSize',14,'FontWeight','bold')
    axis off
    V_vector = calculate_energy(A(:,i)',N,bcrits(i)*ones(N,1)',0.1*pi*ones(N,1)',N_modes);
    title(sprintf("%.6f",sum(V_vector)))
 

    figure(9899); hold on
    scatter3(A(1,i),A(2,i),bcrits(i),200,"MarkerFaceColor",'c',"MarkerEdgeColor","k")
    text(A(1,i),A(2,i),bcrits(i),num2str(lbl),'HorizontalAlignment','center','FontSize',14,'FontWeight','bold')
end

end