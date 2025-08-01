function percent_trans = snap_through_question(t,A,data)
%% Plot
    N = data.N;
    N_modes = data.N_modes;
    
    font_name = 'Helvetica';
    
    f = figure(76); clf;
    hold on;
    tl = tiledlayout(f, 2,2,'Padding', 'loose', 'TileSpacing', 'compact');
    ax = nexttile(tl);
    
    sgtitle(data.shape_name + ", b = " + num2str(data.b) + ", beta = " + num2str(data.beta) + ", NumCells = "+ num2str(data.N_cells) + ", t = "+num2str(data.t)+", sigma = " +num2str(data.sigma), 'FontName', font_name);
    
    imagesc(A(:,1:N*N_modes))
    title("A", 'FontName', font_name);
    
    xlim(ax,[1 N*N_modes]);
    ylim(ax,[1 length(t)]);
    xlabel(ax, "Modes");
    ylabel(ax, "Time");
    colorbar(ax);
    set(gca, 'FontName', font_name);
    
    x_num = 25;
    x = linspace(0,pi,x_num);
    modes = (1:N_modes)';
    
    w_matrix = zeros(length(t),N*x_num);
    w_matrix_norm = w_matrix;
    w_matrix_sum = zeros(length(t),N);
    
    for i = 1:N
        cols_of_interest = N_modes*(i-1)+1:N_modes*i;
    
        as = A(:,cols_of_interest);
    
        w = as*sin(modes.*x);
        w_norm = abs(w - w(1,:))/x_num;
        w_sum = sum(w_norm,2);
    
        %w_matrix(:,x_num*(i-1)+1:x_num*i) = w;
        w_matrix_norm(:,x_num*(i-1)+1:x_num*i) = w_norm;
        w_matrix_sum(:,i) = w_sum;
    end
    clear A;
    w_size = size(w_matrix,2);
    
    
    ax = nexttile(tl); hold on;
    imagesc(w_matrix)
    xlim(ax,[1 w_size]);
    ylim(ax,[1 length(t)]);
    xlabel(ax, "All x's");
    ylabel(ax, "Time");
    colorbar(ax);
    title("Arch position vs x (all arches)", 'FontName', font_name);
    set(gca, 'FontName', font_name);
    
    
    %%
    % for i = 1:N
    %     plot(x_num*i*[1 1],ylim,":","Color","r",'LineWidth',1)
    % end
    
    %%
    w_matrix_norm = abs(w_matrix - w_matrix(1,:));
    clear w_matrix;
    
    w_matrix_norm(isnan(w_matrix_norm)) = data.b;
    
    ax = nexttile(tl); hold on;
    imagesc(w_matrix_norm)
    xlim(ax,[1 w_size]);
    ylim(ax,[1 length(t)]);
    colorbar(ax);
    xlabel(ax, "All x's");
    ylabel(ax, "Time");
    title("\Delta Arch Position vs x (all arches)", 'FontName', font_name);
    
    set(gca, 'FontName', font_name);
    clear w_matrix_norm;
    
    
    
    %%
    ax = nexttile(tl); hold on;
    imagesc(w_matrix_sum)
    xlim(ax,[1 N]);
    ylim(ax,[1 length(t)]);
    colorbar(ax);
    xlabel(ax, "Arch");
    ylabel(ax, "Time");
    title("\Delta Arch Position vs Arch", 'FontName', font_name);
    imagesc(w_matrix_sum)
    set(gca, 'FontName', font_name);
    

    percent_trans = find(w_matrix_sum(end,:) > data.b/4, 1, 'last' )/N*100;
    if isempty(percent_trans)
        percent_trans = 0;
    end
    clear w_matrix_sum;
    disp(num2str(percent_trans) +"%" )
    plot(percent_trans*N/100*[1 1],ylim,":","Color","r",'LineWidth',1)
%% Save
if ~exist("Videos\", 'dir')
    mkdir("Videos\");
end

file_name = "Videos\"+data.file_name+", sigma = " +num2str(data.sigma);

%axes(f,'Position',[0 0 999 999],'Visible','off');

%exportgraphics(f,file_name + " - Transition.png","Resolution",600)
set(f, 'PaperPositionMode', 'auto');
print(f, file_name + " - Transition.png", '-dpng', '-r600');
clear f ax tl;

end