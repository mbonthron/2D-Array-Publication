function percent_trans = snap_through_question(t,A,data)
N = data.N;
N_modes = data.N_modes;

font_name = 'Helvetica';

figure(76); clf; 
subplot(2,2,1);
hold on
sgtitle(data.shape_name + ", b = " + num2str(data.b) + ", beta = " + num2str(data.beta) + ", NumCells = "+ num2str(data.N_cells) + ", t = "+num2str(data.t), 'FontName', font_name);
title("A", 'FontName', font_name);
imagesc(A(:,1:N*N_modes))
set(gca, 'FontName', font_name);

x_num = 100;
x = linspace(0,pi,x_num);
modes = (1:N_modes)';

w_matrix = zeros(length(t),N*x_num);
w_matrix_norm = w_matrix;
w_matrix_sum = zeros(length(t),N);

for i = 1:N
    cols_of_interest = N_modes*(i-1)+1:N_modes*i;

    as = A(:,cols_of_interest);

    w = as*sin(modes.*x);
    w_norm = abs(w - w(1,:))/100;
    w_sum = sum(w_norm,2);

    w_matrix(:,x_num*(i-1)+1:x_num*i) = w;
    w_matrix_norm(:,x_num*(i-1)+1:x_num*i) = w_norm;
    w_matrix_sum(:,i) = w_sum;
end

subplot(2,2,2); hold on;
imagesc(w_matrix)
title("w vs x for all arches", 'FontName', font_name);
set(gca, 'FontName', font_name);


%%
for i = 1:N
    plot(x_num*i*[1 1],ylim,":","Color","r",'LineWidth',1)
end

%%
w_matrix_norm = abs(w_matrix - w_matrix(1,:));

w_matrix_norm(isnan(w_matrix_norm)) = data.b;

subplot(2,2,3); hold on;
title("w norm (last w - first w) vs x for all arches", 'FontName', font_name);
imagesc(w_matrix_norm)
set(gca, 'FontName', font_name);



%%
subplot(2,2,4); hold on;
title("w norm summed for all arches", 'FontName', font_name);
imagesc(w_matrix_sum)
set(gca, 'FontName', font_name);


percent_trans = find(w_matrix_sum(end,:) > data.b/4, 1, 'last' )/N*100;
if isempty(percent_trans)
    percent_trans = 0;
end
disp(num2str(percent_trans) +"%" )
plot(percent_trans*N/100*[1 1],ylim,":","Color","r",'LineWidth',1)

if ~exist("Videos\", 'dir')
    mkdir("Videos\");
end

file_name = "Videos\"+data.file_name;

exportgraphics(gcf,file_name + " - Transition.png","Resolution",600)

end