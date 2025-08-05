load("Analytical Prediction.mat");
Q_anal  = Q_N;
T_anal  = T_sec;
aN_anal = aN_mm;

load("LabView Data.mat")
Q_lab = F;
T_lab = t_sec;
x_lab = x;

load("Tracked Video.mat")
fps_vid = 13.98;
U_vid = node_height_pixels*mm_per_pixel;
nodexN = 1 - x_position_UL;

% Include the node x locations here in the form of a 1xN vector
number_of_nodes = size(nodexN,2);


aN_vid = zeros(length(U_vid),3);
for i = 1:length(aN_vid)
    % Create the coefficient matrix at each time
    coeff_matrix = zeros(number_of_nodes,3);
    for j = 1:number_of_nodes
        coeff_matrix(j,:) = [sin(pi*nodexN(i,j)) sin(2*pi*nodexN(i,j)) sin(3*pi*nodexN(i,j))];
    end
    
    % Use linear algebra to determine what each mode is
    aN_vid(i,:) = coeff_matrix \ -U_vid(i,:).';
end

T_vid = (1:length(U_vid))/fps_vid;
%%
color1 = [0 114 189]/256;  color1L = [151 214 255]/256;
color2 = [217 83 25]/256;  color2L = [242 172 142]/256;
color3 = [237 177 32]/256; color3L = [249 227 176]/256;
colors = [color1;color2;color3];
colorsL = [color1L;color2L;color3L];


t_shift_video =  5.4;
t_shift_lab_view = 5.7;

%%
figure(1)
clf, hold on
for i = 1:3
    plot(T_vid,aN_vid(:,i),'linewidth',3,"color",colors(i,:),"DisplayName","a_"+string(i))
end
title("Experimental Data")
legend("location","best")
xlabel("Time - [sec]")
set(gca,'FontSize',14)
grid()


figure(2)
clf, hold on
for i = 1:3
    plot(T_anal,aN_anal(:,i),":",'linewidth',3,"color",colorsL(i,:),"DisplayName","a_"+string(i))
end
title("Numerical Data")
xlabel("Time - [sec]")
set(gca,'FontSize',14)
legend()
grid()

figure(3)
clf, hold on
for i = 1:3
    plot(T_vid,aN_vid(:,i),'linewidth',3,"color",colors(i,:),"DisplayName","a_"+string(i)+" - Experimental")
    plot(T_anal+t_shift_video,aN_anal(:,i),":",'linewidth',3,"color",colorsL(i,:),"DisplayName","a_"+string(i)+" - Analytical")
end
plot(t_shift_video*[1 1],ylim(),":","Color","k",'linewidth',4,"HandleVisibility","off")
% legend()
xlabel("Time - [sec]")
set(gca,'FontSize',14)
grid()


figure(4)
clf, hold on
plot(T_lab,F,"linewidth",3,"DisplayName","LabView Data")
plot(T_anal+t_shift_lab_view,Q_anal,"linewidth",3,"DisplayName","Analytical Data")
legend()
ylabel("Force - [N]")
xlabel("Time - [sec]")
set(gca,'FontSize',14)
grid()
