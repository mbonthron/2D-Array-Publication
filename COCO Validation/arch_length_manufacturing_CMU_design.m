%% 10-11-2023 MRB
%  Calculating Arch Dimensions
%  SEE ATTACHED PDF FOR DIAGRAM OF RELEVANT DIMENSIONS
syms x
%% Specify Constraints / Desired Dimensions
b = 13.7;      % [mm] Dimensional Arch Rise
L = 100;     % [mm] Dimensional Arch Length
t = 5;   % [mm] Dimensional Distance from centerline of axis to start of tab

%% Calculate ``ideal'' intiial length (if directly to centerline of axis)
L0_prime = int(sqrt(1+(pi/L*b*cos(pi*x/L))^2),x,0,L);

%% Account for the offset from the centerline of the axis
L0 = L0_prime - 2*t;

%% Print the Results:
fprintf("b = %.0f [mm]\n",b)
fprintf("L = %.0f [mm]\n",L)
fprintf("t = %.1f [mm]\n\n",t)
fprintf("L0_prime = %.2f [mm]\n",L0_prime)
fprintf("L0       = %.2f [mm]\n",L0)

%% Plot Shape of Arch
x = linspace(0,L,100);
plot(x,b*sin(pi*x/L),'linewidth',2);
axis equal

%% Plot What the Shape of the Arch should be
notch_length = 3;
notch_height = 1.52;
total_height = 11.25;
figure(2); clf; hold on; axis equal

plot(notch_length+[0 L0],total_height/2*[1 1],"Color","k")
plot(notch_length+[0 L0],-total_height/2*[1 1],"Color","k")

plot([0 0],[0 total_height/2-notch_height],"Color","k")
plot([0 0],[0 -total_height/2+notch_height],"Color","k")

plot([0 notch_length notch_length],[total_height/2-notch_height total_height/2-notch_height total_height/2],"Color","k")
plot([0 notch_length notch_length],-1*[total_height/2-notch_height total_height/2-notch_height total_height/2],"Color","k")

plot((2*notch_length+L0)*[1 1],[0 total_height/2-notch_height],"Color","k")
plot((2*notch_length+L0)*[1 1],[0 -total_height/2+notch_height],"Color","k")

plot((2*notch_length+L0)+[0 -notch_length -notch_length],[total_height/2-notch_height total_height/2-notch_height total_height/2],"Color","k")
plot((2*notch_length+L0)+[0 -notch_length -notch_length],-1*[total_height/2-notch_height total_height/2-notch_height total_height/2],"Color","k")