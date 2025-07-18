%% Dimensional Parameters
L = 100 / 1000;             % Length [m]
t = .3048 / 1000;    % Thickness [m]
w = 11.53 / 1000;     % Width [m]
rise = 16 / 1000;       % Initial Rise [m]
EE = 210e9;                 % Young's Modulus [N/m^2]
% EE = 180e9;

rho = 7930;                 % Volumetric Density [kg/m^3]

indentor_speed_mms = 0.2;      % Indentor Speed mm/s

% beta = 0.65;
beta = 0.4;

% Derived parameters
II = 1/12 * w * t^3;
AA = t*w;
r = sqrt(II/AA);
T_end_sec = (1.5*rise*1000)/indentor_speed_mms;


%% Integrate the system
% opts = odeset('RelTol',1e-3,'AbsTol',1e-3);
% opts = odeset('MaxStep',100,'MinStep',10)

b = rise/r;
T_end = T_end_sec * (pi^2/L^2*sqrt(EE*II/AA/rho));
indentor_speed = (indentor_speed_mms/1000)/r / (pi^2/L^2*sqrt(EE*II/AA/rho));

U = zeros(4,1);
U(1)=+b;		    %a1a
U(2)=1e-6;			%a1b
U(3)=0.0; 			%a1adot
U(4)=0.0; 			%a1bdot

%% Integrate when indentor in contact
tic
[T_contact,U_contact] = ode45(@(t,x0) in_contact(t,x0,b,beta,indentor_speed),linspace(0,T_end,5000),U);
toc

Rxn = zeros(size(T_contact));
a3 = zeros(size(T_contact));
a3dot = zeros(size(T_contact));

for i = 1:1:length(T_contact)
    [Q1,a1cVAL,a1cdotVAL] = in_contact2(T_contact(i),U_contact(i,:),b,beta,indentor_speed);
    Rxn(i) = Q1;
    a3(i) = a1cVAL;
    a3dot(i) = a1cdotVAL;
end

% Find the first point where reaction force is less than zero
snap_through_index = find(Rxn>0,1,'last');

T_contact = T_contact(1:snap_through_index);
Rxn      = Rxn(1:snap_through_index);
U_contact = [U_contact(1:snap_through_index,1:2) a3(1:snap_through_index) U_contact(1:snap_through_index,3:4) a3dot(1:snap_through_index)] ;

% Switch to the free integration
[T_free,U_free] = ode45(@(t,x0) snap_through(t,x0,b,beta),linspace(T_contact(end),T_end,500),U_contact(end,:));

%% Dimensionless Plots
T_combined = [T_contact; T_free(2:end,:)];
U_combined = [U_contact; U_free(2:end,:)];

figure(1)
clf, hold on
plot(T_combined,U_combined(:,1),"DisplayName","a1")
plot(T_combined,U_combined(:,2),"DisplayName","a2")
plot(T_combined,U_combined(:,3),"DisplayName","a3")
legend()

figure(2)
plot(T_combined,U_combined(:,1)-U_combined(:,3),"DisplayName","Midpoint Position")
legend()

figure(3)
plot(T_contact,Rxn(1:snap_through_index),"DisplayName","Reaction Force")
legend()

Rxn_dimensional = Rxn*pi^4*EE*II*r / (2*L^3);
disp_dimensional = (T_contact*indentor_speed)*r;

energy = cumtrapz(Rxn_dimensional,disp_dimensional);

%%
figure(4); hold on
plot(disp_dimensional*1000,Rxn_dimensional(1:snap_through_index),":","LineWidth",3)
xlabel('Displacement')
ylabel('Force')
legend()

save("EE = "+EE/(1e9)+"e9 beta= "+beta+".mat",'disp_dimensional','Rxn_dimensional','snap_through_index')

%% Dimensional Plots
sample_time = linspace(0,T_combined(end),50000);

Q_N = interp1(T_contact,Rxn,sample_time)/(2*L^3/pi^4/EE/II/r);
T_sec = sample_time/(pi^2/L^2*sqrt(EE*II/AA/rho));
aN_mm = interp1(T_combined,U_combined,sample_time)*r*1000;
% save("Analytical Prediction.mat","T_sec","Q_N","aN_mm")