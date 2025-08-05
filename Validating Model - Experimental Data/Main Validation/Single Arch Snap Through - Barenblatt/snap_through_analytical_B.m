function [disp_dimensional, Rxn_dimensional_return] = snap_through_analytical_B(data)

%% Dimensional Parameters
% L = 100 / 1000;             % Length [m]
% thickness = 0.75 / 1000;    % Thickness [m]
% w = 14.18 / 1000;     % Width [m]
% rise = 10 / 1000;       % Initial Rise [m]
% EE = 2.3e9;                 % Young's Modulus [N/m^2]
% 
% rho = 1270;                 % Volumetric Density [kg/m^3]
% 
% indentor_speed_mms = 0.2;      % Indentor Speed mm/s

L = data.L;             % Length [m]
thickness = data.t;    % Thickness [m]
w = data.w;     % Width [m]
rise = data.rise;       % Initial Rise [m]
EE = data.EE;                 % Young's Modulus [N/m^2]
rho = data.rho;                 % Volumetric Density [kg/m^3]
indentor_speed_mms = data.indentor_speed_mms;      % Indentor Speed mm/s
beta_PRL = data.beta_PRL;   
eta = data.eta;   

% Derived parameters
II = 1/12 * w * thickness^3;
AA = thickness*w;
r = sqrt(II/AA);

% beta_PRL = 0.25;
beta_dimensional = beta_PRL * pi^2*sqrt(rho*AA*EE*II)/L^2;
beta = beta_dimensional*pi/(L*sqrt(rho*EE));

% eta = 0.5;


% Nondimensionalization
AAbar = AA * pi^2 / L^2;
IIbar = II * pi^4 / L^4;
tbar  = thickness * pi / L;
b     = rise * pi / L;

indentor_speed = (indentor_speed_mms/1000)*sqrt(rho/EE);

% Time Bounds
T_end_sec = (1.5*rise*1000)/indentor_speed_mms;
T_end  = T_end_sec/(L*sqrt(rho/EE)/pi);

% U2_PRL = 1e-6;
U2_PRL = data.U2_PRL;
U2_dimensional = U2_PRL*r;
U2_BB = U2_dimensional*pi/L;


%% Integrate the system
U = zeros(4,1);
U(1)=+sqrt(12*(b/2)^2-tbar^2)/sqrt(3);		    %a1a
U(2)=U2_BB;			%a1b
U(3)=0.0; 			%a1adot
U(4)=0.0; 			%a1bdot

%% Integrate when indentor in contact
options = odeset(RelTol=1e-11,AbsTol=1e-12);

tic
[T_contact,U_contact] = ode45(@(t,x0) in_contact_B(t,x0,eta,b,tbar,AAbar,beta,indentor_speed),linspace(0,T_end,50000) ,U,options);
toc

Rxn = zeros(size(T_contact));
a3 = zeros(size(T_contact));
a3dot = zeros(size(T_contact));

for i = 1:1:length(T_contact)
    [Q1,a1cVAL,a1cdotVAL] = in_contact2_B(T_contact(i),U_contact(i,:),eta,b,tbar,AAbar,beta,indentor_speed);
    Rxn(i) = Q1;
    a3(i) = a1cVAL;
    a3dot(i) = a1cdotVAL;
end

% Find the first point where reaction force is less than zero
snap_through_index = find(Rxn>0,1,'last')+1;
% snap_through_index = find(Rxn<0,1,'first')+1;

T_contact = T_contact(1:snap_through_index);
Rxn      = Rxn(1:snap_through_index);
U_contact = [U_contact(1:snap_through_index,1:2) a3(1:snap_through_index) U_contact(1:snap_through_index,3:4) a3dot(1:snap_through_index)] ;

% Switch to the free integration
[T_free,U_free] = ode45(@(t,x0) snap_through_B(t,x0,b,tbar,AAbar,beta),linspace(T_contact(end),T_end,500),U_contact(end,:),options);

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

Rxn_dimensional = Rxn*EE*L^2/pi^2;
disp_dimensional = (T_contact*indentor_speed)*L/pi;

energy = cumtrapz(Rxn_dimensional,disp_dimensional);

%%
T_combined_dimensional = T_combined*(L*sqrt(rho/EE)/pi);
U_combined_dimensional = U_combined*L/pi;

figure(10)
clf, hold on
plot(T_combined_dimensional,U_combined_dimensional(:,1)*1000,"DisplayName","a1")
plot(T_combined_dimensional,U_combined_dimensional(:,2)*1000,"DisplayName","a2")
plot(T_combined_dimensional,U_combined_dimensional(:,3)*1000,"DisplayName","a3")
legend()


%%
figure(4); hold on
plot(disp_dimensional*1000,Rxn_dimensional(1:snap_through_index),"LineWidth",3,'DisplayName',"BB eta = " + eta+ " t = "+ thickness*1000+ " b = "+ rise*1000)
xlabel('Displacement')
ylabel('Force')
legend()

Rxn_dimensional_return = Rxn_dimensional(1:snap_through_index);

save("b = "+data.rise*1000+"mm L =" +data.L*1000 +"mm EE = "+EE/(1e9)+"e9 beta = "+beta_PRL+" eta =" + eta + ".mat",'disp_dimensional','Rxn_dimensional','snap_through_index','T_combined_dimensional','U_combined_dimensional')

%% Dimensional Plots
sample_time = linspace(0,T_combined(end),50000);

Q_N = interp1(T_contact,Rxn,sample_time)/(2*L^3/pi^4/EE/II/r);
T_sec = sample_time/(pi^2/L^2*sqrt(EE*II/AA/rho));
aN_mm = interp1(T_combined,U_combined,sample_time)*r*1000;
% save("Analytical Prediction.mat","T_sec","Q_N","aN_mm")