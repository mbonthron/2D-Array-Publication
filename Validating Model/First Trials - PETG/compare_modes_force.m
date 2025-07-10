%%
L = 100 / 1000;    % Length [m]
t = 0.75 / 1000;   % Thickness [m]
w = 14.18 / 1000;  % Width [m]
rise = 10 / 1000;  % Initial Rise [m]
EE = 4.3e9;        % Young's Modulus [N/m^2]

rho = 1270;        % Volumetric Density [kg/m^3]

indentor_speed_mms = 0.2;      % Indentor Speed mm/s
indentor_speed = (indentor_speed_mms/1000);

b = rise;

eta = 0.5;

II = 1/12 * w * t^3;
AA = t*w;
r = sqrt(II/AA);
T_end_sec = (1.5*rise*1000)/indentor_speed_mms;

beta_PRL = 0.25;
beta = beta_PRL * pi^2*sqrt(rho*AA*EE*II)/L^2;

data.rho = rho;
data.AA = AA;
data.II = II;
data.EE = EE;
data.L = L;

%% ============== Plot PRL
figure(1); clf; hold on

cd('PRL')
matfiles = dir('*.mat');
for i = 1:length(matfiles)
    load(matfiles(i).name)

    plot(T_combined_dimensional,U_combined_dimensional(:,1:3)*1000,'Color','r','LineWidth',3)

end
cd ..

tPRL = T_combined_dimensional;
for i = 1:length(T_combined_dimensional)
    fPRL(i) = in_contact2(T_combined_dimensional(i),U_combined_dimensional(i,:),eta,b,beta,indentor_speed,data);
end
%% ============== Plot Barenblatt
cd('Barenblatt')
matfiles = dir('*.mat');
for i = 1:length(matfiles)
    load(matfiles(i).name)

    plot(T_combined_dimensional,U_combined_dimensional(:,1:3)*1000,"--",'Color','g','LineWidth',2)
end
cd ..
tBB = T_combined_dimensional;

for i = 1:length(T_combined_dimensional)
    fBB(i) = in_contact2(T_combined_dimensional(i),U_combined_dimensional(i,:),eta,b,beta,indentor_speed,data);
end

%% ============== Plot Dimensional
cd('Dimensional')
matfiles = dir('*.mat');
for i = 1:length(matfiles)
    load(matfiles(i).name)

    plot(T_combined_dimensional,U_combined_dimensional(:,1:3)*1000,'Color','c','LineWidth',1)
end
cd ..
tdim = T_combined_dimensional;
for i = 1:length(T_combined_dimensional)
    fdim(i) = in_contact2(T_combined_dimensional(i),U_combined_dimensional(i,:),eta,b,beta,indentor_speed,data);
end
%% ========================================================================
%  Compare the force
figure(2); clf; hold on
plot(tdim,fdim,'LineWidth',2)
plot(tBB,fBB,'LineWidth',2)
plot(tPRL,fPRL,'LineWidth',2)
ylim([0 5])