
function snap_through_analytical(material)
%material = 1; % Spring steel 1, PETG 2


%% Dimensional Parameters
L = 50 / 1000;             % Length [m]

if material == 1
    disp("####### Spring Steel #######")
    t = 0.35 / 1000;    % Thickness [m]
    w = 11.53 / 1000;     % Width [m]
    EE = 210e9;                 % Young's Modulus [N/m^2]
    rho = 7850;                 % Volumetric Density [kg/m^3]
    rise = 16.5 / 1000;       % Initial Rise [m]
    beta = .4;
elseif material == 2
    disp("####### PETG #######")
    t = 0.88 / 1000;    % Thickness [m]
    w = 14.18 / 1000;     % Width [m]
    EE = 8.9e9;                 % Young's Modulus [N/m^2]
    rho = 1270;                 % Volumetric Density [kg/m^3]
    rise = 13.7 / 1000;       % Initial Rise [m]
    beta = 1.4;

else 
    error("No valid material specified")
end


indentor_speed_mms = 0.2;      % Indentor Speed mm/s


% Derived parameters
II = 1/12 * w * t^3;
AA = t*w;
r = sqrt(II/AA);
T_end_sec = (1.5*rise*1000)/indentor_speed_mms;


%% Integrate the system
opts = odeset('RelTol',1e-5,'AbsTol',1e-5);


b = rise/r;
T_end = T_end_sec * (pi^2/L^2*sqrt(EE*II/AA/rho));
indentor_speed = (indentor_speed_mms/1000)/r / (pi^2/L^2*sqrt(EE*II/AA/rho));

U = zeros(4,1);
U(1)=+b;		    %a1a
U(2)=1e-3;			%a1b
U(3)=0.0; 			%a1adot
U(4)=0.0; 			%a1bdot

%% Integrate when indentor in contact
[T_contact,U_contact] = ode45(@(t,x0) in_contact(t,x0,b,beta,indentor_speed),[0 T_end],U,opts);

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
snap_through_index = find(Rxn<0,20);
snap_through_index = snap_through_index(end);

T_contact = T_contact(1:snap_through_index);
Rxn      = Rxn(1:snap_through_index);
U_contact = [U_contact(1:snap_through_index,1:2) a3(1:snap_through_index) U_contact(1:snap_through_index,3:4) a3dot(1:snap_through_index)] ;

% Switch to the free integration
[T_free,U_free] = ode45(@(t,x0) snap_through(t,x0,b,beta),[T_contact(end) T_end],U_contact(end,:),opts);

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

%% Dimensional Plots
sample_time = linspace(0,T_combined(end),50000);

Q_N = interp1(T_contact,Rxn,sample_time)/(2*L^3/pi^4/EE/II/r);
T_sec = sample_time/(pi^2/L^2*sqrt(EE*II/AA/rho));
aN_mm = interp1(T_combined,U_combined,sample_time)*r*1000;
save("Analytical Prediction "+num2str(material) +".mat","T_sec","Q_N","aN_mm")

end