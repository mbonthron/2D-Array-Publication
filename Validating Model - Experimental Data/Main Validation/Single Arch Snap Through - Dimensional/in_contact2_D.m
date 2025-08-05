function [Q1,a3,a3dot] = in_contact2_D(T,U,eta,b,beta,indentor_speed,data)
F = zeros(4,length(T));

a1 = U(1);
a2 = U(2);
a1dot = U(3);
a2dot = U(4);

%% === Load in from data
rho = data.rho;
AA  = data.AA;
II  = data.II;
EE  = data.EE;
L   = data.L;

%% === Values specific to arch
deltaL1=(b*pi/2.0/L).^2;
    
a3   = (1/sin(3*eta*pi))*(b*sin(eta*pi) - T*indentor_speed - a1*sin(eta*pi) - a2*sin(2*eta*pi));
a3dot= (1/sin(3*eta*pi))*(              -   indentor_speed - a1dot*sin(eta*pi) - a2dot*sin(2*eta*pi));
QQQ1 = 0;

bending = EE*II*pi^4/L^4;
stretch = AA*EE*pi^2/L^2;


dVa1=beta.*a1dot+bending*    a1-stretch*   (deltaL1-0.25*pi^2/L^2.*(a1.^2+4.*a2.^2+9.*a3.^2)).*a1;
dVa2=beta.*a2dot+bending*16.*a2-stretch*4.*(deltaL1-0.25*pi^2/L^2.*(a1.^2+4.*a2.^2+9.*a3.^2)).*a2;
dVa3=beta.*a3dot+bending*81.*a3-stretch*9.*(deltaL1-0.25*pi^2/L^2.*(a1.^2+4.*a2.^2+9.*a3.^2)).*a3;


%% First Order System
if eta == .5
    m = AA*rho;
    LHS = [m 0 0 2/L ;
        0 m 0 0;
        0 0 m -2/L;
        1 0 -1 0];
else
    LHS = [1 0 0 2*sin(eta*pi)/L; 0 1 0  2*sin(2*eta*pi)/L; 0 0 1 2*sin(3*eta*pi)/L; sin((1:3)*eta*pi) 0];
end

sol = LHS \ (-1*[dVa1 ; dVa2 ; dVa3 ; 0]);
F(1,:) = a1dot;
F(2,:) = a2dot;
F(3,:) = sol(1);
F(4,:) = sol(2);

Q1 = sol(4);


end

