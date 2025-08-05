function [F] = snap_through(T,U,b,beta,data)
F = zeros(6,1);

a1 = U(1);
a2 = U(2);
a3 = U(3);
a1dot = U(4);
a2dot = U(5);
a3dot = U(6);

%% === Load in from data
rho = data.rho;
AA  = data.AA;
II  = data.II;
EE  = data.EE;
L   = data.L;

%% === Values specific to arch
deltaL1=(b*pi/2.0/L).^2;

bending = EE*II*pi^4/L^4;
stretch = AA*EE*pi^2/L^2;

dVa1=beta.*a1dot+bending*    a1-stretch*   (deltaL1-0.25*pi^2/L^2.*(a1.^2+4.*a2.^2+9.*a3.^2)).*a1;
dVa2=beta.*a2dot+bending*16.*a2-stretch*4.*(deltaL1-0.25*pi^2/L^2.*(a1.^2+4.*a2.^2+9.*a3.^2)).*a2;
dVa3=beta.*a3dot+bending*81.*a3-stretch*9.*(deltaL1-0.25*pi^2/L^2.*(a1.^2+4.*a2.^2+9.*a3.^2)).*a3;


%% First Order System
F(1)=a1dot;
F(2)=a2dot;
F(3)=a3dot;
F(4)=-dVa1;
F(5)=-dVa2;
F(6)=-dVa3;

end

