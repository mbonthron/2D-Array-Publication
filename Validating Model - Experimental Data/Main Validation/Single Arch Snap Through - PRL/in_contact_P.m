function [F] = in_contact_P(T,U,eta,b,beta,indentor_speed)
F = zeros(4,length(T));

a1 = U(1);
a2 = U(2);
a1dot = U(3);
a2dot = U(4);

%% === Values specific to arch
deltaL1=1+(+b/2.0).^2;
    
a3   = (1/sin(3*eta*pi))*(b*sin(eta*pi) - T*indentor_speed - a1*sin(eta*pi) - a2*sin(2*eta*pi));
a3dot= (1/sin(3*eta*pi))*(              -   indentor_speed - a1dot*sin(eta*pi) - a2dot*sin(2*eta*pi));
QQQ1 = 0;


dVa1=beta.*a1dot+    a1-   (deltaL1-0.25.*(a1.^2+4.*a2.^2+9.*a3.^2)).*a1;
dVa2=beta.*a2dot+16.*a2-4.*(deltaL1-0.25.*(a1.^2+4.*a2.^2+9.*a3.^2)).*a2;
dVa3=beta.*a3dot+81.*a3-9.*(deltaL1-0.25.*(a1.^2+4.*a2.^2+9.*a3.^2)).*a3;


%% First Order System
if eta == .5
    LHS = [1 0 0 1 ; 0 1 0 0 ; 0 0 1 -1; 1 0 -1 0];
else
    LHS = [1 0 0 sin(eta*pi); 0 1 0  sin(2*eta*pi); 0 0 1 sin(3*eta*pi); sin((1:3)*eta*pi) 0];
end

sol = LHS \ (-1*[dVa1 ; dVa2 ; dVa3 ; 0]);
F(1,:)=a1dot;
F(2,:)=a2dot;
F(3,:)=sol(1);
F(4,:)=sol(2);

Q1 = sol(4);

end

