function [F] = snap_through(t,U,b,tbar,AAbar,beta)
F = zeros(6,1);


a1 = U(1);
a2 = U(2);
a3 = U(3);
a1dot = U(4);
a2dot = U(5);
a3dot = U(6);

%% === Values specific to arch
deltaL1=(b/2.0)^2;

a1a=(beta/AAbar).*a1dot+(1/12*tbar^2)*a1    -1.*(deltaL1-0.25.*(a1.^2+4.*a2.^2+9.*a3.^2)).*a1;
a1b=(beta/AAbar).*a2dot+(1/12*tbar^2)*16.*a2-4.*(deltaL1-0.25.*(a1.^2+4.*a2.^2+9.*a3.^2)).*a2;
a1c=(beta/AAbar).*a3dot+(1/12*tbar^2)*81.*a3-9.*(deltaL1-0.25.*(a1.^2+4.*a2.^2+9.*a3.^2)).*a3;

%% First Order System
F(1)=U(4);
F(2)=U(5);
F(3)=U(6);
F(4)=-a1a;
F(5)=-a1b;
F(6)=-a1c;

end

