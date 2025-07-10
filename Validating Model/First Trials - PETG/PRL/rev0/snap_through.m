function [F] = snap_through(T,U,b,beta)
F = zeros(6,1);

%% === Values specific to arch
deltaL1=1+(+b/2.0)^2;
a1a=beta*U(4)+U(1)-(deltaL1-0.25*(U(1)^2+4*U(2)^2+9*U(3)^2))*U(1);
a1b=beta*U(5)+16*U(2)-4*(deltaL1-0.25*(U(1)^2+4*U(2)^2+9*U(3)^2))*U(2);
a1c=beta*U(6)+81*U(3)-9*(deltaL1-0.25*(U(1)^2+4*U(2)^2+9*U(3)^2))*U(3);


%% First Order System
F(1)=U(4);
F(2)=U(5);
F(3)=U(6);
F(4)=-a1a;
F(5)=-a1b;
F(6)=-a1c;

end

