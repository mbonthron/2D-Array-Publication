function [Q1,a3,a3dot] = in_contact2(t,U,b,tbar,AAbar,beta,indentor_speed)
F = zeros(4,length(t));

%% === Values specific to arch
deltaL1=(b/2.0).^2;

a1 = U(1);
a2 = U(2);
a1dot = U(3);
a2dot = U(4);

    
starting_height = sqrt(12*(b/2)^2-tbar^2)/sqrt(3);

a3=a1-starting_height+t*indentor_speed;
a3dot=a1dot+indentor_speed;
QQQ1=0;


Va1=(beta/AAbar).*a1dot+(1/12*tbar^2)*a1    -1.*(deltaL1-0.25.*(a1.^2+4.*a2.^2+9.*a3.^2)).*a1;
Va2=(beta/AAbar).*a2dot+(1/12*tbar^2)*16.*a2-4.*(deltaL1-0.25.*(a1.^2+4.*a2.^2+9.*a3.^2)).*a2;
Va3=(beta/AAbar).*a3dot+(1/12*tbar^2)*81.*a3-9.*(deltaL1-0.25.*(a1.^2+4.*a2.^2+9.*a3.^2)).*a3;


%% First Order System
LHS = [1 0 0 1/AAbar ;
       0 1 0 0;
       0 0 1 -1/AAbar  ;
       1 0 -1 0];

RHS = -1*[Va1 ; Va2 ; Va3 ; 0];

sol = LHS \ RHS;

F(1,:)=a1dot;
F(2,:)=a2dot;
F(3,:)=sol(1);
F(4,:)=sol(2);

Q1 = sol(4);

end