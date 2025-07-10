function [F] = in_contact(T,U,b,beta,indentor_speed)
F = zeros(4,length(T));

%% === Values specific to arch
deltaL1=1+(+b/2.0).^2;
% a1cVAL=U(1)-b.*cos(omn01.*T);
% a1cdotVAL=U(3)+b.*omn01.*sin(omn01.*T);
% QQQ1=+b.*omn01.^2.*cos(omn01.*T);
    
a1cVAL=U(1)-b+T*indentor_speed;
a1cdotVAL=U(3)+indentor_speed;
QQQ1=0;


a1a=beta.*U(3)+U(1)-(deltaL1-0.25.*(U(1).^2+4.*U(2).^2+9.*a1cVAL.^2)).*U(1);
a1b=beta.*U(4)+16.*U(2)-4.*(deltaL1-0.25.*(U(1).^2+4.*U(2).^2+9.*a1cVAL.^2)).*U(2);
a1c=beta.*a1cdotVAL+81.*a1cVAL-9.*(deltaL1-0.25.*(U(1).^2+4.*U(2).^2+9.*a1cVAL.^2)).*a1cVAL;


%% First Order System
F(1,:)=U(3);
F(2,:)=U(4);
F(3,:)=-0.5.*QQQ1-0.5.*a1a-0.5.*a1c;
F(4,:)=-a1b;

Q1 = (0.5.*QQQ1 - 0.5.*a1a + 0.5.*a1c);

end

