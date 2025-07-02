t_vector = [0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.11 0.12]*pi;
b = 0.1*pi;
eta  = 0.5;

energy = {};

figure(1); clf
figure(2); clf

for i = 1:length(t_vector)
    thick = t_vector(i);
    file_string = "b = "+b/pi + "pi t = "+thick/pi+"pi eta = "+eta+".mat";
    load(file_string)

    QQ = Q(1:last_positive);
    xx = data.alpha*t1(1:last_positive);
    tt = thick*ones(last_positive,1);
    
    energy{i} = cumtrapz(xx,QQ);

    figure(1); hold on  
    plot3(tt,xx,QQ,'linewidth',2,"DisplayName","t = "+data.t/pi+"\pi")

    figure(2); hold on
    plot3(tt,xx,energy{i},'linewidth',2,"DisplayName","t = "+data.t/pi+"\pi")

end

%%
figure(1)
xlabel('t'); ylabel('Displacement'); zlabel('F')
legend(); grid()
view([45 45])


figure(2)
xlabel('t'); ylabel('Displacement'); zlabel('V')
legend(); grid()
view([45 45])