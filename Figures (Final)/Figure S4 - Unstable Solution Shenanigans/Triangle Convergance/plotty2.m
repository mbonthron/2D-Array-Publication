V = zeros(3,7);
V(1,1) = nan;
V(2,1) = nan;
V(3,1) = nan;

for i = 2:7
    load("n="+string(i)+".mat","A","data")
    data.b_vector = 0.1*pi*[1 1 1]';
    V(:,i) = calculate_energy(A',data);
end

%%
error1 = ((V(1,:) - V(1,end))./V(1,end))*100;
error2 = ((V(2,:) - V(2,end))./V(2,end))*100;
error3 = ((V(3,:) - V(3,end))./V(3,end))*100;

f = figure(1); clf; hold on
f.Units = "inches";
f.Position(3:4) = [2.5 1.5];
set(gca,"FontSize",10);
xlim([1 7]); xticks(1:7)
ylim([-10 1])

plot(1:7,error1,"o-","LineWidth",1,"DisplayName","$V(w_1)$","MarkerSize",7)
plot(1:7,error2,"o-","LineWidth",1,"DisplayName","$V(w_2)$","MarkerSize",5)
plot(1:7,error3,"o-","LineWidth",1,"DisplayName","$V(w_3)$")
plot(xlim(),-5*[1 1],"k-","LineWidth",1,"HandleVisibility","off")


legend("Location","southeast")
xlabel("N")
ylabel("Error \%")