figure(1); clf; hold on
for i = 1:5
    plot(position_t_05{i},force_t_05{i})
end
set(gca,'FontSize',16)
xlabel('Displacement - [mm]')
ylabel('Force - [N]')

figure(2); clf; hold on
for i = 1:5
    plot(position_t_10{i},force_t_10{i})
end
set(gca,'FontSize',16)
xlabel('Displacement - [mm]')
ylabel('Force - [N]')

figure(3); clf; hold on
for i = 1:5
    plot(position_t_15{i},force_t_15{i})
end
set(gca,'FontSize',16)
xlabel('Displacement - [mm]')
ylabel('Force - [N]')

%%
N = 3;
color1 = [230 65 0]/255;
color2 = [11 22 42]/255;
colormat = color1 + (color2 - color1).*linspace(0,1,N).';

figure(4); clf; hold on
plot(energy_t_05{1},'LineWidth',2,'Color',colormat(1,:),"DisplayName","t = 0.05 \pi")
plot(energy_t_10{1},'LineWidth',2,'Color',colormat(2,:),"DisplayName","t = 0.10 \pi")
plot(energy_t_15{1},'LineWidth',2,'Color',colormat(3,:),"DisplayName","t = 0.15 \pi")

for k = 2:5
    plot(energy_t_05{k},'LineWidth',2,'Color',colormat(1,:),"HandleVisibility","off")
    plot(energy_t_10{k},'LineWidth',2,'Color',colormat(2,:),"HandleVisibility","off")
    plot(energy_t_15{k},'LineWidth',2,'Color',colormat(3,:),"HandleVisibility","off")
end

legend('Location','northwest')
set(gca,'FontSize',16)
xlabel('index')
ylabel('V - [Nmm]')

%%
figure(5); clf; hold on
for k = 1:5
    plot(energy_t_10{k},'LineWidth',2)
end