x = linspace(0,4*pi,100);
y = sin(x);

figure(1); clf; hold on
plot(x,y)
scatter(0:pi:4*pi,zeros(1,5),100,"filled","MarkerFaceColor","k")

axis off