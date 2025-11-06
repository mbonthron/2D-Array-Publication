x = linspace(0,pi,500);

load("n=2.mat")
[w2,x] = w_from_A(A,2);

load("n=3.mat")
[w3,x] = w_from_A(A,3);

load("n=4.mat")
[w4,x] = w_from_A(A,4);

load("n=5.mat")
[w5,x] = w_from_A(A,5);

load("n=6.mat")
[w6,x] = w_from_A(A,6);

load("n=7.mat")
[w7,x] = w_from_A(A,7);

%%
error = zeros(1,7); error(1) = nan;
error(2) = calc_error(x,w2,w7);
error(3) = calc_error(x,w3,w7);
error(4) = calc_error(x,w4,w7);
error(5) = calc_error(x,w5,w7);
error(6) = calc_error(x,w6,w7);
error(7) = calc_error(x,w7,w7);



%%
i = 3;
f = figure(1); clf; hold on
f.Units = "inches";
f.Position(3:4) = [2 0.65];
set(gca,"FontSize",10);
xticks([0 pi]); xticklabels(["0" "$\pi$"])

plot([0 pi],[0 0],"k-")
plot(x,1*((1/pi)*w2(i,:)))
plot(x,1*((1/pi)*w7(i,:)))
ylim([-.175 .175])
%%
f = figure(2); clf; hold on
f.Units = "inches";
f.Position(3:4) = [2.5 1.5];
set(gca,"FontSize",10);
xlim([1 7]); xticks(1:7)
% ylim([0 10])

plot(1:7,error*100,"o-")
plot(xlim(),5*[1 1],"k-","LineWidth",1,"HandleVisibility","off")

xlabel("N")
ylabel("Error \%")