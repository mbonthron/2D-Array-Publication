a1 = zeros(1,6);
a2 = zeros(1,6);
a3 = zeros(1,6);
a4 = zeros(1,6);
a5 = zeros(1,6);
a6 = zeros(1,6);
a7 = zeros(1,6);


load("n=2.mat","A")
a1(1) = A(1);
a2(1) = A(2);
a3(1) = nan;
a4(1) = nan;
a5(1) = nan;
a6(1) = nan;
a7(1) = nan;

load("n=3.mat","A")
a1(2) = A(1);
a2(2) = A(2);
a3(2) = A(3);
a4(2) = nan;
a5(2) = nan;
a6(2) = nan;
a7(2) = nan;

load("n=4.mat","A")
a1(3) = A(1);
a2(3) = A(2);
a3(3) = A(3);
a4(3) = A(4);
a5(3) = nan;
a6(3) = nan;
a7(3) = nan;

load("n=5.mat","A")
a1(4) = A(1);
a2(4) = A(2);
a3(4) = A(3);
a4(4) = A(4);
a5(4) = A(5);
a6(4) = nan;
a7(4) = nan;

load("n=6.mat","A")
a1(5) = A(1);
a2(5) = A(2);
a3(5) = A(3);
a4(5) = A(4);
a5(5) = A(5);
a6(5) = A(6);
a7(5) = nan;

load("n=7.mat","A")
a1(6) = A(1);
a2(6) = A(2);
a3(6) = A(3);
a4(6) = A(4);
a5(6) = A(5);
a6(6) = A(6);
a7(6) = A(7);

%%
f = figure(1); clf; hold on
f.Units = "inches";
f.Position(3:4) = [2 1];
set(gca,'FontSize',10)
xlim([1 7])

n = 2:7;


mode = a2;
error = 100*((mode-mode(end))./mode(end));
plot(n,error)

mode = a4;
error = 100*((mode-mode(end))./mode(end));
plot(n,error)


ylim([0 10])



