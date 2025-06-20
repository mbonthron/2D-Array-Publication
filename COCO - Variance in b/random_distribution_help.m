%% Generate the Distribution
mu = 1;         % Mean of distribution
sigma = 0.05;   % Standard deviation

rng('default')  % For reproducibility
R = normrnd(mu,sigma,[1,1000]);

%% Plot the distribution
figure(1); clf;
hist(R,50)
xlim([0.5 1.5])