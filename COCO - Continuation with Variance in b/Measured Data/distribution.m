% Load the data
L0_vector = readmatrix('b = 10 mm L = 100 mm.csv');

% remove the four largest values (they seem like outliers)
[~, idx] = maxk(L0_vector, 4);

% Remove those elements
L0_vector(idx) = [];

%
b_vector = zeros(size(L0_vector));

% Plot the histogram
figure(1); clf
histogram(L0_vector,15)
xlabel('Length [mm]')
% Physical values from design
t = 10;   % [mm] Dimensional Distance from centerline of axis to start of tab
L = 100;  % [mm] Centerline distance between the hinges

%%
for i = 1:length(L0_vector)
    L0 = L0_vector(i);
    L0_prime = L0 + 2*t;
    if L0_prime > L
        b = fzero(@(b) integral(@(x) sqrt(1+(pi/L*b*cos(pi*x/L)).^2),0,L)-L0_prime,1);
        b_vector(i) = abs(b);
    end
end

%%
figure(2); clf; hold on
histogram(b_vector,15)
xlabel('Rise [mm]')

%%
mu = mean(b_vector)
sigma = std(b_vector)

%% Overlay the fit
x = linspace(8.5,10.5,200);
y = 1/(sigma*sqrt(2*pi))*exp(-(x-mu).^2/(2*sigma^2));

plot(x,5.5*y,'LineWidth',2)