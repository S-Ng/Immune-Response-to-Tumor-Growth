% Immune Response to Tumor Growth
% 20/11/2019
% Get populations applying 1 model variation through entire time

clear all;
clc;

%P(1) = T : tumor cells
%P(2) = N : natural killer cells
%P(3) = L : CD8+ T cells
% variable w in MATLAB is variable p in dN/dt in paper

model = 0; % mouse = 0, human = 1
variation = 4; % nn = 1, nl = 2, ln = 3, ll = 4
[a, b, sigma, f, h, w, m, k ,q, r, c, g, d, lambda, s, j] = getParameters(model, variation);

D = @(t, p) d*(p(3)/p(1))^lambda / (s + (p(3)/p(1))^lambda) * p(1); % functional form for (CD8+ T)-tumor kill term
longfunctionname = @(t, p) [a*p(1)*(1-b*p(1)) - c * p(2)*p(1) - D(t,p) ; % dT/dt
    sigma - f*p(2) + (g*p(1)^2 / (h+p(1)^2)) * p(2) - w*p(2)*p(1) ; % dN/dt
    -m*p(3) + j*D(t,p)^2 / (k+D(t,p)^2) * p(3) - q*p(3)*p(1) + r*p(2)*p(1) ]; % dL/dt

% Initial conditions
T0 = 5*10^3;
N0 = 10^3;
L0 = 1;

% Solve first regime
time = linspace(0, 35, 1000);
[t, sol] = ode45(longfunctionname, time, [T0, N0, L0]);

plottype = 1; % Separate plots: 0 , Overlaid plots: 1
population_plot(t, sol, plottype, true, 'BasicModel_ll')