% Immune Response to Tumor Growth with CAR T cell modifications
% 09/12/2019

clear all; clc;

%P(1) = T : tumor cells
%P(2) = N : natural killer cells
%P(3) = L : CD8+ T cells
% variable w in MATLAB is variable p in dN/dt in paper

model = 0; % mouse = 0, human = 1
variation = 1; % nn = 1, nl = 2, ln = 3, ll = 4
ttl = 'Tumor response to immunotherapy (nn)';
[a, b, sigma, f, h, w, m, k ,q, r, c, g, d, lambda, s, j] = getParameters(model, variation);

% Base equations
D = @(t, p) d*((p(3)/p(1))^lambda) / (s + (p(3)/p(1))^lambda) * p(1); % functional form for (CD8+ T)-tumor kill term
ODEs = @(t, p) [a*p(1)*(1-b*p(1)) - c * p(2)*p(1) - D(t,p) ; % dT/dt
    sigma - f*p(2) + (g*p(1)^2 / (h+p(1)^2)) * p(2) - w*p(2)*p(1) ; % dN/dt
    -m*p(3) + j*D(t,p)^2 / (k+D(t,p)^2) * p(3) - q*p(3)*p(1) + r*p(2)*p(1) ]; % dL/dt

% Initial conditions
T0 = 5*10^3;
N0 = 5*10^2;
L0 = 10;

time = linspace(0, 10, 100); % initial time for 10 days
[t1, sol1] = ode45(ODEs, time, [T0, N0, L0]); % get solution for first time span

% adjust relevant parameters to CAR T cells (only in D)
d = 5; % 7.17 maximum originally in ll 

% Redefine base equations for updated d value
D = @(t, p) d*((p(3)/p(1))^lambda) / (s + (p(3)/p(1))^lambda) * p(1); % functional form for (CD8+ T)-tumor kill term
ODEs2 = @(t, p) [a*p(1)*(1-b*p(1)) - c * p(2)*p(1) - D(t,p) ; % dT/dt
    sigma - f*p(2) + (g*p(1)^2 / (h+p(1)^2)) * p(2) - w*p(2)*p(1) ; % dN/dt
    -m*p(3) + j*D(t,p)^2 / (k+D(t,p)^2) * p(3) - q*p(3)*p(1) + r*p(2)*p(1) ]; % dL/dt

%Get final values from first regime to be initial conditions for second
T02 = sol1(length(time),1);
N02 = sol1(length(time),2);
L02 = sol1(length(time),3);
CART_input = 5*10^4; % CAR T cell injection
Tchallenge = 0; % tumor challenge

time2 = linspace(10, 35, 250); % day 10-35
[t2, sol2] = ode45(ODEs2, time2, [T02+Tchallenge, N02, L02+CART_input]); % get solution for second time span

% Concatenate regimes
sol = [sol1 ; sol2];
t = [t1 ; t2];
for i = 1:length(t)
    if sol(i, 1) < 1
        sol(i, 1) = 1;
    end
end

% Plot
plottype = 1;
population_plot(t, sol, plottype, ttl)
ylim([10^2,10^9]);