close all; clear all; clc;

%% Graph Topology (Adjacency & Laplacian Matrices)              
A = [ 0, 1, 1, 1, 0; ...
      1, 0, 0, 1, 1; ...
      1, 0, 0, 0, 1; ...
      1, 1, 0, 0, 1; ...
      0, 1, 1, 1, 0];

D = diag(sum(A, 2));  % Degree matrix
L = D - A;            % Laplacian matrix
N = size(L, 1);       % Number of agents

B = [1, 0, 0, 0, 0];  % leader connection
B = diag(B);
L_bar = L + B; 

%% Parameters
a1 = -0.1;
a2 = -0.2;

k1 = 20;
k2 = 10;
k3 = 5;
uB = 5;  % Saturation limit

x_star = 10;

%% Initial Conditions
xF = randi([-60, 40], N, 1); % Initial positions of followers
vF = zeros(N, 1);            % Initial velocities
nF = zeros(N, 1);            % Integral term for PI control    
eF = xF - x_star;
uF = 0;

d  = a1*x_star*ones(N, 1);

%% Simulation 
ft = 50;             % Total simulation time
dt = 0.001;          % Discretization step size
t_rec = 0:dt:ft;
num_steps = length(t_rec);

for k = 1:num_steps
    % Compute consensus errors using Laplacian & Integral term update
    vF = vF + dt*(a1*eF + a2*vF + uF + d);
    eF = eF + dt*vF;
    nF = -L_bar*eF - L_bar*vF; 
    
    % Compute control input with saturation
    uF = -k1*L_bar*eF - k2*L_bar*vF + k3*nF;
    uF = max(min(uF, uB), -uB);  % Saturation function

    % System dynamics update
    xF = xF + dt*vF;
    
    % Data logging
    xF_rec(:, k) = xF;
    vF_rec(:, k) = vF;
    uF_rec(:, k) = uF;
end

%% Plot Results
fig1 = figure(1);
hold on;
p0 = plot(t_rec, xF_rec(1, :), 'r', 'LineWidth', 1);
p1 = plot(t_rec, xF_rec(2, :), 'g', 'LineWidth', 1);
p2 = plot(t_rec, xF_rec(3, :), 'b', 'LineWidth', 1);
p3 = plot(t_rec, xF_rec(4, :), 'm', 'LineWidth', 1);
p4 = plot(t_rec, xF_rec(5, :), 'c', 'LineWidth', 1);
yline(10, '--k', 'LineWidth',1); 
grid on; 
xlabel('Time (s)');
ylabel('x_i,i=1,...5');
title('Position of Agents');
legend('Agent 1', 'Agent 2', 'Agent 3', 'Agent 4', 'Agent 5');

%%
fig2 = figure(2);
hold on;
a0 = plot(t_rec, vF_rec(1, :), 'r', 'LineWidth', 1);
a1 = plot(t_rec, vF_rec(2, :), 'g', 'LineWidth', 1);
a2 = plot(t_rec, vF_rec(3, :), 'b', 'LineWidth', 1);
a3 = plot(t_rec, vF_rec(4, :), 'm', 'LineWidth', 1);
a4 = plot(t_rec, vF_rec(5, :), 'c', 'LineWidth', 1);
grid on; 
xlabel('time (s)')
ylabel('v_i, i= = 1,...,5');
title('Velocity of Agents');
legend('Agent 1', 'Agent 2', 'Agent 3', 'Agent 4', 'Agent 5');

%%
fig3 = figure(3);
hold on;
b0 = plot(t_rec, uF_rec(1, :), 'r', 'LineWidth', 1);
b1 = plot(t_rec, uF_rec(2, :), 'g', 'LineWidth', 1);
b2 = plot(t_rec, uF_rec(3, :), 'b', 'LineWidth', 1);
b3 = plot(t_rec, uF_rec(4, :), 'm', 'LineWidth', 1);
b4 = plot(t_rec, uF_rec(5, :), 'c', 'LineWidth', 1);
yline(-5, '--k', 'LineWidth',1); 
yline(5, '--k', 'LineWidth', 1);  
grid on; 
xlabel('time (s)')
ylabel('sat(u_i), i = 1,...,5');
title('Control Inputs with Saturation');
legend('Agent 1', 'Agent 2', 'Agent 3', 'Agent 4', 'Agent 5');
ylim([-6, 6]);