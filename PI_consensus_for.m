clear all; close all; clc; 

%% Graph Topology (Adjacency & Laplacian Matrices) 
Aij = [ 0, 1, 1, 1, 0; ...
        1, 0, 0, 1, 1; ...
        1, 0, 0, 0, 1; ...
        1, 1, 0, 0, 1; ...
        0, 1, 1, 1, 0];

Dij = diag(sum(Aij, 2));  % Degree matrix
Lij = Dij - Aij;          % Laplacian matrix
N = size(Lij, 1);         % Number of agents

Bij = [ 1, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0; ...
        0, 0, 0, 0, 0];

%% Parameters
a1 = -0.1;
a2 = -0.2;

k1 = 20;
k2 = 10;
k3 = 5;
uB = 5;  % Saturation limit

x_star = 10;

%% Initial Conditions
xL = 10;  % Leader's position
vL = 0;   % Leader's velocity

% xF = 5 * randn(1, N);  % Initial positions of followers
xF = randi([-60, 40], 1, N); % Initial positions of followers
vF = zeros(1, N);      % Initial velocities
uF = zeros(1, N);      % Control inputs

nF = zeros(N, 1);      % PI integral term

%% Simulation
ft = 50;             % Total simulation time
dt = 0.001;          % Euler discretization step
t_rec = 0:dt:ft;
num_steps = length(t_rec);

for index = 1:num_steps
    t = t_rec(index);
    
    z1 = zeros(1, N);
    z2 = zeros(1, N);
    
    for i = 1:N
        for j = 1:N
            z1(i) = z1(i) + Aij(i, j) * (xF(j) - xF(i)) + Bij(i, j) * (xL - xF(i));
            z2(i) = z2(i) + Aij(i, j) * (vF(j) - vF(i)) + Bij(i, j) * (vL - vF(i));
        end
        nF(i) = nF(i) + dt * (z1(i) + z2(i));  % Integral term update
        
        % Compute control input with saturation
        uF(i) = k1 * z1(i) + k2 * z2(i) + k3 * nF(i);
        uF(i) = max(min(uF(i), uB), -uB);  % Saturation function

        % System dynamics update
        vF(i) = vF(i) + dt * (a1 * xF(i) + a2 * vF(i) + uF(i));
        xF(i) = xF(i) + dt * vF(i);
    end
    
    % Store data
    xF_rec(:, index) = xF';
    vF_rec(:, index) = vF';
    u_rec(:, index)  = uF';
end

%% Plotting Results
figure;
% subplot(2,1,1);
plot(t_rec, xF_rec(1,:), 'r', 'LineWidth', 1); hold on;
plot(t_rec, xF_rec(2,:), 'g', 'LineWidth', 1);
plot(t_rec, xF_rec(3,:), 'b', 'LineWidth', 1);
plot(t_rec, xF_rec(4,:), 'm', 'LineWidth', 1);
plot(t_rec, xF_rec(5,:), 'c', 'LineWidth', 1);
grid on; title('Position of Agents');
legend('Agent 1', 'Agent 2', 'Agent 3', 'Agent 4', 'Agent 5');

figure;
% subplot(2,1,2);
plot(t_rec, vF_rec(1,:), 'r', 'LineWidth', 1); hold on;
plot(t_rec, vF_rec(2,:), 'g', 'LineWidth', 1);
plot(t_rec, vF_rec(3,:), 'b', 'LineWidth', 1);
plot(t_rec, vF_rec(4,:), 'm', 'LineWidth', 1);
plot(t_rec, vF_rec(5,:), 'c', 'LineWidth', 1);
grid on; title('Velocity of Agents');
legend('Agent 1', 'Agent 2', 'Agent 3', 'Agent 4', 'Agent 5');

figure;
plot(t_rec, u_rec(1,:), 'r', 'LineWidth', 1); hold on;
plot(t_rec, u_rec(2,:), 'g', 'LineWidth', 1);
plot(t_rec, u_rec(3,:), 'b', 'LineWidth', 1);
plot(t_rec, u_rec(4,:), 'm', 'LineWidth', 1);
plot(t_rec, u_rec(5,:), 'c', 'LineWidth', 1);
grid on; title('Control Inputs with Saturation');
legend('Agent 1', 'Agent 2', 'Agent 3', 'Agent 4', 'Agent 5');
