close all; clear all; clc;

%% Graph Topology (Adjacency & Laplacian Matrices) 
A = [0 1 0 0 0 0 0 0 0 1;
     1 0 1 0 0 0 0 0 0 0;
     0 1 0 1 0 0 0 0 0 0;
     0 0 1 0 1 0 0 0 0 0;
     0 0 0 1 0 1 0 0 0 0;
     0 0 0 0 1 0 1 0 0 0;
     0 0 0 0 0 1 0 1 0 0;
     0 0 0 0 0 0 1 0 1 0;
     0 0 0 0 0 0 0 1 0 1;
     1 0 0 0 0 0 0 0 1 0];

D = diag(sum(A, 2));   % Degree matrix 
L = D - A;             % Laplacian matrix
N = size(L, 1);        % Number of agents

a0 = [1 1 1 1 1 0 0 0 0 0];  % leader connection
d0 = diag(a0);
L_bar = L + d0; 

%% Parameters
alpha = 5;
beta  = 5;

k = 5 * a0;
K = diag(k);

v_star = 1;

%% Nonlinear function
f = @(z) sign(z) .* min(abs(z), 1);

%% Initial Conditions
xL = randi([-50, 50], N, 1);   
vL = zeros(N, 1);     

xF = 100;
vF = 10;

eL = vL - v_star*ones(N, 1);
zL = xL - xF;
eF = vF - v_star;

%% Simulation
sim_time = 300;
dt = 0.001;

idx = 1;
for k = 0:dt:sim_time

    fz = f(zL);
    eL = eL + dt*(-K*eL - beta*L*eL - fz);
    zL = zL + dt*(-alpha*L*fz + eL - eF*ones(N, 1));
    vL = eL + v_star*ones(N, 1);
    uL = -alpha*L*fz + vL;

    xL = xL + dt*(uL);


    eF = vF - v_star;
    vF = vF + dt*(sum(fz));

    xF = xF + dt*(vF);

    % vL = eL + v_star*ones(N, 1);  
    % xL = zL + xF;
    % 
    % eF = eF + dt*(sum(fz));
    % 
    % vF = eF + v_star;
    % xF = xF + dt*vF;


    % Data recording
    t_rec(:, idx)  = k;
    xL_rec(:, idx) = xL;
    vL_rec(:, idx) = vL;
    xF_rec(:, idx) = xF;
    vF_rec(:, idx) = vF;
    fz_rec(:, idx) = fz;

    idx = idx + 1;
end

%% Plot the results of simulation
% 1. Position trajectory of leaders (solid line) and follower (dashed line)
fig1 = figure(1);
set(fig1, 'OuterPosition', [100, 700, 500, 300])
a1 = plot(t_rec, xL_rec(1,:), '-k', 'LineWidth', 1.5); hold on
a2 = plot(t_rec, xL_rec(2,:), '-k', 'LineWidth', 1.5); hold on
a3 = plot(t_rec, xL_rec(3,:), '-k', 'LineWidth', 1.5); hold on
a4 = plot(t_rec, xL_rec(4,:), '-k', 'LineWidth', 1.5); hold on
a5 = plot(t_rec, xL_rec(5,:), '-k', 'LineWidth', 1.5); hold on
a6 = plot(t_rec, xL_rec(6,:), '-k', 'LineWidth', 1.5); hold on
a7 = plot(t_rec, xL_rec(7,:), '-k', 'LineWidth', 1.5); hold on
a8 = plot(t_rec, xL_rec(8,:), '-k', 'LineWidth', 1.5); hold on
a9 = plot(t_rec, xL_rec(9,:), '-k', 'LineWidth', 1.5); hold on
a10 = plot(t_rec, xL_rec(10,:), '-k', 'LineWidth', 1.5); hold on
a11 = plot(t_rec, xF_rec(1,:), '-. r', 'LineWidth', 1.5); hold on
grid on
xlabel('Time (s)')
ylabel('$x_i$', 'Interpreter', 'latex')



% 2. Trajectory of agent group state v
fig2 = figure(2);
set(fig2, 'OuterPosition', [600, 700, 500, 300])
a1 = plot(t_rec, vL_rec(1,:), '-k', 'LineWidth', 1.5); hold on
a2 = plot(t_rec, vL_rec(2,:), '-k', 'LineWidth', 1.5); hold on
a3 = plot(t_rec, vL_rec(3,:), '-k', 'LineWidth', 1.5); hold on
a4 = plot(t_rec, vL_rec(4,:), '-k', 'LineWidth', 1.5); hold on
a5 = plot(t_rec, vL_rec(5,:), '-k', 'LineWidth', 1.5); hold on
a6 = plot(t_rec, vL_rec(6,:), '-k', 'LineWidth', 1.5); hold on
a7 = plot(t_rec, vL_rec(7,:), '-k', 'LineWidth', 1.5); hold on
a8 = plot(t_rec, vL_rec(8,:), '-k', 'LineWidth', 1.5); hold on
a9 = plot(t_rec, vL_rec(9,:), '-k', 'LineWidth', 1.5); hold on
a10 = plot(t_rec, vL_rec(10,:), '-k', 'LineWidth', 1.5); hold on
a11 = plot(t_rec, vF_rec(1,:), '-. r', 'LineWidth', 1.5); hold on
grid on
xlabel('Time (s)')
ylabel('$w_i, v_(N+1)$', 'Interpreter', 'latex')

% 3. Trajectory of agent group input u
fig3 = figure(3);
set(fig3, 'OuterPosition', [1100, 700, 500, 300])
a1 = plot(t_rec, fz_rec(1,:), '-k', 'LineWidth', 1.5); hold on
a2 = plot(t_rec, fz_rec(2,:), '-k', 'LineWidth', 1.5); hold on
a3 = plot(t_rec, fz_rec(3,:), '-k', 'LineWidth', 1.5); hold on
a4 = plot(t_rec, fz_rec(4,:), '-k', 'LineWidth', 1.5); hold on
a5 = plot(t_rec, fz_rec(5,:), '-k', 'LineWidth', 1.5); hold on
a6 = plot(t_rec, fz_rec(6,:), '-k', 'LineWidth', 1.5); hold on
a7 = plot(t_rec, fz_rec(7,:), '-k', 'LineWidth', 1.5); hold on
a8 = plot(t_rec, fz_rec(8,:), '-k', 'LineWidth', 1.5); hold on
a9 = plot(t_rec, fz_rec(9,:), '-k', 'LineWidth', 1.5); hold on
a10 = plot(t_rec, fz_rec(10,:), '-k', 'LineWidth', 1.5); hold on
grid on
xlabel('Time (s)')
ylabel('$f(z_i)$', 'Interpreter', 'latex')