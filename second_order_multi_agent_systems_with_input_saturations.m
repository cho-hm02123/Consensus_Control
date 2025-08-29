%% Display initialization
close all
clear
clc

%% Parameters
% Input staturation
N = 5;

% Topology
A = [0 1 1 1 0;
     1 0 0 1 1;
     1 0 0 0 1;
     1 1 0 0 1;
     0 1 1 1 0
    ];

% Known target
Alpa0 = [1 0 0 0 0];

% Target
x_star = 10*ones(N,1);

% System parameter
a1 = -0.1;
a2 = -0.2;
u_bar = 5;

% Control gain
k1 = 20;
k2 = 10;
k3 = 5;

%% Laplacian Matrix
D = sum(A,2);
D = diag(D);
L = D - A;                                                                 % Laplacian Matrix
L_bar = L + diag(Alpa0);

%% Initial Vector
x = randi([-60,40], N, 1);
e = zeros(N,1);
v = zeros(N,1);
eta = zeros(N,1);

%% Simulation
% Simulation Time and Rate
dt = 0.01;

idx = 1;
for k = 0:dt:50
    e = x - x_star;
    e_dot = v;
    eta = -L_bar*e - L_bar*v;

    u = -k1*L_bar*e -k2*L_bar*v +k3*eta;
    sat = max(min(u, u_bar), -u_bar);

    v_dot = a1*e + a2*v + sat;
    v = v + dt*v_dot;
    x = x + dt*v;                                                          % x_dot = v

    % Data recording
    x_rec(:,idx) = x;
    e_rec(:,idx) = e;
    v_rec(:,idx) = v;
    u_rec(:,idx)= sat;
    t_rec(:,idx) = k;

    idx = idx + 1;
end

%% Plot the results of simulation
% 1. Trajectory of agent group state x
fig1 = figure(1);
set(fig1, 'OuterPosition', [100, 700, 500, 300])
a1 = plot(t_rec, x_rec(1,:), '-r', 'LineWidth', 1.5); hold on
a2 = plot(t_rec, x_rec(2,:), '-g', 'LineWidth', 1.5); hold on
a3 = plot(t_rec, x_rec(3,:), '-b', 'LineWidth', 1.5); hold on
a4 = plot(t_rec, x_rec(4,:), '-k', 'LineWidth', 1.5); hold on
a5 = plot(t_rec, x_rec(5,:), '-y', 'LineWidth', 1.5); hold on
grid on
xlabel('Time (s)')
ylabel('$x_i, i = 1, ..., 5$', 'Interpreter', 'latex')
legend([a1, a2, a3, a4, a5], {'Agent1', 'Agent2', 'Agent3', 'Agent4', ...
    'Agent5'}, 'Interpreter', 'latex')


% 2. Trajectory of agent group state v
fig2 = figure(2);
set(fig2, 'OuterPosition', [600, 700, 500, 300])
a1 = plot(t_rec, v_rec(1,:), '-r', 'LineWidth', 1.5); hold on
a2 = plot(t_rec, v_rec(2,:), '-g', 'LineWidth', 1.5); hold on
a3 = plot(t_rec, v_rec(3,:), '-b', 'LineWidth', 1.5); hold on
a4 = plot(t_rec, v_rec(4,:), '-k', 'LineWidth', 1.5); hold on
a5 = plot(t_rec, v_rec(5,:), '-y', 'LineWidth', 1.5); hold on
grid on
xlabel('Time (s)')
ylabel('$v_i, i = 1, ..., 5$', 'Interpreter', 'latex')
legend([a1, a2, a3, a4, a5], {'Agent1', 'Agent2', 'Agent3', 'Agent4', ...
    'Agent5'}, 'Interpreter', 'latex')

% 3. Trajectory of agent group input u
fig3 = figure(3);
set(fig3, 'OuterPosition', [1100, 700, 500, 300])
a1 = plot(t_rec, u_rec(1,:), '-r', 'LineWidth', 1.5); hold on
a2 = plot(t_rec, u_rec(2,:), '-g', 'LineWidth', 1.5); hold on
a3 = plot(t_rec, u_rec(3,:), '-b', 'LineWidth', 1.5); hold on
a4 = plot(t_rec, u_rec(4,:), '-k', 'LineWidth', 1.5); hold on
a5 = plot(t_rec, u_rec(5,:), '-y', 'LineWidth', 1.5); hold on
grid on
xlabel('Time (s)')
ylabel('$sat(u_i), i = 1, ..., 5$', 'Interpreter', 'latex')
legend([a1, a2, a3, a4, a5], {'Agent1', 'Agent2', 'Agent3', 'Agent4', ...
    'Agent5'}, 'Interpreter', 'latex')