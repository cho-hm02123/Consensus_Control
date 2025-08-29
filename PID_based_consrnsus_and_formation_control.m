%% Display initialization
close all
clear
clc

%% Run
% If you want Fig2, Set Fig1 = 1, Velo = 0
% If you want Fig3 and Fig 4, Set Fig1 = 1, Velo = 1
% If you want Fig6 and Fig 7, Set Fig1 = 0, Velo = 1
Fig1 = 0;
Velo = 1;

%% Parameters
% Input staturation
N = 5;

% Topology
if Fig1
    Ax = [0 1 0 1 0;
          1 0 1 0 1;
          0 1 0 0 1;
          1 0 0 0 0;
          0 1 1 0 0];

    Av = [0 1 0 1 0;
          1 0 0 0 0;
          0 0 0 0 1;
          1 0 0 0 0;
          0 0 1 0 0];

    Gx = [1 0 1 0 0];

    Gv = [0 1 0 0 0];
end

if ~Fig1
    Ax = [0 1 0 1 0;
          1 0 1 0 1;
          0 1 0 0 1;
          1 0 0 0 0;
          0 1 1 0 0];

    Av = [0 1 0 1 0;
          1 0 0 1 1;
          0 0 0 0 1;
          1 1 0 0 1;
          0 1 1 1 0];

    Gx = [1 0 1 0 0];

    Gv = [0 1 1 0 0];
end

% Target
% Leader와 연결된 노드는 x*_i = x*_ip+1 + leader_x 를 만족한다.
% 또한 x*_ij = x*_ik + x*_kj를 만족하므로 이를 이용해 x*_ip+1을 구한다.
% 따라서 ex = x_i - (x_star_i + leader_x)를 만족한다.
x_star(:,1) = [0 -1]';
x_star(:,2) = [2 -0]';
x_star(:,3) = [1 1]';
x_star(:,4) = [-1 1]';
x_star(:,5) = [-2 0]';

% System parameter
gamma = -2;
us = 10;
d1 = [5 2];
d2 = [7 -3];
d3 = [-3 3];
d4 = [6 -2];
d5 = [-1 -4];
d = [d1; d2; d3; d4; d5]';

% Control gain
k1 = 5;
k2 = 2;
k3 = 5;
k4 = 2;

%% Laplacian Matrix
% Position
Dx = sum(Ax,2);
Dx = diag(Dx);

Gx = diag(Gx);

Lx = Dx - Ax + Gx;

% Velocity
Dv = sum(Av,2);
Dv = diag(Dv);

Gv = diag(Gv);

Lv = Dv - Av + Gv;

%% Saturation function
sat = @(u) max(min(u, us), -us);

%% Initial
x = zeros(2,6);
v = zeros(2,6);
z = zeros(5,2);

% Set Leader's position
x(:,6) = [5 5]';

% Set leader's velocity
if Velo
    v(:,6) = [0.5 0.2]';
end

%% Simulation
% Simulation Time and Rate
dt = 0.001;

idx = 1;
for k = 0:dt:30
    follower_x = x(:,1:5);
    follower_v = v(:,1:5);

    leader_x = x(:,6);
    leader_v = v(:,6);

    d_bar = gamma*leader_v + d;

    ev = follower_v(:,1:5)' - leader_v';
    ex = follower_x(:,1:5)' - (x_star(:,1:5)' + leader_x');
    
    z_dot = -Lx*ex - Lv*ev;
    z = z + dt*z_dot;

    q = -Lx*ev;

    u = -k1*Lx*ex - k2*Lv*ev + k3*z + k4*q;

    ev_dot = gamma*ev +sat(u) + d_bar';
    ev = ev + dt*ev_dot;

    ex = ex + dt*ev;                                                       % ex_dot = ev
    leader_x = leader_x + dt*leader_v;
    follower_v = ev' + leader_v;
    follower_x = ex' + x_star(:,1:5) + leader_x;
    
    x(:,1:5) = follower_x;
    x(:,6) = leader_x;
    v(:,1:5) = follower_v;
    
    % Data recording
    t_rec(:,idx) = k;
    v_rec(:,:,idx) = v;
    x_rec(:,:,idx) = x;

    idx = idx + 1;
end

%% Plot the results of simulation
% Position trajectory of follower agent(dashed line) and leader agent
fx = squeeze(x_rec(:, 1:5, :));
lx = squeeze(x_rec(:,6,:));

follower_start_position = fx(:,:,1);
follower_end_position = fx(:,:,end);

fig1 = figure(1);
hold on;

for i = 1:5
    plot(squeeze(fx(1,i,:)), squeeze(fx(2,i,:)), ...
        '--', 'LineWidth', 1.5);
end

plot(lx(1,:), lx(2,:), '-. k', 'LineWidth', 1.5, ...
    'DisplayName', 'Leader agent');

plot(follower_start_position(1,:), follower_start_position(2,:), ...
    'DisplayName', 'Initial Positions');

plot(follower_end_position(1,:), follower_end_position(2,:), ...
    's', 'MarkerSize', 8, 'LineWidth', 1.5, ...
    'DisplayName', 'Final Positions');

plot(lx(1,1), lx(2,1), '*', ...
    'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Leader agent');

xlabel('x-axis $x_{1,i}$, $i = 1,\dots,6$', 'Interpreter', 'latex');
ylabel('y-axis $x_{2,i}$, $i = 1,\dots,6$', 'Interpreter', 'latex');
title('Position trajectories of follower agents and leader agent', ...
    'Interpreter', 'latex');
legend('show', 'Location', 'best');
grid on;
hold off;

% Velocity trajectories of follower agent(solid line) and leader agent(dashed line)
if Velo
    % (a)
    fig2 = figure(2);
    set(fig2, 'OuterPosition', [100, 700, 500, 300])
    a1 = plot(t_rec, squeeze(v_rec(1,1,:)), '-k', 'LineWidth', 1.5); hold on
    a2 = plot(t_rec, squeeze(v_rec(1,2,:)), '-k', 'LineWidth', 1.5); hold on
    a3 = plot(t_rec, squeeze(v_rec(1,3,:)), '-k', 'LineWidth', 1.5); hold on
    a4 = plot(t_rec, squeeze(v_rec(1,4,:)), '-k', 'LineWidth', 1.5); hold on
    a5 = plot(t_rec, squeeze(v_rec(1,5,:)), '-k', 'LineWidth', 1.5); hold on
    a6 = plot(t_rec, squeeze(v_rec(1,6,:)), '-. r', 'LineWidth', 1.5); hold on
    grid on
    xlabel('Time (s)')
    ylabel('$v_{1,i}$, $i = 1,\dots,6$', 'Interpreter', 'latex')

    % b
    fig3 = figure(3);
    set(fig3, 'OuterPosition', [600, 700, 500, 300])
    a1 = plot(t_rec, squeeze(v_rec(2,1,:)), '-k', 'LineWidth', 1.5); hold on
    a2 = plot(t_rec, squeeze(v_rec(2,2,:)), '-k', 'LineWidth', 1.5); hold on
    a3 = plot(t_rec, squeeze(v_rec(2,3,:)), '-k', 'LineWidth', 1.5); hold on
    a4 = plot(t_rec, squeeze(v_rec(2,4,:)), '-k', 'LineWidth', 1.5); hold on
    a5 = plot(t_rec, squeeze(v_rec(2,5,:)), '-k', 'LineWidth', 1.5); hold on
    a6 = plot(t_rec, squeeze(v_rec(2,6,:)), '-. r', 'LineWidth', 1.5); hold on
    grid on
    xlabel('Time (s)')
    ylabel('$v_{2,i}$, $i = 1,\dots,6$', 'Interpreter', 'latex')
end