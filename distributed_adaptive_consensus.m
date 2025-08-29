%% Display initialization
close all
clear
clc

%% Run
% If you want to see graph for leaderless, Set LEADER = 0
% If you want to see graph for leader, Set LEADER = 1
LEADER = 0;

%% Parameters
% Input staturation
N = 5;

% System parameter
A = [0 1;
     0 0];

B = [0 1]';

% Positive definite matrix
P = [2.7881 -0.8250;
     -0.8250 0.9886];
P_inv = inv(P);

% Gain
k_i = 0.001;                                                               % For i = 1,2,3,4,5
H = P_inv*(B*B')*P_inv;
K = -B'*P_inv;

% Topology
% Graph of Follower (Adjacency Matrix)
A1 = [0 1 1 1 0;
      1 0 0 1 1;
      1 0 0 0 1;
      1 1 0 0 1;
      0 1 1 1 0];

% Graph of Leader-Follower Connections
A2 = [0 1 0 0 0 0;
      1 0 1 1 1 0;
      0 1 0 0 1 1;
      0 1 0 0 0 1;
      0 1 1 0 0 1;
      0 0 1 1 1 0];

%% Laplacian Matrix
% Laplacian of Follower node
D1 = diag(sum(A1,2));
L1 = D1 - A1;                                                              % Laplacian Matrix

% Laplacian of Leader node
D2 = diag(sum(A2,2));
L_bar = D2 - A2; 

%% Initial
x0 = [22; 2];                                                              % Leader state
x1 = [15; 3];   
x2 = [10; -2];  
x3 = [5; -5];  
x4 = [-5; 2];  
x5 = [-10; -1];
x = [x1 x2 x3 x4 x5];                                                      % Followers state
                                                             
c_dot = ones(N,N);
c = zeros(N,N);
rho = zeros(N,N);

%% Input Gain Nonlinear Function
f = @(u) atan(u)/2;

%% Simulation
% Simulation Time and Rate
dt = 0.001;

idx = 1;

for k=1:dt:250

    if ~LEADER
        x = reshape(x,10,1);
        z = kron(L1,eye(2))*x;
        z = reshape (z,2,N);

        c_dot(:,1:N) = z(:,1:N)'*H*z(:,1:N);
        c = c + dt*c_dot;
        
        rho(:,1:N) = k_i*z(:,1:N)'*P_inv*z(:,1:N);
        c = diag(c);
        C = diag(c);
        rho = diag(rho);
        Rho = diag(rho);

        u = kron(C + Rho, K)*reshape(z,10,1);
        u = reshape(u,1,N);

        x = reshape(x, 2, N);
        x_dot = A*x(:,1:N) + B*(u + f(u));
        x = x + dt*x_dot;
    end

    if LEADER
        X = [x0 x];
        X = reshape (X,12,1);
        z = kron(L_bar, eye(2)) * (X - kron(ones(N+1,1), x0));
        z = reshape (z,2,N+1);

        c_dot(:,1:N) = z(:,2:N+1)'*H*z(:,2:N+1);
        c = c + dt*c_dot;

        rho(:,1:N) = k_i*z(:,2:N+1)'*P_inv*z(:,2:N+1);
        c = diag(c);
        C = diag(c);
        rho = diag(rho);
        Rho = diag(rho);

        u = kron(C + Rho, K)*reshape(z(:,2:N+1),10,1);
        u = reshape(u,1,N);

        X = reshape(X, 2, N+1);
        x0_dot = A*X(:,1);
        x0 = x0 + dt*x0_dot;

        x_dot = A*X(:,2:N+1) + B*(u + f(u));
        X(:,2:N+1) = X(:,2:N+1) + dt*x_dot;
        x = X(:,2:N+1);
    end

    t_rec(:,idx) = k;
    x_rec(:,:,idx) = x;
    x0_rec(:,:,idx) = x0;
    c_rec(:,idx) = c;
    
    idx = idx + 1;
end

%% Plot the results of simulation

if ~LEADER
    % 1. Trajectory of States xi,1 and xi,2 in leaderless case.
    fig1 = figure(1);
    set(fig1, 'OuterPosition', [100, 700, 500, 300])
    a1 = plot(t_rec, squeeze(x_rec(1,1,:)), '-k', 'LineWidth', 1.5); hold on
    a2 = plot(t_rec, squeeze(x_rec(1,2,:)), '-k', 'LineWidth', 1.5); hold on
    a3 = plot(t_rec, squeeze(x_rec(1,3,:)), '-k', 'LineWidth', 1.5); hold on
    a4 = plot(t_rec, squeeze(x_rec(1,4,:)), '-k', 'LineWidth', 1.5); hold on
    a5 = plot(t_rec, squeeze(x_rec(1,5,:)), '-k', 'LineWidth', 1.5); hold on
    a6 = plot(t_rec, squeeze(x_rec(2,1,:)), '-r', 'LineWidth', 1.5); hold on
    a7 = plot(t_rec, squeeze(x_rec(2,2,:)), '-r', 'LineWidth', 1.5); hold on
    a8 = plot(t_rec, squeeze(x_rec(2,3,:)), '-r', 'LineWidth', 1.5); hold on
    a9 = plot(t_rec, squeeze(x_rec(2,4,:)), '-r', 'LineWidth', 1.5); hold on
    a10 = plot(t_rec, squeeze(x_rec(2,5,:)), '-r', 'LineWidth', 1.5); hold on
    grid on
    xlabel('t')
    ylabel('$x_i$', 'Interpreter', 'latex')
    legend([a1,a10], ...
        {'$xi,1$', '$xi,2$'}, 'Interpreter', 'latex')


    % 2. Trajectory of coupling gain ci in leaderless case.
    fig2 = figure(2);
    set(fig2, 'OuterPosition', [600, 700, 500, 300])
    a1 = plot(t_rec, c_rec(1,:), '-k', 'LineWidth', 1.5); hold on
    a2 = plot(t_rec, c_rec(2,:), '-k', 'LineWidth', 1.5); hold on
    a3 = plot(t_rec, c_rec(3,:), '-k', 'LineWidth', 1.5); hold on
    a4 = plot(t_rec, c_rec(4,:), '-k', 'LineWidth', 1.5); hold on
    a5 = plot(t_rec, c_rec(5,:), '-k', 'LineWidth', 1.5); hold on
    grid on
    xlabel('t')
    ylabel('$ci$', 'Interpreter', 'latex')
end

if LEADER
    % 1. Trajectory of States xi,1 and xi,2 in leaderless case.
    fig1 = figure(1);
    set(fig1, 'OuterPosition', [100, 700, 500, 300])
    a1 = plot(t_rec, squeeze(x_rec(1,1,:)), '-k', 'LineWidth', 1.5); hold on
    a2 = plot(t_rec, squeeze(x_rec(1,2,:)), '-k', 'LineWidth', 1.5); hold on
    a3 = plot(t_rec, squeeze(x_rec(1,3,:)), '-k', 'LineWidth', 1.5); hold on
    a4 = plot(t_rec, squeeze(x_rec(1,4,:)), '-k', 'LineWidth', 1.5); hold on
    a5 = plot(t_rec, squeeze(x_rec(1,5,:)), '-k', 'LineWidth', 1.5); hold on
    a6 = plot(t_rec, squeeze(x_rec(2,1,:)), '-r', 'LineWidth', 1.5); hold on
    a7 = plot(t_rec, squeeze(x_rec(2,2,:)), '-r', 'LineWidth', 1.5); hold on
    a8 = plot(t_rec, squeeze(x_rec(2,3,:)), '-r', 'LineWidth', 1.5); hold on
    a9 = plot(t_rec, squeeze(x_rec(2,4,:)), '-r', 'LineWidth', 1.5); hold on
    a10 = plot(t_rec, squeeze(x_rec(2,5,:)), '-r', 'LineWidth', 1.5); hold on
    a11 = plot(t_rec, squeeze(x0_rec(1,:,:)), '-. b', 'LineWidth', 1.5); hold on
    a12 = plot(t_rec, squeeze(x0_rec(2,:,:)), '-. y', 'LineWidth', 1.5); hold on
    grid on
    xlabel('t')
    ylabel('$x_i$', 'Interpreter', 'latex')
    legend([a1,a10, a11, a12], ...
        {'$xi,1$', '$xi,2$', '$x0,1$', '$x0,2$'}, 'Interpreter', 'latex')


    % 2. Trajectory of coupling gain ci in leaderless case.
    fig2 = figure(2);
    set(fig2, 'OuterPosition', [600, 700, 500, 300])
    a1 = plot(t_rec, c_rec(1,:), '-k', 'LineWidth', 1.5); hold on
    a2 = plot(t_rec, c_rec(2,:), '-k', 'LineWidth', 1.5); hold on
    a3 = plot(t_rec, c_rec(3,:), '-k', 'LineWidth', 1.5); hold on
    a4 = plot(t_rec, c_rec(4,:), '-k', 'LineWidth', 1.5); hold on
    a5 = plot(t_rec, c_rec(5,:), '-k', 'LineWidth', 1.5); hold on
    grid on
    xlabel('t')
    ylabel('$ci$', 'Interpreter', 'latex')
end