%% Display initialization
close all
clear
clc

%% Parameters
% Input staturation
N = 10;

% Topology
A = [0 1 0 0 0 0 0 0 0 1;
     1 0 1 0 0 0 0 0 0 0;
     0 1 0 1 0 0 0 0 0 0;
     0 0 1 0 1 0 0 0 0 0;
     0 0 0 1 0 1 0 0 0 0;
     0 0 0 0 1 0 1 0 0 0;
     0 0 0 0 0 1 0 1 0 0;
     0 0 0 0 0 0 1 0 1 0;
     0 0 0 0 0 0 0 1 0 1;
     1 0 0 0 0 0 0 0 1 0
    ];

% Known target
Alpha0 = [1 1 1 1 1 0 0 0 0 0];

% Target
v_star = ones(N,1);

% Control gain
alpha = 5;
beta = 5;
K = 5*Alpha0;
K = diag(K);

%% Laplacian Matrix
D = sum(A,2);
D = diag(D);
L = D - A;                                                                 % Laplacian Matrix

%% Interaction function (Saturation)
f = @(z) sign(z) .* min(abs(z), 1);

%% Initial Vector
% Leader
x_leader = randi([-50,50], N, 1);
v_leader = zeros(N,1);                                                     % w

% Follower
x_follower = 100;
v_follower = 10;

%% Simulation
% Simulation Time and Rate
dt = 0.001;

idx = 1;
for k = 0:dt:300
    zi = x_leader - x_follower;
    ei = v_leader - v_star;
    fz = f(zi);

    en = v_follower - v_star;
    en_dot = sum(fz);                                                      % By Eq. (8)
    en = en + dt*en_dot;
    v_follower = en + v_star;
    x_follower = x_follower + dt*v_follower;
       
    z_dot = -alpha*L*fz + ei - en;
    e_dot = -K*ei - beta*L*ei - fz;

    zi = zi + dt*z_dot;
    ei = ei + dt*e_dot;
        
    v_leader = ei + v_star;  
    x_leader = zi + x_follower;

    % Data recording
    t_rec(:,idx) = k;
    x_leader_rec(:,idx) = x_leader;
    v_leader_rec(:,idx) = v_leader;
    x_follower_rec(:,idx) = x_follower;
    v_follower_rec(:,idx) = v_follower;
    fz_rec(:,idx) = fz;

    idx = idx + 1;
end

%% Plot the results of simulation
% 1. Position trajectory of leaders (solid line) and follower (dashed line)
fig1 = figure(1);
set(fig1, 'OuterPosition', [100, 700, 500, 300])
a1 = plot(t_rec, x_leader_rec(1,:), '-k', 'LineWidth', 1.5); hold on
a2 = plot(t_rec, x_leader_rec(2,:), '-k', 'LineWidth', 1.5); hold on
a3 = plot(t_rec, x_leader_rec(3,:), '-k', 'LineWidth', 1.5); hold on
a4 = plot(t_rec, x_leader_rec(4,:), '-k', 'LineWidth', 1.5); hold on
a5 = plot(t_rec, x_leader_rec(5,:), '-k', 'LineWidth', 1.5); hold on
a6 = plot(t_rec, x_leader_rec(6,:), '-k', 'LineWidth', 1.5); hold on
a7 = plot(t_rec, x_leader_rec(7,:), '-k', 'LineWidth', 1.5); hold on
a8 = plot(t_rec, x_leader_rec(8,:), '-k', 'LineWidth', 1.5); hold on
a9 = plot(t_rec, x_leader_rec(9,:), '-k', 'LineWidth', 1.5); hold on
a10 = plot(t_rec, x_leader_rec(10,:), '-k', 'LineWidth', 1.5); hold on
a11 = plot(t_rec, x_follower_rec(1,:), '-. r', 'LineWidth', 1.5); hold on
grid on
xlabel('Time (s)')
ylabel('$x_i$', 'Interpreter', 'latex')



% 2. Trajectory of agent group state v
fig2 = figure(2);
set(fig2, 'OuterPosition', [600, 700, 500, 300])
a1 = plot(t_rec, v_leader_rec(1,:), '-k', 'LineWidth', 1.5); hold on
a2 = plot(t_rec, v_leader_rec(2,:), '-k', 'LineWidth', 1.5); hold on
a3 = plot(t_rec, v_leader_rec(3,:), '-k', 'LineWidth', 1.5); hold on
a4 = plot(t_rec, v_leader_rec(4,:), '-k', 'LineWidth', 1.5); hold on
a5 = plot(t_rec, v_leader_rec(5,:), '-k', 'LineWidth', 1.5); hold on
a6 = plot(t_rec, v_leader_rec(6,:), '-k', 'LineWidth', 1.5); hold on
a7 = plot(t_rec, v_leader_rec(7,:), '-k', 'LineWidth', 1.5); hold on
a8 = plot(t_rec, v_leader_rec(8,:), '-k', 'LineWidth', 1.5); hold on
a9 = plot(t_rec, v_leader_rec(9,:), '-k', 'LineWidth', 1.5); hold on
a10 = plot(t_rec, v_leader_rec(10,:), '-k', 'LineWidth', 1.5); hold on
a11 = plot(t_rec, v_follower_rec(1,:), '-. r', 'LineWidth', 1.5); hold on
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

% 환경에 따라 시스템의 특정 상태를 알 수 없는 경우가 있으므로
% 시스템의 상태변수가 아닌 상호작용에 대한 정보만을 이용.
% 다중 에이전트들의 협력 제어 모델을 기반으로 다중 리더와 단일 추종자를
% 고려한 일치 문제를 연구. 즉, 다중 리더는 제어 가능하고 추종자는 리더와의
% 상호작용을 통해 제어. 이때 리더의 상호작용만을 측정할 수 있다.
% 리더의 위치는 입력에 의해 제어되고, 추종의 가속도는 모든 리더들의 
% 상호작용에 의해 제어된다. 
% 즉, v_follower_dot = sigma[f(x_leader - x_follwer)]
% 이때 f는 상호작용 함수.
% 논문의 목표: 리더들의 협력을 통해 전체 그룹의 일치를 달성.