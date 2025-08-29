clc; clear; close all;

%% Run
%% If you want Fig2, Set Fig1 = 1, Velo = 0
%% If you want Fig3 and Fig 4, Set Fig1 = 1, Velo = 1
%% If you want Fig6 and Fig 7, Set Fig1 = 0, Velo = 1

Fig1 = 1;
Velo = 1;

%% Graph Topology (Adjacency Matrices)
if Fig1
    Ax = [ 0, 1, 0, 1, 0; 
           1, 0, 1, 0, 1; 
           0, 1, 0, 0, 1; 
           1, 0, 0, 0, 0; 
           0, 1, 1, 0, 0];

    Av = [ 0, 1, 0, 1, 0; 
           1, 0, 0, 0, 0; 
           0, 0, 0, 0, 1; 
           1, 0, 0, 0, 0; 
           0, 0, 1, 0, 0];

    Gx = [1, 0, 1, 0, 0]; % leader와 연결
    Gv = [0, 1, 0, 0, 0]; % leader와 연결

else
    Ax = [0, 1, 0, 1, 0;
          1, 0, 1, 0, 1;
          0, 1, 0, 0, 1;
          1, 0, 0, 0, 0;
          0, 1, 1, 0, 0];

    Av = [0, 1, 0, 1, 0;
          1, 0, 0, 1, 1;
          0, 0, 0, 0, 1;
          1, 1, 0, 0, 1;
          0, 1, 1, 1, 0];

    Gx = [1, 0, 1, 0, 0]; % leader와 연결
    Gv = [0, 1, 1, 0, 0]; % leader와 연결

end

%% Graph Topology Laplacian Matrix
N = size(Ax, 1);    % Number of agents
Dx = diag(sum(Ax, 2)); 
Dv = diag(sum(Av, 2)); 

Gx = diag(Gx);      % Leader Connection 
Gv = diag(Gv);      % Leader Connection 

Lx = Dx - Ax + Gx;  % Position Laplacian
Lv = Dv - Av + Gv;  % Velocity Laplacian

%% parameter
gamma = -2;
uB = 10;
d1 = [ 5  2];
d2 = [ 7 -3];
d3 = [-3  3];
d4 = [ 6 -2];
d5 = [-1 -4];
d  = [d1; d2; d3; d4; d5]';

%% Control gain
%% gamma*k1 + k3 + (gamma*gamma/4)*k4 < 0; gamma*k2 + k3 > 0 (PID based Consensus)
k1 = 5; k2 = 2; k3 = 5; k4 = 2;

%% Saturation function
sat = @(u) max(min(u, uB), -uB);

%% x_star from relative distance 
%% x*_ij = x*_i - x*_j
%% x*_16 = x*_1 - x*_6
%% x*_36 = x*_3 - x*_6
x_star(:,1) = [ 0 -1]';
x_star(:,2) = [ 2  0]';
x_star(:,3) = [ 1  1]';
x_star(:,4) = [-1  1]';
x_star(:,5) = [-2  0]';

%% Initial
%% 1:5 -> follower, 6 -> leader
xF = zeros(N, 2);
vF = zeros(N, 2);
z  = zeros(N, 2);
u  = zeros(N, 2);

xL = [5, 5];
if Velo
    vL = [0.5, 0.2];
else
    vL = [0, 0];
end

ev = vF - vL;
ex = xF - (x_star' + xL);

d_bar = gamma*vL' + d;

%% Simulation
sim_time = 20; 
dt       = 0.05;
iter     = sim_time / dt;

for k = 1:iter
    t = k * dt;
    
    % Update leader dynamics
    xL   = xL + dt*vL;

    % Control law
    ev = ev + dt*(gamma*ev +sat(u) + d_bar');
    ex = ex + dt*ev;     

    z_dot = -Lx*ex - Lv*ev;
    z = z + dt*z_dot;

    q = -Lx*ev;

    u = -k1*Lx*ex - k2*Lv*ev + k3*z + k4*q;
    
    % Update follower dynamics
    % vF = ev + vL;
    % xF = ex + x_star' + xL;

    vF = vF + dt*(gamma*vF + sat(u) + d');
    xF = xF + dt*(vF);

    % Data recording
    t_rec(k) = t;
    xL_rec(k, :) = xL;
    xF_rec(k, :, :) = xF;
    vF_rec(k, :, :) = vF;

    % real-time plot
    clf;
    hold on; grid on;
    plot(xL(1), xL(2), 'ro');       % leader 
    plot(xF(:, 1), xF(:, 2), 'bo'); % follower

    for i = 1:N
        for j = 1:N
            if Lx(i, j) ~= 0
                plot([xF(i, 1), xF(j, 1)], [xF(i, 2), xF(j, 2)], 'k--');
            end
        end
    end
    axis([0 20 0 20]);
    drawnow;
end

%% Plot the results of simulation
figure(1);
hold on; grid on; axis equal;

% leader trajectory
plot(xL_rec(:,1), xL_rec(:,2), 'k', 'LineWidth', 2, 'DisplayName', 'Leader');

% follow trajectory
colors = {'r', 'g', 'b', 'm', 'c'};
line_styles = {'-.', '--', ':', '-.', '--'};

for i = 1:N
    plot(xF_rec(:, i, 1), xF_rec(:, i, 2), line_styles{i}, 'Color', colors{i}, 'LineWidth', 1.5, ...
        'DisplayName', ['Follower ' num2str(i)]);
    scatter(xF_rec(end, i, 1), xF_rec(end, i, 2), 80, colors{i}, 'filled');
end

% legend show;
title('Position Trajectories of Leader and Followers');
xlabel('X Position');
ylabel('Y Position');

if Velo
    % Velocity trajectories of follower agent(solid line) and leader agent(dashed line)
    fig2 = figure(2);
    set(fig2, 'OuterPosition', [100, 700, 500, 300])
    a1 = plot(t_rec, (vF_rec(:,1,1)), '-k', 'LineWidth', 1.5); hold on
    a2 = plot(t_rec, (vF_rec(:,2,1)), '-k', 'LineWidth', 1.5); hold on
    a3 = plot(t_rec, (vF_rec(:,3,1)), '-k', 'LineWidth', 1.5); hold on
    a4 = plot(t_rec, (vF_rec(:,4,1)), '-k', 'LineWidth', 1.5); hold on
    a5 = plot(t_rec, (vF_rec(:,5,1)), '-k', 'LineWidth', 1.5); hold on
    grid on
    xlabel('Time (s)')
    ylabel('$v_{1,i}$, $i = 1,\dots,6$', 'Interpreter', 'latex')

    % Velocity trajectories of follower agent(solid line) and leader agent(dashed line)
    fig3 = figure(3);
    set(fig3, 'OuterPosition', [600, 700, 500, 300])
    a1 = plot(t_rec, (vF_rec(:,1,2)), '-k', 'LineWidth', 1.5); hold on
    a2 = plot(t_rec, (vF_rec(:,2,2)), '-k', 'LineWidth', 1.5); hold on
    a3 = plot(t_rec, (vF_rec(:,3,2)), '-k', 'LineWidth', 1.5); hold on
    a4 = plot(t_rec, (vF_rec(:,4,2)), '-k', 'LineWidth', 1.5); hold on
    a5 = plot(t_rec, (vF_rec(:,5,2)), '-k', 'LineWidth', 1.5); hold on
    grid on
    xlabel('Time (s)')
    ylabel('$v_{2,i}$, $i = 1,\dots,6$', 'Interpreter', 'latex')

end

% figure(3)
% plot(xL_rec(:, 1),  xL_rec(:, 2), 'k','LineWidth',1); hold on; grid on;
% for i = 1:5
%     plot(xF_rec(:, i, 1), xF_rec(:, i, 2),'LineWidth',1); hold on; grid on;
%     text(xF_rec(end, i, 1), xF_rec(end, i, 2), ['x=' num2str(xF_rec(end, i, 1)) ', y=' num2str(xF_rec(end, i, 2))])
% end
% text(xL_rec(end, 1), xL_rec(end, 2), ['x=' num2str(xL_rec(end, 1)) ', y=' num2str(xL_rec(end, 2))])
% for i = 1:4
%     head = scatter(xL_rec(i*100, 1), xL_rec(i*100, 2), 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
%     for j = 1:5
%         head = scatter(xF_rec(i*100, j, 1), xF_rec(i*100, j, 2), 'filled');
%     end
% end