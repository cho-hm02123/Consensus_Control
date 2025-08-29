clc; 
clear; 
close all;

%% System Params
A = [0 1; 0 0];
B = [0 1]';
N =5;

Adj_f = [0 1 0 0 0 0;
    0 0 1 1 1 0;
    0 1 0 0 1 1;
    0 1 0 0 0 1;
    0 1 1 0 0 1;
    0 0 1 1 1 0];

f = @(u) atan(u)/2;

D_column = sum(Adj_f,2);
D = diag(D_column);   
L_f = D - Adj_f;
Adj = Adj_f(2:end,2:end);
L = L_f(2:end, 2:end);


P = [2.7881 -0.8250;
    -0.8250 0.9886];

K = -B' * inv(P);
ki = 0.001;
H = inv(P) * B * B' * inv(P);

x1 = [15; 3];   
x2 = [10; -2];  
x3 = [5; -5];  
x4 = [-5; 2];  
x5 = [-10; -1]; 
x = [x1; x2; x3; x4; x5];







c_dot = ones(1,N);
c = zeros(1,N);
p_i = zeros(1,N);


%% Simulate
T = 250;
dt = 0.001;
index = 1;


for t =  0:dt:T
    z = kron(L, eye(2)) * x;

    for i = 1: 5
        z_i = z(2 * i -1: 2 * i);
        c_dot(i) = z_i' * H * z_i;
        p_i(i) = ki*(z_i' * inv(P) * z_i);
        c(i) = c(i) + dt * c_dot(i);
    end

    C = diag(c);
    p_diag = diag(p_i);

    u = (kron(C+p_diag,K) )* z;
    f_u = f(u);

    % z_dot = (kron(eye(N), A) + kron(L, B * K)) * z + kron(L, B) * f_u;
    % z = z + dt * z_dot;
    x_dot = zeros(2*N,1);
    for i = 1: 5
        x_i = x(2 * i -1: 2 * i);
        x_dot(2 * i -1: 2 * i) = A * x_i + B * (u(i) + f_u(i));
    end

    x = x + dt * x_dot;
    x_rec(:,index) = x;
    t_rec(:,index) = t;
    c_rec(:,index) = c;

    index = index +1;
    
    
end

%%

fig1 = figure(1);
h1 = plot(t_rec, x_rec(1,:), 'k', 'LineWidth', 1.5);
hold on;
plot(t_rec, x_rec(3,:), 'k', 'LineWidth', 1.5);
plot(t_rec, x_rec(5,:), 'k', 'LineWidth', 1.5);
plot(t_rec, x_rec(7,:), 'k', 'LineWidth', 1.5);
plot(t_rec, x_rec(9,:), 'k', 'LineWidth', 1.5);


h2 = plot(t_rec, x_rec(2,:), 'r--', 'LineWidth', 1.5); 
plot(t_rec, x_rec(4,:), 'r--', 'LineWidth', 1.5);
plot(t_rec, x_rec(6,:), 'r--', 'LineWidth', 1.5);
plot(t_rec, x_rec(8,:), 'r--', 'LineWidth', 1.5);
plot(t_rec, x_rec(10,:), 'r--', 'LineWidth', 1.5);


xlabel('Time (t)', 'FontSize', 12);
ylabel('State Variables', 'FontSize', 12);
title('Trajectories of $x_{i,1}$ and $x_{i,2}$ in Leaderless Case', 'Interpreter', 'latex', 'FontSize', 14);
ylim([-inf, 20]);
grid on;

legend([h1, h2], {'$x_{i,1}$ (Position)', '$x_{i,2}$ (Velocity)'}, ...
       'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');

hold off;
%%
fig2 = figure(2);
plot(t_rec, c_rec', 'k', 'LineWidth', 1.5); 
xlabel('t', 'FontSize', 12);
ylabel('$c_i$', 'Interpreter', 'latex', 'FontSize', 12);
title('Trajectories of Coupling Gains $c_i$ in leaderless case', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12); 
ylim([0,20]);
%%

clc; 
clear; 
%%
A = [0 1; 0 0];
B = [0 1]';
N_total = 6;
N_follower = 5;

Adj_f = [0 1 0 0 0 0;
    1 0 1 1 1 0;
    0 1 0 0 1 1;
    0 1 0 0 0 1;
    0 1 1 0 0 1;
    0 0 1 1 1 0];

f = @(u) atan(u)/2;

D_column = sum(Adj_f,2);
D = diag(D_column);   
L1 = D - Adj_f;
Adj = Adj_f(2:end,2:end);
L = L1(2:end, 2:end);

P = [2.7881 -0.8250;
    -0.8250 0.9886];

K = -B' * inv(P);
ki = 0.001;
H = inv(P) * B * B' * inv(P);

x0 = [22;2];
x1 = [15; 3];   
x2 = [10; -2];  
x3 = [5; -5];   
x4 = [-5; 2];   
x5 = [-10; -1];
x = [x0;x1; x2; x3; x4; x5];

c_dot = ones(1,N_follower);
c = zeros(1,N_follower);
p_i = zeros(1,N_follower);

%%
T = 250;
dt = 0.001;
index = 1;

for t =  0:dt:T
    
    x_leader = x(1:2);       
    x_followers = x(3:end);   
  
    z = kron(L1, eye(2)) * (x - kron(ones(6,1), x_leader));

    for i = 2: 6
        j = i - 1;
        z_i = z(2*i-1:2*i);
        c_dot(j) = z_i' * H * z_i;
        p_i(j) = ki * (z_i' * inv(P) * z_i);
        c(j) = c(j) + dt * c_dot(j);
    end



    C = diag(c);
    p_diag = diag(p_i);

    for i = 2:N_total
        j = i - 1;
        u(j) = (c(j) + p_i(j)) * K * z(2*i-1:2*i);
    end
    f_u = f(u);

    % z_dot = (kron(eye(N), A) + kron(L, B * K)) * z + kron(L, B) * f_u;
    % z = z + dt * z_dot;
    x_dot = zeros(10,1);
    for i = 1: 5
        x_i = x_followers(2 * i -1: 2 * i);
        x_dot(2 * i -1: 2 * i) = A * x_i + B * (u(i) + f_u(i));
    end
    x_followers = x_followers + dt * x_dot;
    x_leader = x_leader + dt* A*x_leader;
    x =[x_leader;x_followers];
    x_rec(:,index) = x;
    t_rec(:,index) = t;
    c_rec(:,index) = c;

    index = index +1;
    
    
end


%%
fig3 = figure(3);


plot(t_rec, x_rec(3,:), 'k', 'LineWidth', 1); hold on;
plot(t_rec, x_rec(5,:), 'k', 'LineWidth', 1);
plot(t_rec, x_rec(7,:), 'k', 'LineWidth', 1);
plot(t_rec, x_rec(9,:), 'k', 'LineWidth', 1);
plot(t_rec, x_rec(11,:), 'k', 'LineWidth', 1);


plot(t_rec, x_rec(4,:), 'r', 'LineWidth', 1);
plot(t_rec, x_rec(6,:), 'r', 'LineWidth', 1);
plot(t_rec, x_rec(8,:), 'r', 'LineWidth', 1);
plot(t_rec, x_rec(10,:), 'r', 'LineWidth', 1);
plot(t_rec, x_rec(12,:), 'r', 'LineWidth', 1);


plot(t_rec, x_rec(1,:), 'b', 'LineWidth', 1);
plot(t_rec, x_rec(2,:), 'b--', 'LineWidth', 1);


xlabel('Time (t)', 'FontSize', 12);
ylabel('x_i', 'FontSize', 12);
title('Trajectories of $x_{i,1}$ and $x_{i,2}$ in Leader-following Case', 'Interpreter', 'latex', 'FontSize', 14);


legend('$x_{i,1}$', '$x_{i,2}$', '$x_{0,1}$', '$x_{0,2}$', ...
       'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');

grid on;

hold off;
%%
fig4 = figure(4);
plot(t_rec, c_rec', 'k', 'LineWidth', 1.5);
xlabel('t', 'FontSize', 12);
ylabel('$c_i$', 'Interpreter', 'latex', 'FontSize', 12);
title('Trajectories of Coupling Gains $c_i$ in leader-following case', 'Interpreter', 'latex', 'FontSize', 14);
grid on;
ylim([0, 20]);
set(gca, 'FontSize', 12);