% MATLAB脚本：动态模拟大理石台阶磨损量随时间变化的模型

% 时间参数
T = [0, 200, 400, 600, 800]; % 分段时间点
dt = 0.1;                   % 时间步长
Time = T(1):dt:T(end);      % 总时间范围

% 初始条件
W = zeros(size(Time));       % 磨损量初始化
H0 = 100;                   % 初始硬度 (单位：任意)
H_min = 20;                 % 硬度最小值
mu0 = 0.7;                  % 初始摩擦系数
mu_min = 0.3;               % 摩擦系数最小值
K0 = 1e-3;                  % 初始磨损系数
K_max = 1e-2;               % 最大磨损系数
F0 = 50;                    % 初始正压力 (单位：任意)

% 分段游客频率函数
v_low = 10;                 % 阶段1游客频率
k2 = 0.05;                  % 阶段2游客频率增长率
v_med = 20;                 % 阶段3初始游客频率
k3 = 0.01;                  % 阶段3游客指数增长率
v_high = 100;               % 阶段4游客饱和频率
A = 10;                     % 阶段4游客波动幅度
omega = 2*pi/50;            % 阶段4波动频率

% 定义时间动态函数
v_t = @(t) (t < T(2)) .* v_low + ...                           % 阶段1
            ((t >= T(2)) & (t < T(3))) .* (v_low + k2*(t-T(2))) + ... % 阶段2
            ((t >= T(3)) & (t < T(4))) .* (v_med .* exp(k3*(t-T(3)))) + ... % 阶段3
            (t >= T(4)) .* (v_high + A*sin(omega*(t-T(4))));  % 阶段4

% 模型参数动态变化
H_t = @(W) max(H_min, H0 - 0.01*W);     % 硬度随磨损变化
mu_t = @(t) max(mu_min, mu0*exp(-0.001*t)); % 摩擦系数随时间变化
K_t = @(W) min(K_max, K0*(1 + 0.005*W));    % 磨损系数随磨损变化

% 动态模拟
for i = 2:length(Time)
    t = Time(i); % 当前时间
    v = v_t(t);  % 当前游客频率
    H = H_t(W(i-1)); % 当前硬度
    mu = mu_t(t); % 当前摩擦系数
    K = K_t(W(i-1)); % 当前磨损系数

    % 磨损速率计算
    dW_dt = K * mu * (F0 * v) / H;

    % 更新磨损量
    W(i) = W(i-1) + dW_dt * dt;
end

% 可视化动态变化
figure;
subplot(4, 1, 1);
plot(Time, W, 'b', 'LineWidth', 2);
xlabel('Time (years)');
ylabel('Cumulative Wear Volume');
title('Cumulative Wear over Time');
grid on;

H_dynamic = H_t(W);
subplot(4, 1, 2);
plot(Time, H_dynamic, 'g', 'LineWidth', 2);
xlabel('Time (years)');
ylabel('Hardness');
title('Material Hardness over Time');
grid on;

mu_dynamic = mu_t(Time);
subplot(4, 1, 3);
plot(Time, mu_dynamic, 'r', 'LineWidth', 2);
xlabel('Time (years)');
ylabel('Friction Coefficient');
title('Friction Coefficient over Time');
grid on;

visitor_dynamic = arrayfun(v_t, Time);
subplot(4, 1, 4);
plot(Time, visitor_dynamic, 'orange', 'LineWidth', 2);
xlabel('Time (years)');
ylabel('Visitor Frequency');
title('Visitor Frequency over Time');
grid on;

sgtitle('Dynamic Changes in Wear Model Parameters');
