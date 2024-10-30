addpath('common/'); % 添加路径 'common/'
rng(1); % 设置随机数种子

% 检查是否存在 'generated_data.mat' 文件
if isfile('generated_data.mat')
    load('generated_data.mat'); % 加载已有数据
else
    %%
    N = 1; % 设置迭代次数
    settings = getSettings(); % 获取设置参数
    save('settings.mat', 'settings'); % 保存设置参数到 'settings.mat'
    numStates = settings.numStates; % 获取状态数量
    stateMasks = settings.stateMask; % 获取状态掩码
    numErrorStates = settings.numErrorStates; % 获取误差状态数量
    errorStateMasks = settings.errorStateMask; % 获取误差状态掩码

    numSamples = settings.numSamples; % 获取样本数量
    dT = settings.dT; % 获取时间间隔
    P0 = settings.P; % 获取初始协方差矩阵

    timeVector = 0:dT:(settings.duration-dT); % 生成时间向量

    % 生成轨迹数据
    [position,orientation,velocity,acceleration,angularVelocity] = trajectory_gen(settings);
    % 生成传感器数据
    [ImuMag_data, ImuMag_bias, theta_cell, aux] = sensor_data_gen(settings, position, orientation, acceleration, angularVelocity, N);

    % 保存生成的数据到 'generated_data.mat'
    save('generated_data.mat', 'position', 'orientation', 'velocity', 'acceleration', 'angularVelocity', ...
        'ImuMag_data', 'ImuMag_bias', 'theta_cell', 'aux', 'settings', 'timeVector');
end

% 创建地面真实数据
xs = cell(N, 1);
for i = 1 : N
    xs{i} = [position velocity compact(orientation) ImuMag_bias(i).IMU ImuMag_bias(i).MAG theta_cell{i}];
end

XData = cell(N, 1); % 初始化状态数据
PData = cell(N, 1); % 初始化协方差数据

t = now; % 获取当前时间
d = datetime(t,'ConvertFrom','datenum'); % 转换为日期时间格式

processNoiseCov =  settings.Q; % 获取过程噪声协方差矩阵
measurementNoiseCov = settings.R; % 获取测量噪声协方差矩阵

NEES_metric = zeros(N, numSamples - 1); % 初始化NEES指标
Phi = settings.Phi; % 获取Phi矩阵
tic % 开始计时

for iter = 1 : N
    if(mod(iter, 10) == 1)
        fprintf('%d iter\n', iter); % 每10次迭代打印一次迭代次数
    end
    
    magnetometerReadings = ImuMag_data(iter).MAG; % 获取磁力计读数
    accelerometerReadings = ImuMag_data(iter).IMU(:, 1:3); % 获取加速度计读数
    gyroReadings = ImuMag_data(iter).IMU(:, 4:6); % 获取陀螺仪读数
    
    % 获取多项式系数的初始值
    coeff = inv(Phi.'*Phi) * Phi.'*magnetometerReadings(1, :).';

    % 初始化状态向量
    x0 = [position(1,:) velocity(1,:) compact(orientation(1,:)) settings.init_acc_bias settings.init_gyro_bias settings.init_mag_bias coeff.'].'; 
    X = zeros(numSamples, numStates); % 初始化状态矩阵
    Ps = zeros(numSamples, numStates - 1); % 初始化协方差矩阵
    
    X(1, :) = x0; % 设置初始状态
    Ps(1, :) = diag(P0); % 设置初始协方差

    x = x0; % 当前状态
    P = P0; % 当前协方差
    for i = 1 : numSamples - 1
            % 【预测（Prediction）】：预测步骤用于估计下一时刻的状态x和协方差P。
            acc_m = accelerometerReadings(i, :)'; % 当前加速度计读数
            omega_m = gyroReadings(i, :)'; % 当前陀螺仪读数
            u = [acc_m; omega_m]; % 输入向量u
            [xh, F, Q] =  Nav_eq(x, u, dT, processNoiseCov, settings); % 预测步骤
            % 协方差传播
            P = F * P * F' + Q;

            % 【更新（Update）】：使用观测模型得到校正值 $\text{bel}(\delta \mathbf{x}_{t})$，并计算卡尔曼增益K和协方差矩阵P：
            H = settings.H; % 获取测量矩阵
            
            % 0~20秒：**位置**和磁场测量均可用
            % 20~40秒：仅磁场测量可用
            if i < 2000
                % 如果当前索引 i 小于 2000，执行以下操作
                % 构建观测矩阵 H，添加新的观测模型
                H = [settings.H; [eye(3) zeros(3, numErrorStates - 3)]];
                % 构建测量噪声协方差矩阵，增加新的噪声项
                measurementNoiseCov_ =  blkdiag(measurementNoiseCov, 0.01^2*eye(3));
                % 计算卡尔曼增益 K
                K = P * H' / (H * P * H' + measurementNoiseCov_);
                % 计算测量残差 delta_z
                delta_z = [magnetometerReadings(i + 1, :).' - (Phi * xh(stateMasks.theta) + [xh(stateMasks.mag_bias); zeros(3, 1)]);...
                           (position(i+1, :).'+ 0.01*randn(3,1))-xh(stateMasks.pos)]; % 位置测量
                % 更新协方差矩阵 P
                P = (eye(size(P)) - K * H) * P * (eye(size(P)) - K * H)' + K * measurementNoiseCov_ * K';
                % 计算状态增量 delta_x
                delta_x = K * delta_z;
            elseif i >= 2000 && i < 4000
                % 计算卡尔曼增益 K
                K = P * H' / (H * P * H' + measurementNoiseCov);   
                % 计算测量残差 delta_z
                delta_z = magnetometerReadings(i + 1, :).' - (Phi * xh(stateMasks.theta) + [xh(stateMasks.mag_bias); zeros(3, 1)]);
                % 更新协方差矩阵 P
                P = (eye(size(P)) - K * H) * P * (eye(size(P)) - K * H)' + K * measurementNoiseCov * K';
                % 计算状态增量 delta_x
                delta_x = K * delta_z;
            else

            end

            % 更新状态向量
            xh(stateMasks.pos) = xh(stateMasks.pos) +  delta_x(errorStateMasks.pos);
            xh(stateMasks.vel) = xh(stateMasks.vel) +  delta_x(errorStateMasks.vel);
            xh(stateMasks.q_nb) = quatmultiply(xh(stateMasks.q_nb).',  [1 1/2*delta_x(7:9).']);
            xh(stateMasks.q_nb) = xh(stateMasks.q_nb).' / norm(xh(stateMasks.q_nb)); 
            xh(stateMasks.acc_bias) = xh(stateMasks.acc_bias) +  delta_x(errorStateMasks.acc_bias);
            xh(stateMasks.gyro_bias) = xh(stateMasks.gyro_bias) +  delta_x(errorStateMasks.gyro_bias);
            xh(stateMasks.mag_bias) = xh(stateMasks.mag_bias) +  delta_x(errorStateMasks.mag_bias);
            xh(stateMasks.theta) = xh(stateMasks.theta) +  delta_x(errorStateMasks.theta);
             
            x = xh; % 更新当前状态
            X(i + 1, :) = x; % 保存状态
            Ps(i + 1, :) = diag(P); % 保存协方差
            G = blkdiag(eye(6), eye(3) - 1 / 2 * vect2skew(delta_x(7:9)) ,eye(3 + 3 + (settings.numSensors - 1) * 3 + 15));
            P = G * P * G.'; % 更新协方差
            
            % 计算估计误差
            est_error = xh- xs{iter}(i + 1, :).';
            est_error = est_error(~stateMasks.mag_bias);
            est_error = est_error([1:6, 8:length(est_error)]);
            NEES_metric(iter, i) = est_error.' * (P(~errorStateMasks.mag_bias, ~errorStateMasks.mag_bias) \ est_error);
    end
    
    XData{iter} = X; % 保存状态数据
    PData{iter} = Ps; % 保存协方差数据
end
toc % 结束计时

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             trajectory plot             %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot3(position(:, 1), position(:, 2), position(:, 3), 'k');
hold on;
plot3(position(end, 1), position(end, 2), position(end, 3), '.', 'markerSize', 15);
title("Trajectory",'FontSize',12,'FontName','Times New Roman');
for n=1:100:size(position, 1)
    u = rotatepoint(orientation(n), [0.1 0 0]);
    quiver3(position(n,1),position(n,2),position(n,3),u(1),u(2),u(3),'r');
    u = rotatepoint(orientation(n), [0 0.1 0]);
    quiver3(position(n,1),position(n,2),position(n,3),u(1),u(2),u(3),'b')
    u = rotatepoint(orientation(n), [0 0 0.1]);
    quiver3(position(n,1),position(n,2),position(n,3),u(1),u(2),u(3),'g')
end
axis equal;
grid minor;
xlabel('x [m]','FontSize',12,'FontName','Times New Roman')
ylabel('y [m]','FontSize',12,'FontName','Times New Roman')
zlabel('z [m]','FontSize',12,'FontName','Times New Roman')
saveas(gcf, 'figures/groundtruth', 'epsc');

figure;
plot3(X(:, stateMasks.pos(1)), X(:, stateMasks.pos(2)), X(:, stateMasks.pos(3)), 'b');
hold on;
plot3(X(end, stateMasks.pos(1)), X(end, stateMasks.pos(2)), X(end, stateMasks.pos(3)), '.', 'markerSize', 15);
title("ESKF Estimated Trajectory",'FontSize',12,'FontName','Times New Roman');
for n=1:100:size(X, 1)
    u = rotatepoint(quaternion(X(n, stateMasks.q_nb)), [0.1 0 0]);
    quiver3(X(n, stateMasks.pos(1)), X(n, stateMasks.pos(2)), X(n, stateMasks.pos(3)), u(1), u(2), u(3), 'r');
    u = rotatepoint(quaternion(X(n, stateMasks.q_nb)), [0 0.1 0]);
    quiver3(X(n, stateMasks.pos(1)), X(n, stateMasks.pos(2)), X(n, stateMasks.pos(3)), u(1), u(2), u(3), 'b');
    u = rotatepoint(quaternion(X(n, stateMasks.q_nb)), [0 0 0.1]);
    quiver3(X(n, stateMasks.pos(1)), X(n, stateMasks.pos(2)), X(n, stateMasks.pos(3)), u(1), u(2), u(3), 'g');
end
axis equal;
grid minor;
xlabel('x [m]','FontSize',12,'FontName','Times New Roman');
ylabel('y [m]','FontSize',12,'FontName','Times New Roman');
zlabel('z [m]','FontSize',12,'FontName','Times New Roman');
saveas(gcf, 'figures/eskf_estimated_trajectory', 'epsc');