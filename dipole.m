% 计算磁偶极子位于`pref`时，在点`p`处的磁场向量`h`的函数。
%
% 输入:
%   p    - 3x1 向量，表示计算磁场的位置。
%   pref - 3x1 向量，表示磁偶极子的位置。
%   m    - 3x1 向量，表示磁偶极矩。
%
% 输出:
%   h    - 3x1 向量，表示点`p`处的磁场。
%
% 磁场`h`的计算公式为:
%   h = (3 * (r * r') / norm(r)^5 - eye(3) / norm(r)^3) * m
% 其中`r`是从`pref`到`p`的向量（即 r = p - pref）。
% 该公式来源于磁偶极子的理论磁场表达式。
function h=dipole(p,pref,m)
    % p: 观测点的位置向量 (1x3 vector)
    % pref: 磁偶极子的位置向量 (1x3 vector)
    % m: 磁偶极子的磁矩向量 (1x3 vector)

    % 常数 k = mu_0 / (4 * pi)
    k = 1e-7;  % mu_0 / (4 * pi) 的常数值

    % 计算 r 向量和其模长 norm(r)
    r=p-pref; % 从偶极子位置指向观测点的向量
    r_norm = norm(r); % r 的模长

    if r_norm == 0
        h = [0, 0, 0]; % 避免在偶极子位置处除以零
        return;
    end

    % 单位方向向量 r_hat
    r_hat = r / r_norm;

    % 磁场计算公式
    h = k * (3 * (dot(m, r) * r) / r_norm^5 - m / r_norm^3);
    % h = k * (3. * (r*r')./norm(r)^5-eye(3)./norm(r)^3)*m;
end
