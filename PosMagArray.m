function r=PosMagArray()
% % spatial locations of magnetometer sensors
%   r = NaN(3, 30);
%   dx = 0.064;
%   dy = 0.055;
%   kk = 0;

%   for jj=1:5
%     for ii=1:6
%       kk=kk+1;

%       r(1, kk) = (ii-3.5)*dx;
%       r(2, kk) = -(jj-3)*dy;
%       r(3, kk) = 0;
%     end
%   end

    % 修改为一个传感器
    r = NaN(3, 1);

    % 设置传感器位置
    r(:, 1) = [0; 0; 0];  % 传感器位置 (x, y, z)
end
