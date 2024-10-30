load('generated_data.mat');
settings = getSettings();
timeVector = 0:1/settings.fs:(settings.duration-1/settings.fs);

% %% generate magnetic field

% [ImuMag_data, ImuMag_bias, theta_cell, aux] = sensor_data_gen(settings, position, orientation, acceleration, angularVelocity, 4);
% save('data/generated_data.mat', 'position', 'orientation', 'velocity', 'acceleration', 'angularVelocity', ...
%     'ImuMag_data', 'ImuMag_bias', 'theta_cell', 'aux', 'settings', 'timeVector');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             Magnetic field              %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
tmp = load(settings.model);

m.moments = tmp.mm(:, 1:end-2);
m.pos_dipoles = tmp.MM;
m.f_earth = tmp.mm(:, end-1);
clear tmp
[Xq,Yq,Zq]=meshgrid(-1:0.4:1,-1:0.4:1,0:0.2:1);
m.pos_grid=[reshape(Xq,1,numel(Xq)); reshape(Yq,1,numel(Yq)); reshape(Zq,1,numel(Zq))];

% Allocate memory
m.f_grid=m.f_earth*ones(1,size(m.pos_grid,2));

% Location loop
for ii=1:size(m.pos_grid,2)
    % Dipole loop
    for kk=1:size(m.moments,2)
        m.f_grid(:,ii)= m.f_grid(:,ii)+dipole(m.pos_grid(:,ii),m.pos_dipoles(:,kk),m.moments(:,kk));
    end
end

quiverC3D(m.pos_grid(1,:)',m.pos_grid(2,:)',m.pos_grid(3,:)',m.f_grid(1,:)',m.f_grid(2,:)',m.f_grid(3,:)','Colorbar',true,'LineWidth',0.5);

title('Generated magnetic-field','FontSize',12,'FontName','Times New Roman')
xlabel('x [m]','FontSize',12,'FontName','Times New Roman')
ylabel('y [m]','FontSize',12,'FontName','Times New Roman')
zlabel('z [m]','FontSize',12,'FontName','Times New Roman')
axis tight;

print('figures/magnetic-field', '-depsc', '-r600', '-painters' );

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%             trajectory plot             %% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% plot3(position(:, 1), position(:, 2), position(:, 3), 'k');
% hold on;
% plot3(position(end, 1), position(end, 2), position(end, 3), '.', 'markerSize', 15);
% title("Trajectory",'FontSize',12,'FontName','Times New Roman');
% for n=1:100:size(position, 1)
%     u = rotatepoint(orientation(n), [0.1 0 0]);
%     quiver3(position(n,1),position(n,2),position(n,3),u(1),u(2),u(3),'r');
%     u = rotatepoint(orientation(n), [0 0.1 0]);
%     quiver3(position(n,1),position(n,2),position(n,3),u(1),u(2),u(3),'b')
%     u = rotatepoint(orientation(n), [0 0 0.1]);
%     quiver3(position(n,1),position(n,2),position(n,3),u(1),u(2),u(3),'g')
% end
% axis equal;
% grid minor;
% xlabel('x [m]','FontSize',12,'FontName','Times New Roman')
% ylabel('y [m]','FontSize',12,'FontName','Times New Roman')
% zlabel('z [m]','FontSize',12,'FontName','Times New Roman')
% saveas(gcf, 'figures/groundtruth', 'epsc');