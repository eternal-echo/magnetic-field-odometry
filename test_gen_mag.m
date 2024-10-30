% Define grid for visualization
[x, y, z] = meshgrid(-5:0.1:5, -5:0.1:5, -5:0.1:5);
grid_points = [x(:), y(:), z(:)];  % Flatten grid for vectorized calculation

% Define dipole positions and moments
dipole_positions = [0, 0, 0]; % four dipoles in a square arrangement
dipole_moments = [0, 0, 1];  % moments aligned along z-axis

% Initialize magnetic field components
Hx = zeros(size(x));
Hy = zeros(size(y));
Hz = zeros(size(z));

% Calculate the magnetic field at each grid point
for i = 1:length(grid_points)
    p = grid_points(i, :);  % 观测点的位置
    H = [0, 0, 0];  % Initialize field at this point
    
    % Sum field contributions from each dipole
    for j = 1:size(dipole_positions, 1)
        pref = dipole_positions(j, :); % 磁偶极子的位置
        m = dipole_moments(j, :)'; % 磁矩
        H = H + dipole(p, pref, m)';  % Using the provided dipole function
    end
    
    % Store the computed field components
    Hx(i) = H(1);
    Hy(i) = H(2);
    Hz(i) = H(3);
end

% Reshape for 3D plotting
Hx = reshape(Hx, size(x));
Hy = reshape(Hy, size(y));
Hz = reshape(Hz, size(z));

% Visualize the field with quiver3 for a 3D vector field plot
figure;
quiver3(x, y, z, Hx, Hy, Hz, 'AutoScale', 'on', 'AutoScaleFactor', 1.5);
% quiverC3D(x(:), y(:), z(:), Hx(:), Hy(:), Hz(:), 'Colorbar', true, 'LineWidth', 0.5);
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Magnetic Field of Four Magnetic Dipoles');
grid on;

% Reshape for 3D plotting
Hx = reshape(Hx, size(x));
Hy = reshape(Hy, size(y));
Hz = reshape(Hz, size(z));

% Define streamline starting points
[startX, startY, startZ] = meshgrid(-1.5:1:1.5, -1.5:1:1.5, -1.5:1:1.5);

% Visualize the field with streamlines in 3D
figure;
streamline(x, y, z, Hx, Hy, Hz, startX, startY, startZ);
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Magnetic Field of Two Magnetic Dipoles (Streamline Visualization)');
grid on;
view(3);