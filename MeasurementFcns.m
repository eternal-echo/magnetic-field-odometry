function [yk] = MeasurementFcns(x, r)
    sensor_locs = [[r/2; r; 0] [-r/2; r; 0] [r/2; 0; 0] [-r/2; 0; 0] [r/2; -r; 0] [-r/2; -r; 0]];
    persistent H
    % global invAB;
    if isempty(H)
            
        H = [calcAB(sensor_locs(:, 1)); ...
         calcAB(sensor_locs(:, 2)); ...
         calcAB(sensor_locs(:, 3)); ...
         calcAB(sensor_locs(:, 4)); ...
         calcAB(sensor_locs(:, 5)); ...
         calcAB(sensor_locs(:, 6))];
    end
    

    yk = H * x(end-14:end) + [x(10:24); zeros(3, 1)];
end





