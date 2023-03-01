function [y,param] = eyeVelocity(x,param)

% calculate the euklidian distance (square or abs?)
% Abs is needed as pitch angle ia in the range of [-90,90]
dx = [abs(diff(x,1)); 0];

% calculate the velocity (resampled fsr if resampled)
y = dx*param.fsr;

% replace first value of the velocity by avg
y(1) = mean(y,1);

% save pitch angle
param.pitch = x;
param.velocity = y; 

end

% eof