function [param] = eyeAngle(x,param)


% calculate yaw
yaw = atan2(x(:,1),x(:,3));
yaw = yaw-mean(yaw,1);
param.yaw = yaw*180/pi;

% calculate pitch
x4= sqrt(x(:,1).^2+x(:,3).^2);
pitch = atan2(x(:,2),x4);
% pitch = x(:,2);
pitch = pitch-mean(pitch,1);
param.pitch = pitch*180/pi;
% param.pitch = pitch;

% calculate roll
param.roll = [];

if param.debug == 1
    figure;
    histogram(param.pitch,30);
end

end

% eof