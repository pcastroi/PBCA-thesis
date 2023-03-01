function fixation = fix_duration(Gaze,MaxAngle)
%   FIXATION_DURATION - Obtain the fixation duration from 3D gaze data
% calculate yaw
yaw = atan2(Gaze(:,1),Gaze(:,3));
yaw = yaw-mean(yaw,'omitnan');
angle = yaw*180/pi;

fixation = 0;
end

