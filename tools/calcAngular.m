function [az, el] = calcAngular(x,y,z)
[az,el] = cart2sph(x,z,y);   % matlab's Z and Y are reversed compared to tobiis
    az  =  az*180/pi-90;    % I have checked sign and offset of azi and ele so that things match the gaze position on the scene video in the data file (gp)
    el  = el*180/pi;
end