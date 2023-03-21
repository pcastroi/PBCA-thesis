function fixation = fix_duration(Gaze,MaxVisualAngle,Fs)
%FIXATION_DURATION - Count the number of samples before and after each gaze sample that fall within a sphere of a given radius.
%
% Inputs:
%   - Gaze: a Nx3 matrix containing the x, y, and z coordinates of each gaze sample
%   - MaxAngle: the maximum visual angle (in degrees) for the sphere radius
%   - Fs: sampling frequency of Gaze vector
%
% Output:
%   - fixation_duration: number of samples within the sphere for each gaze sample

% Initialize the output variable
fixation = zeros(size(Gaze, 1), 1);
radius = zeros(size(Gaze, 1), 1);


% Loop over each gaze sample
for i = 1:size(Gaze, 1)
    % Get the current gaze sample
    current_gaze_xy = Gaze(i, 1:2);
    current_gaze_z = Gaze(i, 3);
    
    % If current_gaze = NaN -> fixation = NaN
    if isnan(current_gaze_xy)
       fixation(i) = NaN;
       continue;
    end
    
    % Radius calculated from MaxVisualAngle and distance
    b = current_gaze_z; % from 0,0,0 to GazeX,GazeY,GazeZ
    c = b/cos(deg2rad(MaxVisualAngle)/2); % from 0,0,0 to GazeX + radius,GazeY,GazeZ
    radius(i) = 2*sqrt(c^2-b^2);

%     radius(i)=4;
            
    % Initialize counters for pre- and post-samples
    pre_samples = 0;
    post_samples = 0;

    % Loop over all previous samples
    for j = i-1:-1:1
        % Calculate the distance between the current gaze sample and the previous sample
        distance = norm(current_gaze_xy - Gaze(j, 1:2));
        
        % If a data value of any pre- or post sample within the 0.5° visual angle location had been coded NaN (i.e., was missing), the fixation duration of the corresponding time point was set to NaN
        if isnan(distance)
            fixation(i) = NaN;
            break;
        end
        
        % If the distance is greater than the sphere radius, stop counting pre-samples
        if distance > radius(i)
            break;
        end

        % If the distance is within the sphere radius, increment the pre-sample counter
        pre_samples = pre_samples + 1;
    end

    % Loop over all subsequent samples
    for j = i+1:size(Gaze, 1)
        % Calculate the distance between the current gaze sample and the subsequent sample
        distance = norm(current_gaze_xy - Gaze(j, 1:2));

        % If a data value of any pre- or post sample within the 0.5° visual angle location had been coded NaN (i.e., was missing), the fixation duration of the corresponding time point was set to NaN
        if isnan(distance)
            fixation(i) = NaN;
            break;
        end
        
        % If the distance is greater than the sphere radius, stop counting post-samples
        if distance > radius(i)
            break;
        end

        % If the distance is within the sphere radius, increment the post-sample counter
        post_samples = post_samples + 1;
    end
    
    % Count the number of continuous pre- and post-samples within the sphere
    if ~isnan(fixation(i))
        if pre_samples + post_samples ~= 0
            fixation(i) = (pre_samples + post_samples)/Fs;
        else
            fixation(i) = NaN;
        end
    end
end
% % 3D Plot
% tic
% figure
% plot3(Gaze(:,1),Gaze(:,2),Gaze(:,3),'b*')
% hold on
% plot3(0,0,0,'r*')
% grid on
% for j = 1:length(radius)
%     [x,y,z] = sphere(10);
%     x = x * radius(j);
%     y = y * radius(j);
%     z = z * radius(j);
%     surf(x+Gaze(j,1),y+Gaze(j,2),z+Gaze(j,3),'FaceColor','none','EdgeColor','r')
% end
% toc
% 
% % 2D Plot
% figure
% plot(Gaze(:,1),Gaze(:,2),'b*')
% hold on
% plot(0,0,'r*')
% grid on
% viscircles(Gaze(:, 1:2),radius);

end

