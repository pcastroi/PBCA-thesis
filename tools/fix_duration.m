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

% Convert the maximum visual angle to radians
MaxAngle_rad = deg2rad(MaxVisualAngle);

% Initialize the output variable
fixation = zeros(size(Gaze, 1), 1);

% Loop over each gaze sample
for i = 1:size(Gaze, 1)
    % Get the current gaze sample
    current_gaze = Gaze(i, :);
    
    % If current_gaze = NaN -> fixation = NaN
    if isnan(current_gaze)
       fixation(i) = NaN;
       continue;
    end
    
    % From deg to mm - Sphere radius
    Sph_r=current_gaze(3)*tan(deg2rad(MaxVisualAngle)/2);
    
    % Initialize counters for pre- and post-samples
    pre_samples = 0;
    post_samples = 0;

    % Loop over all previous samples
    for j = i-1:-1:1
        % Calculate the distance between the current gaze sample and the previous sample
        distance = norm(current_gaze - Gaze(j, :));

        % If the distance is greater than the sphere radius, stop counting pre-samples
        if distance > Sph_r || isnan(distance)
            break;
        end

        % If the distance is within the sphere radius, increment the pre-sample counter
        pre_samples = pre_samples + 1;
    end

    % Loop over all subsequent samples
    for j = i+1:size(Gaze, 1)
        % Calculate the distance between the current gaze sample and the subsequent sample
        distance = norm(current_gaze - Gaze(j, :));

        % If the distance is greater than the sphere radius, stop counting post-samples
        if distance > Sph_r || isnan(distance)
            break;
        end

        % If the distance is within the sphere radius, increment the post-sample counter
        post_samples = post_samples + 1;
    end

    % Count the number of continuous pre- and post-samples within the sphere
    fixation(i) = (pre_samples + post_samples)/Fs;
end

end

