%% Function to merge windows whenever the gap in between the end of a window and the start of the next one is < TimeMerge
function Merge = merge_windows(Raw, Fs, TimeMerge)
% Initialize variables
num_windows = size(Raw, 1);
merged_windows = [];

% Loop through the windows and merge adjacent windows if necessary
i = 1;
while i <= num_windows
    window_start = Raw(i, 2);
    window_end = Raw(i, 3);
    window_duration = (window_end - window_start) / Fs;
    
    % Find the closest start time of a window within the TimeMerge
    closest_start = Inf;
    for j = i+1:num_windows
        start_time = Raw(j, 2);
        if start_time - window_end <= TimeMerge*Fs && start_time > window_end
            closest_start = start_time;
            break;
        end
    end
    
    % If there is a window within the TimeMerge and the gap between the end of this
    % window and the start of the next one is less than TimeMerge, merge them
    if closest_start ~= Inf
        gap_duration = (closest_start - window_end) / Fs;
        if gap_duration <= TimeMerge
            merged_start = window_start;
            idx_merged = find(Raw(:,2) == closest_start, 1);
            merged_end = Raw(idx_merged, 3);
            merged_duration = (merged_end - merged_start) / Fs;
            merged_windows = [merged_windows; merged_duration, merged_start, merged_end];
            i = idx_merged + 1;
            
            % Check if the last merged window can be merged with the next one
            while i <= num_windows
                next_start = Raw(i, 2);
                next_end = Raw(i, 3);
                next_duration = Raw(i, 1);
                if next_start - merged_end <= TimeMerge*Fs
                    merged_end = next_end;
                    merged_duration = (merged_end - merged_start) / Fs;
                    merged_windows(end, :) = [merged_duration, merged_start, merged_end];
                    i = i + 1;
                else
                    break;
                end
            end
        else
            % Otherwise, add the current window to the list of merged windows
            merged_windows = [merged_windows; [window_duration, window_start, window_end]];
            i = i + 1;
        end
    else
        % Otherwise, add the current window to the list of merged windows
        merged_windows = [merged_windows; [window_duration, window_start, window_end]];
        i = i + 1;
    end
end

% Output the list of merged windows
Merge = merged_windows;

end
