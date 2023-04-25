%% Function to fix overlap between windows by: Initial window END = Secondary window START
function [SpeakOut, ListenOut] = overlap_windows(SpeakIn, ListenIn, Fs)
% Initialize output arrays
SpeakOut = SpeakIn;
ListenOut = ListenIn;

% Loop over all pairs of speaking and listening windows
for i = 1:size(SpeakIn,1)
    for j = 1:size(ListenIn,1)
        if SpeakIn(i,2) <= ListenIn(j,2) && SpeakIn(i,2) <= ListenIn(j,3) && SpeakIn(i,3) >= ListenIn(j,2) && SpeakIn(i,3) >= ListenIn(j,3)
            % Speaking window entirely contains the listening window
            SpeakOut(i,3) = ListenIn(j,2);
            
            % Add the new speaking window to the end of output array
            SpeakOut(end+1,:) = [(SpeakIn(i,3)-ListenIn(j,3))/Fs, ListenIn(j,3), SpeakIn(i,3)];
            
        elseif ListenIn(j,2) <= SpeakIn(i,2) && ListenIn(j,2) <= SpeakIn(i,3) && ListenIn(j,3) >= SpeakIn(i,2) && ListenIn(j,3) >= SpeakIn(i,3)
            % Speaking window entirely contains the listening window
            ListenOut(i,3) = SpeakIn(i,2);
            
            % Add the new speaking window to the end of output array
            ListenOut(end+1,:) = [(ListenIn(j,3)-SpeakIn(i,3))/Fs, SpeakIn(i,3), ListenIn(j,3)];
            
        elseif SpeakIn(i,3) >= ListenIn(j,2) && SpeakIn(i,3) <= ListenIn(j,3)
            % If there is an overlap, set the end of the speaking window to the start of the listening window
            SpeakOut(i,3) = ListenIn(j,2);
            
        elseif ListenIn(j,3) >= SpeakIn(i,2) && ListenIn(j,3) <= SpeakIn(i,3)
            % If there is an overlap, set the end of the listening window to the start of the speaking window
            ListenOut(i,3) = SpeakIn(i,2); 
            
        end
    end
end


% Ensure that start times are not greater than end times
SpeakOut(SpeakOut(:,2) > SpeakOut(:,3),2:3) = SpeakOut(SpeakOut(:,2) > SpeakOut(:,3),[3 2]);
ListenOut(ListenOut(:,2) > ListenOut(:,3),2:3) = ListenOut(ListenOut(:,2) > ListenOut(:,3),[3 2]);

% Sort arrays
SpeakOut=sortrows(SpeakOut,[2 3]);
ListenOut=sortrows(ListenOut,[2 3]);

% Change durations accordingly
SpeakOut(:,1) = (SpeakOut(:,3)-SpeakOut(:,2))/Fs;
ListenOut(:,1) = (ListenOut(:,3)-ListenOut(:,2))/Fs;

end

