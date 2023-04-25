function [SpeakOut, ListenOut] = overlap_windows2(SpeakIn, ListenIn, Fs)

% Initialize output arrays
SpeakOut = SpeakIn;
ListenOut = ListenIn;

% Create all possible pairs of speaking and listening windows
pairList = [];
for i = 1:size(SpeakIn,1)
    for j = 1:size(ListenIn,1)
        pairList(end+1,:) = [i j];
    end
end

% Loop over all pairs of speaking and listening windows
for p = 1:size(pairList,1)
    i = pairList(p,1);
    j = pairList(p,2);
    
    if SpeakIn(i,2) <= ListenIn(j,2) && SpeakIn(i,3) >= ListenIn(j,3)
        % Speaking window entirely contains the listening window
        % Check if there's any future windows inside
        if size(SpeakIn,1) > i && SpeakIn(i+1,2) <= ListenIn(j,3)
            SpeakOut(i,3) = SpeakIn(i+1,2);
        else
            SpeakOut(i,3) = ListenIn(j,2);
        end
        
        % Add the new speaking window to the end of output array
        SpeakOut(end+1,:) = [(SpeakIn(i,3)-SpeakOut(i,3))/Fs, SpeakOut(i,3), SpeakIn(i,3)];
        
        % Check if there's any future windows inside
        if size(ListenIn,1) > j && ListenIn(j+1,2) <= SpeakOut(end,3)
            ListenOut(end,3) = ListenIn(j+1,2);
        else
            ListenOut(end,3) = SpeakOut(end,3);
        end
        
    elseif ListenIn(j,2) <= SpeakIn(i,2) && ListenIn(j,3) >= SpeakIn(i,3)
        % Listening window entirely contains the speaking window
        % Check if there's any future windows inside
        if size(ListenIn,1) > j && ListenIn(j+1,2) <= SpeakIn(i,3)
            ListenOut(j,3) = ListenIn(j+1,2);
        else
            ListenOut(j,3) = SpeakIn(i,2);
        end
        
        % Add the new speaking window to the end of output array
        ListenOut(end+1,:) = [(ListenIn(j,3)-ListenOut(j,3))/Fs, ListenOut(j,3), ListenIn(j,3)];
        
        % Check if there's any future windows inside
        if size(SpeakIn,1) > i && SpeakIn(i+1,2) <= ListenOut(end,3)
            SpeakOut(end,3) = SpeakIn(i+1,2);
        else
            SpeakOut(end,3) = ListenOut(end,3);
        end
        
    elseif SpeakIn(i,3) >= ListenIn(j,2) && SpeakIn(i,3) <= ListenIn(j,3)
        % If there is an overlap, set the end of the speaking window to the start of the listening window
        SpeakOut(i,3) = ListenIn(j,2);
        
        % Check if there's
    end