%% Function to fix overlap between windows by: Initial window END = Secondary window START
function [SpeakOut, ListenOut] = overlap_windows(SpeakIn, ListenIn, Fs)
% Remove duplicate windows (Speaking inside Speaking and viceversa)

% Initialize output arrays
SpeakOut = SpeakIn;
ListenOut = ListenIn;

NumSW = size(SpeakIn,1);
NumLW = size(ListenIn,1);
% Loop over all pairs of speaking and listening windows
i=1;
while i <= NumSW
    j=1;
    while j <= NumLW
        if i <= NumSW && j <= NumLW
            % CASE 1: Speaking window contains at least 1 listening window
            if SpeakIn(i,2) <= ListenIn(j,2) && SpeakIn(i,2) <= ListenIn(j,3) && SpeakIn(i,3) >= ListenIn(j,2) && SpeakIn(i,3) >= ListenIn(j,3)      
                SWCont=sum(SpeakIn(i,2) <= ListenIn(:,2) & SpeakIn(i,3) >= ListenIn(:,3));
                SWCFind=find(SpeakIn(i,2) <= ListenIn(:,2) & SpeakIn(i,3) >= ListenIn(:,3));
                for k = 1:SWCont
                    if k == 1
                        SpeakOut(i,3) = ListenIn(SWCFind(k),2);
                    end

%                     if k == SWCont
%                         SpeakOut(end+1,:) = [(SpeakIn(i,3)-ListenIn(SWCFind(k),3))/Fs, ListenIn(SWCFind(k),3), SpeakIn(i,3)];
%                     else
                        ClosestW = ListenIn(SWCFind(k),3) + min(abs([SpeakIn(:,2:3);ListenIn(1:end ~= SWCFind(k),2:3)]-ListenIn(SWCFind(k),3)),[],1:2);
                        SpeakOut(end+1,:) = [(ClosestW-ListenIn(SWCFind(k),3))/Fs, ListenIn(SWCFind(k),3), ClosestW];
%                     end
                end
                j=j+k;
            % CASE 2: Listening window contains at least 1 speaking window
            elseif ListenIn(j,2) <= SpeakIn(i,2) && ListenIn(j,2) <= SpeakIn(i,3) && ListenIn(j,3) >= SpeakIn(i,2) && ListenIn(j,3) >= SpeakIn(i,3)
                LWCont=sum(ListenIn(j,2) <= SpeakIn(:,2) & ListenIn(j,3) >= SpeakIn(:,3));
                LWCFind=find(ListenIn(j,2) <= SpeakIn(:,2) & ListenIn(j,3) >= SpeakIn(:,3));
                for k = 1:LWCont
                    if k == 1
                        ListenOut(j,3) = SpeakIn(LWCFind(k),2);
                    end

%                     if k == LWCont
%                         ListenOut(end+1,:) = [(ListenIn(j,3)-SpeakIn(LWCFind(k),3))/Fs, SpeakIn(LWCFind(k),3), ListenIn(j,3)];
%                     else
%                         [ClosestWRow,ClosestWCol] = find((ListenIn(:,2:3)>=SpeakIn(LWCFind(k),3)));
                        ClosestW = min(ListenIn(ListenIn>=SpeakIn(LWCFind(k),3)));
%                         ClosestW = SpeakIn(LWCFind(k),3) + min(abs([ListenIn(:,2:3);SpeakIn(1:end ~= LWCFind(k),2:3)]-SpeakIn(LWCFind(k),3)),[],1:2);
                        ListenOut(end+1,:) = [(ClosestW-SpeakIn(LWCFind(k),3))/Fs, SpeakIn(LWCFind(k),3), ClosestW];
%                     end
                end
                i=i+k;
            % CASE 3: Overlap, first Speaking, then Listening
            elseif SpeakIn(i,3) >= ListenIn(j,2) && SpeakIn(i,3) <= ListenIn(j,3) && SpeakIn(i,2) <= ListenIn(j,2) && SpeakIn(i,2) <= ListenIn(j,3)
                SpeakOut(i,3) = ListenIn(j,2);
            % CASE 4: Overlap, first Listening, then Speaking
            elseif ListenIn(j,3) >= SpeakIn(i,2) && ListenIn(j,3) <= SpeakIn(i,3) && ListenIn(j,2) <= SpeakIn(i,2) && ListenIn(j,2) <= SpeakIn(i,3)
                ListenOut(j,3) = SpeakIn(i,2);
            end
        end
        j=j+1;
    end
        i=i+1;
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

