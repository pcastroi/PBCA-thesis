function [Merged] = MergeWin(Raw, Fs, TimeMergeGap)
PreWin = 0;
    for j=1:size(Raw,1)
        if j < size(Raw,1)
            if (Raw(j+1,2)-Raw(j,3))/Fs <= TimeMergeGap
                Merged(j,2) = Raw(j,2);
                if j < size(Raw,1)-1
                    if (Raw(j+2,2)-Raw(j+1,3))/Fs <= TimeMergeGap
                        Merged(j,3) = Raw(j+2,3);
                        Merged(j,1) = (Merged(j,3)-Merged(j,2))/Fs;
                        PreWin = j+2;
                    else
                        Merged(j,3) = Raw(j+1,3);
                        Merged(j,1) = (Merged(j,3)-Merged(j,2))/Fs;
                        PreWin = j+1;
                    end
                else
                    Merged(j,3) = Raw(j+1,3);
                    Merged(j,1) = (Merged(j,3)-Merged(j,2))/Fs;
                    PreWin = j+1;
                end
            elseif j ~= PreWin
                Merged(j,:) = Raw(j,:);
            end
        elseif j > 1
            if Merged(end,2)~=Raw(j-1,2) && Merged(end,3)~=Raw(j,3)
                Merged(j,:) = Raw(j,:);
            end
        else
            Merged(j,:) = Raw(j,:);
        end
    end

    Merged(~any(Merged,2),:)=[];
    % Merge all
    Merged(Merged(1:end-1,3)>Merged(2:end,2),3) = Merged(find(Merged(1:end-1,3)>Merged(2:end,2))+1,3);
    Merged(find(Merged(1:end-1,3)>Merged(2:end,2))+1,:)=[];
end
