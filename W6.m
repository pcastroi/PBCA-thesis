%% PBCA-Thesis - Week 6 & Week 7 - Plot histograms with lenght of gaps, turns, pauses, overlap-between and overlap-within - Multiple Durations
% Pathing
clear all; clc; close all;
BPath = strsplit(pwd,'PBCA-thesis');
addpath('tools')
addpath([BPath{1} 'Pupil-preprocessing-tools\tools']) % For preprocessing

% Colors
SpeakColor = [53, 155, 67]./255;
ListenColor = [204, 36, 0]./255;
GapColor=[126, 206, 253]./255;
OLWColor=[35, 36, 139]./255;
OLBColor=[139, 86, 35]./255;
TurnColor=[234, 111, 250]./255;
PauseColor=[77, 100, 111]./255;

[subDirs_I] = GetSubDirsFirstLevelOnly('data\AMEND_I');
FileNames_I={'P1_Quiet_B1.mat','P1_Quiet_B2.mat','P1_SHL_B1.mat','P1_SHL_B2.mat','P1_Noise60_B1.mat','P1_Noise60_B2.mat','P1_Noise70_B1.mat','P1_Noise70_B2.mat','P2_Quiet_B1.mat','P2_Quiet_B2.mat','P2_SHL_B1.mat','P2_SHL_B2.mat','P2_Noise60_B1.mat','P2_Noise60_B2.mat','P2_Noise70_B1.mat','P2_Noise70_B2.mat'};
LoadUtt_I=load('data\AMEND_I\utterances1110.mat');
LoadDelays_I=load('data\AMEND_I\delays1110.mat');
LoadTPsOrder_I=load('data\AMEND_I\TPsOrder_I.mat');

[subDirs_II] = GetSubDirsFirstLevelOnly('data\AMEND_II');
subDirs_II(contains(subDirs_II,{'Pilot 1','Pair01','Pair14'}))=[]; % Only using data from Pair01 to Pair13, others removed
FileNames_II={'UNHI_N0.mat','UNHI_N60.mat','UNHI_N70.mat','AAHI_N0.mat','AAHI_N60.mat','AAHI_N70.mat','ABHI_N0.mat','ABHI_N60.mat','ABHI_N70.mat'};
LoadUtt_II=load('data\AMEND_II\Utterances0805.mat');
LoadTPsOrder_II=load('data\AMEND_II\TPsOrder_II.mat');

TimeMinWin = 0.5; % [s], Minimum time of a window
TimeInitialMerge = 0.3; % [s], Time threshold for merging windows initially
TimeMerge = 2; % [s], Time threshold for merging windows after rejecting small windows
Param.Fs = 250; % Sampling frequency of Utt data
NPs = 2; % N of TPs per Pair
NCond_I = numel(FileNames_I)/NPs; % N of conditions

% Preallocate durations
GapDurRaw_I=zeros(200,1000);
GapDur_I=GapDurRaw_I;
OLWDurRaw_I=GapDurRaw_I;
OLWDur_I=GapDurRaw_I;
OLBDurRaw_I=GapDurRaw_I;
OLBDur_I=GapDurRaw_I;
TurnDurRaw_I=GapDurRaw_I;
TurnDur_I=GapDurRaw_I;
TurnDiff_I=GapDurRaw_I;
TurnDurQuiet_I=GapDurRaw_I;
TurnDurSHL_I=GapDurRaw_I;
TurnDurN60_I=GapDurRaw_I;
TurnDurN70_I=GapDurRaw_I;
PauseDurRaw_I=GapDurRaw_I;
PauseDur_I=GapDurRaw_I;
SpeakDurRaw_I=GapDurRaw_I;
ListenDurRaw_I=GapDurRaw_I;
SpeakDur_I=GapDurRaw_I;
ListenDur_I=GapDurRaw_I;

GapDurRaw_II=zeros(200,1000);
GapDur_II=GapDurRaw_II;
OLWDurRaw_II=GapDurRaw_II;
OLWDur_II=GapDurRaw_II;
OLBDurRaw_II=GapDurRaw_II;
OLBDur_II=GapDurRaw_II;
TurnDurRaw_II=GapDurRaw_II;
TurnDur_II=GapDurRaw_II;
TurnDiff_II=GapDurRaw_II;
TurnDurQuiet_II=GapDurRaw_II;
TurnDurSHL_II=GapDurRaw_II;
TurnDurN60_II=GapDurRaw_II;
TurnDurN70_II=GapDurRaw_II;
PauseDurRaw_II=GapDurRaw_II;
PauseDur_II=GapDurRaw_II;
SpeakDurRaw_II=GapDurRaw_II;
ListenDurRaw_II=GapDurRaw_II;
SpeakDur_II=GapDurRaw_II;
ListenDur_II=GapDurRaw_II;

nbins=100;
x_I=1;
x_II=1;

% AMEND I
for q=1:numel(subDirs_I)
    PairIn_I = q;
    PairFiles_I=dir(['data\AMEND_I\Main',sprintf('%d',PairIn_I),'\*.mat']);
    PairUtt_I=LoadUtt_I.Utterances(PairIn_I,:);
    PairDelay_I=LoadDelays_I.TobAudDelay(PairIn_I,:);
    
    for i=1:numel(FileNames_I)
        if contains(cell2mat(FileNames_I(i)),'P2')
            if LoadTPsOrder_I.TPsOrder(2*q,i-NCond_I) ~= x_I
                disp(['Warning: File ',PairFiles_I(1).folder, '\', cell2mat(FileNames_I(i)),' was rejected in AMEND I (P2) analysis.']);
                continue
            end
        elseif contains(cell2mat(FileNames_I(i)),'P1')
            if LoadTPsOrder_I.TPsOrder(2*q-1,i) ~= x_I
                disp(['Warning: File ',PairFiles_I(1).folder, '\', cell2mat(FileNames_I(i)),' was rejected in AMEND I (P1) analysis.']);
                continue
            end
        end
        % Retrieve Utterances
        if contains(cell2mat(FileNames_I(i)),'P2')
            SpeakKey = 'utteranceCH1';
            ListenKey = 'utteranceCH2';
            GapKey = 'gapCH1';
            OLWKey = 'overlapWCH1';
            OLBKey = 'overlapBCH1';
            TurnKey = 'turnCH1';
            PauseKey = 'pauseCH1';
        elseif contains(cell2mat(FileNames_I(i)),'P1')
            SpeakKey = 'utteranceCH2';
            ListenKey = 'utteranceCH1';
            GapKey = 'gapCH2';
            OLWKey = 'overlapWCH2';
            OLBKey = 'overlapBCH2';
            TurnKey = 'turnCH2';
            PauseKey = 'pauseCH2';
        end

        if contains(cell2mat(FileNames_I(i)),'B1')
            SpeB = 0;
        elseif contains(cell2mat(FileNames_I(i)),'B2')
            SpeB = 1;
        end

        if contains(cell2mat(FileNames_I(i)),'Quiet')
            SpeCond = SpeB + 1;
            TurnDurQuiet_I(x_I,1:size(PairUtt_I{1,SpeCond}.(TurnKey)(PairUtt_I{1,SpeCond}.(TurnKey)(:,1)'>TimeMinWin),2)) = PairUtt_I{1,SpeCond}.(TurnKey)(PairUtt_I{1,SpeCond}.(TurnKey)(:,1)'>TimeMinWin);
        elseif contains(cell2mat(FileNames_I(i)),'SHL')
            SpeCond = SpeB + 3;
            TurnDurSHL_I(x_I,1:size(PairUtt_I{1,SpeCond}.(TurnKey)(PairUtt_I{1,SpeCond}.(TurnKey)(:,1)'>TimeMinWin),2)) = PairUtt_I{1,SpeCond}.(TurnKey)(PairUtt_I{1,SpeCond}.(TurnKey)(:,1)'>TimeMinWin);
        elseif contains(cell2mat(FileNames_I(i)),'Noise60')
            SpeCond = SpeB + 5;
            TurnDurN60_I(x_I,1:size(PairUtt_I{1,SpeCond}.(TurnKey)(PairUtt_I{1,SpeCond}.(TurnKey)(:,1)'>TimeMinWin),2)) = PairUtt_I{1,SpeCond}.(TurnKey)(PairUtt_I{1,SpeCond}.(TurnKey)(:,1)'>TimeMinWin);
        elseif contains(cell2mat(FileNames_I(i)),'Noise70')
            SpeCond = SpeB + 7;
            TurnDurN70_I(x_I,1:size(PairUtt_I{1,SpeCond}.(TurnKey)(PairUtt_I{1,SpeCond}.(TurnKey)(:,1)'>TimeMinWin),2)) = PairUtt_I{1,SpeCond}.(TurnKey)(PairUtt_I{1,SpeCond}.(TurnKey)(:,1)'>TimeMinWin);
        end

        SpeakRaw = PairUtt_I{1,SpeCond}.(SpeakKey);
        ListenRaw = PairUtt_I{1,SpeCond}.(ListenKey);
        GapRaw = PairUtt_I{1,SpeCond}.(GapKey);
        OLWRaw = PairUtt_I{1,SpeCond}.(OLWKey);
        OLBRaw = PairUtt_I{1,SpeCond}.(OLBKey);
        TurnRaw = PairUtt_I{1,SpeCond}.(TurnKey);
        PauseRaw = PairUtt_I{1,SpeCond}.(PauseKey);
        binResUtt = PairUtt_I{1,SpeCond}.binRes;

        Speak = SpeakRaw;
        Listen = ListenRaw;
        Gap = GapRaw;
        OLW = OLWRaw;
        OLB = OLBRaw;
        Turn = TurnRaw;
        Pause = PauseRaw;

        % SAME PROCESSING AS IN W1.m -> No delay, No start at 20 s
%         [Speak,Listen] = overlap_windows(Speak,Listen,Param.Fs);

        % MI = {SpeakMI,ListenMI,GapMI,OLWMI,OLBMI,TurnMI,PauseMI};
        MI = {};
        % Merge windows if duration between windows <= TimeInitialMerge (300 ms)
        MI{1} = merge_windows(Speak, Param.Fs, TimeInitialMerge);
        MI{2} = merge_windows(Listen, Param.Fs, TimeInitialMerge);
        MI{3} = merge_windows(Gap, Param.Fs, TimeInitialMerge);
        MI{4} = merge_windows(OLW, Param.Fs, TimeInitialMerge);
        MI{5} = merge_windows(OLB, Param.Fs, TimeInitialMerge);
        MI{6} = merge_windows(Turn, Param.Fs, TimeInitialMerge);
        MI{7} = merge_windows(Pause, Param.Fs, TimeInitialMerge);
        
        % Check empty variables, assign nan
        if any(cellfun('isempty',MI))
            [MI{cellfun('isempty',MI)}]=deal(nan);
        end
        
        % Discard windows if duration is < TimeMinWin (500 ms)
        DI{1} = MI{1}(MI{1}(:,1)>TimeMinWin,:);
        DI{2} = MI{2}(MI{2}(:,1)>TimeMinWin,:);
        DI{3} = MI{3}(MI{3}(:,1)>TimeMinWin,:);
        DI{4} = MI{4}(MI{4}(:,1)>TimeMinWin,:);
        DI{5} = MI{5}(MI{5}(:,1)>TimeMinWin,:);
        DI{6} = MI{6}(MI{6}(:,1)>TimeMinWin,:);
        DI{7} = MI{7}(MI{7}(:,1)>TimeMinWin,:);
        
        % Check empty variables, assign nan
        if any(cellfun('isempty',DI))
            [DI{cellfun('isempty',DI)}]=deal(nan);
        end
        
        % Merge again if duration between windows <= TimeMerge (2 s)
        M{1} = merge_windows(DI{1}, Param.Fs, TimeMerge);
        M{2} = merge_windows(DI{2}, Param.Fs, TimeMerge);
        M{3} = merge_windows(DI{3}, Param.Fs, TimeMerge);
        M{4} = merge_windows(DI{4}, Param.Fs, TimeMerge);
        M{5} = merge_windows(DI{5}, Param.Fs, TimeMerge);
        M{6} = merge_windows(DI{6}, Param.Fs, TimeMerge);
        M{7} = merge_windows(DI{7}, Param.Fs, TimeMerge);
        
        % Check empty variables, assign nan
        if any(cellfun('isempty',M))
            [M{cellfun('isempty',M)}]=deal(nan);
        end
        
        % Discard windows if duration is < 2*TimeMinWin (1 s)
        D{1} = M{1}(M{1}(:,1)>2*TimeMinWin,:);
        D{2} = M{2}(M{2}(:,1)>2*TimeMinWin,:);
        D{3} = M{3}(M{3}(:,1)>2*TimeMinWin,:);
        D{4} = M{4}(M{4}(:,1)>2*TimeMinWin,:);
        D{5} = M{5}(M{5}(:,1)>2*TimeMinWin,:);
        D{6} = M{6}(M{6}(:,1)>2*TimeMinWin,:);
        D{7} = M{7}(M{7}(:,1)>2*TimeMinWin,:);
        % Check empty variables, assign nan
        if any(cellfun('isempty',D))
            [D{cellfun('isempty',D)}]=deal(nan);
        end
        
        % Overlap Speaking/Listening
        [OM{1},OM{2}] = overlap_windows(D{1},D{2},Param.Fs);

        GapDurRaw_I(x_I,1:size(GapRaw,1))=GapRaw(:,1)';
        OLWDurRaw_I(x_I,1:size(OLWRaw,1))=OLWRaw(:,1)';
        OLBDurRaw_I(x_I,1:size(OLBRaw,1))=OLBRaw(:,1)';
        TurnDurRaw_I(x_I,1:size(TurnRaw,1))=TurnRaw(:,1)';
        PauseDurRaw_I(x_I,1:size(PauseRaw,1))=PauseRaw(:,1)';
        SpeakDurRaw_I(x_I,1:size(SpeakRaw,1))=SpeakRaw(:,1)';
        ListenDurRaw_I(x_I,1:size(ListenRaw,1))=ListenRaw(:,1)';
        
        GapDur_I(x_I,1:size(D{3},1))=D{3}(:,1)';
        OLWDur_I(x_I,1:size(D{4},1))=D{4}(:,1)';
        OLBDur_I(x_I,1:size(D{5},1))=D{5}(:,1)';
        TurnDur_I(x_I,1:size(D{6},1))=D{6}(:,1)';
        PauseDur_I(x_I,1:size(D{7},1))=D{7}(:,1)';
        SpeakDur_I(x_I,1:size(OM{1},1))=OM{1}(:,1)';
        ListenDur_I(x_I,1:size(OM{2},1))=OM{2}(:,1)';
        
        % Find time diff between turns (same talker)
        for j=1:size(Turn,1)-1
            TurnDiff_I(x_I,j)=(Turn(j+1,2)-Turn(j,3))*binResUtt;
        end
        
        % Median
        TurnMed=[median(nonzeros(TurnDurQuiet_I));median(nonzeros(TurnDurSHL_I));median(nonzeros(TurnDurN60_I));median(nonzeros(TurnDurN70_I))];
        
        x_I=x_I+1;

        %Plots
%         figure;
%         startStop = SpeakRaw(:,2:3)*binResUtt;width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),0,width(i),1],'EdgeColor', 'none', 'FaceColor', [SpeakColor 0.3]), 1:size(startStop,1));
%         startStop = ListenRaw(:,2:3)*binResUtt;width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),0,width(i),1],'EdgeColor', 'none', 'FaceColor', [ListenColor 0.3]), 1:size(startStop,1));
%         startStop = MI{1}(:,2:3)*binResUtt;width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),1,width(i),1],'EdgeColor', 'none', 'FaceColor', [SpeakColor 0.3]), 1:size(startStop,1));
%         startStop = MI{2}(:,2:3)*binResUtt;width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),1,width(i),1],'EdgeColor', 'none', 'FaceColor', [ListenColor 0.3]), 1:size(startStop,1));
%         startStop = DI{1}(:,2:3)*binResUtt;width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),2,width(i),1],'EdgeColor', 'none', 'FaceColor', [SpeakColor 0.3]), 1:size(startStop,1));
%         startStop = DI{2}(:,2:3)*binResUtt;width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),2,width(i),1],'EdgeColor', 'none', 'FaceColor', [ListenColor 0.3]), 1:size(startStop,1));
%         startStop = M{1}(:,2:3)*binResUtt;width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),3,width(i),1],'EdgeColor', 'none', 'FaceColor', [SpeakColor 0.3]), 1:size(startStop,1));
%         startStop = M{2}(:,2:3)*binResUtt;width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),3,width(i),1],'EdgeColor', 'none', 'FaceColor', [ListenColor 0.3]), 1:size(startStop,1));
%         startStop = D{1}(:,2:3)*binResUtt;width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),4,width(i),1],'EdgeColor', 'none', 'FaceColor', [SpeakColor 0.3]), 1:size(startStop,1));
%         startStop = D{2}(:,2:3)*binResUtt;width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),4,width(i),1],'EdgeColor', 'none', 'FaceColor', [ListenColor 0.3]), 1:size(startStop,1));
%         startStop = OM{1}(:,2:3)*binResUtt;width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),5,width(i),1],'EdgeColor', 'none', 'FaceColor', [SpeakColor 0.3]), 1:size(startStop,1));
%         startStop = OM{2}(:,2:3)*binResUtt;width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),5,width(i),1],'EdgeColor', 'none', 'FaceColor', [ListenColor 0.3]), 1:size(startStop,1));
%         startStop = GapRaw(:,2:3)*binResUtt;width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),6,width(i),1],'EdgeColor', 'none', 'FaceColor', [GapColor 0.3]), 1:size(startStop,1))
%         startStop = OLWRaw(:,2:3)*binResUtt;width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),6,width(i),1],'EdgeColor', 'none', 'FaceColor', [OLWColor 0.3]), 1:size(startStop,1))
%         startStop = OLBRaw(:,2:3)*binResUtt;width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),6,width(i),1],'EdgeColor', 'none', 'FaceColor', [OLBColor 0.3]), 1:size(startStop,1))
%         startStop = TurnRaw(:,2:3)*binResUtt;width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),6,width(i),1],'EdgeColor', 'none', 'FaceColor', [TurnColor 0.3]), 1:size(startStop,1))
%         startStop = PauseRaw(:,2:3)*binResUtt;width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),6,width(i),1],'EdgeColor', 'none', 'FaceColor', [PauseColor 0.3]), 1:size(startStop,1))
%         yline(1,"--",'HandleVisibility','off')
%         yline(2,"--",'HandleVisibility','off')
%         yline(3,"--",'HandleVisibility','off')
%         yline(4,"--",'HandleVisibility','off')
%         yline(5,"--",'HandleVisibility','off')
%         yline(6,"--",'HandleVisibility','off')
%         yline(7,"--",'HandleVisibility','off')
%         xlabel('Time [s]')
%         yticks(0.5:6.5)
%         yticklabels({'Raw','Merge < 0.3','Discard < 0.5','Merge < 2','Discard < 1','Split overlaps','Other windows'})
%         title(strrep(strrep([PairFiles_I(1).folder(25:end),'\',cell2mat(FileNames_I(i))],'_','-'),'\','\\'))
    end
end

% AMEND II
for q=1:numel(subDirs_II)
    PairIn_II = q;
    for ChosenFolder = {'\NH\','\HI\'}
        PairFolder_II=[pwd,'\data\AMEND_II\',cell2mat(subDirs_II(q)),cell2mat(ChosenFolder)]; % Folder naming changed
        PairFiles_II=dir(PairFolder_II); % Folder naming changed
        PairUtt_II=LoadUtt_II.Utterances(PairIn_II,:);
        
        for i=1:numel(FileNames_II)
            if contains(ChosenFolder,'HI')
                if LoadTPsOrder_II.TPsOrder_II(2*q,i) ~= x_II
                    disp(['Warning: File ',PairFiles_II(1).folder, '\', cell2mat(FileNames_II(i)),' was rejected in AMEND II (P2) analysis.']);
                    continue
                end
            elseif contains(ChosenFolder,'NH')
                if LoadTPsOrder_II.TPsOrder_II(2*q-1,i) ~= x_II
                    disp(['Warning: File ',PairFiles_II(1).folder, '\', cell2mat(FileNames_II(i)),' was rejected in AMEND II (P1) analysis.']);
                    continue
                end
            end
            % Retrieve Utterances
            if contains(ChosenFolder,'HI')
                SpeakKey = 'utteranceCH1';
                ListenKey = 'utteranceCH2';
                GapKey = 'gapCH1';
                OLWKey = 'overlapWCH1';
                OLBKey = 'overlapBCH1';
                TurnKey = 'turnCH1';
                PauseKey = 'pauseCH1';
            elseif contains(ChosenFolder,'NH')
                SpeakKey = 'utteranceCH2';
                ListenKey = 'utteranceCH1';
                GapKey = 'gapCH2';
                OLWKey = 'overlapWCH2';
                OLBKey = 'overlapBCH2';
                TurnKey = 'turnCH2';
                PauseKey = 'pauseCH2';
            end

            if contains(cell2mat(FileNames_II(i)),'UN')
                SpeAid = 0;
            elseif contains(cell2mat(FileNames_II(i)),'AA')
                SpeAid = 3;
            elseif contains(cell2mat(FileNames_II(i)),'AB')
                SpeAid = 6;
            end

            if contains(cell2mat(FileNames_II(i)),'N0')
                SpeCond = SpeAid + 1;
                TurnDurQuiet_II(x_II,1:size(PairUtt_II{1,SpeCond}.(TurnKey)(PairUtt_II{1,SpeCond}.(TurnKey)(:,1)'>TimeMinWin),2)) = PairUtt_II{1,SpeCond}.(TurnKey)(PairUtt_II{1,SpeCond}.(TurnKey)(:,1)'>TimeMinWin);
            elseif contains(cell2mat(FileNames_II(i)),'N60')
                SpeCond = SpeAid + 2;
                TurnDurN60_II(x_II,1:size(PairUtt_II{1,SpeCond}.(TurnKey)(PairUtt_II{1,SpeCond}.(TurnKey)(:,1)'>TimeMinWin),2)) = PairUtt_II{1,SpeCond}.(TurnKey)(PairUtt_II{1,SpeCond}.(TurnKey)(:,1)'>TimeMinWin);
            elseif contains(cell2mat(FileNames_II(i)),'N70')
                SpeCond = SpeAid + 3;
                TurnDurN70_II(x_II,1:size(PairUtt_II{1,SpeCond}.(TurnKey)(PairUtt_II{1,SpeCond}.(TurnKey)(:,1)'>TimeMinWin),2)) = PairUtt_II{1,SpeCond}.(TurnKey)(PairUtt_II{1,SpeCond}.(TurnKey)(:,1)'>TimeMinWin);
            end

            SpeakRaw = PairUtt_II{1,SpeCond}.(SpeakKey);
            ListenRaw = PairUtt_II{1,SpeCond}.(ListenKey);
            GapRaw = PairUtt_II{1,SpeCond}.(GapKey);
            OLWRaw = PairUtt_II{1,SpeCond}.(OLWKey);
            OLBRaw = PairUtt_II{1,SpeCond}.(OLBKey);
            TurnRaw = PairUtt_II{1,SpeCond}.(TurnKey);
            PauseRaw = PairUtt_II{1,SpeCond}.(PauseKey);
            binResUtt = PairUtt_II{1,SpeCond}.binRes;

            Speak = SpeakRaw;
            Listen = ListenRaw;
            Gap = GapRaw;
            OLW = OLWRaw;
            OLB = OLBRaw;
            Turn = TurnRaw;
            Pause = PauseRaw;

            % SAME PROCESSING AS IN W1.m -> No delay, No start at 20 s
            % MI = {SpeakMI,ListenMI,GapMI,OLWMI,OLBMI,TurnMI,PauseMI};
            MI = {};
            % Merge windows if duration between windows <= TimeInitialMerge (300 ms)
            MI{1} = merge_windows(Speak, Param.Fs, TimeInitialMerge);
            MI{2} = merge_windows(Listen, Param.Fs, TimeInitialMerge);
            MI{3} = merge_windows(Gap, Param.Fs, TimeInitialMerge);
            MI{4} = merge_windows(OLW, Param.Fs, TimeInitialMerge);
            MI{5} = merge_windows(OLB, Param.Fs, TimeInitialMerge);
            MI{6} = merge_windows(Turn, Param.Fs, TimeInitialMerge);
            MI{7} = merge_windows(Pause, Param.Fs, TimeInitialMerge);
            % Check empty variables, assign nan
            if any(cellfun('isempty',MI))
                [MI{cellfun('isempty',MI)}]=deal(nan);
            end

            % Discard windows if duration is < TimeMinWin (500 ms)
            DI{1} = MI{1}(MI{1}(:,1)>TimeMinWin,:);
            DI{2} = MI{2}(MI{2}(:,1)>TimeMinWin,:);
            DI{3} = MI{3}(MI{3}(:,1)>TimeMinWin,:);
            DI{4} = MI{4}(MI{4}(:,1)>TimeMinWin,:);
            DI{5} = MI{5}(MI{5}(:,1)>TimeMinWin,:);
            DI{6} = MI{6}(MI{6}(:,1)>TimeMinWin,:);
            DI{7} = MI{7}(MI{7}(:,1)>TimeMinWin,:);
            % Check empty variables, assign nan
            if any(cellfun('isempty',DI))
                [DI{cellfun('isempty',DI)}]=deal(nan);
            end

            % Merge again if duration between windows <= TimeMerge (2 s)
            M{1} = merge_windows(DI{1}, Param.Fs, TimeMerge);
            M{2} = merge_windows(DI{2}, Param.Fs, TimeMerge);
            M{3} = merge_windows(DI{3}, Param.Fs, TimeMerge);
            M{4} = merge_windows(DI{4}, Param.Fs, TimeMerge);
            M{5} = merge_windows(DI{5}, Param.Fs, TimeMerge);
            M{6} = merge_windows(DI{6}, Param.Fs, TimeMerge);
            M{7} = merge_windows(DI{7}, Param.Fs, TimeMerge);
            % Check empty variables, assign nan
            if any(cellfun('isempty',M))
                [M{cellfun('isempty',M)}]=deal(nan);
            end

            % Discard windows if duration is < 2*TimeMinWin (1 s)
            D{1} = M{1}(M{1}(:,1)>2*TimeMinWin,:);
            D{2} = M{2}(M{2}(:,1)>2*TimeMinWin,:);
            D{3} = M{3}(M{3}(:,1)>2*TimeMinWin,:);
            D{4} = M{4}(M{4}(:,1)>2*TimeMinWin,:);
            D{5} = M{5}(M{5}(:,1)>2*TimeMinWin,:);
            D{6} = M{6}(M{6}(:,1)>2*TimeMinWin,:);
            D{7} = M{7}(M{7}(:,1)>2*TimeMinWin,:);
            % Check empty variables, assign nan
            if any(cellfun('isempty',D))
                [D{cellfun('isempty',D)}]=deal(nan);
            end

            % Overlap Speaking/Listening
            [D{1},D{2}] = overlap_windows(D{1},D{2},Param.Fs);

            GapDurRaw_II(x_II,1:size(GapRaw,1))=GapRaw(:,1)';
            OLWDurRaw_II(x_II,1:size(OLWRaw,1))=OLWRaw(:,1)';
            OLBDurRaw_II(x_II,1:size(OLBRaw,1))=OLBRaw(:,1)';
            TurnDurRaw_II(x_II,1:size(TurnRaw,1))=TurnRaw(:,1)';
            PauseDurRaw_II(x_II,1:size(PauseRaw,1))=PauseRaw(:,1)';
            SpeakDurRaw_II(x_II,1:size(SpeakRaw,1))=SpeakRaw(:,1)';
            ListenDurRaw_II(x_II,1:size(ListenRaw,1))=ListenRaw(:,1)';

            GapDur_II(x_II,1:size(D{3},1))=D{3}(:,1)';
            OLWDur_II(x_II,1:size(D{4},1))=D{4}(:,1)';
            OLBDur_II(x_II,1:size(D{5},1))=D{5}(:,1)';
            TurnDur_II(x_II,1:size(D{6},1))=D{6}(:,1)';
            PauseDur_II(x_II,1:size(D{7},1))=D{7}(:,1)';
            SpeakDur_II(x_II,1:size(D{1},1))=D{1}(:,1)';
            ListenDur_II(x_II,1:size(D{2},1))=D{2}(:,1)';

            x_II=x_II+1;
        end
    end
end
%% Out-of-loop: Plots
%             GapDur1(~any(GapDur1,2),:)=[];
%             OLWDur1(~any(OLWDur1,2),:)=[];
%             OLBDur1(~any(OLBDur1,2),:)=[];
%             TurnDur1(~any(TurnDur1,2),:)=[];
%             PauseDur1(~any(PauseDur1,2),:)=[];

GroupRaw_I = struct('SpeakDur',SpeakDurRaw_I,'ListenDur',ListenDurRaw_I,'GapDur',GapDurRaw_I,'PauseDur',PauseDurRaw_I,'OLWDur',OLWDurRaw_I,'OLBDur',OLBDurRaw_I,'TurnDur',TurnDurRaw_I);
GroupPro_I = struct('SpeakDur',SpeakDur_I,'ListenDur',ListenDur_I,'GapDur',GapDur_I,'PauseDur',PauseDur_I,'OLWDur',OLWDur_I,'OLBDur',OLBDur_I,'TurnDur',TurnDur_I);
GroupRaw_II = struct('SpeakDur',SpeakDurRaw_II,'ListenDur',ListenDurRaw_II,'GapDur',GapDurRaw_II,'PauseDur',PauseDurRaw_II,'OLWDur',OLWDurRaw_II,'OLBDur',OLBDurRaw_II,'TurnDur',TurnDurRaw_II);
GroupPro_II = struct('SpeakDur',SpeakDur_II,'ListenDur',ListenDur_II,'GapDur',GapDur_II,'PauseDur',PauseDur_II,'OLWDur',OLWDur_II,'OLBDur',OLBDur_II,'TurnDur',TurnDur_II);
GroupColor = struct('SpeakColor',SpeakColor,'ListenColor',ListenColor,'GapColor',GapColor,'PauseColor',PauseColor,'OLWColor',OLWColor,'OLBColor',OLBColor,'TurnColor',TurnColor);

histplot(GroupRaw_I,GroupColor,nbins);
% sgtitle('Raw')
histplot(GroupPro_I,GroupColor,nbins);
% sgtitle('Processed')

histplot(GroupRaw_II,GroupColor,nbins);
% sgtitle('Raw')
histplot(GroupPro_II,GroupColor,nbins);
% sgtitle('Processed')

flgd=figure;axlgd=gca();hold(axlgd,'on')
plot(axlgd,nan,color=[SpeakColor 0.3],linewidth=10) % plot nans to show color in legend
plot(axlgd,nan,color=[ListenColor 0.3],linewidth=10) % plot nans to show color in legend
% plot(axlgd,nan,color=[GapColor 0.3],linewidth=10) % plot nans to show color in legend
% plot(axlgd,nan,color=[PauseColor 0.3],linewidth=10) % plot nans to show color in legend
% plot(axlgd,nan,color=[OLWColor 0.3],linewidth=10) % plot nans to show color in legend
% plot(axlgd,nan,color=[OLBColor 0.3],linewidth=10) % plot nans to show color in legend
% plot(axlgd,nan,color=[TurnColor 0.3],linewidth=10) % plot nans to show color in legend
lgd=legend(axlgd,'Speak','Listen','Gap','Pause','O.Within','O.Between','Turn','Location','southeastoutside');
lgd.Title.String = 'Types of windows:';

saveLegendToImage(flgd, lgd, 'lgd', 'eps');
exportgraphics(flgd, ['output_' num2str(1) '.jpg']);

% figure;histogram(nonzeros(TurnDiff_I),nbins,'FaceColor',TurnColor(1:3));title(['Time between turns - Avg: ',sprintf('%0.5f',mean(nonzeros(TurnDiff_I))),' s']);
% xlabel('Time [s]'),ylabel('Bin Count');grid on;
% 
% figure;tiledlayout(1,4);ax1=nexttile;ax2=nexttile;ax3=nexttile;ax4=nexttile;histogram(ax1,nonzeros(TurnDurQuiet_I),nbins,'FaceColor',TurnColor(1:3)),histogram(ax2,nonzeros(TurnDurSHL_I),nbins,'FaceColor',TurnColor(1:3));
% histogram(ax3,nonzeros(TurnDurN60_I),nbins,'FaceColor',TurnColor(1:3));histogram(ax4,nonzeros(TurnDurN70_I),nbins,'FaceColor',TurnColor(1:3))
% title(ax1,['Quiet turn duration - Avg: ',sprintf('%0.5f',mean(nonzeros(TurnDurQuiet_I),'omitnan')),' s']);title(ax2,['SHL turn duration - Avg: ',sprintf('%0.2f',mean(nonzeros(TurnDurSHL_I),'omitnan')),' s']);
% title(ax3,['Noise60 turn duration - Avg: ',sprintf('%0.2f',mean(nonzeros(TurnDurN60_I),'omitnan')),' s']);title(ax4,['Noise70 turn duration - Avg: ',sprintf('%0.2f',mean(nonzeros(TurnDurN70_I),'omitnan')),' s']);
% xlabel([ax1 ax2 ax3 ax4],'Time [s]');ylabel([ax1 ax2 ax3 ax4],'Bin Count');grid([ax1 ax2 ax3 ax4],'on');

% figure;h11=boxplotGroup({nonzeros(SpeakDurRaw_I(~isnan(SpeakDurRaw_I))),nonzeros(ListenDurRaw_I(~isnan(ListenDurRaw_I))),nonzeros(GapDurRaw_I(~isnan(GapDurRaw_I))),nonzeros(PauseDurRaw_I(~isnan(PauseDurRaw_I))),nonzeros(OLWDurRaw_I(~isnan(OLWDurRaw_I))),nonzeros(OLBDurRaw_I(~isnan(OLBDurRaw_I))),nonzeros(TurnDurRaw_I(~isnan(TurnDurRaw_I)))},...
% 'PrimaryLabels', {'Speak','Listen','Gap','Pause','OLW','OLB','Turn'}, ...
% 'SecondaryLabels', {'Durations'}, ...
% 'interGroupSpace',2,'groupLabelType','Vertical', ...
% 'PlotStyle','Compact','BoxStyle','filled',...
% 'Colors',[SpeakColor;ListenColor;GapColor;PauseColor;OLWColor;OLBColor;TurnColor],'GroupType','betweenGroups');


function histplot(G,C,nbins)
    figure;tl=tiledlayout(2,4);
    ax1=nexttile;ax2=nexttile;ax3=nexttile;ax4=nexttile;ax5=nexttile;ax6=nexttile;ax7=nexttile;
    h1=histogram(ax1,nonzeros(G.SpeakDur(~isnan(G.SpeakDur))),nbins,'Normalization','PDF','FaceColor',C.SpeakColor,'FaceAlpha',0.3,'EdgeColor', 'none');
    h2=histogram(ax2,nonzeros(G.ListenDur(~isnan(G.ListenDur))),nbins,'Normalization','PDF','FaceColor',C.ListenColor,'FaceAlpha',0.3,'EdgeColor', 'none');
    h3=histogram(ax3,nonzeros(G.GapDur(~isnan(G.GapDur))),nbins,'Normalization','PDF','FaceColor',C.GapColor,'FaceAlpha',0.3,'EdgeColor', 'none');
    h4=histogram(ax4,nonzeros(G.PauseDur(~isnan(G.PauseDur))),nbins,'Normalization','PDF','FaceColor',C.PauseColor,'FaceAlpha',0.3,'EdgeColor', 'none');
    h5=histogram(ax5,nonzeros(G.OLWDur(~isnan(G.OLWDur))),nbins,'Normalization','PDF','FaceColor',C.OLWColor,'FaceAlpha',0.3,'EdgeColor', 'none');
    h6=histogram(ax6,nonzeros(G.OLBDur(~isnan(G.OLBDur))),nbins,'Normalization','PDF','FaceColor',C.OLBColor,'FaceAlpha',0.3,'EdgeColor', 'none');
    h7=histogram(ax7,nonzeros(G.TurnDur(~isnan(G.TurnDur))),10*nbins,'Normalization','PDF','FaceColor',C.TurnColor,'FaceAlpha',0.3,'EdgeColor', 'none');
%     title(ax1,['Speak duration - Avg: ',sprintf('%0.2f',mean(nonzeros(G.SpeakDur),'omitnan')),' s']);
%     title(ax2,['Listen duration - Avg: ',sprintf('%0.2f',mean(nonzeros(G.ListenDur),'omitnan')),' s']);
%     title(ax3,['Gap duration - Avg: ',sprintf('%0.2f',mean(nonzeros(G.GapDur),'omitnan')),' s']);
%     title(ax4,['Pause duration - Avg: ',sprintf('%0.2f',mean(nonzeros(G.PauseDur),'omitnan')),' s']);
%     title(ax5,['Overlap-Within duration - Avg: ',sprintf('%0.2f',mean(nonzeros(G.OLWDur),'omitnan')),' s']);
%     title(ax6,['Overlap-Between duration - Avg: ',sprintf('%0.2f',mean(nonzeros(G.OLBDur),'omitnan')),' s']);
%     title(ax7,['Turn duration - Avg: ',sprintf('%0.2f',mean(nonzeros(G.TurnDur),'omitnan')),' s']);
    xlabel([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'Time [s]')
    ylabel([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'Bin Count')
    xlim([ax1,ax2,ax3,ax4,ax5,ax6,ax7],[0 4])
    xticks([ax1,ax2,ax3,ax4,ax5,ax6,ax7],0:2:4)
    set([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'ytick',[])
    grid([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'on')
    grid([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'minor')
    lgd=legend([h1,h2,h3,h4,h5,h6,h7],'Speak','Listen','Gap','Pause','O.Within','O.Between','Turn');
    lgd.Layout.Tile = 8;
%     exportgraphics(tl, ['output_' num2str(round(100*rand(1))) '.pdf'])
end
