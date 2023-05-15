%% PBCA-Thesis - Week 1, 2, 3 - Processing existing pupil data
% Pathing
clear all; clc; close all;
BPath = strsplit(pwd,'PBCA-thesis');
addpath('tools')
addpath([BPath{1} 'Pupil-preprocessing-tools\tools']) % For preprocessing

% Ask for pair number (pair of test subjects)
% while true
%     PairIn = str2double(input('Enter participant pair (from 1 to 12): ','s'));
%     if fix(PairIn) == PairIn && PairIn >= 1 && PairIn <= 12
%         break;
%     end
%     disp('Number must be an integer between 1 and 12.');
% end

[subDirs] = GetSubDirsFirstLevelOnly('data\AMEND_I');
LoadUtt=load('data\AMEND_I\utterances1110.mat');
LoadDelays=load('data\AMEND_I\delays1110.mat');

% Parameters for processing
Param.Fs = 50; % Sampling frequency of pupil data
Param.RemoveBeforeAndAfter = [35 100]*1e-3; % Samples within the time range before and after NaNs will set NaNs as well.
Param.MinLengthNaNRepair = 5; % Drop values (i.e., change to NaN) before and after NaNs only for contiguous NaNs of at least __ samples. 
LPWinSize = 1; % [s]: Window size of hamming-window for low-pass filtering
LPWindow = hamming(round(LPWinSize*Param.Fs));
LPWindow = LPWindow/sum(LPWindow); % Hamming-window
TimeStartW = 0.5; % [s], time before Utt/Lis starts
TimeEndW = 0; % [s], time after Utt/Lis starts
TimeStart = 20; % [s], time at which simultaneous recording started
TimeBL = [10,15]; % [s], time chosen for baseline
TimeMinWin = 0.5; % [s], Minimum time of a window
TimeInitialMerge = 0.3; % [s], Time threshold for merging windows initially
TimeMerge = 2; % [s], Time threshold for merging windows after rejecting small windows
AdapBL = 0.3; % [s], Duration of baseline prior to event
FilterWidth = round((LPWinSize*Param.Fs)/2); % [samples]: Width of hamming filter used for fixation duration
RejectRatio = 0.4; % Rejection threshold based on the ratio of NaNs in data
RejectDelay = 0.5; % [s], Rejection threshold based on delay between timestamps and n-samples
NFilesMax = 16; % Max number of files

% Initialize variables
AvgLPupilSize = zeros(numel(subDirs),NFilesMax); % [mm], avg Listening pupil size (per file)
AvgLPupilSlope = zeros(numel(subDirs),NFilesMax); % [mm/s], avg Listening baselined pupil slope (per file)
AvgSPupilSize = zeros(numel(subDirs),NFilesMax); % [mm], avg Speaking pupil size (per file)
AvgSPupilSlope = zeros(numel(subDirs),NFilesMax); % [mm/s], avg Speaking baselined pupil slope (per file)
AvgTimeSpe = zeros(numel(subDirs),NFilesMax); % [s], avg time of Utterance (per file)
AvgTimeLis = zeros(numel(subDirs),NFilesMax); % [s], avg time of Listening (per file)

% Run through all subfolders once
for q=1:numel(subDirs)
    PairIn = q;
    
    % Files and Utterances: different conditions
    PairFiles=dir(['data\AMEND_I\Main',sprintf('%d',PairIn),'\*.mat']);
    PairUtt=LoadUtt.Utterances(PairIn,:);
    PairDelay=LoadDelays.TobAudDelay(PairIn,:);
    
    % Assign NaNs - Number of files in folder <= NFilesMax
    if numel(PairFiles) < NFilesMax
        AvgSPupilSize(q,numel(PairFiles):end)=NaN;
        AvgSPupilSlope(q,numel(PairFiles):end)=NaN;
        AvgLPupilSize(q,numel(PairFiles):end)=NaN;
        AvgLPupilSlope(q,numel(PairFiles):end)=NaN;
        AvgTimeSpe(q,numel(PairFiles):end)=NaN;
        AvgTimeLis(q,numel(PairFiles):end)=NaN;
    end
    
    % Run once per file
    for i=1:numel(PairFiles)
        alldata = load([PairFiles(i).folder, '\', PairFiles(i).name]);
        alldata_mat = cell2mat(alldata.data);
        
        % Replace blanks '[]' for 'NaN' in fields diameterLeft and diameterRight
        [alldata_mat(cellfun(@isempty,{alldata_mat.diameterLeft})).diameterLeft] = deal(NaN);
        [alldata_mat(cellfun(@isempty,{alldata_mat.diameterRight})).diameterRight] = deal(NaN);        
        
        % Extract delay [s] (Duration_timeStamps - Duration_nsamples)
        EyeAudDelay=alldata_mat(end).timeStamp-alldata_mat(1).timeStamp-length(alldata_mat)/Param.Fs;
        
        % Skip file if the difference in duration from the number of
        % samples and the duration given by timestamps is bigger than 0.5 s
        if EyeAudDelay > RejectDelay
            disp(['Warning: File ',PairFiles(i).folder, '\', PairFiles(i).name, ' was rejected, too much delay (',sprintf('%0.2f',EyeAudDelay),'s).']);
            AvgSPupilSize(q,i)=NaN;
            AvgSPupilSlope(q,i)=NaN;
            AvgLPupilSize(q,i)=NaN;
            AvgLPupilSlope(q,i)=NaN;
            AvgTimeSpe(q,i)=NaN;
            AvgTimeLis(q,i)=NaN;
            continue
        end
        
        LDiamRaw = [alldata_mat.diameterLeft];
        RDiamRaw = [alldata_mat.diameterRight];
        
         % Preprocessing - Setting outliers as NaNs (remove artifacts)
%             LThreshOut = [mean(LDiamRaw,'omitnan')-2*std(LDiamRaw,'omitnan'),mean(LDiamRaw,'omitnan')+2*std(LDiamRaw,'omitnan')];
%             RThreshOut = [mean(RDiamRaw,'omitnan')-2*std(RDiamRaw,'omitnan'),mean(RDiamRaw,'omitnan')+2*std(RDiamRaw,'omitnan')];
%             for s=1:length(alldata_mat)
%                 if LDiamRaw(1,s) < LThreshOut(1) || LDiamRaw(1,s) > LThreshOut(2)
%                     LDiamRaw(1,s)=NaN;
%                 end
%                 if RDiamRaw(1,s) < RThreshOut(1) || RDiamRaw(1,s) > RThreshOut(2)
%                     RDiamRaw(1,s)=NaN;
%                 end
%             end

        % New artifact-removal method
        standardRawSettings = rawDataFilter();
        [LvalOut,LspeedFiltData,LdevFiltData] = rawDataFilter(linspace(0,length(LDiamRaw)./Param.Fs,length(LDiamRaw))',LDiamRaw',standardRawSettings);
        [RvalOut,RspeedFiltData,RdevFiltData] = rawDataFilter(linspace(0,length(RDiamRaw)./Param.Fs,length(RDiamRaw))',RDiamRaw',standardRawSettings);

        LDiamRaw(~LvalOut)=NaN;
        RDiamRaw(~RvalOut)=NaN;
        
        % Processing - Interpolating NaNs
        [LDiam,LMetadata] = preprocpupil(LDiamRaw,Param);
        [RDiam,RMetadata] = preprocpupil(RDiamRaw,Param);

        % Low-Pass Filtering
        LDiamConv = conv(LDiam,LPWindow,'same'); 
        RDiamConv = conv(RDiam,LPWindow,'same');
        
        % Remove start/end artifacts (peak/dip) originated from Low-Pass Filtering
        LDiamConv(1:round(length(LPWindow)/2-1)) = mean(LDiam(1:round(length(LPWindow)/2-1)));
        LDiamConv(end-round(length(LPWindow)/2-1):end) = mean(LDiam(end-round(length(LPWindow)/2-1):end));
        RDiamConv(1:round(length(LPWindow)/2-1)) = mean(RDiam(1:round(length(LPWindow)/2-1)));
        RDiamConv(end-round(length(LPWindow)/2-1):end) = mean(RDiam(end-round(length(LPWindow)/2-1):end));
        
        LDiam = LDiamConv;
        RDiam = RDiamConv;
        
        % Decide 'better' eye results
        [Min idx_decision] = min([sum(LMetadata.Isnan) sum(RMetadata.Isnan)]);
        if idx_decision == 1
            Diameter = LDiam;
            DiameterRaw = LDiamRaw;
            eyeChosen = 'Left';
            DiamNaN = sum(LMetadata.Isnan);
        elseif idx_decision == 2
            Diameter = RDiam;
            DiameterRaw = RDiamRaw;
            eyeChosen = 'Right';
            DiamNaN = sum(RMetadata.Isnan);
        end
        
        % Reject a file if the chosen eye data has too many NaNs
        if DiamNaN/length(Diameter) >= RejectRatio
            disp(['Warning: File ',PairFiles(i).folder, '\', PairFiles(i).name, ' was rejected because it contains too many NaNs (',sprintf('%0.2f',100*DiamNaN/length(Diameter)),'%).'])
            AvgSPupilSize(q,i)=NaN;
            AvgSPupilSlope(q,i)=NaN;
            AvgLPupilSize(q,i)=NaN;
            AvgLPupilSlope(q,i)=NaN;
            AvgTimeSpe(q,i)=NaN;
            AvgTimeLis(q,i)=NaN;
            continue
        end
        % Baseline the Diameter (substractively)
        BLDiam = Diameter-mean(Diameter(TimeBL(1)*Param.Fs:TimeBL(2)*Param.Fs-1));

        % Retrieve Utterances
        if contains(PairFiles(i).name,'P2')
            SpeakKey = 'utteranceCH1';
            ListenKey = 'utteranceCH2';
            SDelayKey = 'delayCH1';
            LDelayKey = 'delayCH2';
        elseif contains(PairFiles(i).name,'P1')
            SpeakKey = 'utteranceCH2';
            ListenKey = 'utteranceCH1';
            SDelayKey = 'delayCH2';
            LDelayKey = 'delayCH1';
        end

        if contains(PairFiles(i).name,'B1')
            SpeB = 0;
        elseif contains(PairFiles(i).name,'B2')
            SpeB = 1;
        end

        if contains(PairFiles(i).name,'Quiet')
            SpeCond = SpeB + 1;
        elseif contains(PairFiles(i).name,'SHL')
            SpeCond = SpeB + 3;
        elseif contains(PairFiles(i).name,'Noise60')
            SpeCond = SpeB + 5;
        elseif contains(PairFiles(i).name,'Noise70')
            SpeCond = SpeB + 7;
        end

        SpeakRaw = PairUtt{1,SpeCond}.(SpeakKey);
        ListenRaw = PairUtt{1,SpeCond}.(ListenKey);
        binResUtt = PairUtt{1,SpeCond}.binRes;
        try
            SDelayRaw = PairDelay{1,SpeCond}.(SDelayKey);
            LDelayRaw = PairDelay{1,SpeCond}.(LDelayKey);
        catch ME
            disp(['Warning: No associated Delay data for file ', PairFiles(i).folder, '\', PairFiles(i).name, '.']);
            AvgSPupilSize(q,i)=NaN;
            AvgSPupilSlope(q,i)=NaN;
            AvgLPupilSize(q,i)=NaN;
            AvgLPupilSlope(q,i)=NaN;
            AvgTimeSpe(q,i)=NaN;
            AvgTimeLis(q,i)=NaN;
            continue
        end
        
        if or(SDelayRaw < 0,LDelayRaw < 0)
            SDelayRaw=[0,0];
            LDelayRaw=[0,0];
        end
        
        binResDel = PairDelay{1,SpeCond}.binRes;

        if isempty(SpeakRaw) && isempty(ListenRaw)
            disp(['Warning: No associated Utterance/Listening data for file ', PairFiles(i).folder, '\', PairFiles(i).name, '.']);
            AvgSPupilSize(q,i)=NaN;
            AvgSPupilSlope(q,i)=NaN;
            AvgLPupilSize(q,i)=NaN;
            AvgLPupilSlope(q,i)=NaN;
            AvgTimeSpe(q,i)=NaN;
            AvgTimeLis(q,i)=NaN;
            continue
        end

        % Downsample (rounding) Utt from 250 Hz (1/binRes) to 50 Hz, shift
        % in time to account for the time at which the audio recording
        % started (from 0 to 20 s only eye data) plus delay
        SpeakRaw(:,2:3)=round((SpeakRaw(:,2:3)*binResUtt+TimeStart)*Param.Fs+SDelayRaw(1)/2);
        ListenRaw(:,2:3)=round((ListenRaw(:,2:3)*binResUtt+TimeStart)*Param.Fs+LDelayRaw(1)/2);

        % Merge windows if duration between windows <= TimeInitialMerge (300 ms)
        SpeakMI = merge_windows(SpeakRaw, Param.Fs, TimeInitialMerge);
        ListenMI = merge_windows(ListenRaw, Param.Fs, TimeInitialMerge);
        
        % Discard windows if duration is < TimeMinWin (500 ms)
        SpeakD = SpeakMI(SpeakMI(:,1)>TimeMinWin,:);
        ListenD = ListenMI(ListenMI(:,1)>TimeMinWin,:);
        
        % Merge again if duration between windows <= TimeMerge (2 s)
        SpeakM = merge_windows(SpeakD, Param.Fs, TimeMerge);
        ListenM = merge_windows(ListenD, Param.Fs, TimeMerge);
        
        % Discard windows if duration is < 2*TimeMinWin (1 s)
        Speak = SpeakM(SpeakM(:,1)>2*TimeMinWin,:);
        Listen = ListenM(ListenM(:,1)>2*TimeMinWin,:);
        
        % figure;t_Diam = linspace(1,length(Diameter)./Param.Fs,length(Diameter));startStopS = t_Diam(SpeakRaw(:,2:3));yl=ylim();widthS = startStopS(:,2)-startStopS(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStopS(i,1),yl(1),widthS(i),1],'EdgeColor', 'none', 'FaceColor', [1 0 0 .2]), 1:size(startStopS,1));startStopS2 = t_Diam(SpeakMI(:,2:3));widthS = startStopS2(:,2)-startStopS2(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStopS2(i,1),1,widthS(i),1],'EdgeColor', 'none', 'FaceColor', [0 0 1 .2]), 1:size(startStopS2,1));startStopS3 = t_Diam(SpeakD(:,2:3));widthS = startStopS3(:,2)-startStopS3(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStopS3(i,1),2,widthS(i),1],'EdgeColor', 'none', 'FaceColor', [0 1 1 .2]), 1:size(startStopS3,1));startStopS4 = t_Diam(SpeakM(:,2:3));widthS = startStopS4(:,2)-startStopS4(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStopS5(i,1),3,widthS(i),1],'EdgeColor', 'none', 'FaceColor', [1 0 1 .2]), 1:size(startStopS4,1));startStopS5 = t_Diam(Speak(:,2:3));widthS = startStopS5(:,2)-startStopS5(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStopS5(i,1),4,widthS(i),1],'EdgeColor', 'none', 'FaceColor', [1 1 0 .2]), 1:size(startStopS5,1));grid on        
        
        % Time-locked indexes (based on Start or End of events)
        SWSpeakIdx=[Speak(:,2)-TimeStartW*Param.Fs,Speak(:,2),Speak(:,2)+TimeEndW*Param.Fs];
        SWListenIdx=[Listen(:,2)-TimeStartW*Param.Fs,Listen(:,2),Listen(:,2)+TimeEndW*Param.Fs];
        EWSpeakIdx=[Speak(:,3)-TimeStartW*Param.Fs,Speak(:,3),Speak(:,3)+TimeEndW*Param.Fs];
        EWListenIdx=[Listen(:,3)-TimeStartW*Param.Fs,Listen(:,3),Listen(:,3)+TimeEndW*Param.Fs];

        % Time vectors and idx for plotting
        t_Diam = linspace(0,length(BLDiam)./Param.Fs,length(BLDiam));
        startStopS = t_Diam(Speak(:,2:3)); 
        widthS = startStopS(:,2)-startStopS(:,1);
        startStopL = t_Diam(Listen(:,2:3)); 
        widthL = startStopL(:,2)-startStopL(:,1);

        % Features (each window): mean slope baselined, mean pupil raw
        % Speaking & Listening
        for j=1:size(Speak,1)
            AvgSPupilSize(q,i)=AvgSPupilSize(q,i)+mean(Diameter(Speak(j,2):Speak(j,3)));
            SSlope(j,:)=polyfit(t_Diam(Speak(j,2):Speak(j,3)),BLDiam(Speak(j,2):Speak(j,3)),1);
            AvgSPupilSlope(q,i)=AvgSPupilSlope(q,i)+SSlope(1);
        end
        AvgSPupilSize(q,i)=AvgSPupilSize(q,i)./j;
        AvgSPupilSlope(q,i)=AvgSPupilSlope(q,i)./j;

        for j=1:size(Listen,1)
            AvgLPupilSize(q,i)=AvgLPupilSize(q,i)+mean(Diameter(Listen(j,2):Listen(j,3)));
            LSlope(j,:)=polyfit(t_Diam(Listen(j,2):Listen(j,3)),BLDiam(Listen(j,2):Listen(j,3)),1);
            AvgLPupilSlope(q,i)=AvgLPupilSlope(q,i)+LSlope(j,1);
        end
        AvgLPupilSize(q,i)=AvgLPupilSize(q,i)./j;
        AvgLPupilSlope(q,i)=AvgLPupilSlope(q,i)./j;

        % Average duration of speaking/listening per file
        AvgTimeSpe(q,i) = mean(Speak(:,1)); % Possible NaNs
        AvgTimeLis(q,i)= mean(Listen(:,1)); % Possible NaNs
        
        % Plots
        figure
        subplot(2,2,[1 2])
        plot(t_Diam,Diameter,color='black');
        hold on
        grid on
        xline(TimeStart,"--",'HandleVisibility','off')
        yl=ylim();
        ylim(yl);
        % Plot rectangles (Utterance and listening time windows)
        arrayfun(@(i)rectangle('Position', [startStopS(i,1),yl(1),widthS(i),range(yl)],'EdgeColor', 'none', 'FaceColor', [0 1 0 .2]), 1:size(startStopS,1))
        arrayfun(@(i)rectangle('Position', [startStopL(i,1),yl(1),widthL(i),range(yl)],'EdgeColor', 'none', 'FaceColor', [1 0 1 .2]), 1:size(startStopL,1))
        eline1=line(NaN,NaN,'LineWidth',2,'Color',[0 1 0 .2]);
        eline2=line(NaN,NaN,'LineWidth',2,'Color',[1 0 1 .2]);
        legend(['Baselined diameter (', eyeChosen,' Eye)'],'Speaking windows','Listening windows')
        sgtitle(strrep(PairFiles(i).name,'_','-'))
%         xticks([0:TimeStart:t_Diam(end)])
        xlabel('Time [s]')
        ylabel('Pupil diameter [mm]');

        % Speak-Diam-Window Adaptive Baselined
        subplot(2,2,3)
        hold on
        SDW = zeros(size(Speak,1),ceil(Param.Fs*(TimeStartW+max(Speak(:,1))))+1);
        for j=1:size(Speak,1)
            if SWSpeakIdx(j,3)-1 <= length(Diameter)
                % Add nan-padding when necessary
                SDW(j,:)=[Diameter(SWSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SDW)-length(SWSpeakIdx(j,1):Speak(j,3)))'];
                plot(linspace(-TimeStartW,Speak(j,1),length(SDW)),SDW(j,:),color=[0 1 0 .2],LineWidth=0.3)
            end
        end
        plot(linspace(-TimeStartW,max(Speak(:,1)),length(SDW)),ndnanfilter(mean(SDW,1,'omitnan'),'hamming',FilterWidth),color=[1 0 0 0.8],LineWidth=1.5)
        xline(0,"--")
        grid on
%         xticks([-TimeStartW:TimeStartW:max(Speak(:,1))])
        title('Baselined Speaking Windows')
        xlabel('Time [s]')
        ylabel('Pupil baseline difference [mm]')
        
        % Listen-Diam-Window Adaptive Baselined
        subplot(2,2,4)
        hold on
        LDW = zeros(size(Listen,1),ceil(Param.Fs*(TimeStartW+max(Listen(:,1))))+1);
        for j=1:size(Listen,1)
            if SWListenIdx(j,3)-1 <= length(Diameter)
                % Add nan-padding when necessary
                LDW(j,:)=[Diameter(SWListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LDW)-length(SWListenIdx(j,1):Listen(j,3)))'];
                plot(linspace(-TimeStartW,Listen(j,1),length(LDW)),LDW(j,:),color=[0 1 0 .2],LineWidth=0.3)
            end
        end
        plot(linspace(-TimeStartW,max(Listen(:,1)),length(LDW)),ndnanfilter(mean(LDW,1,'omitnan'),'hamming',FilterWidth),color=[1 0 0 0.8],LineWidth=1.5)
        xline(0,"--")
        grid on
%         xticks([-TimeStartW:TimeStartW:max(Listen(:,1))])
        title('Baselined Listening Windows')
        xlabel('Time [s]')
        ylabel('Pupil baseline difference [mm]')
    end
end

% Global parameters from Utterance windows 
GlobalAvgLPupilSize = mean(mean(AvgLPupilSize,'omitnan')); % Avg Listening Diameter
GlobalAvgLPupilSlope = mean(mean(AvgLPupilSlope,'omitnan')); % Avg Listening Slope
GlobalAvgSPupilSize = mean(mean(AvgSPupilSize,'omitnan')); % Avg Speaking Diameter
GlobalAvgSPupilSlope = mean(mean(AvgSPupilSlope,'omitnan')); % Avg Speaking Slope
GlobalAvgTimeSpe = mean(mean(AvgTimeSpe,'omitnan')); % Avg duration of Utt
GlobalAvgTimeLis = mean(mean(AvgTimeLis,'omitnan')); % Avg duration of Lis
