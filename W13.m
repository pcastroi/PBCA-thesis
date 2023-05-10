%% PBCA-Thesis - Week 13,14,15 - Different types of data processing/normalization methods - Pupillometry
% Pathing
clear all; clc; close all;
BPath = strsplit(pwd,'PBCA-thesis');
addpath('tools')
addpath([BPath{1} 'Pupil-preprocessing-tools\tools']) % For preprocessing

[subDirs] = GetSubDirsFirstLevelOnly('data\AMEND_I');
FileNames={'P1_Quiet_B1.mat','P1_Quiet_B2.mat','P1_SHL_B1.mat','P1_SHL_B2.mat','P1_Noise60_B1.mat','P1_Noise60_B2.mat','P1_Noise70_B1.mat','P1_Noise70_B2.mat','P2_Quiet_B1.mat','P2_Quiet_B2.mat','P2_SHL_B1.mat','P2_SHL_B2.mat','P2_Noise60_B1.mat','P2_Noise60_B2.mat','P2_Noise70_B1.mat','P2_Noise70_B2.mat'};
LoadUtt=load('data\AMEND_I\utterances1110.mat');
LoadDelays=load('data\AMEND_I\delays1110.mat');

% Parameters for processing
Param.Fs = 50; % Sampling frequency of pupil data
Param.RemoveBeforeAndAfter = [35 100]*1e-3; % Samples within the time range before and after NaNs will set NaNs as well.
Param.MinLengthNaNRepair = 5; % Drop values (i.e., change to NaN) before and after NaNs only for contiguous NaNs of at least __ samples. 
LPWinSize = 0.5; % [s]: Window size of hamming-window for low-pass filtering
FilterWidth = round((LPWinSize*Param.Fs)/2); % [samples]: Width of hamming filter used for fixation duration
LPWindow = hamming(round(LPWinSize*Param.Fs));
LPWindow = LPWindow/sum(LPWindow); % Hamming-window
AudFs = 48000; % [Hz], sampling frequency of the audio files
AdapBL = 0.3; % [s], Duration of baseline prior to event
TimeStart = 20; % [s], Start of trial (disregard data prior to 'TimeStart')
TimeMinWin = 0.5; % [s], Minimum time of a window
TimeInitialMerge = 0.3; % [s], Time threshold for merging windows initially
TimeMerge = 2; % [s], Time threshold for merging windows after rejecting small windows
RejectRatio = 0.4; % Rejection threshold based on the ratio of NaNs in data
RejectDelay = 0.5; % [s], Rejection threshold based on delay between timestamps and n-samples
TimeStartW = 0.5; % [s], time before Utt/Lis starts
TimeEndW = 0; % [s], time after Utt/Lis starts
x = 1; % idx to store global values
NPs = 2; % N of TPs per Pair
NCond = numel(FileNames)/NPs; % N of conditions
NTPs = numel(subDirs)*numel(FileNames)/NCond; % Total N of TPs
TPsOrder = zeros(NTPs,NCond); % Vector that will contain the indexes of trials/file for each TP

GSW_MinMax = zeros(NTPs,2); % Vector with Min and Max pupil size (Speaking) for each TP
GLW_MinMax = zeros(NTPs,2); % Vector with Min and Max pupil size (Listening) for each TP
GSW_B_MinMax = zeros(NTPs,2); % Vector with Min and Max pupil size (Speaking Baselined) for each TP
GLW_B_MinMax = zeros(NTPs,2); % Vector with Min and Max pupil size (Listening Baselined) for each TP
GSW_MPDSTD = zeros(NTPs,2); % Vector with MPD and STD (Speaking) for each TP
GLW_MPDSTD = zeros(NTPs,2); % Vector with MPD and STD (Listening) for each TP

SW_Quiet_MinMax = zeros(NTPs,2); % Vector with Min and Max pupil size (Speaking) for each TP (Quiet)
LW_Quiet_MinMax = zeros(NTPs,2); % Vector with Min and Max pupil size (Listening) for each TP (Quiet)
SW_B_Quiet_MinMax = zeros(NTPs,2); % Vector with Min and Max pupil size (Speaking Baselined) for each TP (Quiet)
LW_B_Quiet_MinMax = zeros(NTPs,2); % Vector with Min and Max pupil size (Listening Baselined) for each TP (Quiet)
SW_Quiet_MPDSTD = zeros(NTPs,2); % Vector with MPD and STD (Speaking) for each TP (Quiet)
LW_Quiet_MPDSTD = zeros(NTPs,2); % Vector with MPD and STD (Listening) for each TP (Quiet)

SW_SHL_MinMax = zeros(NTPs,2); % Vector with Min and Max pupil size (Speaking) for each TP (SHL)
LW_SHL_MinMax = zeros(NTPs,2); % Vector with Min and Max pupil size (Listening) for each TP (SHL)
SW_B_SHL_MinMax = zeros(NTPs,2); % Vector with Min and Max pupil size (Speaking Baselined) for each TP (SHL)
LW_B_SHL_MinMax = zeros(NTPs,2); % Vector with Min and Max pupil size (Listening Baselined) for each TP (SHL)
SW_SHL_MPDSTD = zeros(NTPs,2); % Vector with MPD and STD (Speaking) for each TP (SHL)
LW_SHL_MPDSTD = zeros(NTPs,2); % Vector with MPD and STD (Listening) for each TP (SHL)

SW_N60_MinMax = zeros(NTPs,2); % Vector with Min and Max pupil size (Speaking) for each TP (N60)
LW_N60_MinMax = zeros(NTPs,2); % Vector with Min and Max pupil size (Listening) for each TP (N60)
SW_B_N60_MinMax = zeros(NTPs,2); % Vector with Min and Max pupil size (Speaking Baselined) for each TP (N60)
LW_B_N60_MinMax = zeros(NTPs,2); % Vector with Min and Max pupil size (Listening Baselined) for each TP (N60)
SW_N60_MPDSTD = zeros(NTPs,2); % Vector with MPD and STD (Speaking) for each TP (N60)
LW_N60_MPDSTD = zeros(NTPs,2); % Vector with MPD and STD (Listening) for each TP (N60)

SW_N70_MinMax = zeros(NTPs,2); % Vector with Min and Max pupil size (Speaking) for each TP (N70)
LW_N70_MinMax = zeros(NTPs,2); % Vector with Min and Max pupil size (Listening) for each TP (N70)
SW_B_N70_MinMax = zeros(NTPs,2); % Vector with Min and Max pupil size (Speaking Baselined) for each TP (N70)
LW_B_N70_MinMax = zeros(NTPs,2); % Vector with Min and Max pupil size (Listening Baselined) for each TP (N70)
SW_N70_MPDSTD = zeros(NTPs,2); % Vector with MPD and STD (Speaking) for each TP (N70)
LW_N70_MPDSTD = zeros(NTPs,2); % Vector with MPD and STD (Listening) for each TP (N70)

% Colors
SColor = [53, 155, 67]./255;
LColor = [204, 36, 0]./255;
QuietColor = [204, 152, 0]./255;
SHLColor = [123, 31, 162]./255;
N60Color = [0, 196, 215]./255;
N70Color = [2, 36, 223]./255;

% Variables
NCols=60*Param.Fs; % Duration (samples) of each window
NRows=100; % Number of windows per trial
NLayers=numel(FileNames)*numel(subDirs); % Number of trials

GSW = zeros(NLayers,NRows,NCols); % Global Speaking Windows
GLW = zeros(NLayers,NRows,NCols); % Global Listening Windows
SW_Quiet = zeros(NLayers,NRows,NCols); % Quiet Speaking Windows
LW_Quiet = zeros(NLayers,NRows,NCols); % Quiet Listening Windows
SW_SHL = zeros(NLayers,NRows,NCols); % SHL Speaking Windows
LW_SHL = zeros(NLayers,NRows,NCols); % SHL Listening Windows
SW_N60 = zeros(NLayers,NRows,NCols); % N60 Speaking Windows
LW_N60 = zeros(NLayers,NRows,NCols); % N60 Listening Windows
SW_N70 = zeros(NLayers,NRows,NCols); % N70 Speaking Windows
LW_N70 = zeros(NLayers,NRows,NCols); % N70 Listening Windows

GSW_B = zeros(NLayers,NRows,NCols); % Global Speaking Windows Adaptive-Baseline corrected
GLW_B = zeros(NLayers,NRows,NCols); % Global Listening Windows Adaptive-Baseline corrected
SW_B_Quiet = zeros(NLayers,NRows,NCols); % Quiet Speaking Windows Adaptive-Baseline corrected
LW_B_Quiet = zeros(NLayers,NRows,NCols); % Quiet Listening Windows Adaptive-Baseline corrected
SW_B_SHL = zeros(NLayers,NRows,NCols); % SHL Speaking Windows Adaptive-Baseline corrected
LW_B_SHL = zeros(NLayers,NRows,NCols); % SHL Listening Windows Adaptive-Baseline corrected
SW_B_N60 = zeros(NLayers,NRows,NCols); % N60 Speaking Windows Adaptive-Baseline corrected
LW_B_N60 = zeros(NLayers,NRows,NCols); % N60 Listening Windows Adaptive-Baseline corrected
SW_B_N70 = zeros(NLayers,NRows,NCols); % N70 Speaking Windows Adaptive-Baseline corrected
LW_B_N70 = zeros(NLayers,NRows,NCols); % N70 Listening Windows Adaptive-Baseline corrected

GSDur = zeros(NLayers,NRows); % Global Speaking Windows Duration
GLDur = zeros(NLayers,NRows); % Global Listening Windows Duration
SDur_Quiet = zeros(NLayers,NRows); % Quiet Speaking Windows Duration
LDur_Quiet = zeros(NLayers,NRows); % Quiet Listening Windows Duration
SDur_SHL = zeros(NLayers,NRows); % SHL Speaking Windows Duration
LDur_SHL = zeros(NLayers,NRows); % SHL Listening Windows Duration
SDur_N60 = zeros(NLayers,NRows); % N60 Speaking Windows Duration
LDur_N60 = zeros(NLayers,NRows); % N60 Listening Windows Duration
SDur_N70 = zeros(NLayers,NRows); % N70 Speaking Windows Duration
LDur_N70 = zeros(NLayers,NRows); % N70 Listening Windows Duration

figure;tiledlayout(1,2);ax1 = nexttile;ax2 = nexttile;
figure;tiledlayout(1,2);ax3 = nexttile;ax4 = nexttile;
figure;tiledlayout(1,2);ax5 = nexttile;ax6 = nexttile;
figure;tiledlayout(1,2);ax7 = nexttile;ax8 = nexttile;
figure;tiledlayout(1,2);ax9 = nexttile;ax10 = nexttile;

hold([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10],'on')
grid([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10],'on')

for q=1:numel(subDirs)
    PairIn = q;
    PairFiles=dir(['data\AMEND_I\Main',sprintf('%d',PairIn),'\*.mat']);
    PairUtt=LoadUtt.Utterances(PairIn,:);
    PairDelay=LoadDelays.TobAudDelay(PairIn,:);
    
    for i=1:numel(FileNames)        
        try
            alldata = load([PairFiles(1).folder, '\', cell2mat(FileNames(i))]);
        catch ME
            disp(['Warning: File ', PairFiles(1).folder, '\', cell2mat(FileNames(i)), ' not found (no Gaze data).']);
            continue
        end
        alldata_mat = cell2mat(alldata.data);
        
        % Replace blanks '[]' for 'NaN' in fields diameterLeft and diameterRight
        [alldata_mat(cellfun(@isempty,{alldata_mat.diameterLeft})).diameterLeft] = deal(NaN);
        [alldata_mat(cellfun(@isempty,{alldata_mat.diameterRight})).diameterRight] = deal(NaN);
        
        % Extract delay [s] (Duration_timeStamps - Duration_nsamples)
        EyeAudDelay=alldata_mat(end).timeStamp-alldata_mat(1).timeStamp-length(alldata_mat)/Param.Fs;
        
        % Skip file if the difference in duration from the number of
        % samples and the duration given by timestamps is bigger than 0.5 s
        if EyeAudDelay > RejectDelay
            disp(['Warning: File ',PairFiles(i).folder, '\', cell2mat(FileNames(i)), ' was rejected, too much delay (',sprintf('%0.2f',EyeAudDelay),'s).']);
            continue
        end
        
        LDiamRaw = [alldata_mat.diameterLeft];
        RDiamRaw = [alldata_mat.diameterRight];
        
        % Preprocessing - Setting outliers as NaNs (remove artifacts)
        LThreshOut = [mean(LDiamRaw,'omitnan')-2*std(LDiamRaw,'omitnan'),mean(LDiamRaw,'omitnan')+2*std(LDiamRaw,'omitnan')];
        RThreshOut = [mean(RDiamRaw,'omitnan')-2*std(RDiamRaw,'omitnan'),mean(RDiamRaw,'omitnan')+2*std(RDiamRaw,'omitnan')];
        for s=1:length(alldata_mat)
            if LDiamRaw(1,s) < LThreshOut(1) || LDiamRaw(1,s) > LThreshOut(2)
                LDiamRaw(1,s)=NaN;
            end
            if RDiamRaw(1,s) < RThreshOut(1) || RDiamRaw(1,s) > RThreshOut(2)
                RDiamRaw(1,s)=NaN;
            end
        end
        
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
        
        % Comment -> No LP filtering
%         LDiam = LDiamConv;
%         RDiam = RDiamConv;
        
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
        
        if DiamNaN/length(Diameter) >= RejectRatio
            disp(['Warning: File ',PairFiles(1).folder, '\', cell2mat(FileNames(i)), ' was rejected because it contains too many NaNs (',sprintf('%0.2f',100*DiamNaN/length(Diameter)),'%).'])
            continue
        end
        
        % Retrieve Utterances
        if contains(cell2mat(FileNames(i)),'P2')
            SpeakKey = 'utteranceCH1';
            ListenKey = 'utteranceCH2';
            SDelayKey = 'delayCH1';
            LDelayKey = 'delayCH2';
        elseif contains(cell2mat(FileNames(i)),'P1')
            SpeakKey = 'utteranceCH2';
            ListenKey = 'utteranceCH1';
            SDelayKey = 'delayCH2';
            LDelayKey = 'delayCH1';
        end

        if contains(cell2mat(FileNames(i)),'B1')
            SpeB = 0;
        elseif contains(cell2mat(FileNames(i)),'B2')
            SpeB = 1;
        end

        if contains(cell2mat(FileNames(i)),'Quiet')
            SpeCond = SpeB + 1;
        elseif contains(cell2mat(FileNames(i)),'SHL')
            SpeCond = SpeB + 3;
        elseif contains(cell2mat(FileNames(i)),'Noise60')
            SpeCond = SpeB + 5;
        elseif contains(cell2mat(FileNames(i)),'Noise70')
            SpeCond = SpeB + 7;
        end

        SpeakRaw = PairUtt{1,SpeCond}.(SpeakKey);
        ListenRaw = PairUtt{1,SpeCond}.(ListenKey);
        binResUtt = PairUtt{1,SpeCond}.binRes;
        try
            SDelayRaw = PairDelay{1,SpeCond}.(SDelayKey);
            LDelayRaw = PairDelay{1,SpeCond}.(LDelayKey);
        catch ME
            disp(['Warning: File ',PairFiles(1).folder, '\', cell2mat(FileNames(i)),' was rejected for not having an associated delay.']);
            continue
        end
        
        if ((SDelayRaw <= 0 | SDelayRaw >= 1) | (LDelayRaw <= 0 | LDelayRaw >= 1))
            disp(['Warning: File ', PairFiles(1).folder, '\', cell2mat(FileNames(i)),' was rejected because the delay between Tobii and Utterances is out of proportions (Speaking: ',sprintf('%0.3f',SDelayRaw(2)),' [s], Listening: ',sprintf('%0.3f',LDelayRaw(2)),' [s]).']);
            continue
        end
        
        if isempty(SpeakRaw) && isempty(ListenRaw)
            disp(['Warning: File ',PairFiles(1).folder, '\', cell2mat(FileNames(i)),' was rejected for not having associated Utterance windows.']);
            continue
        end
        
        % SAME PROCESSING AS IN W1.m
        % Downsample (rounding) Utt from 250 Hz (1/binRes) to 50 Hz, shift
        % in time to account for the time at which the audio recording
        % started (from 0 to 20 s only eye data) plus delay
        SpeakRaw(:,2:3)=round((SpeakRaw(:,2:3)*binResUtt+TimeStart)*Param.Fs+SDelayRaw(1)/2);
        ListenRaw(:,2:3)=round((ListenRaw(:,2:3)*binResUtt+TimeStart)*Param.Fs+LDelayRaw(1)/2);

        % Merge windows if duration between windows <= TimeInitialMerge (300 ms)
        SpeakMI = merge_windows(SpeakRaw, Param.Fs, TimeInitialMerge);
        ListenMI = merge_windows(ListenRaw, Param.Fs, TimeInitialMerge);
        
        % Discard windows if duration is < TimeMinWin (500 ms)
        SpeakDI = SpeakMI(SpeakMI(:,1)>TimeMinWin,:);
        ListenDI = ListenMI(ListenMI(:,1)>TimeMinWin,:);
        
        % Merge again if duration between windows <= TimeMerge (2 s)
        SpeakM = merge_windows(SpeakDI, Param.Fs, TimeMerge);
        ListenM = merge_windows(ListenDI, Param.Fs, TimeMerge);
        
        % Discard windows if duration is < 2*TimeMinWin (1 s)
        SpeakD = SpeakM(SpeakM(:,1)>2*TimeMinWin,:);
        ListenD = ListenM(ListenM(:,1)>2*TimeMinWin,:);
        
        % Added from W6.m -> Cut/Split overlaps
        [Speak,Listen] = overlap_windows(SpeakD,ListenD,Param.Fs);
        
%         Speak = SpeakRaw;
%         Listen = ListenRaw;
                
        % Time-locked indexes (based on Start or End of events)
        WSpeakIdx=[Speak(:,2)-TimeStartW*Param.Fs,Speak(:,2),Speak(:,3)+TimeEndW*Param.Fs];
        WListenIdx=[Listen(:,2)-TimeStartW*Param.Fs,Listen(:,2),Listen(:,2)+TimeEndW*Param.Fs];
%         EWSpeakIdx=[Speak(:,3)-TimeStartW*Param.Fs,Speak(:,3),Speak(:,3)+TimeEndW*Param.Fs];
%         EWListenIdx=[Listen(:,3)-TimeStartW*Param.Fs,Listen(:,3),Listen(:,3)+TimeEndW*Param.Fs];
        
%         t_Diam = linspace(0,length(Diameter)./Param.Fs,length(Diameter));
        
        % Storing Speaking/Listening by conditions
        if contains(cell2mat(FileNames(i)),'Quiet')
            for j=1:size(Speak,1)
                SW_Quiet(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_Quiet)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                SW_B_Quiet(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_B_Quiet)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                SDur_Quiet(i,j) = Speak(j,1);
            end
            for j=1:size(Listen,1)
                LW_Quiet(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_Quiet)-length(WListenIdx(j,1):Listen(j,3)))'];
                LW_B_Quiet(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_B_Quiet)-length(WListenIdx(j,1):Listen(j,3)))'];
                LDur_Quiet(i,j) = Listen(j,1);
            end
        elseif contains(cell2mat(FileNames(i)),'SHL')
            for j=1:size(Speak,1)
                SW_SHL(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_SHL)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                SW_B_SHL(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_B_SHL)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                SDur_SHL(i,j) = Speak(j,1);
            end
            for j=1:size(Listen,1)
                LW_SHL(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_SHL)-length(WListenIdx(j,1):Listen(j,3)))'];
                LW_B_SHL(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_B_SHL)-length(WListenIdx(j,1):Listen(j,3)))'];
                LDur_SHL(i,j) = Listen(j,1);
            end
        elseif contains(cell2mat(FileNames(i)),'Noise60')
            for j=1:size(Speak,1)
                SW_N60(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                SW_B_N60(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_B_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                SDur_N60(i,j) = Speak(j,1);
            end
            for j=1:size(Listen,1)
                LW_N60(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                LW_B_N60(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_B_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                LDur_N60(i,j) = Listen(j,1);
            end
        elseif contains(cell2mat(FileNames(i)),'Noise70')
            for j=1:size(Speak,1)
                SW_N70(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                SW_B_N70(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_B_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                SDur_N70(i,j) = Speak(j,1);
            end
            for j=1:size(Listen,1)
                LW_N70(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                LW_B_N70(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_B_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                LDur_N70(x,j) = Listen(j,1);
            end
        end
        
        % Storing Global Speaking/Listening pupil sizes
        for j=1:size(Speak,1)
            % Add nan-padding when necessary
            GSW(x,j,:) = [Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(GSW)-length(WSpeakIdx(j,1):Speak(j,3)))'];
            GSW_B(x,j,:) = [Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(GSW_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
            GSDur(x,j) = Speak(j,1);
        end
        
        for j=1:size(Listen,1)
            % Add nan-padding when necessary
            GLW(x,j,:) = [Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(GLW)-length(WListenIdx(j,1):Listen(j,3)))'];
            GLW_B(x,j,:) = [Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(GLW_B)-length(WListenIdx(j,1):Listen(j,3)))'];
            GLDur(x,j) = Listen(j,1);
        end
        % Store idx of non-rejected files associated with TPs
        if contains(cell2mat(FileNames(i)),'P2')
            TPsOrder(2*q,i-NCond) = x;
        elseif contains(cell2mat(FileNames(i)),'P1')
            TPsOrder(2*q-1,i) = x;
        end
        
        % Increase index of num of files used
        x=x+1;
        
    end
end
%% Out-of-files-loop calculations
% Clean empty rows and layers, set 0's to NaN
GSW(~any(GSW,[2 3]),:,:)=[];GSW(:,~any(GSW,[1 3]),:)=[];GSW(GSW==0)=NaN;
GSW_B(~any(GSW_B,[2 3]),:,:)=[];GSW_B(:,~any(GSW_B,[1 3]),:)=[];GSW_B(GSW_B==0)=NaN;
GLW(~any(GLW,[2 3]),:,:)=[];GLW(:,~any(GLW,[1 3]),:)=[];GLW(GLW==0)=NaN;
GLW_B(~any(GLW_B,[2 3]),:,:)=[];GLW_B(:,~any(GLW_B,[1 3]),:)=[];GLW_B(GLW_B==0)=NaN;

SW_Quiet(~any(SW_Quiet,[2 3]),:,:)=[];SW_Quiet(:,~any(SW_Quiet,[1 3]),:)=[];SW_Quiet(SW_Quiet==0)=NaN;
SW_B_Quiet(~any(SW_B_Quiet,[2 3]),:,:)=[];SW_B_Quiet(:,~any(SW_B_Quiet,[1 3]),:)=[];SW_B_Quiet(SW_B_Quiet==0)=NaN;
LW_Quiet(~any(LW_Quiet,[2 3]),:,:)=[];LW_Quiet(:,~any(LW_Quiet,[1 3]),:)=[];LW_Quiet(LW_Quiet==0)=NaN;
LW_B_Quiet(~any(LW_B_Quiet,[2 3]),:,:)=[];LW_B_Quiet(:,~any(LW_B_Quiet,[1 3]),:)=[];LW_B_Quiet(LW_B_Quiet==0)=NaN;

SW_SHL(~any(SW_SHL,[2 3]),:,:)=[];SW_SHL(:,~any(SW_SHL,[1 3]),:)=[];SW_SHL(SW_SHL==0)=NaN;
SW_B_SHL(~any(SW_B_SHL,[2 3]),:,:)=[];SW_B_SHL(:,~any(SW_B_SHL,[1 3]),:)=[];SW_B_SHL(SW_B_SHL==0)=NaN;
LW_SHL(~any(LW_SHL,[2 3]),:,:)=[];LW_SHL(:,~any(LW_SHL,[1 3]),:)=[];LW_SHL(LW_SHL==0)=NaN;
LW_B_SHL(~any(LW_B_SHL,[2 3]),:,:)=[];LW_B_SHL(:,~any(LW_B_SHL,[1 3]),:)=[];LW_B_SHL(LW_B_SHL==0)=NaN;

SW_N60(~any(SW_N60,[2 3]),:,:)=[];SW_N60(:,~any(SW_N60,[1 3]),:)=[];SW_N60(SW_N60==0)=NaN;
SW_B_N60(~any(SW_B_N60,[2 3]),:,:)=[];SW_B_N60(:,~any(SW_B_N60,[1 3]),:)=[];SW_B_N60(SW_B_N60==0)=NaN;
LW_N60(~any(LW_N60,[2 3]),:,:)=[];LW_N60(:,~any(LW_N60,[1 3]),:)=[];LW_N60(LW_N60==0)=NaN;
LW_B_N60(~any(LW_B_N60,[2 3]),:,:)=[];LW_B_N60(:,~any(LW_B_N60,[1 3]),:)=[];LW_B_N60(LW_B_N60==0)=NaN;

SW_N70(~any(SW_N70,[2 3]),:,:)=[];SW_N70(:,~any(SW_N70,[1 3]),:)=[];SW_N70(SW_N70==0)=NaN;
SW_B_N70(~any(SW_B_N70,[2 3]),:,:)=[];SW_B_N70(:,~any(SW_B_N70,[1 3]),:)=[];SW_B_N70(SW_B_N70==0)=NaN;
LW_N70(~any(LW_N70,[2 3]),:,:)=[];LW_N70(:,~any(LW_N70,[1 3]),:)=[];LW_N70(LW_N70==0)=NaN;
LW_B_N70(~any(LW_B_N70,[2 3]),:,:)=[];LW_B_N70(:,~any(LW_B_N70,[1 3]),:)=[];LW_B_N70(LW_B_N70==0)=NaN;

% Retrieve Min and Max from Global Windows for each participant
for i=1:NTPs    
    if ~isempty(nonzeros(TPsOrder(i,:)))
        GSW_MinMax(i,1) = min(reshape(GSW(min(nonzeros(TPsOrder(i,:))):max(nonzeros(TPsOrder(i,:))),:,:),[],1));
        GLW_MinMax(i,1) = min(reshape(GLW(min(nonzeros(TPsOrder(i,:))):max(nonzeros(TPsOrder(i,:))),:,:),[],1));
        GSW_MinMax(i,2) = max(reshape(GSW(min(nonzeros(TPsOrder(i,:))):max(nonzeros(TPsOrder(i,:))),:,:),[],1));
        GLW_MinMax(i,2) = max(reshape(GLW(min(nonzeros(TPsOrder(i,:))):max(nonzeros(TPsOrder(i,:))),:,:),[],1));
        GSW_B_MinMax(i,1) = min(reshape(GSW_B(min(nonzeros(TPsOrder(i,:))):max(nonzeros(TPsOrder(i,:))),:,:),[],1));
        GLW_B_MinMax(i,1) = min(reshape(GLW_B(min(nonzeros(TPsOrder(i,:))):max(nonzeros(TPsOrder(i,:))),:,:),[],1));
        GSW_B_MinMax(i,2) = max(reshape(GSW_B(min(nonzeros(TPsOrder(i,:))):max(nonzeros(TPsOrder(i,:))),:,:),[],1));
        GLW_B_MinMax(i,2) = max(reshape(GLW_B(min(nonzeros(TPsOrder(i,:))):max(nonzeros(TPsOrder(i,:))),:,:),[],1));
        GSW_MPDSTD(i,1) = mean(reshape(GSW(min(nonzeros(TPsOrder(i,:))):max(nonzeros(TPsOrder(i,:))),:,:),[],1),[1 2],'omitnan');
        GLW_MPDSTD(i,1) = mean(reshape(GLW(min(nonzeros(TPsOrder(i,:))):max(nonzeros(TPsOrder(i,:))),:,:),[],1),[1 2],'omitnan');
        GSW_MPDSTD(i,2) = std(reshape(GSW(min(nonzeros(TPsOrder(i,:))):max(nonzeros(TPsOrder(i,:))),:,:),[],1),[],'omitnan');
        GLW_MPDSTD(i,2) = std(reshape(GLW(min(nonzeros(TPsOrder(i,:))):max(nonzeros(TPsOrder(i,:))),:,:),[],1),[],'omitnan');
        
        % Quiet condition
        if ~isempty(nonzeros(TPsOrder(i,1:2)))
            TPCondIdx = nonzeros(TPsOrder(i,1:2));
            if length(TPCondIdx)==1
                TPCondIdx = TPCondIdx*ones(2,1);
            end
            SW_Quiet_MinMax(i,1) = min(reshape(GSW(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            LW_Quiet_MinMax(i,1) = min(reshape(GLW(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            SW_Quiet_MinMax(i,2) = max(reshape(GSW(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            LW_Quiet_MinMax(i,2) = max(reshape(GLW(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            SW_B_Quiet_MinMax(i,1) = min(reshape(GSW_B(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            LW_B_Quiet_MinMax(i,1) = min(reshape(GLW_B(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            SW_B_Quiet_MinMax(i,2) = max(reshape(GSW_B(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            LW_B_Quiet_MinMax(i,2) = max(reshape(GLW_B(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            SW_Quiet_MPDSTD(i,1) = mean(reshape(GSW(TPCondIdx(1):TPCondIdx(2),:,:),[],1),[1 2],'omitnan');
            LW_Quiet_MPDSTD(i,1) = mean(reshape(GLW(TPCondIdx(1):TPCondIdx(2),:,:),[],1),[1 2],'omitnan');
            SW_Quiet_MPDSTD(i,2) = std(reshape(GSW(TPCondIdx(1):TPCondIdx(2),:,:),[],1),[],'omitnan');
            LW_Quiet_MPDSTD(i,2) = std(reshape(GLW(TPCondIdx(1):TPCondIdx(2),:,:),[],1),[],'omitnan');
        end
        
        % SHL condition
        if ~isempty(nonzeros(TPsOrder(i,3:4)))
            TPCondIdx = nonzeros(TPsOrder(i,3:4));
            if length(TPCondIdx)==1
                TPCondIdx = TPCondIdx*ones(2,1);
            end
            SW_SHL_MinMax(i,1) = min(reshape(GSW(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            LW_SHL_MinMax(i,1) = min(reshape(GLW(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            SW_SHL_MinMax(i,2) = max(reshape(GSW(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            LW_SHL_MinMax(i,2) = max(reshape(GLW(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            SW_B_SHL_MinMax(i,1) = min(reshape(GSW_B(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            LW_B_SHL_MinMax(i,1) = min(reshape(GLW_B(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            SW_B_SHL_MinMax(i,2) = max(reshape(GSW_B(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            LW_B_SHL_MinMax(i,2) = max(reshape(GLW_B(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            SW_SHL_MPDSTD(i,1) = mean(reshape(GSW(TPCondIdx(1):TPCondIdx(2),:,:),[],1),[1 2],'omitnan');
            LW_SHL_MPDSTD(i,1) = mean(reshape(GLW(TPCondIdx(1):TPCondIdx(2),:,:),[],1),[1 2],'omitnan');
            SW_SHL_MPDSTD(i,2) = std(reshape(GSW(TPCondIdx(1):TPCondIdx(2),:,:),[],1),[],'omitnan');
            LW_SHL_MPDSTD(i,2) = std(reshape(GLW(TPCondIdx(1):TPCondIdx(2),:,:),[],1),[],'omitnan');
        end
        
        % N60 condition
        if ~isempty(nonzeros(TPsOrder(i,5:6)))
            TPCondIdx = nonzeros(TPsOrder(i,5:6));
            if length(TPCondIdx)==1
                TPCondIdx = TPCondIdx*ones(2,1);
            end
            SW_N60_MinMax(i,1) = min(reshape(GSW(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            LW_N60_MinMax(i,1) = min(reshape(GLW(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            SW_N60_MinMax(i,2) = max(reshape(GSW(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            LW_N60_MinMax(i,2) = max(reshape(GLW(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            SW_B_N60_MinMax(i,1) = min(reshape(GSW_B(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            LW_B_N60_MinMax(i,1) = min(reshape(GLW_B(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            SW_B_N60_MinMax(i,2) = max(reshape(GSW_B(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            LW_B_N60_MinMax(i,2) = max(reshape(GLW_B(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            SW_N60_MPDSTD(i,1) = mean(reshape(GSW(TPCondIdx(1):TPCondIdx(2),:,:),[],1),[1 2],'omitnan');
            LW_N60_MPDSTD(i,1) = mean(reshape(GLW(TPCondIdx(1):TPCondIdx(2),:,:),[],1),[1 2],'omitnan');
            SW_N60_MPDSTD(i,2) = std(reshape(GSW(TPCondIdx(1):TPCondIdx(2),:,:),[],1),[],'omitnan');
            LW_N60_MPDSTD(i,2) = std(reshape(GLW(TPCondIdx(1):TPCondIdx(2),:,:),[],1),[],'omitnan');
        end
        
        % N70 condition
        if ~isempty(nonzeros(TPsOrder(i,7:8)))
            TPCondIdx = nonzeros(TPsOrder(i,7:8));
            if length(TPCondIdx)==1
                TPCondIdx = TPCondIdx*ones(2,1);
            end
            SW_N70_MinMax(i,1) = min(reshape(GSW(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            LW_N70_MinMax(i,1) = min(reshape(GLW(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            SW_N70_MinMax(i,2) = max(reshape(GSW(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            LW_N70_MinMax(i,2) = max(reshape(GLW(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            SW_B_N70_MinMax(i,1) = min(reshape(GSW_B(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            LW_B_N70_MinMax(i,1) = min(reshape(GLW_B(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            SW_B_N70_MinMax(i,2) = max(reshape(GSW_B(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            LW_B_N70_MinMax(i,2) = max(reshape(GLW_B(TPCondIdx(1):TPCondIdx(2),:,:),[],1));
            SW_N70_MPDSTD(i,1) = mean(reshape(GSW(TPCondIdx(1):TPCondIdx(2),:,:),[],1),[1 2],'omitnan');
            LW_N70_MPDSTD(i,1) = mean(reshape(GLW(TPCondIdx(1):TPCondIdx(2),:,:),[],1),[1 2],'omitnan');
            SW_N70_MPDSTD(i,2) = std(reshape(GSW(TPCondIdx(1):TPCondIdx(2),:,:),[],1),[],'omitnan');
            LW_N70_MPDSTD(i,2) = std(reshape(GLW(TPCondIdx(1):TPCondIdx(2),:,:),[],1),[],'omitnan');
        end
    end
end

% Initialize X_R, X_RB and X_Z
GSW_R = zeros(size(GSW)); % Global Speaking Windows Range normalized
GLW_R = zeros(size(GLW)); % Global Listening Windows Range normalized
SW_R_Quiet = zeros(size(SW_Quiet)); % Quiet Speaking Windows Range normalized
LW_R_Quiet = zeros(size(LW_Quiet)); % Quiet Listening Windows Range normalized
SW_R_SHL = zeros(size(SW_SHL)); % SHL Speaking Windows Range normalized
LW_R_SHL = zeros(size(LW_SHL)); % SHL Listening Windows Range normalized
SW_R_N60 = zeros(size(SW_N60)); % N60 Speaking Windows Range normalized
LW_R_N60 = zeros(size(LW_N60)); % N60 Listening Windows Range normalized
SW_R_N70 = zeros(size(SW_N70)); % N70 Speaking Windows Range normalized
LW_R_N70 = zeros(size(LW_N70)); % N70 Listening Windows Range normalized

GSW_RB = zeros(size(GSW)); % Global Speaking Windows Baselined Range normalized
GLW_RB = zeros(size(GLW)); % Global Listening Windows Baselined Range normalized
SW_RB_Quiet = zeros(size(SW_Quiet)); % Quiet Speaking Windows Baselined Range normalized
LW_RB_Quiet = zeros(size(LW_Quiet)); % Quiet Listening Windows Baselined Range normalized
SW_RB_SHL = zeros(size(SW_SHL)); % SHL Speaking Windows Baselined Range normalized
LW_RB_SHL = zeros(size(LW_SHL)); % SHL Listening Windows Baselined Range normalized
SW_RB_N60 = zeros(size(SW_N60)); % N60 Speaking Windows Baselined Range normalized
LW_RB_N60 = zeros(size(LW_N60)); % N60 Listening Windows Baselined Range normalized
SW_RB_N70 = zeros(size(SW_N70)); % N70 Speaking Windows Baselined Range normalized
LW_RB_N70 = zeros(size(LW_N70)); % N70 Listening Windows Baselined Range normalized

GSW_Z = zeros(size(GSW)); % Global Speaking Windows Z-score normalized
GLW_Z = zeros(size(GLW)); % Global Listening Windows Z-score normalized
SW_Z_Quiet = zeros(size(SW_Quiet)); % Quiet Speaking Windows Z-score normalized
LW_Z_Quiet = zeros(size(LW_Quiet)); % Quiet Listening Windows Z-score normalized
SW_Z_SHL = zeros(size(SW_SHL)); % SHL Speaking Windows Z-score normalized
LW_Z_SHL = zeros(size(LW_SHL)); % SHL Listening Windows Z-score normalized
SW_Z_N60 = zeros(size(SW_N60)); % N60 Speaking Windows Z-score normalized
LW_Z_N60 = zeros(size(LW_N60)); % N60 Listening Windows Z-score normalized
SW_Z_N70 = zeros(size(SW_N70)); % N70 Speaking Windows Z-score normalized
LW_Z_N70 = zeros(size(LW_N70)); % N70 Listening Windows Z-score normalized

% Calculate X_R, X_RB and X_Z
for j = 1:size(GSW,1)
    TPRow = ismember(TPsOrder,j);
    GSW_R(j,:,:) = (GSW(j,:,:)-GSW_MinMax(mod(find(TPRow,1),size(TPRow,1)),1))./(GSW_MinMax(mod(find(TPRow,1),size(TPRow,1)),2)-GSW_MinMax(mod(find(TPRow,1),size(TPRow,1)),1));
    GSW_RB(j,:,:) = (GSW_B(j,:,:)-GSW_B_MinMax(mod(find(TPRow,1),size(TPRow,1)),1))./(GSW_B_MinMax(mod(find(TPRow,1),size(TPRow,1)),2)-GSW_B_MinMax(mod(find(TPRow,1),size(TPRow,1)),1));
    GSW_Z(j,:,:) = (GSW(j,:,:)-GSW_MPDSTD(mod(find(TPRow,1),size(TPRow,1)),1))./GSW_MPDSTD(mod(find(TPRow,1),size(TPRow,1)),2);
    GLW_R(j,:,:) = (GLW(j,:,:)-GLW_MinMax(mod(find(TPRow,1),size(TPRow,1)),1))./(GLW_MinMax(mod(find(TPRow,1),size(TPRow,1)),2)-GLW_MinMax(mod(find(TPRow,1),size(TPRow,1)),1));
    GLW_RB(j,:,:) = (GLW_B(j,:,:)-GLW_B_MinMax(mod(find(TPRow,1),size(TPRow,1)),1))./(GLW_B_MinMax(mod(find(TPRow,1),size(TPRow,1)),2)-GLW_B_MinMax(mod(find(TPRow,1),size(TPRow,1)),1));
    GLW_Z(j,:,:) = (GLW(j,:,:)-GLW_MPDSTD(mod(find(TPRow,1),size(TPRow,1)),1))./GLW_MPDSTD(mod(find(TPRow,1),size(TPRow,1)),2);
end

TPsQuiet = nonzeros(TPsOrder(:,1:2));
for j = 1:size(SW_Quiet,1)
    TPRow = ismember(TPsOrder,TPsQuiet(j));
    if ~isempty(nonzeros(TPRow))
        SW_R_Quiet(j,:,:) = (SW_Quiet(j,:,:)-SW_Quiet_MinMax(mod(find(TPRow,1),size(TPRow,1)),1))./(SW_Quiet_MinMax(mod(find(TPRow,1),size(TPRow,1)),2)-SW_Quiet_MinMax(mod(find(TPRow,1),size(TPRow,1)),1));
        SW_RB_Quiet(j,:,:) = (SW_B_Quiet(j,:,:)-SW_B_Quiet_MinMax(mod(find(TPRow,1),size(TPRow,1)),1))./(SW_B_Quiet_MinMax(mod(find(TPRow,1),size(TPRow,1)),2)-SW_B_Quiet_MinMax(mod(find(TPRow,1),size(TPRow,1)),1));
        SW_Z_Quiet(j,:,:) = (SW_Quiet(j,:,:)-SW_Quiet_MPDSTD(mod(find(TPRow,1),size(TPRow,1)),1))./SW_Quiet_MPDSTD(mod(find(TPRow,1),size(TPRow,1)),2);
        LW_R_Quiet(j,:,:) = (LW_Quiet(j,:,:)-LW_Quiet_MinMax(mod(find(TPRow,1),size(TPRow,1)),1))./(LW_Quiet_MinMax(mod(find(TPRow,1),size(TPRow,1)),2)-LW_Quiet_MinMax(mod(find(TPRow,1),size(TPRow,1)),1));
        LW_RB_Quiet(j,:,:) = (LW_B_Quiet(j,:,:)-LW_B_Quiet_MinMax(mod(find(TPRow,1),size(TPRow,1)),1))./(LW_B_Quiet_MinMax(mod(find(TPRow,1),size(TPRow,1)),2)-LW_B_Quiet_MinMax(mod(find(TPRow,1),size(TPRow,1)),1));
        LW_Z_Quiet(j,:,:) = (LW_Quiet(j,:,:)-LW_Quiet_MPDSTD(mod(find(TPRow,1),size(TPRow,1)),1))./LW_Quiet_MPDSTD(mod(find(TPRow,1),size(TPRow,1)),2);
    end
end

TPsSHL = nonzeros(TPsOrder(:,3:4));
for j = 1:size(SW_SHL,1)
    TPRow = ismember(TPsOrder,TPsSHL(j));
    if ~isempty(nonzeros(TPRow))
        SW_R_SHL(j,:,:) = (SW_SHL(j,:,:)-SW_SHL_MinMax(mod(find(TPRow,1),size(TPRow,1)),1))./(SW_SHL_MinMax(mod(find(TPRow,1),size(TPRow,1)),2)-SW_SHL_MinMax(mod(find(TPRow,1),size(TPRow,1)),1));
        SW_RB_SHL(j,:,:) = (SW_B_SHL(j,:,:)-SW_B_SHL_MinMax(mod(find(TPRow,1),size(TPRow,1)),1))./(SW_B_SHL_MinMax(mod(find(TPRow,1),size(TPRow,1)),2)-SW_B_SHL_MinMax(mod(find(TPRow,1),size(TPRow,1)),1));
        SW_Z_SHL(j,:,:) = (SW_SHL(j,:,:)-SW_SHL_MPDSTD(mod(find(TPRow,1),size(TPRow,1)),1))./SW_SHL_MPDSTD(mod(find(TPRow,1),size(TPRow,1)),2);
        LW_R_SHL(j,:,:) = (LW_SHL(j,:,:)-LW_SHL_MinMax(mod(find(TPRow,1),size(TPRow,1)),1))./(LW_SHL_MinMax(mod(find(TPRow,1),size(TPRow,1)),2)-LW_SHL_MinMax(mod(find(TPRow,1),size(TPRow,1)),1));
        LW_RB_SHL(j,:,:) = (LW_B_SHL(j,:,:)-LW_B_SHL_MinMax(mod(find(TPRow,1),size(TPRow,1)),1))./(LW_B_SHL_MinMax(mod(find(TPRow,1),size(TPRow,1)),2)-LW_B_SHL_MinMax(mod(find(TPRow,1),size(TPRow,1)),1));
        LW_Z_SHL(j,:,:) = (LW_SHL(j,:,:)-LW_SHL_MPDSTD(mod(find(TPRow,1),size(TPRow,1)),1))./LW_SHL_MPDSTD(mod(find(TPRow,1),size(TPRow,1)),2);
    end
end

TPsN60 = nonzeros(TPsOrder(:,5:6));
for j = 1:size(SW_N60,1)
    TPRow = ismember(TPsOrder,TPsN60(j));
    if ~isempty(nonzeros(TPRow))
        SW_R_N60(j,:,:) = (SW_N60(j,:,:)-SW_N60_MinMax(mod(find(TPRow,1),size(TPRow,1)),1))./(SW_N60_MinMax(mod(find(TPRow,1),size(TPRow,1)),2)-SW_N60_MinMax(mod(find(TPRow,1),size(TPRow,1)),1));
        SW_RB_N60(j,:,:) = (SW_B_N60(j,:,:)-SW_B_N60_MinMax(mod(find(TPRow,1),size(TPRow,1)),1))./(SW_B_N60_MinMax(mod(find(TPRow,1),size(TPRow,1)),2)-SW_B_N60_MinMax(mod(find(TPRow,1),size(TPRow,1)),1));
        SW_Z_N60(j,:,:) = (SW_N60(j,:,:)-SW_N60_MPDSTD(mod(find(TPRow,1),size(TPRow,1)),1))./SW_N60_MPDSTD(mod(find(TPRow,1),size(TPRow,1)),2);
        LW_R_N60(j,:,:) = (LW_N60(j,:,:)-LW_N60_MinMax(mod(find(TPRow,1),size(TPRow,1)),1))./(LW_N60_MinMax(mod(find(TPRow,1),size(TPRow,1)),2)-LW_N60_MinMax(mod(find(TPRow,1),size(TPRow,1)),1));
        LW_RB_N60(j,:,:) = (LW_B_N60(j,:,:)-LW_B_N60_MinMax(mod(find(TPRow,1),size(TPRow,1)),1))./(LW_B_N60_MinMax(mod(find(TPRow,1),size(TPRow,1)),2)-LW_B_N60_MinMax(mod(find(TPRow,1),size(TPRow,1)),1));
        LW_Z_N60(j,:,:) = (LW_N60(j,:,:)-LW_N60_MPDSTD(mod(find(TPRow,1),size(TPRow,1)),1))./LW_N60_MPDSTD(mod(find(TPRow,1),size(TPRow,1)),2);
    end
end

TPsN70 = nonzeros(TPsOrder(:,7:8));
for j = 1:size(SW_N70,1)
    TPRow = ismember(TPsOrder,TPsN70(j));
    if ~isempty(nonzeros(TPRow))
        SW_R_N70(j,:,:) = (SW_N70(j,:,:)-SW_N70_MinMax(mod(find(TPRow,1),size(TPRow,1)),1))./(SW_N70_MinMax(mod(find(TPRow,1),size(TPRow,1)),2)-SW_N70_MinMax(mod(find(TPRow,1),size(TPRow,1)),1));
        SW_RB_N70(j,:,:) = (SW_B_N70(j,:,:)-SW_B_N70_MinMax(mod(find(TPRow,1),size(TPRow,1)),1))./(SW_B_N70_MinMax(mod(find(TPRow,1),size(TPRow,1)),2)-SW_B_N70_MinMax(mod(find(TPRow,1),size(TPRow,1)),1));
        SW_Z_N70(j,:,:) = (SW_N70(j,:,:)-SW_N70_MPDSTD(mod(find(TPRow,1),size(TPRow,1)),1))./SW_N70_MPDSTD(mod(find(TPRow,1),size(TPRow,1)),2);
        LW_R_N70(j,:,:) = (LW_N70(j,:,:)-LW_N70_MinMax(mod(find(TPRow,1),size(TPRow,1)),1))./(LW_N70_MinMax(mod(find(TPRow,1),size(TPRow,1)),2)-LW_N70_MinMax(mod(find(TPRow,1),size(TPRow,1)),1));
        LW_RB_N70(j,:,:) = (LW_B_N70(j,:,:)-LW_B_N70_MinMax(mod(find(TPRow,1),size(TPRow,1)),1))./(LW_B_N70_MinMax(mod(find(TPRow,1),size(TPRow,1)),2)-LW_B_N70_MinMax(mod(find(TPRow,1),size(TPRow,1)),1));
        LW_Z_N70(j,:,:) = (LW_N70(j,:,:)-LW_N70_MPDSTD(mod(find(TPRow,1),size(TPRow,1)),1))./LW_N70_MPDSTD(mod(find(TPRow,1),size(TPRow,1)),2);
    end
end

% Calculate means omitting NaNs, LP filtering with mean-padding at start
% end of each group of Speaking/Listening windows
GSW_Mean = ndnanfilter(reshape(mean(GSW,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
GSW_B_Mean = ndnanfilter(reshape(mean(GSW_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
GSW_R_Mean = ndnanfilter(reshape(mean(GSW_R,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
GSW_RB_Mean = ndnanfilter(reshape(mean(GSW_RB,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
GSW_Z_Mean = ndnanfilter(reshape(mean(GSW_Z,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
GLW_Mean = ndnanfilter(reshape(mean(GLW,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
GLW_B_Mean = ndnanfilter(reshape(mean(GLW_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
GLW_R_Mean = ndnanfilter(reshape(mean(GLW_R,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
GLW_RB_Mean = ndnanfilter(reshape(mean(GLW_RB,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
GLW_Z_Mean = ndnanfilter(reshape(mean(GLW_Z,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

SW_Quiet_Mean = ndnanfilter(reshape(mean(SW_Quiet,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_B_Quiet_Mean = ndnanfilter(reshape(mean(SW_B_Quiet,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_R_Quiet_Mean = ndnanfilter(reshape(mean(SW_R_Quiet,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_RB_Quiet_Mean = ndnanfilter(reshape(mean(SW_RB_Quiet,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_Z_Quiet_Mean = ndnanfilter(reshape(mean(SW_Z_Quiet,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_Quiet_Mean = ndnanfilter(reshape(mean(LW_Quiet,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_B_Quiet_Mean = ndnanfilter(reshape(mean(LW_B_Quiet,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_R_Quiet_Mean = ndnanfilter(reshape(mean(LW_R_Quiet,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_RB_Quiet_Mean = ndnanfilter(reshape(mean(LW_RB_Quiet,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_Z_Quiet_Mean = ndnanfilter(reshape(mean(LW_Z_Quiet,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

SW_SHL_Mean = ndnanfilter(reshape(mean(SW_SHL,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_B_SHL_Mean = ndnanfilter(reshape(mean(SW_B_SHL,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_R_SHL_Mean = ndnanfilter(reshape(mean(SW_R_SHL,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_RB_SHL_Mean = ndnanfilter(reshape(mean(SW_RB_SHL,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_Z_SHL_Mean = ndnanfilter(reshape(mean(SW_Z_SHL,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_SHL_Mean = ndnanfilter(reshape(mean(LW_SHL,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_B_SHL_Mean = ndnanfilter(reshape(mean(LW_B_SHL,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_R_SHL_Mean = ndnanfilter(reshape(mean(LW_R_SHL,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_RB_SHL_Mean = ndnanfilter(reshape(mean(LW_RB_SHL,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_Z_SHL_Mean = ndnanfilter(reshape(mean(LW_Z_SHL,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

SW_N60_Mean = ndnanfilter(reshape(mean(SW_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_B_N60_Mean = ndnanfilter(reshape(mean(SW_B_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_R_N60_Mean = ndnanfilter(reshape(mean(SW_R_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_RB_N60_Mean = ndnanfilter(reshape(mean(SW_RB_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_Z_N60_Mean = ndnanfilter(reshape(mean(SW_Z_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_N60_Mean = ndnanfilter(reshape(mean(LW_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_B_N60_Mean = ndnanfilter(reshape(mean(LW_B_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_R_N60_Mean = ndnanfilter(reshape(mean(LW_R_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_RB_N60_Mean = ndnanfilter(reshape(mean(LW_RB_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_Z_N60_Mean = ndnanfilter(reshape(mean(LW_Z_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

SW_N70_Mean = ndnanfilter(reshape(mean(SW_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_B_N70_Mean = ndnanfilter(reshape(mean(SW_B_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_R_N70_Mean = ndnanfilter(reshape(mean(SW_R_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_RB_N70_Mean = ndnanfilter(reshape(mean(SW_RB_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_Z_N70_Mean = ndnanfilter(reshape(mean(SW_Z_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_N70_Mean = ndnanfilter(reshape(mean(LW_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_B_N70_Mean = ndnanfilter(reshape(mean(LW_B_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_R_N70_Mean = ndnanfilter(reshape(mean(LW_R_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_RB_N70_Mean = ndnanfilter(reshape(mean(LW_RB_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_Z_N70_Mean = ndnanfilter(reshape(mean(LW_Z_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

% Calculate SEM as: std(X)/sqrt(squeeze(sum(~isnan(X),[1 2]))))
GSW_SEM = (squeeze(std(GSW,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(GSW),[1 2]))))';
GSW_B_SEM = (squeeze(std(GSW_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(GSW_B),[1 2]))))';
GSW_R_SEM = (squeeze(std(GSW_R,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(GSW_R),[1 2]))))';
GSW_RB_SEM = (squeeze(std(GSW_RB,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(GSW_RB),[1 2]))))';
GSW_Z_SEM = (squeeze(std(GSW_Z,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(GSW_Z),[1 2]))))';
GLW_SEM = (squeeze(std(GLW,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(GLW),[1 2]))))';
GLW_B_SEM = (squeeze(std(GLW_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(GLW_B),[1 2]))))';
GLW_R_SEM = (squeeze(std(GLW_R,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(GLW_R),[1 2]))))';
GLW_RB_SEM = (squeeze(std(GLW_RB,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(GLW_RB),[1 2]))))';
GLW_Z_SEM = (squeeze(std(GLW_Z,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(GLW_Z),[1 2]))))';

SW_Quiet_SEM = (squeeze(std(SW_Quiet,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_Quiet),[1 2]))))';
SW_B_Quiet_SEM = (squeeze(std(SW_B_Quiet,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_B_Quiet),[1 2]))))';
SW_R_Quiet_SEM = (squeeze(std(SW_R_Quiet,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_R_Quiet),[1 2]))))';
SW_RB_Quiet_SEM = (squeeze(std(SW_RB_Quiet,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_RB_Quiet),[1 2]))))';
SW_Z_Quiet_SEM = (squeeze(std(SW_Z_Quiet,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_Z_Quiet),[1 2]))))';
LW_Quiet_SEM = (squeeze(std(LW_Quiet,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_Quiet),[1 2]))))';
LW_B_Quiet_SEM = (squeeze(std(LW_B_Quiet,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_B_Quiet),[1 2]))))';
LW_R_Quiet_SEM = (squeeze(std(LW_R_Quiet,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_R_Quiet),[1 2]))))';
LW_RB_Quiet_SEM = (squeeze(std(LW_RB_Quiet,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_RB_Quiet),[1 2]))))';
LW_Z_Quiet_SEM = (squeeze(std(LW_Z_Quiet,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_Z_Quiet),[1 2]))))';

SW_SHL_SEM = (squeeze(std(SW_SHL,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_SHL),[1 2]))))';
SW_B_SHL_SEM = (squeeze(std(SW_B_SHL,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_B_SHL),[1 2]))))';
SW_R_SHL_SEM = (squeeze(std(SW_R_SHL,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_R_SHL),[1 2]))))';
SW_RB_SHL_SEM = (squeeze(std(SW_RB_SHL,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_RB_SHL),[1 2]))))';
SW_Z_SHL_SEM = (squeeze(std(SW_Z_SHL,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_Z_SHL),[1 2]))))';
LW_SHL_SEM = (squeeze(std(LW_SHL,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_SHL),[1 2]))))';
LW_B_SHL_SEM = (squeeze(std(LW_B_SHL,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_B_SHL),[1 2]))))';
LW_R_SHL_SEM = (squeeze(std(LW_R_SHL,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_R_SHL),[1 2]))))';
LW_RB_SHL_SEM = (squeeze(std(LW_RB_SHL,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_RB_SHL),[1 2]))))';
LW_Z_SHL_SEM = (squeeze(std(LW_Z_SHL,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_Z_SHL),[1 2]))))';

SW_N60_SEM = (squeeze(std(SW_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_N60),[1 2]))))';
SW_B_N60_SEM = (squeeze(std(SW_B_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_B_N60),[1 2]))))';
SW_R_N60_SEM = (squeeze(std(SW_R_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_R_N60),[1 2]))))';
SW_RB_N60_SEM = (squeeze(std(SW_RB_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_RB_N60),[1 2]))))';
SW_Z_N60_SEM = (squeeze(std(SW_Z_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_Z_N60),[1 2]))))';
LW_N60_SEM = (squeeze(std(LW_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_N60),[1 2]))))';
LW_B_N60_SEM = (squeeze(std(LW_B_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_B_N60),[1 2]))))';
LW_R_N60_SEM = (squeeze(std(LW_R_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_R_N60),[1 2]))))';
LW_RB_N60_SEM = (squeeze(std(LW_RB_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_RB_N60),[1 2]))))';
LW_Z_N60_SEM = (squeeze(std(LW_Z_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_Z_N60),[1 2]))))';

SW_N70_SEM = (squeeze(std(SW_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_N70),[1 2]))))';
SW_B_N70_SEM = (squeeze(std(SW_B_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_B_N70),[1 2]))))';
SW_R_N70_SEM = (squeeze(std(SW_R_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_R_N70),[1 2]))))';
SW_RB_N70_SEM = (squeeze(std(SW_RB_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_RB_N70),[1 2]))))';
SW_Z_N70_SEM = (squeeze(std(SW_Z_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_Z_N70),[1 2]))))';
LW_N70_SEM = (squeeze(std(LW_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_N70),[1 2]))))';
LW_B_N70_SEM = (squeeze(std(LW_B_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_B_N70),[1 2]))))';
LW_R_N70_SEM = (squeeze(std(LW_R_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_R_N70),[1 2]))))';
LW_RB_N70_SEM = (squeeze(std(LW_RB_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_RB_N70),[1 2]))))';
LW_Z_N70_SEM = (squeeze(std(LW_Z_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_Z_N70),[1 2]))))';

%% Global Plots
% Plot event onset and baseline markers
xline(ax1,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax2,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax3,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax4,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax5,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax6,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax7,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax8,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax9,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax10,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')

% xline(ax1,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
% xline(ax2,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
xline(ax3,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
xline(ax4,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
% xline(ax5,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
% xline(ax6,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
xline(ax7,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
xline(ax8,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
% xline(ax9,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
% xline(ax10,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')

plot(ax1,linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3)),GSW_Mean,color=SColor,linewidth=2)
plot(ax1,linspace(-TimeStartW,size(SW_Quiet,3)/Param.Fs,size(SW_Quiet,3)),SW_Quiet_Mean,color=QuietColor)
plot(ax1,linspace(-TimeStartW,size(SW_SHL,3)/Param.Fs,size(SW_SHL,3)),SW_SHL_Mean,color=SHLColor)
plot(ax1,linspace(-TimeStartW,size(SW_N60,3)/Param.Fs,size(SW_N60,3)),SW_N60_Mean,color=N60Color)
plot(ax1,linspace(-TimeStartW,size(SW_N70,3)/Param.Fs,size(SW_N70,3)),SW_N70_Mean,color=N70Color)
plot(ax2,nan,color=SColor,linewidth=2) % plot nans to show color in legend
plot(ax2,linspace(-TimeStartW,size(GLW,3)/Param.Fs,size(GLW,3)),GLW_Mean,color=LColor,linewidth=2)
plot(ax2,linspace(-TimeStartW,size(LW_Quiet,3)/Param.Fs,size(LW_Quiet,3)),LW_Quiet_Mean,color=QuietColor)
plot(ax2,linspace(-TimeStartW,size(LW_SHL,3)/Param.Fs,size(LW_SHL,3)),LW_SHL_Mean,color=SHLColor)
plot(ax2,linspace(-TimeStartW,size(LW_N60,3)/Param.Fs,size(LW_N60,3)),LW_N60_Mean,color=N60Color)
plot(ax2,linspace(-TimeStartW,size(LW_N70,3)/Param.Fs,size(LW_N70,3)),LW_N70_Mean,color=N70Color)

% plot(ax1,linspace(-TimeStartW,size(SW_SHL,3)/Param.Fs,size(SW_SHL,3)),mean([SW_Quiet_Mean;SW_SHL_Mean;SW_N60_Mean;SW_N70_Mean],1,'omitnan'),'k--')
% plot(ax2,linspace(-TimeStartW,size(LW_SHL,3)/Param.Fs,size(LW_SHL,3)),mean([LW_Quiet_Mean;LW_SHL_Mean;LW_N60_Mean;LW_N70_Mean],1,'omitnan'),'k--')

plot(ax3,linspace(-TimeStartW,size(GSW_B,3)/Param.Fs,size(GSW_B,3)),GSW_B_Mean,color=SColor,linewidth=2)
plot(ax3,linspace(-TimeStartW,size(SW_B_Quiet,3)/Param.Fs,size(SW_B_Quiet,3)),SW_B_Quiet_Mean,color=QuietColor)
plot(ax3,linspace(-TimeStartW,size(SW_B_SHL,3)/Param.Fs,size(SW_B_SHL,3)),SW_B_SHL_Mean,color=SHLColor)
plot(ax3,linspace(-TimeStartW,size(SW_B_N60,3)/Param.Fs,size(SW_B_N60,3)),SW_B_N60_Mean,color=N60Color)
plot(ax3,linspace(-TimeStartW,size(SW_B_N70,3)/Param.Fs,size(SW_B_N70,3)),SW_B_N70_Mean,color=N70Color)
plot(ax4,nan,color=SColor,linewidth=2) % plot nans to show color in legend
plot(ax4,linspace(-TimeStartW,size(GLW_B,3)/Param.Fs,size(GLW_B,3)),GLW_B_Mean,color=LColor,linewidth=2)
plot(ax4,linspace(-TimeStartW,size(LW_B_Quiet,3)/Param.Fs,size(LW_B_Quiet,3)),LW_B_Quiet_Mean,color=QuietColor)
plot(ax4,linspace(-TimeStartW,size(LW_B_SHL,3)/Param.Fs,size(LW_B_SHL,3)),LW_B_SHL_Mean,color=SHLColor)
plot(ax4,linspace(-TimeStartW,size(LW_B_N60,3)/Param.Fs,size(LW_B_N60,3)),LW_B_N60_Mean,color=N60Color)
plot(ax4,linspace(-TimeStartW,size(LW_B_N70,3)/Param.Fs,size(LW_B_N70,3)),LW_B_N70_Mean,color=N70Color)

plot(ax5,linspace(-TimeStartW,size(GSW_R,3)/Param.Fs,size(GSW_R,3)),GSW_R_Mean,color=SColor,linewidth=2)
plot(ax5,linspace(-TimeStartW,size(SW_R_Quiet,3)/Param.Fs,size(SW_R_Quiet,3)),SW_R_Quiet_Mean,color=QuietColor)
plot(ax5,linspace(-TimeStartW,size(SW_R_SHL,3)/Param.Fs,size(SW_R_SHL,3)),SW_R_SHL_Mean,color=SHLColor)
plot(ax5,linspace(-TimeStartW,size(SW_R_N60,3)/Param.Fs,size(SW_R_N60,3)),SW_R_N60_Mean,color=N60Color)
plot(ax5,linspace(-TimeStartW,size(SW_R_N70,3)/Param.Fs,size(SW_R_N70,3)),SW_R_N70_Mean,color=N70Color)
plot(ax6,nan,color=SColor,linewidth=2) % plot nans to show color in legend
plot(ax6,linspace(-TimeStartW,size(GLW_R,3)/Param.Fs,size(GLW_R,3)),GLW_R_Mean,color=LColor,linewidth=2)
plot(ax6,linspace(-TimeStartW,size(LW_R_Quiet,3)/Param.Fs,size(LW_R_Quiet,3)),LW_R_Quiet_Mean,color=QuietColor)
plot(ax6,linspace(-TimeStartW,size(LW_R_SHL,3)/Param.Fs,size(LW_R_SHL,3)),LW_R_SHL_Mean,color=SHLColor)
plot(ax6,linspace(-TimeStartW,size(LW_R_N60,3)/Param.Fs,size(LW_R_N60,3)),LW_R_N60_Mean,color=N60Color)
plot(ax6,linspace(-TimeStartW,size(LW_R_N70,3)/Param.Fs,size(LW_R_N70,3)),LW_R_N70_Mean,color=N70Color)

plot(ax7,linspace(-TimeStartW,size(GSW_RB,3)/Param.Fs,size(GSW_RB,3)),GSW_RB_Mean,color=SColor,linewidth=2)
plot(ax7,linspace(-TimeStartW,size(SW_RB_Quiet,3)/Param.Fs,size(SW_RB_Quiet,3)),SW_RB_Quiet_Mean,color=QuietColor)
plot(ax7,linspace(-TimeStartW,size(SW_RB_SHL,3)/Param.Fs,size(SW_RB_SHL,3)),SW_RB_SHL_Mean,color=SHLColor)
plot(ax7,linspace(-TimeStartW,size(SW_RB_N60,3)/Param.Fs,size(SW_RB_N60,3)),SW_RB_N60_Mean,color=N60Color)
plot(ax7,linspace(-TimeStartW,size(SW_RB_N70,3)/Param.Fs,size(SW_RB_N70,3)),SW_RB_N70_Mean,color=N70Color)
plot(ax8,nan,color=SColor,linewidth=2) % plot nans to show color in legend
plot(ax8,linspace(-TimeStartW,size(GLW_RB,3)/Param.Fs,size(GLW_RB,3)),GLW_RB_Mean,color=LColor,linewidth=2)
plot(ax8,linspace(-TimeStartW,size(LW_RB_Quiet,3)/Param.Fs,size(LW_RB_Quiet,3)),LW_RB_Quiet_Mean,color=QuietColor)
plot(ax8,linspace(-TimeStartW,size(LW_RB_SHL,3)/Param.Fs,size(LW_RB_SHL,3)),LW_RB_SHL_Mean,color=SHLColor)
plot(ax8,linspace(-TimeStartW,size(LW_RB_N60,3)/Param.Fs,size(LW_RB_N60,3)),LW_RB_N60_Mean,color=N60Color)
plot(ax8,linspace(-TimeStartW,size(LW_RB_N70,3)/Param.Fs,size(LW_RB_N70,3)),LW_RB_N70_Mean,color=N70Color)

plot(ax9,linspace(-TimeStartW,size(GSW_Z,3)/Param.Fs,size(GSW_Z,3)),GSW_Z_Mean,color=SColor,linewidth=2)
plot(ax9,linspace(-TimeStartW,size(SW_Z_Quiet,3)/Param.Fs,size(SW_Z_Quiet,3)),SW_Z_Quiet_Mean,color=QuietColor)
plot(ax9,linspace(-TimeStartW,size(SW_Z_SHL,3)/Param.Fs,size(SW_Z_SHL,3)),SW_Z_SHL_Mean,color=SHLColor)
plot(ax9,linspace(-TimeStartW,size(SW_Z_N60,3)/Param.Fs,size(SW_Z_N60,3)),SW_Z_N60_Mean,color=N60Color)
plot(ax9,linspace(-TimeStartW,size(SW_Z_N70,3)/Param.Fs,size(SW_Z_N70,3)),SW_Z_N70_Mean,color=N70Color)
plot(ax10,nan,color=SColor,linewidth=2) % plot nans to show color in legend
plot(ax10,linspace(-TimeStartW,size(GLW_Z,3)/Param.Fs,size(GLW_Z,3)),GLW_Z_Mean,color=LColor,linewidth=2)
plot(ax10,linspace(-TimeStartW,size(LW_Z_Quiet,3)/Param.Fs,size(LW_Z_Quiet,3)),LW_Z_Quiet_Mean,color=QuietColor)
plot(ax10,linspace(-TimeStartW,size(LW_Z_SHL,3)/Param.Fs,size(LW_Z_SHL,3)),LW_Z_SHL_Mean,color=SHLColor)
plot(ax10,linspace(-TimeStartW,size(LW_Z_N60,3)/Param.Fs,size(LW_Z_N60,3)),LW_Z_N60_Mean,color=N60Color)
plot(ax10,linspace(-TimeStartW,size(LW_Z_N70,3)/Param.Fs,size(LW_Z_N70,3)),LW_Z_N70_Mean,color=N70Color)

% Mean and SEM: Set NaNs to 0 in order to use "fill" correctly
GSW_Mean(isnan(GSW_Mean))=0;GSW_SEM(isnan(GSW_SEM))=0;
SW_Quiet_Mean(isnan(SW_Quiet_Mean))=0;SW_Quiet_SEM(isnan(SW_Quiet_SEM))=0;
SW_SHL_Mean(isnan(SW_SHL_Mean))=0;SW_SHL_SEM(isnan(SW_SHL_SEM))=0;
SW_N60_Mean(isnan(SW_N60_Mean))=0;SW_N60_SEM(isnan(SW_N60_SEM))=0;
SW_N70_Mean(isnan(SW_N70_Mean))=0;SW_N70_SEM(isnan(SW_N70_SEM))=0;
GLW_Mean(isnan(GLW_Mean))=0;GLW_SEM(isnan(GLW_SEM))=0;
LW_Quiet_Mean(isnan(LW_Quiet_Mean))=0;LW_Quiet_SEM(isnan(LW_Quiet_SEM))=0;
LW_SHL_Mean(isnan(LW_SHL_Mean))=0;LW_SHL_SEM(isnan(LW_SHL_SEM))=0;
LW_N60_Mean(isnan(LW_N60_Mean))=0;LW_N60_SEM(isnan(LW_N60_SEM))=0;
LW_N70_Mean(isnan(LW_N70_Mean))=0;LW_N70_SEM(isnan(LW_N70_SEM))=0;

GSW_B_Mean(isnan(GSW_B_Mean))=0;GSW_B_SEM(isnan(GSW_B_SEM))=0;
SW_B_Quiet_Mean(isnan(SW_B_Quiet_Mean))=0;SW_B_Quiet_SEM(isnan(SW_B_Quiet_SEM))=0;
SW_B_SHL_Mean(isnan(SW_B_SHL_Mean))=0;SW_B_SHL_SEM(isnan(SW_B_SHL_SEM))=0;
SW_B_N60_Mean(isnan(SW_B_N60_Mean))=0;SW_B_N60_SEM(isnan(SW_B_N60_SEM))=0;
SW_B_N70_Mean(isnan(SW_B_N70_Mean))=0;SW_B_N70_SEM(isnan(SW_B_N70_SEM))=0;
GLW_B_Mean(isnan(GLW_B_Mean))=0;GLW_B_SEM(isnan(GLW_B_SEM))=0;
LW_B_Quiet_Mean(isnan(LW_B_Quiet_Mean))=0;LW_B_Quiet_SEM(isnan(LW_B_Quiet_SEM))=0;
LW_B_SHL_Mean(isnan(LW_B_SHL_Mean))=0;LW_B_SHL_SEM(isnan(LW_B_SHL_SEM))=0;
LW_B_N60_Mean(isnan(LW_B_N60_Mean))=0;LW_B_N60_SEM(isnan(LW_B_N60_SEM))=0;
LW_B_N70_Mean(isnan(LW_B_N70_Mean))=0;LW_B_N70_SEM(isnan(LW_B_N70_SEM))=0;

GSW_R_Mean(isnan(GSW_R_Mean))=0;GSW_R_SEM(isnan(GSW_R_SEM))=0;
SW_R_Quiet_Mean(isnan(SW_R_Quiet_Mean))=0;SW_R_Quiet_SEM(isnan(SW_R_Quiet_SEM))=0;
SW_R_SHL_Mean(isnan(SW_R_SHL_Mean))=0;SW_R_SHL_SEM(isnan(SW_R_SHL_SEM))=0;
SW_R_N60_Mean(isnan(SW_R_N60_Mean))=0;SW_R_N60_SEM(isnan(SW_R_N60_SEM))=0;
SW_R_N70_Mean(isnan(SW_R_N70_Mean))=0;SW_R_N70_SEM(isnan(SW_R_N70_SEM))=0;
GLW_R_Mean(isnan(GLW_R_Mean))=0;GLW_R_SEM(isnan(GLW_R_SEM))=0;
LW_R_Quiet_Mean(isnan(LW_R_Quiet_Mean))=0;LW_R_Quiet_SEM(isnan(LW_R_Quiet_SEM))=0;
LW_R_SHL_Mean(isnan(LW_R_SHL_Mean))=0;LW_R_SHL_SEM(isnan(LW_R_SHL_SEM))=0;
LW_R_N60_Mean(isnan(LW_R_N60_Mean))=0;LW_R_N60_SEM(isnan(LW_R_N60_SEM))=0;
LW_R_N70_Mean(isnan(LW_R_N70_Mean))=0;LW_R_N70_SEM(isnan(LW_R_N70_SEM))=0;

GSW_RB_Mean(isnan(GSW_RB_Mean))=0;GSW_RB_SEM(isnan(GSW_RB_SEM))=0;
SW_RB_Quiet_Mean(isnan(SW_RB_Quiet_Mean))=0;SW_RB_Quiet_SEM(isnan(SW_RB_Quiet_SEM))=0;
SW_RB_SHL_Mean(isnan(SW_RB_SHL_Mean))=0;SW_RB_SHL_SEM(isnan(SW_RB_SHL_SEM))=0;
SW_RB_N60_Mean(isnan(SW_RB_N60_Mean))=0;SW_RB_N60_SEM(isnan(SW_RB_N60_SEM))=0;
SW_RB_N70_Mean(isnan(SW_RB_N70_Mean))=0;SW_RB_N70_SEM(isnan(SW_RB_N70_SEM))=0;
GLW_RB_Mean(isnan(GLW_RB_Mean))=0;GLW_RB_SEM(isnan(GLW_RB_SEM))=0;
LW_RB_Quiet_Mean(isnan(LW_RB_Quiet_Mean))=0;LW_RB_Quiet_SEM(isnan(LW_RB_Quiet_SEM))=0;
LW_RB_SHL_Mean(isnan(LW_RB_SHL_Mean))=0;LW_RB_SHL_SEM(isnan(LW_RB_SHL_SEM))=0;
LW_RB_N60_Mean(isnan(LW_RB_N60_Mean))=0;LW_RB_N60_SEM(isnan(LW_RB_N60_SEM))=0;
LW_RB_N70_Mean(isnan(LW_RB_N70_Mean))=0;LW_RB_N70_SEM(isnan(LW_RB_N70_SEM))=0;

GSW_Z_Mean(isnan(GSW_Z_Mean))=0;GSW_Z_SEM(isnan(GSW_Z_SEM))=0;
SW_Z_Quiet_Mean(isnan(SW_Z_Quiet_Mean))=0;SW_Z_Quiet_SEM(isnan(SW_Z_Quiet_SEM))=0;
SW_Z_SHL_Mean(isnan(SW_Z_SHL_Mean))=0;SW_Z_SHL_SEM(isnan(SW_Z_SHL_SEM))=0;
SW_Z_N60_Mean(isnan(SW_Z_N60_Mean))=0;SW_Z_N60_SEM(isnan(SW_Z_N60_SEM))=0;
SW_Z_N70_Mean(isnan(SW_Z_N70_Mean))=0;SW_Z_N70_SEM(isnan(SW_Z_N70_SEM))=0;
GLW_Z_Mean(isnan(GLW_Z_Mean))=0;GLW_Z_SEM(isnan(GLW_Z_SEM))=0;
LW_Z_Quiet_Mean(isnan(LW_Z_Quiet_Mean))=0;LW_Z_Quiet_SEM(isnan(LW_Z_Quiet_SEM))=0;
LW_Z_SHL_Mean(isnan(LW_Z_SHL_Mean))=0;LW_Z_SHL_SEM(isnan(LW_Z_SHL_SEM))=0;
LW_Z_N60_Mean(isnan(LW_Z_N60_Mean))=0;LW_Z_N60_SEM(isnan(LW_Z_N60_SEM))=0;
LW_Z_N70_Mean(isnan(LW_Z_N70_Mean))=0;LW_Z_N70_SEM(isnan(LW_Z_N70_SEM))=0;

fill(ax1,[linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3)), flipud(linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3))')'],[(GSW_Mean+GSW_SEM), flipud((GSW_Mean-GSW_SEM)')'],SColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_Quiet,3)/Param.Fs,size(SW_Quiet,3)), flipud(linspace(-TimeStartW,size(SW_Quiet,3)/Param.Fs,size(SW_Quiet,3))')'],[(SW_Quiet_Mean+SW_Quiet_SEM), flipud((SW_Quiet_Mean-SW_Quiet_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_SHL,3)/Param.Fs,size(SW_SHL,3)), flipud(linspace(-TimeStartW,size(SW_SHL,3)/Param.Fs,size(SW_SHL,3))')'],[(SW_SHL_Mean+SW_SHL_SEM), flipud((SW_SHL_Mean-SW_SHL_SEM)')'],SHLColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_N60,3)/Param.Fs,size(SW_N60,3)), flipud(linspace(-TimeStartW,size(SW_N60,3)/Param.Fs,size(SW_N60,3))')'],[(SW_N60_Mean+SW_N60_SEM), flipud((SW_N60_Mean-SW_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_N70,3)/Param.Fs,size(SW_N70,3)), flipud(linspace(-TimeStartW,size(SW_N70,3)/Param.Fs,size(SW_N70,3))')'],[(SW_N70_Mean+SW_N70_SEM), flipud((SW_N70_Mean-SW_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(GLW,3)/Param.Fs,size(GLW,3)), flipud(linspace(-TimeStartW,size(GLW,3)/Param.Fs,size(GLW,3))')'],[(GLW_Mean+GLW_SEM), flipud((GLW_Mean-GLW_SEM)')'],LColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_Quiet,3)/Param.Fs,size(LW_Quiet,3)), flipud(linspace(-TimeStartW,size(LW_Quiet,3)/Param.Fs,size(LW_Quiet,3))')'],[(LW_Quiet_Mean+LW_Quiet_SEM), flipud((LW_Quiet_Mean-LW_Quiet_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_SHL,3)/Param.Fs,size(LW_SHL,3)), flipud(linspace(-TimeStartW,size(LW_SHL,3)/Param.Fs,size(LW_SHL,3))')'],[(LW_SHL_Mean+LW_SHL_SEM), flipud((LW_SHL_Mean-LW_SHL_SEM)')'],SHLColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_N60,3)/Param.Fs,size(LW_N60,3)), flipud(linspace(-TimeStartW,size(LW_N60,3)/Param.Fs,size(LW_N60,3))')'],[(LW_N60_Mean+LW_N60_SEM), flipud((LW_N60_Mean-LW_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_N70,3)/Param.Fs,size(LW_N70,3)), flipud(linspace(-TimeStartW,size(LW_N70,3)/Param.Fs,size(LW_N70,3))')'],[(LW_N70_Mean+LW_N70_SEM), flipud((LW_N70_Mean-LW_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')

fill(ax3,[linspace(-TimeStartW,size(GSW_B,3)/Param.Fs,size(GSW_B,3)), flipud(linspace(-TimeStartW,size(GSW_B,3)/Param.Fs,size(GSW_B,3))')'],[(GSW_B_Mean+GSW_B_SEM), flipud((GSW_B_Mean-GSW_B_SEM)')'],SColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SW_B_Quiet,3)/Param.Fs,size(SW_B_Quiet,3)), flipud(linspace(-TimeStartW,size(SW_B_Quiet,3)/Param.Fs,size(SW_B_Quiet,3))')'],[(SW_B_Quiet_Mean+SW_B_Quiet_SEM), flipud((SW_B_Quiet_Mean-SW_B_Quiet_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SW_B_SHL,3)/Param.Fs,size(SW_B_SHL,3)), flipud(linspace(-TimeStartW,size(SW_B_SHL,3)/Param.Fs,size(SW_B_SHL,3))')'],[(SW_B_SHL_Mean+SW_B_SHL_SEM), flipud((SW_B_SHL_Mean-SW_B_SHL_SEM)')'],SHLColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SW_B_N60,3)/Param.Fs,size(SW_B_N60,3)), flipud(linspace(-TimeStartW,size(SW_B_N60,3)/Param.Fs,size(SW_B_N60,3))')'],[(SW_B_N60_Mean+SW_B_N60_SEM), flipud((SW_B_N60_Mean-SW_B_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SW_B_N70,3)/Param.Fs,size(SW_B_N70,3)), flipud(linspace(-TimeStartW,size(SW_B_N70,3)/Param.Fs,size(SW_B_N70,3))')'],[(SW_B_N70_Mean+SW_B_N70_SEM), flipud((SW_B_N70_Mean-SW_B_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(GLW_B,3)/Param.Fs,size(GLW_B,3)), flipud(linspace(-TimeStartW,size(GLW_B,3)/Param.Fs,size(GLW_B,3))')'],[(GLW_B_Mean+GLW_B_SEM), flipud((GLW_B_Mean-GLW_B_SEM)')'],LColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LW_B_Quiet,3)/Param.Fs,size(LW_B_Quiet,3)), flipud(linspace(-TimeStartW,size(LW_B_Quiet,3)/Param.Fs,size(LW_B_Quiet,3))')'],[(LW_B_Quiet_Mean+LW_B_Quiet_SEM), flipud((LW_B_Quiet_Mean-LW_B_Quiet_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LW_B_SHL,3)/Param.Fs,size(LW_B_SHL,3)), flipud(linspace(-TimeStartW,size(LW_B_SHL,3)/Param.Fs,size(LW_B_SHL,3))')'],[(LW_B_SHL_Mean+LW_B_SHL_SEM), flipud((LW_B_SHL_Mean-LW_B_SHL_SEM)')'],SHLColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LW_B_N60,3)/Param.Fs,size(LW_B_N60,3)), flipud(linspace(-TimeStartW,size(LW_B_N60,3)/Param.Fs,size(LW_B_N60,3))')'],[(LW_B_N60_Mean+LW_B_N60_SEM), flipud((LW_B_N60_Mean-LW_B_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LW_B_N70,3)/Param.Fs,size(LW_B_N70,3)), flipud(linspace(-TimeStartW,size(LW_B_N70,3)/Param.Fs,size(LW_B_N70,3))')'],[(LW_B_N70_Mean+LW_B_N70_SEM), flipud((LW_B_N70_Mean-LW_B_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')

fill(ax5,[linspace(-TimeStartW,size(GSW_R,3)/Param.Fs,size(GSW_R,3)), flipud(linspace(-TimeStartW,size(GSW_R,3)/Param.Fs,size(GSW_R,3))')'],[(GSW_R_Mean+GSW_R_SEM), flipud((GSW_R_Mean-GSW_R_SEM)')'],SColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax5,[linspace(-TimeStartW,size(SW_R_Quiet,3)/Param.Fs,size(SW_R_Quiet,3)), flipud(linspace(-TimeStartW,size(SW_R_Quiet,3)/Param.Fs,size(SW_R_Quiet,3))')'],[(SW_R_Quiet_Mean+SW_R_Quiet_SEM), flipud((SW_R_Quiet_Mean-SW_R_Quiet_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax5,[linspace(-TimeStartW,size(SW_R_SHL,3)/Param.Fs,size(SW_R_SHL,3)), flipud(linspace(-TimeStartW,size(SW_R_SHL,3)/Param.Fs,size(SW_R_SHL,3))')'],[(SW_R_SHL_Mean+SW_R_SHL_SEM), flipud((SW_R_SHL_Mean-SW_R_SHL_SEM)')'],SHLColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax5,[linspace(-TimeStartW,size(SW_R_N60,3)/Param.Fs,size(SW_R_N60,3)), flipud(linspace(-TimeStartW,size(SW_R_N60,3)/Param.Fs,size(SW_R_N60,3))')'],[(SW_R_N60_Mean+SW_R_N60_SEM), flipud((SW_R_N60_Mean-SW_R_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax5,[linspace(-TimeStartW,size(SW_R_N70,3)/Param.Fs,size(SW_R_N70,3)), flipud(linspace(-TimeStartW,size(SW_R_N70,3)/Param.Fs,size(SW_R_N70,3))')'],[(SW_R_N70_Mean+SW_R_N70_SEM), flipud((SW_R_N70_Mean-SW_R_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax6,[linspace(-TimeStartW,size(GLW_R,3)/Param.Fs,size(GLW_R,3)), flipud(linspace(-TimeStartW,size(GLW_R,3)/Param.Fs,size(GLW_R,3))')'],[(GLW_R_Mean+GLW_R_SEM), flipud((GLW_R_Mean-GLW_R_SEM)')'],LColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax6,[linspace(-TimeStartW,size(LW_R_Quiet,3)/Param.Fs,size(LW_R_Quiet,3)), flipud(linspace(-TimeStartW,size(LW_R_Quiet,3)/Param.Fs,size(LW_R_Quiet,3))')'],[(LW_R_Quiet_Mean+LW_R_Quiet_SEM), flipud((LW_R_Quiet_Mean-LW_R_Quiet_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax6,[linspace(-TimeStartW,size(LW_R_SHL,3)/Param.Fs,size(LW_R_SHL,3)), flipud(linspace(-TimeStartW,size(LW_R_SHL,3)/Param.Fs,size(LW_R_SHL,3))')'],[(LW_R_SHL_Mean+LW_R_SHL_SEM), flipud((LW_R_SHL_Mean-LW_R_SHL_SEM)')'],SHLColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax6,[linspace(-TimeStartW,size(LW_R_N60,3)/Param.Fs,size(LW_R_N60,3)), flipud(linspace(-TimeStartW,size(LW_R_N60,3)/Param.Fs,size(LW_R_N60,3))')'],[(LW_R_N60_Mean+LW_R_N60_SEM), flipud((LW_R_N60_Mean-LW_R_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax6,[linspace(-TimeStartW,size(LW_R_N70,3)/Param.Fs,size(LW_R_N70,3)), flipud(linspace(-TimeStartW,size(LW_R_N70,3)/Param.Fs,size(LW_R_N70,3))')'],[(LW_R_N70_Mean+LW_R_N70_SEM), flipud((LW_R_N70_Mean-LW_R_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')

fill(ax7,[linspace(-TimeStartW,size(GSW_RB,3)/Param.Fs,size(GSW_RB,3)), flipud(linspace(-TimeStartW,size(GSW_RB,3)/Param.Fs,size(GSW_RB,3))')'],[(GSW_RB_Mean+GSW_RB_SEM), flipud((GSW_RB_Mean-GSW_RB_SEM)')'],SColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax7,[linspace(-TimeStartW,size(SW_RB_Quiet,3)/Param.Fs,size(SW_RB_Quiet,3)), flipud(linspace(-TimeStartW,size(SW_RB_Quiet,3)/Param.Fs,size(SW_RB_Quiet,3))')'],[(SW_RB_Quiet_Mean+SW_RB_Quiet_SEM), flipud((SW_RB_Quiet_Mean-SW_RB_Quiet_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax7,[linspace(-TimeStartW,size(SW_RB_SHL,3)/Param.Fs,size(SW_RB_SHL,3)), flipud(linspace(-TimeStartW,size(SW_RB_SHL,3)/Param.Fs,size(SW_RB_SHL,3))')'],[(SW_RB_SHL_Mean+SW_RB_SHL_SEM), flipud((SW_RB_SHL_Mean-SW_RB_SHL_SEM)')'],SHLColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax7,[linspace(-TimeStartW,size(SW_RB_N60,3)/Param.Fs,size(SW_RB_N60,3)), flipud(linspace(-TimeStartW,size(SW_RB_N60,3)/Param.Fs,size(SW_RB_N60,3))')'],[(SW_RB_N60_Mean+SW_RB_N60_SEM), flipud((SW_RB_N60_Mean-SW_RB_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax7,[linspace(-TimeStartW,size(SW_RB_N70,3)/Param.Fs,size(SW_RB_N70,3)), flipud(linspace(-TimeStartW,size(SW_RB_N70,3)/Param.Fs,size(SW_RB_N70,3))')'],[(SW_RB_N70_Mean+SW_RB_N70_SEM), flipud((SW_RB_N70_Mean-SW_RB_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax8,[linspace(-TimeStartW,size(GLW_RB,3)/Param.Fs,size(GLW_RB,3)), flipud(linspace(-TimeStartW,size(GLW_RB,3)/Param.Fs,size(GLW_RB,3))')'],[(GLW_RB_Mean+GLW_RB_SEM), flipud((GLW_RB_Mean-GLW_RB_SEM)')'],LColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax8,[linspace(-TimeStartW,size(LW_RB_Quiet,3)/Param.Fs,size(LW_RB_Quiet,3)), flipud(linspace(-TimeStartW,size(LW_RB_Quiet,3)/Param.Fs,size(LW_RB_Quiet,3))')'],[(LW_RB_Quiet_Mean+LW_RB_Quiet_SEM), flipud((LW_RB_Quiet_Mean-LW_RB_Quiet_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax8,[linspace(-TimeStartW,size(LW_RB_SHL,3)/Param.Fs,size(LW_RB_SHL,3)), flipud(linspace(-TimeStartW,size(LW_RB_SHL,3)/Param.Fs,size(LW_RB_SHL,3))')'],[(LW_RB_SHL_Mean+LW_RB_SHL_SEM), flipud((LW_RB_SHL_Mean-LW_RB_SHL_SEM)')'],SHLColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax8,[linspace(-TimeStartW,size(LW_RB_N60,3)/Param.Fs,size(LW_RB_N60,3)), flipud(linspace(-TimeStartW,size(LW_RB_N60,3)/Param.Fs,size(LW_RB_N60,3))')'],[(LW_RB_N60_Mean+LW_RB_N60_SEM), flipud((LW_RB_N60_Mean-LW_RB_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax8,[linspace(-TimeStartW,size(LW_RB_N70,3)/Param.Fs,size(LW_RB_N70,3)), flipud(linspace(-TimeStartW,size(LW_RB_N70,3)/Param.Fs,size(LW_RB_N70,3))')'],[(LW_RB_N70_Mean+LW_RB_N70_SEM), flipud((LW_RB_N70_Mean-LW_RB_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')

fill(ax9,[linspace(-TimeStartW,size(GSW_Z,3)/Param.Fs,size(GSW_Z,3)), flipud(linspace(-TimeStartW,size(GSW_Z,3)/Param.Fs,size(GSW_Z,3))')'],[(GSW_Z_Mean+GSW_Z_SEM), flipud((GSW_Z_Mean-GSW_Z_SEM)')'],SColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax9,[linspace(-TimeStartW,size(SW_Z_Quiet,3)/Param.Fs,size(SW_Z_Quiet,3)), flipud(linspace(-TimeStartW,size(SW_Z_Quiet,3)/Param.Fs,size(SW_Z_Quiet,3))')'],[(SW_Z_Quiet_Mean+SW_Z_Quiet_SEM), flipud((SW_Z_Quiet_Mean-SW_Z_Quiet_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax9,[linspace(-TimeStartW,size(SW_Z_SHL,3)/Param.Fs,size(SW_Z_SHL,3)), flipud(linspace(-TimeStartW,size(SW_Z_SHL,3)/Param.Fs,size(SW_Z_SHL,3))')'],[(SW_Z_SHL_Mean+SW_Z_SHL_SEM), flipud((SW_Z_SHL_Mean-SW_Z_SHL_SEM)')'],SHLColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax9,[linspace(-TimeStartW,size(SW_Z_N60,3)/Param.Fs,size(SW_Z_N60,3)), flipud(linspace(-TimeStartW,size(SW_Z_N60,3)/Param.Fs,size(SW_Z_N60,3))')'],[(SW_Z_N60_Mean+SW_Z_N60_SEM), flipud((SW_Z_N60_Mean-SW_Z_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax9,[linspace(-TimeStartW,size(SW_Z_N70,3)/Param.Fs,size(SW_Z_N70,3)), flipud(linspace(-TimeStartW,size(SW_Z_N70,3)/Param.Fs,size(SW_Z_N70,3))')'],[(SW_Z_N70_Mean+SW_Z_N70_SEM), flipud((SW_Z_N70_Mean-SW_Z_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax10,[linspace(-TimeStartW,size(GLW_Z,3)/Param.Fs,size(GLW_Z,3)), flipud(linspace(-TimeStartW,size(GLW_Z,3)/Param.Fs,size(GLW_Z,3))')'],[(GLW_Z_Mean+GLW_Z_SEM), flipud((GLW_Z_Mean-GLW_Z_SEM)')'],LColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax10,[linspace(-TimeStartW,size(LW_Z_Quiet,3)/Param.Fs,size(LW_Z_Quiet,3)), flipud(linspace(-TimeStartW,size(LW_Z_Quiet,3)/Param.Fs,size(LW_Z_Quiet,3))')'],[(LW_Z_Quiet_Mean+LW_Z_Quiet_SEM), flipud((LW_Z_Quiet_Mean-LW_Z_Quiet_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax10,[linspace(-TimeStartW,size(LW_Z_SHL,3)/Param.Fs,size(LW_Z_SHL,3)), flipud(linspace(-TimeStartW,size(LW_Z_SHL,3)/Param.Fs,size(LW_Z_SHL,3))')'],[(LW_Z_SHL_Mean+LW_Z_SHL_SEM), flipud((LW_Z_SHL_Mean-LW_Z_SHL_SEM)')'],SHLColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax10,[linspace(-TimeStartW,size(LW_Z_N60,3)/Param.Fs,size(LW_Z_N60,3)), flipud(linspace(-TimeStartW,size(LW_Z_N60,3)/Param.Fs,size(LW_Z_N60,3))')'],[(LW_Z_N60_Mean+LW_Z_N60_SEM), flipud((LW_Z_N60_Mean-LW_Z_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax10,[linspace(-TimeStartW,size(LW_Z_N70,3)/Param.Fs,size(LW_Z_N70,3)), flipud(linspace(-TimeStartW,size(LW_Z_N70,3)/Param.Fs,size(LW_Z_N70,3))')'],[(LW_Z_N70_Mean+LW_Z_N70_SEM), flipud((LW_Z_N70_Mean-LW_Z_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')

xlabel([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10],'Time [s]')
ylabel([ax1 ax2],'Pupil diameter [mm]')
ylabel([ax3 ax4],'Pupil baseline difference [mm]')
ylabel([ax5 ax6],'Range-normalized pupil diameter [%]')
ylabel([ax7 ax8],'Normalized pupil baseline difference [mm]')
ylabel([ax9 ax10],'Pupil z-score difference [mm]')
xlim([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10],[-TimeStartW 3])
lgd2=legend(ax2,'Speaking','Listening','Quiet','SHL','N60','N70','Location','southeastoutside');
lgd2.Title.String = 'Types of windows:';
lgd4=legend(ax4,'Speaking','Listening','Quiet','SHL','N60','N70','Location','southeastoutside');
lgd4.Title.String = 'Types of windows:';
lgd6=legend(ax6,'Speaking','Listening','Quiet','SHL','N60','N70','Location','southeastoutside');
lgd6.Title.String = 'Types of windows:';
lgd8=legend(ax8,'Speaking','Listening','Quiet','SHL','N60','N70','Location','southeastoutside');
lgd8.Title.String = 'Types of windows:';
lgd10=legend(ax10,'Speaking','Listening','Quiet','SHL','N60','N70','Location','southeastoutside');
lgd10.Title.String = 'Types of windows:';
title(ax1,'Global Speaking-evoked')
title(ax2,'Global Listening-evoked')
title(ax3,'Global adaptive baselined Speaking-evoked')
title(ax4,'Global adaptive baselined Listening-evoked')
title(ax5,'Global range normalized Speaking-evoked')
title(ax6,'Global range normalized Listening-evoked')
title(ax7,'Global range normalized baselined Speaking-evoked')
title(ax8,'Global range normalized baselined Listening-evoked')
title(ax9,'Global z-score normalized Speaking-evoked')
title(ax10,'Global z-score normalized Listening-evoked')

% figure;
% plot(linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3)),squeeze(mean(GSW,[1 2],'omitnan')));
% hold on;
% plot(linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3)),squeeze(mean(mean(GSW,'omitnan'),'omitnan')))
% t_vec=linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3));
% for i=1:size(GSW,3)
%    xline(t_vec(i),':')
% end
% for i=1:size(GSW,1)
%     for j=1:size(GSW,2)
%         plot(linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3)),squeeze(GSW(i,j,:)),color=[0 0 0 0.1],linewidth=0.5)
%     end
% end

% figure;hold on;for i=1:size(GSW,1)
% plot(linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3)),squeeze(mean(GSW(i,:,:),'omitnan')),color=[0 0 0 0.1],linewidth=0.5)
% end;plot(linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3)),squeeze(mean(GSW,[1 2],'omitnan')),color=[1 0 0],linewidth=1);plot(linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3)),squeeze(mean(mean(GSW,'omitnan'),'omitnan')),color=[0 0 1],linewidth=1)
