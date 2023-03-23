%% PBCA-Thesis - Week 13 - Different types of data processing - Pupillometry
% Pathing
clear all; clc; close all;
BPath = strsplit(pwd,'PBCA-thesis');
addpath('tools')
addpath([BPath{1} 'Pupil-preprocessing-tools\tools']) % For preprocessing

[subDirs] = GetSubDirsFirstLevelOnly('data');
FileNames={'P1_Quiet_B1.mat','P1_Quiet_B2.mat','P1_SHL_B1.mat','P1_SHL_B2.mat','P1_Noise60_B1.mat','P1_Noise60_B2.mat','P1_Noise70_B1.mat','P1_Noise70_B2.mat','P2_Quiet_B1.mat','P2_Quiet_B2.mat','P2_SHL_B1.mat','P2_SHL_B2.mat','P2_Noise60_B1.mat','P2_Noise60_B2.mat','P2_Noise70_B1.mat','P2_Noise70_B2.mat'};
LoadUtt=load('data\utterances1110.mat');
LoadDelays=load('data\delays1110.mat');

% Parameters for processing
% Files and Utterances: different conditions
Param.Fs = 50; % Sampling frequency of pupil data
Param.RemoveBeforeAndAfter = [35 100]*1e-3; % Samples within the time range before and after NaNs will set NaNs as well.
Param.MinLengthNaNRepair = 5; % Drop values (i.e., change to NaN) before and after NaNs only for contiguous NaNs of at least __ samples. 
LPWinSize = 0.5; % [s]: Window size of hamming-window for low-pass filtering
LPWindow = hamming(round(LPWinSize*Param.Fs));
LPWindow = LPWindow/sum(LPWindow); % Hamming-window
AudFs = 48000; % [Hz], sampling frequency of the audio files
AdapBL = 0.1; % [s], Duration of baseline prior to event
TimeStart = 20; % [s], Start of trial (disregard data prior to 'TimeStart')
TimeMinWin = 0.5; % [s], Minimum time of a window
TimeInitialMerge = 0.3; % [s], Time threshold for merging windows initially
TimeMerge = 2; % [s], Time threshold for merging windows after rejecting small windows
RejectRatio = 0.4; % Rejection threshold based on the ratio of NaNs in data
RejectDelay = 0.5; % [s], Rejection threshold based on delay between timestamps and n-samples
TimeStartW = 0.5; % [s], time before Utt/Lis starts
TimeEndW = 0; % [s], time after Utt/Lis starts
x=1; % idx to store global values

% Colors
SColor = [53, 155, 67]./255;
LColor = [204, 36, 0]./255;
QuietColor = [204, 152, 0]./255;
SHLColor = [123, 31, 162]./255;
N60Color = [0, 196, 215]./255;
N70Color = [2, 36, 223]./255;

% Variables
NCols=3000; % Duration (samples) of each window
NRows=100; % Number of windows per trial
NLayers=16*12; % Number of trials

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

GSW_R = zeros(NLayers,NRows,NCols); % Global Speaking Windows Range normalized
GLW_R = zeros(NLayers,NRows,NCols); % Global Listening Windows Range normalized
SW_R_Quiet = zeros(NLayers,NRows,NCols); % Quiet Speaking Windows Range normalized
LW_R_Quiet = zeros(NLayers,NRows,NCols); % Quiet Listening Windows Range normalized
SW_R_SHL = zeros(NLayers,NRows,NCols); % SHL Speaking Windows Range normalized
LW_R_SHL = zeros(NLayers,NRows,NCols); % SHL Listening Windows Range normalized
SW_R_N60 = zeros(NLayers,NRows,NCols); % N60 Speaking Windows Range normalized
LW_R_N60 = zeros(NLayers,NRows,NCols); % N60 Listening Windows Range normalized
SW_R_N70 = zeros(NLayers,NRows,NCols); % N70 Speaking Windows Range normalized
LW_R_N70 = zeros(NLayers,NRows,NCols); % N70 Listening Windows Range normalized

GSW_RB = zeros(NLayers,NRows,NCols); % Global Speaking Windows Baselined Range normalized
GLW_RB = zeros(NLayers,NRows,NCols); % Global Listening Windows Baselined Range normalized
SW_RB_Quiet = zeros(NLayers,NRows,NCols); % Quiet Speaking Windows Baselined Range normalized
LW_RB_Quiet = zeros(NLayers,NRows,NCols); % Quiet Listening Windows Baselined Range normalized
SW_RB_SHL = zeros(NLayers,NRows,NCols); % SHL Speaking Windows Baselined Range normalized
LW_RB_SHL = zeros(NLayers,NRows,NCols); % SHL Listening Windows Baselined Range normalized
SW_RB_N60 = zeros(NLayers,NRows,NCols); % N60 Speaking Windows Baselined Range normalized
LW_RB_N60 = zeros(NLayers,NRows,NCols); % N60 Listening Windows Baselined Range normalized
SW_RB_N70 = zeros(NLayers,NRows,NCols); % N70 Speaking Windows Baselined Range normalized
LW_RB_N70 = zeros(NLayers,NRows,NCols); % N70 Listening Windows Baselined Range normalized

GSW_Z = zeros(NLayers,NRows,NCols); % Global Speaking Windows Z-score normalized
GLW_Z = zeros(NLayers,NRows,NCols); % Global Listening Windows Z-score normalized
SW_Z_Quiet = zeros(NLayers,NRows,NCols); % Quiet Speaking Windows Z-score normalized
LW_Z_Quiet = zeros(NLayers,NRows,NCols); % Quiet Listening Windows Z-score normalized
SW_Z_SHL = zeros(NLayers,NRows,NCols); % SHL Speaking Windows Z-score normalized
LW_Z_SHL = zeros(NLayers,NRows,NCols); % SHL Listening Windows Z-score normalized
SW_Z_N60 = zeros(NLayers,NRows,NCols); % N60 Speaking Windows Z-score normalized
LW_Z_N60 = zeros(NLayers,NRows,NCols); % N60 Listening Windows Z-score normalized
SW_Z_N70 = zeros(NLayers,NRows,NCols); % N70 Speaking Windows Z-score normalized
LW_Z_N70 = zeros(NLayers,NRows,NCols); % N70 Listening Windows Z-score normalized

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

for q=1:numel(subDirs)
    PairIn = q;
    PairFiles=dir(['data\Main',sprintf('%d',PairIn),'\*.mat']);
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
            disp(['Warning: File ',PairFiles(i).folder, '\', PairFiles(i).name, ' was rejected, too much delay (',sprintf('%0.2f',EyeAudDelay),'s).']);
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
            elseif RDiamRaw(1,s) < RThreshOut(1) || RDiamRaw(1,s) > RThreshOut(2)
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
        
        if DiamNaN/length(Diameter) >= RejectRatio
            disp(['Warning: File ',PairFiles(i).folder, '\', PairFiles(i).name, ' was rejected because it contains too many NaNs (',sprintf('%0.2f',100*DiamNaN/length(Diameter)),'%).'])
            continue
        end
        
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
            disp(['Warning: File ',PairFiles(i).folder, '\', PairFiles(i).name,' was rejected for not having an associated delay.']);
            continue
        end
        
        if (SDelayRaw <= 0 || SDelayRaw >= 1 || LDelayRaw <= 0 || LDelayRaw >= 1)
            disp(['Warning: File ', PairFiles(i).folder, '\', PairFiles(i).name,' was rejected because the delay between Tobii and Utterances is out of proportions (Speaking: ',sprintf('%0.3f',SDelayRaw(2)),' [s], Listening: ',sprintf('%0.3f',LDelayRaw(2)),' [s]).']);
            continue
        end
        
        if isempty(SpeakRaw) && isempty(ListenRaw)
            disp(['Warning: File ',PairFiles(i).folder, '\', PairFiles(i).name,' was rejected for not having associated Utterance windows.']);
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
        SpeakD = SpeakMI(SpeakMI(:,1)>TimeMinWin,:);
        ListenD = ListenMI(ListenMI(:,1)>TimeMinWin,:);
        
        % Merge again if duration between windows <= TimeMerge (2 s)
        SpeakM = merge_windows(SpeakD, Param.Fs, TimeMerge);
        ListenM = merge_windows(ListenD, Param.Fs, TimeMerge);
        
        % Discard windows if duration is < 2*TimeMinWin (1 s)
        Speak = SpeakM(SpeakM(:,1)>2*TimeMinWin,:);
        Listen = ListenM(ListenM(:,1)>2*TimeMinWin,:);
        
        % Time-locked indexes (based on Start or End of events)
        WSpeakIdx=[Speak(:,2)-TimeStartW*Param.Fs,Speak(:,2),Speak(:,2)+TimeEndW*Param.Fs];
        WListenIdx=[Listen(:,2)-TimeStartW*Param.Fs,Listen(:,2),Listen(:,2)+TimeEndW*Param.Fs];
%         EWSpeakIdx=[Speak(:,3)-TimeStartW*Param.Fs,Speak(:,3),Speak(:,3)+TimeEndW*Param.Fs];
%         EWListenIdx=[Listen(:,3)-TimeStartW*Param.Fs,Listen(:,3),Listen(:,3)+TimeEndW*Param.Fs];
        
        t_Diam = linspace(0,length(Diameter)./Param.Fs,length(Diameter));
        
        % Storing Speaking/Listening by conditions
        if contains(cell2mat(FileNames(i)),'Quiet')
            for j=1:size(Speak,1)
                SW_Quiet(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3)-1);NaN*ones(1,length(SW_Quiet)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SW_B_Quiet(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3)-1)-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_B_Quiet)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SW_R_Quiet(x,j,:)=[(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)-min(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)))/(max(Diameter(WSpeakIdx(j,1):Speak(j,3)-1))-min(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)));NaN*ones(1,length(SW_R_Quiet)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SW_RB_Quiet(x,j,:)=[(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)))-min(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)))/(max(Diameter(WSpeakIdx(j,1):Speak(j,3)-1))-min(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)));NaN*ones(1,length(SW_RB_Quiet)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SW_Z_Quiet(x,j,:)=[(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)-mean(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)))/std(mean(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)));NaN*ones(1,length(SW_Z_Quiet)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SDur_Quiet(i,j) = Speak(j,1);
            end
            for j=1:size(Listen,1)
                LW_Quiet(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3)-1);NaN*ones(1,length(LW_Quiet)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LW_B_Quiet(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3)-1)-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_B_Quiet)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LW_R_Quiet(x,j,:)=[(Diameter(WListenIdx(j,1):Listen(j,3)-1)-min(Diameter(WListenIdx(j,1):Listen(j,3)-1)))/(max(Diameter(WListenIdx(j,1):Listen(j,3)-1))-min(Diameter(WListenIdx(j,1):Listen(j,3)-1)));NaN*ones(1,length(LW_R_Quiet)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LW_RB_Quiet(x,j,:)=[(Diameter(WListenIdx(j,1):Listen(j,3)-1)-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)))-min(Diameter(WListenIdx(j,1):Listen(j,3)-1)))/(max(Diameter(WListenIdx(j,1):Listen(j,3)-1))-min(Diameter(WListenIdx(j,1):Listen(j,3)-1)));NaN*ones(1,length(LW_RB_Quiet)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LW_Z_Quiet(x,j,:)=[(Diameter(WListenIdx(j,1):Listen(j,3)-1)-mean(Diameter(WListenIdx(j,1):Listen(j,3)-1)))/std(mean(Diameter(WListenIdx(j,1):Listen(j,3)-1)));NaN*ones(1,length(LW_Z_Quiet)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LDur_Quiet(i,j) = Listen(j,1);
            end
        elseif contains(cell2mat(FileNames(i)),'SHL')
            for j=1:size(Speak,1)
                SW_SHL(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3)-1);NaN*ones(1,length(SW_SHL)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SW_B_SHL(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3)-1)-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_B_SHL)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SW_R_SHL(x,j,:)=[(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)-min(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)))/(max(Diameter(WSpeakIdx(j,1):Speak(j,3)-1))-min(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)));NaN*ones(1,length(SW_R_SHL)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SW_RB_SHL(x,j,:)=[(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)))-min(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)))/(max(Diameter(WSpeakIdx(j,1):Speak(j,3)-1))-min(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)));NaN*ones(1,length(SW_RB_SHL)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SW_Z_SHL(x,j,:)=[(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)-mean(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)))/std(mean(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)));NaN*ones(1,length(SW_Z_SHL)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SDur_SHL(i,j) = Speak(j,1);
            end
            for j=1:size(Listen,1)
                LW_SHL(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3)-1);NaN*ones(1,length(LW_SHL)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LW_B_SHL(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3)-1)-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_B_SHL)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LW_R_SHL(x,j,:)=[(Diameter(WListenIdx(j,1):Listen(j,3)-1)-min(Diameter(WListenIdx(j,1):Listen(j,3)-1)))/(max(Diameter(WListenIdx(j,1):Listen(j,3)-1))-min(Diameter(WListenIdx(j,1):Listen(j,3)-1)));NaN*ones(1,length(LW_R_SHL)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LW_RB_SHL(x,j,:)=[(Diameter(WListenIdx(j,1):Listen(j,3)-1)-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)))-min(Diameter(WListenIdx(j,1):Listen(j,3)-1)))/(max(Diameter(WListenIdx(j,1):Listen(j,3)-1))-min(Diameter(WListenIdx(j,1):Listen(j,3)-1)));NaN*ones(1,length(LW_RB_SHL)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LW_Z_SHL(x,j,:)=[(Diameter(WListenIdx(j,1):Listen(j,3)-1)-mean(Diameter(WListenIdx(j,1):Listen(j,3)-1)))/std(mean(Diameter(WListenIdx(j,1):Listen(j,3)-1)));NaN*ones(1,length(LW_Z_SHL)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LDur_SHL(i,j) = Listen(j,1);
            end
        elseif contains(cell2mat(FileNames(i)),'Noise60')
            for j=1:size(Speak,1)
                SW_N60(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3)-1);NaN*ones(1,length(SW_N60)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SW_B_N60(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3)-1)-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_B_N60)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SW_R_N60(x,j,:)=[(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)-min(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)))/(max(Diameter(WSpeakIdx(j,1):Speak(j,3)-1))-min(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)));NaN*ones(1,length(SW_R_N60)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SW_RB_N60(x,j,:)=[(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)))-min(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)))/(max(Diameter(WSpeakIdx(j,1):Speak(j,3)-1))-min(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)));NaN*ones(1,length(SW_RB_N60)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SW_Z_N60(x,j,:)=[(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)-mean(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)))/std(mean(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)));NaN*ones(1,length(SW_Z_N60)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SDur_N60(i,j) = Speak(j,1);
            end
            for j=1:size(Listen,1)
                LW_N60(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3)-1);NaN*ones(1,length(LW_N60)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LW_B_N60(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3)-1)-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_B_N60)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LW_R_N60(x,j,:)=[(Diameter(WListenIdx(j,1):Listen(j,3)-1)-min(Diameter(WListenIdx(j,1):Listen(j,3)-1)))/(max(Diameter(WListenIdx(j,1):Listen(j,3)-1))-min(Diameter(WListenIdx(j,1):Listen(j,3)-1)));NaN*ones(1,length(LW_R_N60)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LW_RB_N60(x,j,:)=[(Diameter(WListenIdx(j,1):Listen(j,3)-1)-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)))-min(Diameter(WListenIdx(j,1):Listen(j,3)-1)))/(max(Diameter(WListenIdx(j,1):Listen(j,3)-1))-min(Diameter(WListenIdx(j,1):Listen(j,3)-1)));NaN*ones(1,length(LW_RB_N60)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LW_Z_N60(x,j,:)=[(Diameter(WListenIdx(j,1):Listen(j,3)-1)-mean(Diameter(WListenIdx(j,1):Listen(j,3)-1)))/std(mean(Diameter(WListenIdx(j,1):Listen(j,3)-1)));NaN*ones(1,length(LW_Z_N60)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LDur_N60(i,j) = Listen(j,1);
            end
        elseif contains(cell2mat(FileNames(i)),'Noise70')
            for j=1:size(Speak,1)
                SW_N70(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3)-1);NaN*ones(1,length(SW_N70)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SW_B_N70(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3)-1)-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_B_N70)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SW_R_N70(x,j,:)=[(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)-min(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)))/(max(Diameter(WSpeakIdx(j,1):Speak(j,3)-1))-min(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)));NaN*ones(1,length(SW_R_N70)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SW_RB_N70(x,j,:)=[(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)))-min(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)))/(max(Diameter(WSpeakIdx(j,1):Speak(j,3)-1))-min(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)));NaN*ones(1,length(SW_RB_N70)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SW_Z_N70(x,j,:)=[(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)-mean(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)))/std(mean(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)));NaN*ones(1,length(SW_Z_N70)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SDur_N70(i,j) = Speak(j,1);
            end
            for j=1:size(Listen,1)
                LW_N70(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3)-1);NaN*ones(1,length(LW_N70)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LW_B_N70(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3)-1)-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_B_N70)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LW_R_N70(x,j,:)=[(Diameter(WListenIdx(j,1):Listen(j,3)-1)-min(Diameter(WListenIdx(j,1):Listen(j,3)-1)))/(max(Diameter(WListenIdx(j,1):Listen(j,3)-1))-min(Diameter(WListenIdx(j,1):Listen(j,3)-1)));NaN*ones(1,length(LW_R_N70)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LW_RB_N70(x,j,:)=[(Diameter(WListenIdx(j,1):Listen(j,3)-1)-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)))-min(Diameter(WListenIdx(j,1):Listen(j,3)-1)))/(max(Diameter(WListenIdx(j,1):Listen(j,3)-1))-min(Diameter(WListenIdx(j,1):Listen(j,3)-1)));NaN*ones(1,length(LW_RB_N70)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LW_Z_N70(x,j,:)=[(Diameter(WListenIdx(j,1):Listen(j,3)-1)-mean(Diameter(WListenIdx(j,1):Listen(j,3)-1)))/std(mean(Diameter(WListenIdx(j,1):Listen(j,3)-1)));NaN*ones(1,length(LW_Z_N70)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LDur_N70(x,j) = Listen(j,1);
            end
        end
        
        % Storing Global Speaking/Listening pupil sizes
        for j=1:size(Speak,1)
            % Add nan-padding when necessary
            GSW(x,j,:) = [Diameter(WSpeakIdx(j,1):Speak(j,3)-1);NaN*ones(1,length(GSW)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
            GSW_B(x,j,:) = [Diameter(WSpeakIdx(j,1):Speak(j,3)-1)-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(GSW_B)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
            GSW_R(x,j,:)=[(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)-min(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)))/(max(Diameter(WSpeakIdx(j,1):Speak(j,3)-1))-min(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)));NaN*ones(1,length(GSW_R)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
            GSW_RB(x,j,:)=[(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)))-min(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)))/(max(Diameter(WSpeakIdx(j,1):Speak(j,3)-1))-min(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)));NaN*ones(1,length(GSW_RB)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
            GSW_Z(x,j,:)=[(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)-mean(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)))/std(mean(Diameter(WSpeakIdx(j,1):Speak(j,3)-1)));NaN*ones(1,length(GSW_Z)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
            GSDur(x,j) = Speak(j,1);
        end
        
        for j=1:size(Listen,1)
            % Add nan-padding when necessary
            GLW(x,j,:) = [Diameter(WListenIdx(j,1):Listen(j,3)-1);NaN*ones(1,length(GLW)-length(WListenIdx(j,1):Listen(j,3)-1))'];
            GLW_B(x,j,:) = [Diameter(WListenIdx(j,1):Listen(j,3)-1)-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(GLW_B)-length(WListenIdx(j,1):Listen(j,3)-1))'];
            GLW_R(x,j,:)=[(Diameter(WListenIdx(j,1):Listen(j,3)-1)-min(Diameter(WListenIdx(j,1):Listen(j,3)-1)))/(max(Diameter(WListenIdx(j,1):Listen(j,3)-1))-min(Diameter(WListenIdx(j,1):Listen(j,3)-1)));NaN*ones(1,length(GLW_R)-length(WListenIdx(j,1):Listen(j,3)-1))'];
            GLW_RB(x,j,:)=[(Diameter(WListenIdx(j,1):Listen(j,3)-1)-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)))-min(Diameter(WListenIdx(j,1):Listen(j,3)-1)))/(max(Diameter(WListenIdx(j,1):Listen(j,3)-1))-min(Diameter(WListenIdx(j,1):Listen(j,3)-1)));NaN*ones(1,length(GLW_RB)-length(WListenIdx(j,1):Listen(j,3)-1))'];
            GLW_Z(x,j,:)=[(Diameter(WListenIdx(j,1):Listen(j,3)-1)-mean(Diameter(WListenIdx(j,1):Listen(j,3)-1)))/std(mean(Diameter(WListenIdx(j,1):Listen(j,3)-1)));NaN*ones(1,length(GLW_Z)-length(WListenIdx(j,1):Listen(j,3)-1))'];
            GLDur(x,j) = Listen(j,1);
        end
        
        % Increase index of num of files used
        x=x+1;
        
    end
end
