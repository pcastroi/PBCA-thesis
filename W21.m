%% PBCA-Thesis - Week 21 - Pupil Features in 2 phases for AMEND I - Copy of W4.m
% Pathing
clear all; clc; close all;
BPath = strsplit(pwd,'PBCA-thesis');
addpath('tools')
addpath([BPath{1} 'Pupil-preprocessing-tools\tools']) % For preprocessing

% Files and Utterances: different conditions
Param.Fs = 50; % Sampling frequency of pupil data
Param.RemoveBeforeAndAfter = [35 100]*1e-3; % Samples within the time range before and after NaNs will set NaNs as well.
Param.MinLengthNaNRepair = 5; % Drop values (i.e., change to NaN) before and after NaNs only for contiguous NaNs of at least __ samples. 
LPWinSize = 0.5; % [s]: Window size of hamming-window for low-pass filtering
LPWindow = hamming(round(LPWinSize*Param.Fs));
LPWindow = LPWindow/sum(LPWindow); % Hamming-window
AudFs = 48000;
WD = [-0.5,4]; % [s], Duration limits of phase 1 & 2
BLPeriod = [0,20]; % [s]
TimeMinWin = 0.5; % [s], Minimum time of a window
TimeInitialMerge = 0.3; % [s], Time threshold for merging windows initially
TimeMerge = 2; % [s], Time threshold for merging windows after rejecting small windows
RejectRatio = 0.4; % Rejection threshold based on the ratio of NaNs in data
RejectDelay = 0.5; % [s], Rejection threshold based on delay between timestamps and n-samples
Med_S=0; % [mm], Median of Speaking windows

% Colors
SColor = [53, 155, 67]./255;
LColor = [204, 36, 0]./255;
QuietColor = [204, 152, 0]./255;
SHLColor = [123, 31, 162]./255;
N60Color = [0, 196, 215]./255;
N70Color = [2, 36, 223]./255;

% Preallocate features - Phase 1, Phase 2
F1_S_Q=zeros(200,2); % Feature 1 (Quiet): Mean pupil size (Speaking)
F1_L_Q=F1_S_Q; % Feature 1 (Quiet): Mean pupil size (Listening)
F2_S_Q=F1_S_Q; % Feature 2 (Quiet): Mean slope (Speaking)
F2_L_Q=F1_S_Q; % Feature 2 (Quiet): Mean slope (Listening)
F3_S_Q=F1_S_Q; % Feature 3 (Quiet): Mean Peak Pupil Size (Speaking)
F3_L_Q=F1_S_Q; % Feature 3 (Quiet): Mean Peak Pupil Size (Listening)
F1_S_SHL=F1_S_Q; % Feature 1 (SHL): Mean pupil size (Speaking)
F1_L_SHL=F1_S_Q; % Feature 1 (SHL): Mean pupil size (Listening)
F2_S_SHL=F1_S_Q; % Feature 2 (SHL): Mean slope (Speaking)
F2_L_SHL=F1_S_Q; % Feature 2 (SHL): Mean slope (Listening)
F3_S_SHL=F1_S_Q; % Feature 3 (SHL): Mean Peak Pupil Size (Speaking)
F3_L_SHL=F1_S_Q; % Feature 3 (SHL): Mean Peak Pupil Size (Listening)
F1_S_N60=F1_S_Q; % Feature 1 (N60): Mean pupil size (Speaking)
F1_L_N60=F1_S_Q; % Feature 1 (N60): Mean pupil size (Listening)
F2_S_N60=F1_S_Q; % Feature 2 (N60): Mean slope (Speaking)
F2_L_N60=F1_S_Q; % Feature 2 (N60): Mean slope (Listening)
F3_S_N60=F1_S_Q; % Feature 3 (N60): Mean Peak Pupil Size (Speaking)
F3_L_N60=F1_S_Q; % Feature 3 (N60): Mean Peak Pupil Size (Listening)
F1_S_N70=F1_S_Q; % Feature 1 (N70): Mean pupil size (Speaking)
F1_L_N70=F1_S_Q; % Feature 1 (N70): Mean pupil size (Listening)
F2_S_N70=F1_S_Q; % Feature 2 (N70): Mean slope (Speaking)
F2_L_N70=F1_S_Q; % Feature 2 (N70): Mean slope (Listening)
F3_S_N70=F1_S_Q; % Feature 3 (N70): Mean Peak Pupil Size (Speaking)
F3_L_N70=F1_S_Q; % Feature 3 (N70): Mean Peak Pupil Size (Listening)

x=1;

[subDirs] = GetSubDirsFirstLevelOnly('data\AMEND_I');
FileNames={'P1_Quiet_B1.mat','P1_Quiet_B2.mat','P1_SHL_B1.mat','P1_SHL_B2.mat','P1_Noise60_B1.mat','P1_Noise60_B2.mat','P1_Noise70_B1.mat','P1_Noise70_B2.mat','P2_Quiet_B1.mat','P2_Quiet_B2.mat','P2_SHL_B1.mat','P2_SHL_B2.mat','P2_Noise60_B1.mat','P2_Noise60_B2.mat','P2_Noise70_B1.mat','P2_Noise70_B2.mat'};
LoadDelays=load('data\AMEND_I\delays1110.mat');

for q=1:numel(subDirs)
    PairIn = q;
    PairFiles=dir(['data\AMEND_I\Main',sprintf('%d',PairIn),'\*.mat']);
    PairUtt=load('data\AMEND_I\utterances1110.mat');
    PairUtt=PairUtt.Utterances(PairIn,:);
    PairDelay=LoadDelays.TobAudDelay(PairIn,:);
    
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
            continue
        end
        
        LDiamRaw = [alldata_mat.diameterLeft];
        RDiamRaw = [alldata_mat.diameterRight];
        
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
        
        if (SDelayRaw <= 0 | SDelayRaw >= 1 | LDelayRaw <= 0 | LDelayRaw >= 1)
            disp(['Warning: File ', PairFiles(i).folder, '\', PairFiles(i).name,' was rejected because the delay between Tobii and Utterances is out of proportions (Speaking: ',sprintf('%0.3f',SDelayRaw(2)),' [s], Listening: ',sprintf('%0.3f',LDelayRaw(2)),' [s]).']);
            continue
        end
        
        if isempty(SpeakRaw) && isempty(ListenRaw)
            disp(['Warning: File ',PairFiles(i).folder, '\', PairFiles(i).name,' was rejected for not having associated Utterance windows.']);
            continue
        end
        
        SpeakRaw(:,2:3)=round((SpeakRaw(:,2:3)*binResUtt+BLPeriod(2))*Param.Fs+SDelayRaw(1)/2);
        ListenRaw(:,2:3)=round((ListenRaw(:,2:3)*binResUtt+BLPeriod(2))*Param.Fs+LDelayRaw(1)/2);
        
        % SAME PROCESSING AS IN W1.m
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
        
        % Average median, mean of duration of Speaking/Listening
        Med_S=Med_S+median(Speak(:,1));
        
        t_Diam = linspace(0,length(Diameter)./Param.Fs,length(Diameter));
        
        if isempty(Speak)  % Could be empty
            Speak=Listen(1,:)+1;
        elseif isempty(Listen)
            Listen=Speak(1,:)+1;
        end
        
        % Imposed delay makes windows be larger than possible
        if Speak(end,3) > length(Diameter) 
            Speak(end,3) = length(Diameter);
        elseif Listen(end,3) > length(Diameter)
            Listen(end,3) = length(Diameter);
        end
        
        % Extract features for different conditions in different phases
        % Phase 1: from 0 to Pupil Peak
        % Phase 2: from Pupil Peak to 4 s
        temp=zeros(100,2);
        temp2=zeros(100,4);
        temp3=temp;
        if contains(PairFiles(i).name,'Quiet')
            for j=1:size(Speak,1)
                if Speak(j,2)+WD(2)*Param.Fs >= Speak(j,3)
                    PeakLim = Speak(j,3);
                else
                    PeakLim = Speak(j,2)+WD(2)*Param.Fs;
                end
                [~,Peak_Idx] = max(Diameter(Speak(j,2):PeakLim));
                if Peak_Idx == 1
                    Peak_Idx = 2;
                elseif Peak_Idx == length(Diameter(Speak(j,2):PeakLim)) || Peak_Idx == length(Diameter(Speak(j,2):PeakLim))-1
                    Peak_Idx = length(Diameter(Speak(j,2):PeakLim))-2;
                end
                Peak = Peak_Idx+Speak(j,2);
                temp(j,:)=[mean(Diameter(Speak(j,2):Peak)), mean(Diameter(Peak:Speak(j,3)))];
                temp2(j,:)=[polyfit(t_Diam(Speak(j,2):Peak),Diameter(Speak(j,2):Peak),1), polyfit(t_Diam(Peak:Speak(j,3)),Diameter(Peak:Speak(j,3)),1)];
                temp3(j,:)=[max(Diameter(Speak(j,2):Peak)), max(Diameter(Peak:Speak(j,3)))];
            end
            F1_S_Q(x,:)=[mean(nonzeros(temp(:,1))),mean(nonzeros(temp(:,2)))];
            F2_S_Q(x,:)=[mean(nonzeros(temp2(:,2))),mean(nonzeros(temp2(:,4)))];
            F3_S_Q(x,:)=[mean(nonzeros(temp3(:,1))),mean(nonzeros(temp3(:,2)))];
            for j=1:size(Listen,1)
                if Listen(j,2)+WD(2)*Param.Fs >= Listen(j,3)
                    PeakLim = Listen(j,3);
                else
                    PeakLim = Listen(j,2)+WD(2)*Param.Fs;
                end
                [~,Peak_Idx] = max(Diameter(Listen(j,2):PeakLim));
                if Peak_Idx == 1
                    Peak_Idx = 2;
                elseif Peak_Idx == length(Diameter(Listen(j,2):PeakLim)) || Peak_Idx == length(Diameter(Listen(j,2):PeakLim))-1
                    Peak_Idx = length(Diameter(Listen(j,2):PeakLim))-2;
                end
                Peak = Peak_Idx+Listen(j,2);
                temp(j,:)=[mean(Diameter(Listen(j,2):Peak)), mean(Diameter(Peak:Listen(j,3)))];
                temp2(j,:)=[polyfit(t_Diam(Listen(j,2):Peak),Diameter(Listen(j,2):Peak),1), polyfit(t_Diam(Peak:Listen(j,3)),Diameter(Peak:Listen(j,3)),1)];
                temp3(j,:)=[max(Diameter(Listen(j,2):Peak)), max(Diameter(Peak:Listen(j,3)))];
            end
            F1_L_Q(x,:)=[mean(nonzeros(temp(:,1))),mean(nonzeros(temp(:,2)))];
            F2_L_Q(x,:)=[mean(nonzeros(temp2(:,2))),mean(nonzeros(temp2(:,4)))];
            F3_L_Q(x,:)=[mean(nonzeros(temp3(:,1))),mean(nonzeros(temp3(:,2)))];
        elseif contains(PairFiles(i).name,'SHL')
            for j=1:size(Speak,1)
                if Speak(j,2)+WD(2)*Param.Fs >= Speak(j,3)
                    PeakLim = Speak(j,3);
                else
                    PeakLim = Speak(j,2)+WD(2)*Param.Fs;
                end
                [~,Peak_Idx] = max(Diameter(Speak(j,2):PeakLim));
                if Peak_Idx == 1
                    Peak_Idx = 2;
                elseif Peak_Idx == length(Diameter(Speak(j,2):PeakLim)) || Peak_Idx == length(Diameter(Speak(j,2):PeakLim))-1
                    Peak_Idx = length(Diameter(Speak(j,2):PeakLim))-2;
                end
                Peak = Peak_Idx+Speak(j,2);
                temp(j,:)=[mean(Diameter(Speak(j,2):Peak)), mean(Diameter(Peak:Speak(j,3)))];
                temp2(j,:)=[polyfit(t_Diam(Speak(j,2):Peak),Diameter(Speak(j,2):Peak),1), polyfit(t_Diam(Peak:Speak(j,3)),Diameter(Peak:Speak(j,3)),1)];
                temp3(j,:)=[max(Diameter(Speak(j,2):Peak)), max(Diameter(Peak:Speak(j,3)))];
            end
            F1_S_SHL(x,:)=mean(nonzeros(temp));
            F2_S_SHL(x,:)=mean(nonzeros(temp2(:,2)));
            F3_S_SHL(x,:)=mean(nonzeros(temp3));
            for j=1:size(Listen,1)
                if Listen(j,2)+WD(2)*Param.Fs >= Listen(j,3)
                    PeakLim = Listen(j,3);
                else
                    PeakLim = Listen(j,2)+WD(2)*Param.Fs;
                end
                [~,Peak_Idx] = max(Diameter(Listen(j,2):PeakLim));
                if Peak_Idx == 1
                    Peak_Idx = 2;
                elseif Peak_Idx == length(Diameter(Listen(j,2):PeakLim)) || Peak_Idx == length(Diameter(Listen(j,2):PeakLim))-1
                    Peak_Idx = length(Diameter(Listen(j,2):PeakLim))-2;
                end
                Peak = Peak_Idx+Listen(j,2);
                temp(j,:)=[mean(Diameter(Listen(j,2):Peak)), mean(Diameter(Peak:Listen(j,3)))];
                temp2(j,:)=[polyfit(t_Diam(Listen(j,2):Peak),Diameter(Listen(j,2):Peak),1), polyfit(t_Diam(Peak:Listen(j,3)),Diameter(Peak:Listen(j,3)),1)];
                temp3(j,:)=[max(Diameter(Listen(j,2):Peak)), max(Diameter(Peak:Listen(j,3)))];
            end
            F1_L_SHL(x,:)=mean(nonzeros(temp));
            F2_L_SHL(x,:)=mean(nonzeros(temp2(:,2)));
            F3_L_SHL(x,:)=mean(nonzeros(temp3));
        elseif contains(PairFiles(i).name,'Noise60')
            for j=1:size(Speak,1)
                if Speak(j,2)+WD(2)*Param.Fs >= Speak(j,3)
                    PeakLim = Speak(j,3);
                else
                    PeakLim = Speak(j,2)+WD(2)*Param.Fs;
                end
                [~,Peak_Idx] = max(Diameter(Speak(j,2):PeakLim));
                if Peak_Idx == 1
                    Peak_Idx = 2;
                elseif Peak_Idx == length(Diameter(Speak(j,2):PeakLim)) || Peak_Idx == length(Diameter(Speak(j,2):PeakLim))-1
                    Peak_Idx = length(Diameter(Speak(j,2):PeakLim))-2;
                end
                Peak = Peak_Idx+Speak(j,2);
                temp(j,:)=[mean(Diameter(Speak(j,2):Peak)), mean(Diameter(Peak:Speak(j,3)))];
                temp2(j,:)=[polyfit(t_Diam(Speak(j,2):Peak),Diameter(Speak(j,2):Peak),1), polyfit(t_Diam(Peak:Speak(j,3)),Diameter(Peak:Speak(j,3)),1)];
                temp3(j,:)=[max(Diameter(Speak(j,2):Peak)), max(Diameter(Peak:Speak(j,3)))];
            end
            F1_S_N60(x,:)=mean(nonzeros(temp));
            F2_S_N60(x,:)=mean(nonzeros(temp2(:,2)));
            F3_S_N60(x,:)=mean(nonzeros(temp3));
            for j=1:size(Listen,1)
                if Listen(j,2)+WD(2)*Param.Fs >= Listen(j,3)
                    PeakLim = Listen(j,3);
                else
                    PeakLim = Listen(j,2)+WD(2)*Param.Fs;
                end
                [~,Peak_Idx] = max(Diameter(Listen(j,2):PeakLim));
                if Peak_Idx == 1
                    Peak_Idx = 2;
                elseif Peak_Idx == length(Diameter(Listen(j,2):PeakLim)) || Peak_Idx == length(Diameter(Listen(j,2):PeakLim))-1
                    Peak_Idx = length(Diameter(Listen(j,2):PeakLim))-2;
                end
                Peak = Peak_Idx+Listen(j,2);
                temp(j,:)=[mean(Diameter(Listen(j,2):Peak)), mean(Diameter(Peak:Listen(j,3)))];
                temp2(j,:)=[polyfit(t_Diam(Listen(j,2):Peak),Diameter(Listen(j,2):Peak),1), polyfit(t_Diam(Peak:Listen(j,3)),Diameter(Peak:Listen(j,3)),1)];
                temp3(j,:)=[max(Diameter(Listen(j,2):Peak)), max(Diameter(Peak:Listen(j,3)))];
            end
            F1_L_N60(x,:)=mean(nonzeros(temp));
            F2_L_N60(x,:)=mean(nonzeros(temp2(:,2)));
            F3_L_N60(x,:)=mean(nonzeros(temp3));
        elseif contains(PairFiles(i).name,'Noise70')
            for j=1:size(Speak,1)
                if Speak(j,2)+WD(2)*Param.Fs >= Speak(j,3)
                    PeakLim = Speak(j,3);
                else
                    PeakLim = Speak(j,2)+WD(2)*Param.Fs;
                end
                [~,Peak_Idx] = max(Diameter(Speak(j,2):PeakLim));
                if Peak_Idx == 1
                    Peak_Idx = 2;
                elseif Peak_Idx == length(Diameter(Speak(j,2):PeakLim)) || Peak_Idx == length(Diameter(Speak(j,2):PeakLim))-1 || Peak_Idx == length(Diameter(Speak(j,2):PeakLim)) || Peak_Idx == length(Diameter(Speak(j,2):PeakLim))-1-1-1
                    Peak_Idx = length(Diameter(Speak(j,2):PeakLim))-2;
                end
                Peak = Peak_Idx+Speak(j,2);
                temp(j,:)=[mean(Diameter(Speak(j,2):Peak)), mean(Diameter(Peak:Speak(j,3)))];
                temp2(j,:)=[polyfit(t_Diam(Speak(j,2):Peak),Diameter(Speak(j,2):Peak),1), polyfit(t_Diam(Peak:Speak(j,3)),Diameter(Peak:Speak(j,3)),1)];
                temp3(j,:)=[max(Diameter(Speak(j,2):Peak)), max(Diameter(Peak:Speak(j,3)))];
            end
            F1_S_N70(x,:)=mean(nonzeros(temp));
            F2_S_N70(x,:)=mean(nonzeros(temp2(:,2)));
            F3_S_N70(x,:)=mean(nonzeros(temp3));
            for j=1:size(Listen,1)
                if Listen(j,2)+WD(2)*Param.Fs >= Listen(j,3)
                    PeakLim = Listen(j,3);
                else
                    PeakLim = Listen(j,2)+WD(2)*Param.Fs;
                end
                [~,Peak_Idx] = max(Diameter(Listen(j,2):PeakLim));
                if Peak_Idx == 1
                    Peak_Idx = 2;
                elseif Peak_Idx == length(Diameter(Listen(j,2):PeakLim)) || Peak_Idx == length(Diameter(Listen(j,2):PeakLim))-1
                    Peak_Idx = length(Diameter(Listen(j,2):PeakLim))-2;
                end
                Peak = Peak_Idx+Listen(j,2);
                temp(j,:)=[mean(Diameter(Listen(j,2):Peak)), mean(Diameter(Peak:Listen(j,3)))];
                temp2(j,:)=[polyfit(t_Diam(Listen(j,2):Peak),Diameter(Listen(j,2):Peak),1), polyfit(t_Diam(Peak:Listen(j,3)),Diameter(Peak:Listen(j,3)),1)];
                temp3(j,:)=[max(Diameter(Listen(j,2):Peak)), max(Diameter(Peak:Listen(j,3)))];
            end
            F1_L_N70(x,:)=mean(nonzeros(temp));
            F2_L_N70(x,:)=mean(nonzeros(temp2(:,2)));
            F3_L_N70(x,:)=mean(nonzeros(temp3));
        end
        x=x+1;
    end
end
%% Plots
figure;tiledlayout(1,3);
ax1 = nexttile;
ax2 = nexttile;
ax3 = nexttile;
hold([ax1 ax2 ax3],'on')
title(ax1,'Feature 1: Mean Pupil Size')
title(ax2,'Feature 2: Mean Slope')
title(ax3,'Feature 3: Mean Peak Pupil Size')
xlabel([ax1 ax2 ax3],'Conditions')
ylabel([ax1 ax3],'Pupil diameter [mm]')
ylabel(ax2,'Slope [mm/s]')
grid([ax1 ax2 ax3],'on')

% Plot features
data13=[mean(nonzeros(F1_S_Q)),mean(nonzeros(F1_L_Q));mean(nonzeros(F1_S_SHL)),mean(nonzeros(F1_L_SHL));mean(nonzeros(F1_S_N60)),mean(nonzeros(F1_L_N60));mean(nonzeros(F1_S_N70)),mean(nonzeros(F1_L_N70))];
b13=bar(ax1,data13,'grouped');
[ngroups,nbars] = size(data13);
x13 = nan(nbars, ngroups);
for e = 1:nbars
    x13(e,:) = b13(e).XEndPoints;
end
errorbar(ax1,x13',data13,[2*std(nonzeros(F1_S_Q))/sqrt(numel(nonzeros(F1_S_Q))),2*std(nonzeros(F1_L_Q))/sqrt(numel(nonzeros(F1_L_Q)));2*std(nonzeros(F1_S_SHL))/sqrt(numel(nonzeros(F1_S_SHL))),2*std(nonzeros(F1_L_SHL))/sqrt(numel(nonzeros(F1_L_SHL)));2*std(nonzeros(F1_S_N60))/sqrt(numel(nonzeros(F1_S_N60))),2*std(nonzeros(F1_L_N60))/sqrt(numel(nonzeros(F1_L_N60)));2*std(nonzeros(F1_S_N70))/sqrt(numel(nonzeros(F1_S_N70))),2*std(nonzeros(F1_L_N70))/sqrt(numel(nonzeros(F1_L_N70)))],'k','linestyle','none','handlevisibility' ,'off')
data14=[mean(nonzeros(F2_S_Q)),mean(nonzeros(F2_L_Q));mean(nonzeros(F2_S_SHL)),mean(nonzeros(F2_L_SHL));mean(nonzeros(F2_S_N60)),mean(nonzeros(F2_L_N60));mean(nonzeros(F2_S_N70)),mean(nonzeros(F2_L_N70))];
b14=bar(ax2,data14,'grouped');
[ngroups,nbars] = size(data14);
x14 = nan(nbars, ngroups);
for e = 1:nbars
    x14(e,:) = b14(e).XEndPoints;
end
errorbar(ax2,x14',data14,[2*std(nonzeros(F2_S_Q))/sqrt(numel(nonzeros(F2_S_Q))),2*std(nonzeros(F2_L_Q))/sqrt(numel(nonzeros(F2_L_Q)));2*std(nonzeros(F2_S_SHL))/sqrt(numel(nonzeros(F2_S_SHL))),2*std(nonzeros(F2_L_SHL))/sqrt(numel(nonzeros(F2_L_SHL)));2*std(nonzeros(F2_S_N60))/sqrt(numel(nonzeros(F2_S_N60))),2*std(nonzeros(F2_L_N60))/sqrt(numel(nonzeros(F2_L_N60)));2*std(nonzeros(F2_S_N70))/sqrt(numel(nonzeros(F2_S_N70))),2*std(nonzeros(F2_L_N70))/sqrt(numel(nonzeros(F2_L_N70)))],'k','linestyle','none','handlevisibility' ,'off')
data15=[mean(nonzeros(F3_S_Q)),mean(nonzeros(F3_L_Q));mean(nonzeros(F3_S_SHL)),mean(nonzeros(F3_L_SHL));mean(nonzeros(F3_S_N60)),mean(nonzeros(F3_L_N60));mean(nonzeros(F3_S_N70)),mean(nonzeros(F3_L_N70))];
b15=bar(ax3,data15,'grouped');
[ngroups,nbars] = size(data15);
x15 = nan(nbars, ngroups);
for e = 1:nbars
    x15(e,:) = b15(e).XEndPoints;
end
errorbar(ax3,x15',data15,[2*std(nonzeros(F3_S_Q))/sqrt(numel(nonzeros(F3_S_Q))),2*std(nonzeros(F3_L_Q))/sqrt(numel(nonzeros(F3_L_Q)));2*std(nonzeros(F3_S_SHL))/sqrt(numel(nonzeros(F3_S_SHL))),2*std(nonzeros(F3_L_SHL))/sqrt(numel(nonzeros(F3_L_SHL)));2*std(nonzeros(F3_S_N60))/sqrt(numel(nonzeros(F3_S_N60))),2*std(nonzeros(F3_L_N60))/sqrt(numel(nonzeros(F3_L_N60)));2*std(nonzeros(F3_S_N70))/sqrt(numel(nonzeros(F3_S_N70))),2*std(nonzeros(F3_L_N70))/sqrt(numel(nonzeros(F3_L_N70)))],'k','linestyle','none','handlevisibility' ,'off')
lgd15=legend(ax3,'Speaking','Non-Speaking','Location','southeastoutside');
lgd15.Title.String = 'Types of Windows:';
xticks([ax1 ax2 ax3],1:4)
xticklabels([ax1 ax2 ax3],{'Quiet','SHL','N60','N70'})
% xlim([ax2 ax1],[9.5, 10.5]);xlim([ax3 ax4],[-0.5, 0.5]);ylim([ax1 ax2 ax3 ax4],[3.4, 3.7])

% Boxplot
figure;tiledlayout(1,2);
ax11 = nexttile;
ax12 = nexttile;
[names{1:size(nonzeros(F1_S_Q(:,1)),1)}] = deal('Speaking');[names{end+1:end+size(nonzeros(F1_L_Q(:,1)),1)}] = deal('Listening');
boxchart(ax11,reshape([nonzeros(F1_S_Q(:,1)),nonzeros(F1_L_Q(:,1))],[],1),'GroupByColor',names);
boxchart(ax12,reshape([nonzeros(F1_S_Q(:,2)),nonzeros(F1_L_Q(:,2))],[],1),'GroupByColor',names);
xticklabels(ax11,'Phase 1')
xticklabels(ax12,'Phase 2')
legend