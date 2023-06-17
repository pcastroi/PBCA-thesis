%% PBCA-Thesis - Week 21 - Pupil Features for AMEND I and II - Copy of W4.m
% Pathing
clear all; clc; close all;
BPath = strsplit(pwd,'PBCA-thesis');
addpath('tools')
addpath([BPath{1} 'Pupil-preprocessing-tools\tools']) % For preprocessing

% Files and Utterances: different conditions
Param.Fs = 50; % Sampling frequency of pupil data
Param.RemoveBeforeAndAfter = [0 0]*1e-3; % Samples within the time range before and after NaNs will set NaNs as well.
Param.MinLengthNaNRepair = 0; % Drop values (i.e., change to NaN) before and after NaNs only for contiguous NaNs of at least __ samples. 
LPWinSize = 0.5; % [s]: Window size of hamming-window for low-pass filtering
LPWindow = hamming(round(LPWinSize*Param.Fs));
LPWindow = LPWindow/sum(LPWindow); % Hamming-window
AudFs = 48000;
WD = [-0.5,4]; % [s], Duration limits of phase 1 & 2
BLPeriod = [0,20]; % [s]
TimeStartW = 0.5; % [s], time before Utt/Lis starts
TimeMinWin = 0.5; % [s], Minimum time of a window
TimeInitialMerge = 0.3; % [s], Time threshold for merging windows initially
TimeMerge = 2; % [s], Time threshold for merging windows after rejecting small windows
TimeEndW = 3; % [s], time after Utt/Lis starts
RejectRatio = 0.4; % Rejection threshold based on the ratio of NaNs in data
RejectDelay = 0.5; % [s], Rejection threshold based on delay between timestamps and n-samples

% Colors
SColor = [53, 155, 67]./255;
LColor = [204, 36, 0]./255;
QuietColor = [204, 152, 0]./255;
SHLColor = [123, 31, 162]./255;
N60Color = [0, 196, 215]./255;
N70Color = [2, 36, 223]./255;
UNColor = [1, 35, 33]./255;
AAColor = [4, 105, 100]./255;
ABColor = [7, 192, 182]./255;

[subDirs_I] = GetSubDirsFirstLevelOnly('data\AMEND_I');
FileNames_I={'P1_Quiet_B1.mat','P1_Quiet_B2.mat','P1_SHL_B1.mat','P1_SHL_B2.mat','P1_Noise60_B1.mat','P1_Noise60_B2.mat','P1_Noise70_B1.mat','P1_Noise70_B2.mat','P2_Quiet_B1.mat','P2_Quiet_B2.mat','P2_SHL_B1.mat','P2_SHL_B2.mat','P2_Noise60_B1.mat','P2_Noise60_B2.mat','P2_Noise70_B1.mat','P2_Noise70_B2.mat'};
LoadDelays_I=load('data\AMEND_I\delays1110.mat');
LoadTPsOrder_I=load('data\AMEND_I\TPsOrder_I.mat');
LoadUtt_I=load('data\AMEND_I\utterances1110.mat');

[subDirs_II] = GetSubDirsFirstLevelOnly('data\AMEND_II');
subDirs_II(contains(subDirs_II,{'Pilot 1','Pair01','Pair14'}))=[]; % Only using data from Pair01 to Pair13, others removed
FileNames_II={'UNHI_N0.mat','UNHI_N60.mat','UNHI_N70.mat','AAHI_N0.mat','AAHI_N60.mat','AAHI_N70.mat','ABHI_N0.mat','ABHI_N60.mat','ABHI_N70.mat'};
LoadUtt_II=load('data\AMEND_II\Utterances0805.mat');
LoadTPsOrder_II=load('data\AMEND_II\TPsOrder_II.mat');

NPs = 2; % N of TPs per Pair
NCond = numel(FileNames_I)/NPs; % N of conditions
NTPs = numel(subDirs_I)*numel(FileNames_I)/NCond; % Total N of TPs
NCols=200; % N of windows
NRows=numel(FileNames_I)*numel(subDirs_I); % Number of trials

% Preallocate tables
Tab_I = table(0,0,0,0,0,0,'VariableNames',{'Subject','Activity','Noise','MPD','Slope','PPD'});
Tab_II = table(0,0,0,0,0,0,'VariableNames',{'Subject','Activity','Noise','MPD','Slope','PPD'});
Tab_S_II = table(0,0,0,0,0,0,'VariableNames',{'Subject','Activity','Setting','MPD','Slope','PPD'});

% Preallocate features - Phase 1, Phase 2
F1_S_Q_I=nan*ones(NRows,NCols); % Feature 1 (Quiet): Mean pupil size (Speaking)
F1_L_Q_I=F1_S_Q_I; % Feature 1 (Quiet): Mean pupil size (Listening)
F2_S_Q_I=F1_S_Q_I; % Feature 2 (Quiet): Mean slope (Speaking)
F2_L_Q_I=F1_S_Q_I; % Feature 2 (Quiet): Mean slope (Listening)
F3_S_Q_I=F1_S_Q_I; % Feature 3 (Quiet): Mean Peak Pupil Size (Speaking)
F3_L_Q_I=F1_S_Q_I; % Feature 3 (Quiet): Mean Peak Pupil Size (Listening)
F1_S_SHL_I=F1_S_Q_I; % Feature 1 (SHL): Mean pupil size (Speaking)
F1_L_SHL_I=F1_S_Q_I; % Feature 1 (SHL): Mean pupil size (Listening)
F2_S_SHL_I=F1_S_Q_I; % Feature 2 (SHL): Mean slope (Speaking)
F2_L_SHL_I=F1_S_Q_I; % Feature 2 (SHL): Mean slope (Listening)
F3_S_SHL_I=F1_S_Q_I; % Feature 3 (SHL): Mean Peak Pupil Size (Speaking)
F3_L_SHL_I=F1_S_Q_I; % Feature 3 (SHL): Mean Peak Pupil Size (Listening)
F1_S_N60_I=F1_S_Q_I; % Feature 1 (N60): Mean pupil size (Speaking)
F1_L_N60_I=F1_S_Q_I; % Feature 1 (N60): Mean pupil size (Listening)
F2_S_N60_I=F1_S_Q_I; % Feature 2 (N60): Mean slope (Speaking)
F2_L_N60_I=F1_S_Q_I; % Feature 2 (N60): Mean slope (Listening)
F3_S_N60_I=F1_S_Q_I; % Feature 3 (N60): Mean Peak Pupil Size (Speaking)
F3_L_N60_I=F1_S_Q_I; % Feature 3 (N60): Mean Peak Pupil Size (Listening)
F1_S_N70_I=F1_S_Q_I; % Feature 1 (N70): Mean pupil size (Speaking)
F1_L_N70_I=F1_S_Q_I; % Feature 1 (N70): Mean pupil size (Listening)
F2_S_N70_I=F1_S_Q_I; % Feature 2 (N70): Mean slope (Speaking)
F2_L_N70_I=F1_S_Q_I; % Feature 2 (N70): Mean slope (Listening)
F3_S_N70_I=F1_S_Q_I; % Feature 3 (N70): Mean Peak Pupil Size (Speaking)
F3_L_N70_I=F1_S_Q_I; % Feature 3 (N70): Mean Peak Pupil Size (Listening)

F1_S_Q_II=nan*ones(NRows,NCols); % Feature 1 (Quiet): Mean pupil size (Speaking)
F1_L_Q_II=F1_S_Q_II; % Feature 1 (Quiet): Mean pupil size (Listening)
F2_S_Q_II=F1_S_Q_II; % Feature 2 (Quiet): Mean slope (Speaking)
F2_L_Q_II=F1_S_Q_II; % Feature 2 (Quiet): Mean slope (Listening)
F3_S_Q_II=F1_S_Q_II; % Feature 3 (Quiet): Mean Peak Pupil Size (Speaking)
F3_L_Q_II=F1_S_Q_II; % Feature 3 (Quiet): Mean Peak Pupil Size (Listening)
F1_S_N60_II=F1_S_Q_II; % Feature 1 (N60): Mean pupil size (Speaking)
F1_L_N60_II=F1_S_Q_II; % Feature 1 (N60): Mean pupil size (Listening)
F2_S_N60_II=F1_S_Q_II; % Feature 2 (N60): Mean slope (Speaking)
F2_L_N60_II=F1_S_Q_II; % Feature 2 (N60): Mean slope (Listening)
F3_S_N60_II=F1_S_Q_II; % Feature 3 (N60): Mean Peak Pupil Size (Speaking)
F3_L_N60_II=F1_S_Q_II; % Feature 3 (N60): Mean Peak Pupil Size (Listening)
F1_S_N70_II=F1_S_Q_II; % Feature 1 (N70): Mean pupil size (Speaking)
F1_L_N70_II=F1_S_Q_II; % Feature 1 (N70): Mean pupil size (Listening)
F2_S_N70_II=F1_S_Q_II; % Feature 2 (N70): Mean slope (Speaking)
F2_L_N70_II=F1_S_Q_II; % Feature 2 (N70): Mean slope (Listening)
F3_S_N70_II=F1_S_Q_II; % Feature 3 (N70): Mean Peak Pupil Size (Speaking)
F3_L_N70_II=F1_S_Q_II; % Feature 3 (N70): Mean Peak Pupil Size (Listening)

F1_S_UN_II=nan*ones(NRows,NCols); 
F1_L_UN_II=F1_S_UN_II; 
F2_S_UN_II=F1_S_UN_II; 
F2_L_UN_II=F1_S_UN_II; 
F3_S_UN_II=F1_S_UN_II; 
F3_L_UN_II=F1_S_UN_II;
F1_S_AA_II=F1_S_UN_II;
F1_L_AA_II=F1_S_UN_II;
F2_S_AA_II=F1_S_UN_II;
F2_L_AA_II=F1_S_UN_II;
F3_S_AA_II=F1_S_UN_II;
F3_L_AA_II=F1_S_UN_II;
F1_S_AB_II=F1_S_UN_II;
F1_L_AB_II=F1_S_UN_II;
F2_S_AB_II=F1_S_UN_II;
F2_L_AB_II=F1_S_UN_II;
F3_S_AB_II=F1_S_UN_II;
F3_L_AB_II=F1_S_UN_II;

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
            TP_I = 2*q;
            if LoadTPsOrder_I.TPsOrder(2*q,i-NCond) ~= x_I
                disp(['Warning: File ',PairFiles_I(1).folder, '\', cell2mat(FileNames_I(i)),' was rejected in AMEND I (P2) analysis.']);
                continue
            end
        elseif contains(cell2mat(FileNames_I(i)),'P1')
            TP_I = 2*q-1;
            if LoadTPsOrder_I.TPsOrder(2*q-1,i) ~= x_I
                disp(['Warning: File ',PairFiles_I(1).folder, '\', cell2mat(FileNames_I(i)),' was rejected in AMEND I (P1) analysis.']);
                continue
            end
        end
        
        alldata = load([PairFiles_I(1).folder, '\', cell2mat(FileNames_I(i))]);
        alldata_mat = cell2mat(alldata.data);
        
        % Replace blanks '[]' for 'NaN' in fields diameterLeft and diameterRight
        [alldata_mat(cellfun(@isempty,{alldata_mat.diameterLeft})).diameterLeft] = deal(NaN);
        [alldata_mat(cellfun(@isempty,{alldata_mat.diameterRight})).diameterRight] = deal(NaN);
        
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
        
        % Retrieve Utterances
        if contains(cell2mat(FileNames_I(i)),'P2')
            SpeakKey = 'utteranceCH1';
            ListenKey = 'utteranceCH2';
            SDelayKey = 'delayCH1';
            LDelayKey = 'delayCH2';
        elseif contains(cell2mat(FileNames_I(i)),'P1')
            SpeakKey = 'utteranceCH2';
            ListenKey = 'utteranceCH1';
            SDelayKey = 'delayCH2';
            LDelayKey = 'delayCH1';
        end

        if contains(cell2mat(FileNames_I(i)),'B1')
            SpeB = 0;
        elseif contains(cell2mat(FileNames_I(i)),'B2')
            SpeB = 1;
        end

        if contains(cell2mat(FileNames_I(i)),'Quiet')
            SpeCond = SpeB + 1;
        elseif contains(cell2mat(FileNames_I(i)),'SHL')
            SpeCond = SpeB + 3;
        elseif contains(cell2mat(FileNames_I(i)),'Noise60')
            SpeCond = SpeB + 5;
        elseif contains(cell2mat(FileNames_I(i)),'Noise70')
            SpeCond = SpeB + 7;
        end

        SpeakRaw = PairUtt_I{1,SpeCond}.(SpeakKey);
        ListenRaw = PairUtt_I{1,SpeCond}.(ListenKey);
        binResUtt = PairUtt_I{1,SpeCond}.binRes;
        SDelayRaw = PairDelay_I{1,SpeCond}.(SDelayKey);
        LDelayRaw = PairDelay_I{1,SpeCond}.(LDelayKey);

        SpeakRaw(:,2:3)=round((SpeakRaw(:,2:3)*binResUtt+BLPeriod(2))*Param.Fs+SDelayRaw(1)/2);
        ListenRaw(:,2:3)=round((ListenRaw(:,2:3)*binResUtt+BLPeriod(2))*Param.Fs+LDelayRaw(1)/2);
        
        % SAME PROCESSING AS IN W1.m
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
        
        % Unwanted start/end means from LP filtering
        Speak = Speak(Speak(:,3)<length(Diameter)-round(length(LPWindow)/2-1)-TimeEndW*Param.Fs,:);
        Listen = Listen(Listen(:,3)<length(Diameter)-round(length(LPWindow)/2-1)-TimeEndW*Param.Fs,:);
        
        % Time-locked indexes (cannot be bigger than diameter itself)
%         SIdxEnd=Speak(:,3);
%         LIdxEnd=Listen(:,3);
        SIdxEnd_I=Speak(:,2)+TimeEndW*Param.Fs;
        SIdxEnd_I(SIdxEnd_I>length(Diameter)) = length(Diameter);
        LIdxEnd_I=Listen(:,2)+TimeEndW*Param.Fs;
        LIdxEnd_I(LIdxEnd_I>length(Diameter)) = length(Diameter);
        
        t_Diam = linspace(0,length(Diameter)./Param.Fs,length(Diameter));
        
        % Extract features for different conditions
        S_temp=zeros(NCols,1);
        S_temp2=zeros(NCols,2);
        S_temp3=S_temp;
        L_temp=S_temp;
        L_temp2=S_temp2;
        L_temp3=S_temp3;
        if contains(cell2mat(FileNames_I(i)),'Quiet')
            for j=1:size(Speak,1)
                S_temp(j,:)=mean(Diameter(Speak(j,2):SIdxEnd_I(j)),'omitnan');
                S_temp2(j,:)=polyfit(t_Diam(Speak(j,2):SIdxEnd_I(j)),Diameter(Speak(j,2):SIdxEnd_I(j)),1);
                S_temp3(j,:)=max(Diameter(Speak(j,2):SIdxEnd_I(j)));
            end
            F1_S_Q_I(x_I,1:size(S_temp(S_temp(:,1)>0,1),1))=S_temp(S_temp(:,1)>0,1);
            F2_S_Q_I(x_I,1:size(S_temp2(S_temp2(:,2)~=0,2),1))=S_temp2(S_temp2(:,2)~=0,2);
            F3_S_Q_I(x_I,1:size(S_temp3(S_temp3(:,1)>0,1),1))=S_temp3(S_temp3(:,1)>0,1);          
            % Store in table for LMM
            Tab_I = [Tab_I;table(TP_I*ones(size(Speak,1),1),...
                    repmat({'Speak'}, 1, size(Speak,1))',...
                    repmat({'Quiet'}, 1, size(Speak,1))',...
                    S_temp(S_temp(:,1)>0,1),...
                    S_temp2(S_temp2(:,2)~=0,2),...
                    S_temp3(S_temp3(:,1)>0,1),...
                    'VariableNames',{'Subject','Activity','Noise','MPD','Slope','PPD'})];
            for j=1:size(Listen,1)
                L_temp(j,:)=mean(Diameter(Listen(j,2):LIdxEnd_I(j)));
                L_temp2(j,:)=polyfit(t_Diam(Listen(j,2):LIdxEnd_I(j)),Diameter(Listen(j,2):LIdxEnd_I(j)),1);
                L_temp3(j,:)=max(Diameter(Listen(j,2):LIdxEnd_I(j)));
            end
            F1_L_Q_I(x_I,1:size(L_temp(L_temp(:,1)>0,1),1))=L_temp(L_temp(:,1)>0,1);
            F2_L_Q_I(x_I,1:size(L_temp2(L_temp2(:,2)~=0,2),1))=L_temp2(L_temp2(:,2)~=0,2);
            F3_L_Q_I(x_I,1:size(L_temp3(L_temp3(:,1)>0,1),1))=L_temp3(L_temp3(:,1)>0,1);
            % Store in table for LMM
            Tab_I = [Tab_I;table(TP_I*ones(size(Listen,1),1),...
                    repmat({'Listen'}, 1, size(Listen,1))',...
                    repmat({'Quiet'}, 1, size(Listen,1))',...
                    L_temp(L_temp(:,1)>0,1),...
                    L_temp2(L_temp2(:,2)~=0,2),...
                    L_temp3(L_temp3(:,1)>0,1),...
                    'VariableNames',{'Subject','Activity','Noise','MPD','Slope','PPD'})];
        elseif contains(cell2mat(FileNames_I(i)),'SHL')
            for j=1:size(Speak,1)
                S_temp(j,:)=mean(Diameter(Speak(j,2):SIdxEnd_I(j)),'omitnan');
                S_temp2(j,:)=polyfit(t_Diam(Speak(j,2):SIdxEnd_I(j)),Diameter(Speak(j,2):SIdxEnd_I(j)),1);
                S_temp3(j,:)=max(Diameter(Speak(j,2):SIdxEnd_I(j)));
            end
            F1_S_SHL_I(x_I,1:size(S_temp(S_temp(:,1)>0,1),1))=S_temp(S_temp(:,1)>0,1);
            F2_S_SHL_I(x_I,1:size(S_temp2(S_temp2(:,2)~=0,2),1))=S_temp2(S_temp2(:,2)~=0,2);
            F3_S_SHL_I(x_I,1:size(S_temp3(S_temp3(:,1)>0,1),1))=S_temp3(S_temp3(:,1)>0,1);
            % Store in table for LMM
            Tab_I = [Tab_I;table(TP_I*ones(size(Speak,1),1),...
                    repmat({'Speak'}, 1, size(Speak,1))',...
                    repmat({'SHL'}, 1, size(Speak,1))',...
                    S_temp(S_temp(:,1)>0,1),...
                    S_temp2(S_temp2(:,2)~=0,2),...
                    S_temp3(S_temp3(:,1)>0,1),...
                    'VariableNames',{'Subject','Activity','Noise','MPD','Slope','PPD'})];
            for j=1:size(Listen,1)
                L_temp(j,:)=mean(Diameter(Listen(j,2):LIdxEnd_I(j)));
                L_temp2(j,:)=polyfit(t_Diam(Listen(j,2):LIdxEnd_I(j)),Diameter(Listen(j,2):LIdxEnd_I(j)),1);
                L_temp3(j,:)=max(Diameter(Listen(j,2):LIdxEnd_I(j)));
            end
            F1_L_SHL_I(x_I,1:size(L_temp(L_temp(:,1)>0,1),1))=L_temp(L_temp(:,1)>0,1);
            F2_L_SHL_I(x_I,1:size(L_temp2(L_temp2(:,2)~=0,2),1))=L_temp2(L_temp2(:,2)~=0,2);
            F3_L_SHL_I(x_I,1:size(L_temp3(L_temp3(:,1)>0,1),1))=L_temp3(L_temp3(:,1)>0,1);
            % Store in table for LMM
            Tab_I = [Tab_I;table(TP_I*ones(size(Listen,1),1),...
                    repmat({'Listen'}, 1, size(Listen,1))',...
                    repmat({'SHL'}, 1, size(Listen,1))',...
                    L_temp(L_temp(:,1)>0,1),...
                    L_temp2(L_temp2(:,2)~=0,2),...
                    L_temp3(L_temp3(:,1)>0,1),...
                    'VariableNames',{'Subject','Activity','Noise','MPD','Slope','PPD'})];
        elseif contains(cell2mat(FileNames_I(i)),'Noise60')
            for j=1:size(Speak,1)
                S_temp(j,:)=mean(Diameter(Speak(j,2):SIdxEnd_I(j)),'omitnan');
                S_temp2(j,:)=polyfit(t_Diam(Speak(j,2):SIdxEnd_I(j)),Diameter(Speak(j,2):SIdxEnd_I(j)),1);
                S_temp3(j,:)=max(Diameter(Speak(j,2):SIdxEnd_I(j)));
            end
            F1_S_N60_I(x_I,1:size(S_temp(S_temp(:,1)>0,1),1))=S_temp(S_temp(:,1)>0,1);
            F2_S_N60_I(x_I,1:size(S_temp2(S_temp2(:,2)~=0,2),1))=S_temp2(S_temp2(:,2)~=0,2);
            F3_S_N60_I(x_I,1:size(S_temp3(S_temp3(:,1)>0,1),1))=S_temp3(S_temp3(:,1)>0,1);
            % Store in table for LMM
            Tab_I = [Tab_I;table(TP_I*ones(size(Speak,1),1),...
                    repmat({'Speak'}, 1, size(Speak,1))',...
                    repmat({'N60'}, 1, size(Speak,1))',...
                    S_temp(S_temp(:,1)>0,1),...
                    S_temp2(S_temp2(:,2)~=0,2),...
                    S_temp3(S_temp3(:,1)>0,1),...
                    'VariableNames',{'Subject','Activity','Noise','MPD','Slope','PPD'})];
            for j=1:size(Listen,1)
                L_temp(j,:)=mean(Diameter(Listen(j,2):LIdxEnd_I(j)));
                L_temp2(j,:)=polyfit(t_Diam(Listen(j,2):LIdxEnd_I(j)),Diameter(Listen(j,2):LIdxEnd_I(j)),1);
                L_temp3(j,:)=max(Diameter(Listen(j,2):LIdxEnd_I(j)));
            end
            F1_L_N60_I(x_I,1:size(L_temp(L_temp(:,1)>0,1),1))=L_temp(L_temp(:,1)>0,1);
            F2_L_N60_I(x_I,1:size(L_temp2(L_temp2(:,2)~=0,2),1))=L_temp2(L_temp2(:,2)~=0,2);
            F3_L_N60_I(x_I,1:size(L_temp3(L_temp3(:,1)>0,1),1))=L_temp3(L_temp3(:,1)>0,1);
            % Store in table for LMM
            Tab_I = [Tab_I;table(TP_I*ones(size(Listen,1),1),...
                    repmat({'Listen'}, 1, size(Listen,1))',...
                    repmat({'N60'}, 1, size(Listen,1))',...
                    L_temp(L_temp(:,1)>0,1),...
                    L_temp2(L_temp2(:,2)~=0,2),...
                    L_temp3(L_temp3(:,1)>0,1),...
                    'VariableNames',{'Subject','Activity','Noise','MPD','Slope','PPD'})];
        elseif contains(cell2mat(FileNames_I(i)),'Noise70')
            for j=1:size(Speak,1)
                S_temp(j,:)=mean(Diameter(Speak(j,2):SIdxEnd_I(j)),'omitnan');
                S_temp2(j,:)=polyfit(t_Diam(Speak(j,2):SIdxEnd_I(j)),Diameter(Speak(j,2):SIdxEnd_I(j)),1);
                S_temp3(j,:)=max(Diameter(Speak(j,2):SIdxEnd_I(j)));
            end
            F1_S_N70_I(x_I,1:size(S_temp(S_temp(:,1)>0,1),1))=S_temp(S_temp(:,1)>0,1);
            F2_S_N70_I(x_I,1:size(S_temp2(S_temp2(:,2)~=0,2),1))=S_temp2(S_temp2(:,2)~=0,2);
            F3_S_N70_I(x_I,1:size(S_temp3(S_temp3(:,1)>0,1),1))=S_temp3(S_temp3(:,1)>0,1);
            % Store in table for LMM
            Tab_I = [Tab_I;table(TP_I*ones(size(Speak,1),1),...
                    repmat({'Speak'}, 1, size(Speak,1))',...
                    repmat({'N70'}, 1, size(Speak,1))',...
                    S_temp(S_temp(:,1)>0,1),...
                    S_temp2(S_temp2(:,2)~=0,2),...
                    S_temp3(S_temp3(:,1)>0,1),...
                    'VariableNames',{'Subject','Activity','Noise','MPD','Slope','PPD'})];
            for j=1:size(Listen,1)
                L_temp(j,:)=mean(Diameter(Listen(j,2):LIdxEnd_I(j)),'omitnan');
                L_temp2(j,:)=polyfit(t_Diam(Listen(j,2):LIdxEnd_I(j)),Diameter(Listen(j,2):LIdxEnd_I(j)),1);
                L_temp3(j,:)=max(Diameter(Listen(j,2):LIdxEnd_I(j)));
            end
            F1_L_N70_I(x_I,1:size(L_temp(L_temp(:,1)>0,1),1))=L_temp(L_temp(:,1)>0,1);
            F2_L_N70_I(x_I,1:size(L_temp2(L_temp2(:,2)~=0,2),1))=L_temp2(L_temp2(:,2)~=0,2);
            F3_L_N70_I(x_I,1:size(L_temp3(L_temp3(:,1)>0,1),1))=L_temp3(L_temp3(:,1)>0,1);
            % Store in table for LMM
            Tab_I = [Tab_I;table(TP_I*ones(size(Listen,1),1),...
                    repmat({'Listen'}, 1, size(Listen,1))',...
                    repmat({'N70'}, 1, size(Listen,1))',...
                    L_temp(L_temp(:,1)>0,1),...
                    L_temp2(L_temp2(:,2)~=0,2),...
                    L_temp3(L_temp3(:,1)>0,1),...
                    'VariableNames',{'Subject','Activity','Noise','MPD','Slope','PPD'})];
        end        
        x_I=x_I+1;
    end
end
% Delete first row of Table
Tab_I(1,:)=[];

% [F1_S_Q_I(~isnan(F1_S_Q_I)),F2_S_Q_I(~isnan(F2_S_Q_I)),F3_S_Q_I(~isnan(F3_S_Q_I));F1_S_SHL_I(~isnan(F1_S_SHL_I)),F2_S_SHL_I(~isnan(F2_S_SHL_I)),F3_S_SHL_I(~isnan(F3_S_SHL_I));F1_S_N60_I(~isnan(F1_S_N60_I)),F2_S_N60_I(~isnan(F2_S_N60_I)),F3_S_N60_I(~isnan(F3_S_N60_I));F1_S_N70_I(~isnan(F1_S_N70_I)),F2_S_N70_I(~isnan(F2_S_N70_I)),F3_S_N70_I(~isnan(F3_S_N70_I))]
% AMEND II
for q=1:numel(subDirs_II)
    PairIn_II = q;
    for ChosenFolder = {'\NH\','\HI\'}
        PairFolder_II=[pwd,'\data\AMEND_II\',cell2mat(subDirs_II(q)),cell2mat(ChosenFolder)]; % Folder naming changed
        PairFiles_II=dir(PairFolder_II); % Folder naming changed
        PairUtt_II=LoadUtt_II.Utterances(PairIn_II,:);
        
        for i=1:numel(FileNames_II)
            if contains(ChosenFolder,'HI')
                TP_II = 2*q;
                if LoadTPsOrder_II.TPsOrder_II(2*q,i) ~= x_II
                    disp(['Warning: File ',PairFiles_II(1).folder, '\', cell2mat(FileNames_II(i)),' was rejected in AMEND II (P2) analysis.']);
                    continue
                end
            elseif contains(ChosenFolder,'NH')
                TP_II = 2*q-1;
                if LoadTPsOrder_II.TPsOrder_II(2*q-1,i) ~= x_II
                    disp(['Warning: File ',PairFiles_II(1).folder, '\', cell2mat(FileNames_II(i)),' was rejected in AMEND II (P1) analysis.']);
                    continue
                end
            end
            alldata = load([PairFiles_II(1).folder, '\', cell2mat(FileNames_II(i))]);
            alldata_mat = cell2mat(alldata.data);

            % Extract Timestamps [s]
            timestamps = vertcat(alldata_mat.timeStamp);
            
            % DIAMETER
            % Replace blanks '[]' for 'NaN' in fields diameterLeft and diameterRight
            [alldata_mat(cellfun(@isempty,{alldata_mat.diameterLeft})).diameterLeft] = deal(NaN);
            [alldata_mat(cellfun(@isempty,{alldata_mat.diameterRight})).diameterRight] = deal(NaN);

            LDiamRaw = [alldata_mat.diameterLeft];
            RDiamRaw = [alldata_mat.diameterRight];

            % Preprocessing - Setting outliers as NaNs (remove artifacts)
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

            % Comment below if No LP filtering wanted
            LDiam = LDiamConv;
            RDiam = RDiamConv;

            % Decide 'better' eye results
            [Min, idx_decision] = min([sum(LMetadata.Isnan) sum(RMetadata.Isnan)]);
            if idx_decision == 1
                Diameter = LDiam;
                DiameterRaw = LDiamRaw';
                eyeChosen = 'Left';
                DiamNaN = sum(LMetadata.Isnan);
            elseif idx_decision == 2
                Diameter = RDiam;
                DiameterRaw = RDiamRaw';
                eyeChosen = 'Right';
                DiamNaN = sum(RMetadata.Isnan);
            end
            
            % SPEAKING/LISTENING WINDOWS
            % W18 - Add nan padding (for valid baseline)
            NaNPadS = 2*TimeStartW*Param.Fs;
            Diameter = [nan*ones(NaNPadS,1);Diameter];
            
            % Retrieve Utterances
            if contains(ChosenFolder,'HI')
                SpeakKey = 'utteranceCH1';
                ListenKey = 'utteranceCH2';
            elseif contains(ChosenFolder,'NH')
                SpeakKey = 'utteranceCH2';
                ListenKey = 'utteranceCH1';
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
            elseif contains(cell2mat(FileNames_II(i)),'N60')
                SpeCond = SpeAid + 2;
            elseif contains(cell2mat(FileNames_II(i)),'N70')
                SpeCond = SpeAid + 3;
            end

            SpeakRaw = PairUtt_II{1,SpeCond}.(SpeakKey);
            ListenRaw = PairUtt_II{1,SpeCond}.(ListenKey);
            binResUtt = PairUtt_II{1,SpeCond}.binRes;

            % DIFFERENT PROCESSING FOR AMEND II (NanPadS)
            % Downsample (rounding) Utt from 250 Hz (1/binRes) to 50 Hz
            SpeakRaw(:,2:3)=round((SpeakRaw(:,2:3)*binResUtt)*Param.Fs)+NaNPadS;
            ListenRaw(:,2:3)=round((ListenRaw(:,2:3)*binResUtt)*Param.Fs)+NaNPadS;

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
            
            % Unwanted start/end means from LP filtering
            Speak = Speak(Speak(:,3)<length(Diameter)-round(length(LPWindow)/2-1)-TimeEndW*Param.Fs,:);
            Listen = Listen(Listen(:,3)<length(Diameter)-round(length(LPWindow)/2-1)-TimeEndW*Param.Fs,:);
            
            % Time-locked indexes (cannot be bigger than diameter itself)
        %         SIdxEnd=Speak(:,3);
        %         LIdxEnd=Listen(:,3);
            SIdxEnd_II=Speak(:,2)+TimeEndW*Param.Fs;
            SIdxEnd_II(SIdxEnd_II>length(Diameter)) = length(Diameter);
            LIdxEnd_II=Listen(:,2)+TimeEndW*Param.Fs;
            LIdxEnd_II(LIdxEnd_II>length(Diameter)) = length(Diameter);

            t_Diam = linspace(0,length(Diameter)./Param.Fs,length(Diameter));

            % Extract features for different conditions
            S_temp=zeros(NCols,1);
            S_temp2=zeros(NCols,2);
            S_temp3=S_temp;
            L_temp=S_temp;
            L_temp2=S_temp2;
            L_temp3=S_temp3;
            % NOISE
            if contains(cell2mat(FileNames_II(i)),'N0')
                for j=1:size(Speak,1)
                    S_temp(j,:)=mean(Diameter(Speak(j,2):SIdxEnd_II(j)));
                    S_temp2(j,:)=polyfit(t_Diam(Speak(j,2):SIdxEnd_II(j)),Diameter(Speak(j,2):SIdxEnd_II(j)),1);
                    S_temp3(j,:)=max(Diameter(Speak(j,2):SIdxEnd_II(j)));
                end
                F1_S_Q_II(x_II,1:size(S_temp(S_temp(:,1)>0,1),1))=S_temp(S_temp(:,1)>0,1);
                F2_S_Q_II(x_II,1:size(S_temp2(S_temp2(:,2)~=0,2),1))=S_temp2(S_temp2(:,2)~=0,2);
                F3_S_Q_II(x_II,1:size(S_temp3(S_temp3(:,1)>0,1),1))=S_temp3(S_temp3(:,1)>0,1);
                % Store in table for LMM
                Tab_II = [Tab_II;table(TP_II*ones(size(Speak,1),1),...
                        repmat({'Speak'}, 1, size(Speak,1))',...
                        repmat({'N0'}, 1, size(Speak,1))',...
                        S_temp(S_temp(:,1)>0,1),...
                        S_temp2(S_temp2(:,2)~=0,2),...
                        S_temp3(S_temp3(:,1)>0,1),...
                        'VariableNames',{'Subject','Activity','Noise','MPD','Slope','PPD'})];
                for j=1:size(Listen,1)
                    L_temp(j,:)=mean(Diameter(Listen(j,2):LIdxEnd_II(j)));
                    L_temp2(j,:)=polyfit(t_Diam(Listen(j,2):LIdxEnd_II(j)),Diameter(Listen(j,2):LIdxEnd_II(j)),1);
                    L_temp3(j,:)=max(Diameter(Listen(j,2):LIdxEnd_II(j)));
                end
                F1_L_Q_II(x_II,1:size(L_temp(L_temp(:,1)>0,1),1))=L_temp(L_temp(:,1)>0,1);
                F2_L_Q_II(x_II,1:size(L_temp2(L_temp2(:,2)~=0,2),1))=L_temp2(L_temp2(:,2)~=0,2);
                F3_L_Q_II(x_II,1:size(L_temp3(L_temp3(:,1)>0,1),1))=L_temp3(L_temp3(:,1)>0,1);
                % Store in table for LMM
                Tab_II = [Tab_II;table(TP_II*ones(size(Listen,1),1),...
                        repmat({'Listen'}, 1, size(Listen,1))',...
                        repmat({'N0'}, 1, size(Listen,1))',...
                        L_temp(L_temp(:,1)>0,1),...
                        L_temp2(L_temp2(:,2)~=0,2),...
                        L_temp3(L_temp3(:,1)>0,1),...
                        'VariableNames',{'Subject','Activity','Noise','MPD','Slope','PPD'})];
            elseif contains(cell2mat(FileNames_II(i)),'N60')
                for j=1:size(Speak,1)
                    S_temp(j,:)=mean(Diameter(Speak(j,2):SIdxEnd_II(j)));
                    S_temp2(j,:)=polyfit(t_Diam(Speak(j,2):SIdxEnd_II(j)),Diameter(Speak(j,2):SIdxEnd_II(j)),1);
                    S_temp3(j,:)=max(Diameter(Speak(j,2):SIdxEnd_II(j)));
                end
                F1_S_N60_II(x_II,1:size(S_temp(S_temp(:,1)>0,1),1))=S_temp(S_temp(:,1)>0,1);
                F2_S_N60_II(x_II,1:size(S_temp2(S_temp2(:,2)~=0,2),1))=S_temp2(S_temp2(:,2)~=0,2);
                F3_S_N60_II(x_II,1:size(S_temp3(S_temp3(:,1)>0,1),1))=S_temp3(S_temp3(:,1)>0,1);
                % Store in table for LMM
                Tab_II = [Tab_II;table(TP_II*ones(size(Speak,1),1),...
                        repmat({'Speak'}, 1, size(Speak,1))',...
                        repmat({'N60'}, 1, size(Speak,1))',...
                        S_temp(S_temp(:,1)>0,1),...
                        S_temp2(S_temp2(:,2)~=0,2),...
                        S_temp3(S_temp3(:,1)>0,1),...
                        'VariableNames',{'Subject','Activity','Noise','MPD','Slope','PPD'})];
                for j=1:size(Listen,1)
                    L_temp(j,:)=mean(Diameter(Listen(j,2):LIdxEnd_II(j)));
                    L_temp2(j,:)=polyfit(t_Diam(Listen(j,2):LIdxEnd_II(j)),Diameter(Listen(j,2):LIdxEnd_II(j)),1);
                    L_temp3(j,:)=max(Diameter(Listen(j,2):LIdxEnd_II(j)));
                end
                F1_L_N60_II(x_II,1:size(L_temp(L_temp(:,1)>0,1),1))=L_temp(L_temp(:,1)>0,1);
                F2_L_N60_II(x_II,1:size(L_temp2(L_temp2(:,2)~=0,2),1))=L_temp2(L_temp2(:,2)~=0,2);
                F3_L_N60_II(x_II,1:size(L_temp3(L_temp3(:,1)>0,1),1))=L_temp3(L_temp3(:,1)>0,1);
                % Store in table for LMM
                Tab_II = [Tab_II;table(TP_II*ones(size(Listen,1),1),...
                        repmat({'Listen'}, 1, size(Listen,1))',...
                        repmat({'N60'}, 1, size(Listen,1))',...
                        L_temp(L_temp(:,1)>0,1),...
                        L_temp2(L_temp2(:,2)~=0,2),...
                        L_temp3(L_temp3(:,1)>0,1),...
                        'VariableNames',{'Subject','Activity','Noise','MPD','Slope','PPD'})];
            elseif contains(cell2mat(FileNames_II(i)),'N70')
                for j=1:size(Speak,1)
                    S_temp(j,:)=mean(Diameter(Speak(j,2):SIdxEnd_II(j)),'omitnan');
                    S_temp2(j,:)=polyfit(t_Diam(Speak(j,2):SIdxEnd_II(j)),Diameter(Speak(j,2):SIdxEnd_II(j)),1);
                    S_temp3(j,:)=max(Diameter(Speak(j,2):SIdxEnd_II(j)));
                end
                F1_S_N70_II(x_II,1:size(S_temp(S_temp(:,1)>0,1),1))=S_temp(S_temp(:,1)>0,1);
                F2_S_N70_II(x_II,1:size(S_temp2(S_temp2(:,2)~=0,2),1))=S_temp2(S_temp2(:,2)~=0,2);
                F3_S_N70_II(x_II,1:size(S_temp3(S_temp3(:,1)>0,1),1))=S_temp3(S_temp3(:,1)>0,1);
                % Store in table for LMM
                Tab_II = [Tab_II;table(TP_II*ones(size(Speak,1),1),...
                        repmat({'Speak'}, 1, size(Speak,1))',...
                        repmat({'N70'}, 1, size(Speak,1))',...
                        S_temp(S_temp(:,1)>0,1),...
                        S_temp2(S_temp2(:,2)~=0,2),...
                        S_temp3(S_temp3(:,1)>0,1),...
                        'VariableNames',{'Subject','Activity','Noise','MPD','Slope','PPD'})];
                for j=1:size(Listen,1)
                    L_temp(j,:)=mean(Diameter(Listen(j,2):LIdxEnd_II(j)));
                    L_temp2(j,:)=polyfit(t_Diam(Listen(j,2):LIdxEnd_II(j)),Diameter(Listen(j,2):LIdxEnd_II(j)),1);
                    L_temp3(j,:)=max(Diameter(Listen(j,2):LIdxEnd_II(j)));
                end
                F1_L_N70_II(x_II,1:size(L_temp(L_temp(:,1)>0,1),1))=L_temp(L_temp(:,1)>0,1);
                F2_L_N70_II(x_II,1:size(L_temp2(L_temp2(:,2)~=0,2),1))=L_temp2(L_temp2(:,2)~=0,2);
                F3_L_N70_II(x_II,1:size(L_temp3(L_temp3(:,1)>0,1),1))=L_temp3(L_temp3(:,1)>0,1);
                % Store in table for LMM
                Tab_II = [Tab_II;table(TP_II*ones(size(Listen,1),1),...
                        repmat({'Listen'}, 1, size(Listen,1))',...
                        repmat({'N70'}, 1, size(Listen,1))',...
                        L_temp(L_temp(:,1)>0,1),...
                        L_temp2(L_temp2(:,2)~=0,2),...
                        L_temp3(L_temp3(:,1)>0,1),...
                        'VariableNames',{'Subject','Activity','Noise','MPD','Slope','PPD'})];
            end
            
            S_temp=zeros(NCols,1);
            S_temp2=zeros(NCols,2);
            S_temp3=S_temp;
            L_temp=S_temp;
            L_temp2=S_temp2;
            L_temp3=S_temp3;
            % HA SETTINGS - ONLY HI
            if contains(ChosenFolder,'HI')
                if contains(cell2mat(FileNames_II(i)),'UN')
                    for j=1:size(Speak,1)
                        S_temp(j,:)=mean(Diameter(Speak(j,2):SIdxEnd_II(j)),'omitnan');
                        S_temp2(j,:)=polyfit(t_Diam(Speak(j,2):SIdxEnd_II(j)),Diameter(Speak(j,2):SIdxEnd_II(j)),1);
                        S_temp3(j,:)=max(Diameter(Speak(j,2):SIdxEnd_II(j)));
                    end
                    F1_S_UN_II(x_II,1:size(S_temp(S_temp(:,1)>0,1),1))=S_temp(S_temp(:,1)>0,1);
                    F2_S_UN_II(x_II,1:size(S_temp2(S_temp2(:,2)~=0,2),1))=S_temp2(S_temp2(:,2)~=0,2);
                    F3_S_UN_II(x_II,1:size(S_temp3(S_temp3(:,1)>0,1),1))=S_temp3(S_temp3(:,1)>0,1);
                    Tab_S_II = [Tab_S_II;table(TP_II*ones(size(Speak,1),1),...
                        repmat({'Speak'}, 1, size(Speak,1))',...
                        repmat({'UN'}, 1, size(Speak,1))',...
                        S_temp(S_temp(:,1)>0,1),...
                        S_temp2(S_temp2(:,2)~=0,2),...
                        S_temp3(S_temp3(:,1)>0,1),...
                        'VariableNames',{'Subject','Activity','Setting','MPD','Slope','PPD'})];
                    for j=1:size(Listen,1)
                        L_temp(j,:)=mean(Diameter(Listen(j,2):LIdxEnd_II(j)));
                        L_temp2(j,:)=polyfit(t_Diam(Listen(j,2):LIdxEnd_II(j)),Diameter(Listen(j,2):LIdxEnd_II(j)),1);
                        L_temp3(j,:)=max(Diameter(Listen(j,2):LIdxEnd_II(j)));
                    end
                    F1_L_UN_II(x_II,1:size(L_temp(L_temp(:,1)>0,1),1))=L_temp(L_temp(:,1)>0,1);
                    F2_L_UN_II(x_II,1:size(L_temp2(L_temp2(:,2)~=0,2),1))=L_temp2(L_temp2(:,2)~=0,2);
                    F3_L_UN_II(x_II,1:size(L_temp3(L_temp3(:,1)>0,1),1))=L_temp3(L_temp3(:,1)>0,1);
                    Tab_S_II = [Tab_S_II;table(TP_II*ones(size(Listen,1),1),...
                        repmat({'Listen'}, 1, size(Listen,1))',...
                        repmat({'UN'}, 1, size(Listen,1))',...
                        L_temp(L_temp(:,1)>0,1),...
                        L_temp2(L_temp2(:,2)~=0,2),...
                        L_temp3(L_temp3(:,1)>0,1),...
                        'VariableNames',{'Subject','Activity','Setting','MPD','Slope','PPD'})];
                elseif contains(cell2mat(FileNames_II(i)),'AA')
                    for j=1:size(Speak,1)
                        S_temp(j,:)=mean(Diameter(Speak(j,2):SIdxEnd_II(j)),'omitnan');
                        S_temp2(j,:)=polyfit(t_Diam(Speak(j,2):SIdxEnd_II(j)),Diameter(Speak(j,2):SIdxEnd_II(j)),1);
                        S_temp3(j,:)=max(Diameter(Speak(j,2):SIdxEnd_II(j)));
                    end
                    F1_S_AA_II(x_II,1:size(S_temp(S_temp(:,1)>0,1),1))=S_temp(S_temp(:,1)>0,1);
                    F2_S_AA_II(x_II,1:size(S_temp2(S_temp2(:,2)~=0,2),1))=S_temp2(S_temp2(:,2)~=0,2);
                    F3_S_AA_II(x_II,1:size(S_temp3(S_temp3(:,1)>0,1),1))=S_temp3(S_temp3(:,1)>0,1);
                    Tab_S_II = [Tab_S_II;table(TP_II*ones(size(Speak,1),1),...
                        repmat({'Speak'}, 1, size(Speak,1))',...
                        repmat({'AA'}, 1, size(Speak,1))',...
                        S_temp(S_temp(:,1)>0,1),...
                        S_temp2(S_temp2(:,2)~=0,2),...
                        S_temp3(S_temp3(:,1)>0,1),...
                        'VariableNames',{'Subject','Activity','Setting','MPD','Slope','PPD'})];
                    for j=1:size(Listen,1)
                        L_temp(j,:)=mean(Diameter(Listen(j,2):LIdxEnd_II(j)));
                        L_temp2(j,:)=polyfit(t_Diam(Listen(j,2):LIdxEnd_II(j)),Diameter(Listen(j,2):LIdxEnd_II(j)),1);
                        L_temp3(j,:)=max(Diameter(Listen(j,2):LIdxEnd_II(j)));
                    end
                    F1_L_AA_II(x_II,1:size(L_temp(L_temp(:,1)>0,1),1))=L_temp(L_temp(:,1)>0,1);
                    F2_L_AA_II(x_II,1:size(L_temp2(L_temp2(:,2)~=0,2),1))=L_temp2(L_temp2(:,2)~=0,2);
                    F3_L_AA_II(x_II,1:size(L_temp3(L_temp3(:,1)>0,1),1))=L_temp3(L_temp3(:,1)>0,1);
                    Tab_S_II = [Tab_S_II;table(TP_II*ones(size(Listen,1),1),...
                        repmat({'Listen'}, 1, size(Listen,1))',...
                        repmat({'AA'}, 1, size(Listen,1))',...
                        L_temp(L_temp(:,1)>0,1),...
                        L_temp2(L_temp2(:,2)~=0,2),...
                        L_temp3(L_temp3(:,1)>0,1),...
                        'VariableNames',{'Subject','Activity','Setting','MPD','Slope','PPD'})];
                elseif contains(cell2mat(FileNames_II(i)),'AB')
                    for j=1:size(Speak,1)
                        S_temp(j,:)=mean(Diameter(Speak(j,2):SIdxEnd_II(j)),'omitnan');
                        S_temp2(j,:)=polyfit(t_Diam(Speak(j,2):SIdxEnd_II(j)),Diameter(Speak(j,2):SIdxEnd_II(j)),1);
                        S_temp3(j,:)=max(Diameter(Speak(j,2):SIdxEnd_II(j)));
                    end
                    F1_S_AB_II(x_II,1:size(S_temp(S_temp(:,1)>0,1),1))=S_temp(S_temp(:,1)>0,1);
                    F2_S_AB_II(x_II,1:size(S_temp2(S_temp2(:,2)~=0,2),1))=S_temp2(S_temp2(:,2)~=0,2);
                    F3_S_AB_II(x_II,1:size(S_temp3(S_temp3(:,1)>0,1),1))=S_temp3(S_temp3(:,1)>0,1);
                    Tab_S_II = [Tab_S_II;table(TP_II*ones(size(Speak,1),1),...
                        repmat({'Speak'}, 1, size(Speak,1))',...
                        repmat({'AB'}, 1, size(Speak,1))',...
                        S_temp(S_temp(:,1)>0,1),...
                        S_temp2(S_temp2(:,2)~=0,2),...
                        S_temp3(S_temp3(:,1)>0,1),...
                        'VariableNames',{'Subject','Activity','Setting','MPD','Slope','PPD'})];
                    for j=1:size(Listen,1)
                        L_temp(j,:)=mean(Diameter(Listen(j,2):LIdxEnd_II(j)),'omitnan');
                        L_temp2(j,:)=polyfit(t_Diam(Listen(j,2):LIdxEnd_II(j)),Diameter(Listen(j,2):LIdxEnd_II(j)),1);
                        L_temp3(j,:)=max(Diameter(Listen(j,2):LIdxEnd_II(j)));
                    end
                    F1_L_AB_II(x_II,1:size(L_temp(L_temp(:,1)>0,1),1))=L_temp(L_temp(:,1)>0,1);
                    F2_L_AB_II(x_II,1:size(L_temp2(L_temp2(:,2)~=0,2),1))=L_temp2(L_temp2(:,2)~=0,2);
                    F3_L_AB_II(x_II,1:size(L_temp3(L_temp3(:,1)>0,1),1))=L_temp3(L_temp3(:,1)>0,1);
                    Tab_S_II = [Tab_S_II;table(TP_II*ones(size(Listen,1),1),...
                        repmat({'Listen'}, 1, size(Listen,1))',...
                        repmat({'AB'}, 1, size(Listen,1))',...
                        L_temp(L_temp(:,1)>0,1),...
                        L_temp2(L_temp2(:,2)~=0,2),...
                        L_temp3(L_temp3(:,1)>0,1),...
                        'VariableNames',{'Subject','Activity','Setting','MPD','Slope','PPD'})];
                end
            end
            x_II=x_II+1;
        end
    end
end
% Delete first row of Table
Tab_II(1,:)=[];
Tab_S_II(1,:)=[];
%% Linear Mixed Model approach - show statistical differences between conditions/settings
% p = 0.05 significacnt difference -> discussion: i dont'a avg across
% comparisons (not divinding by the total number of comparisons - Bonferroni correction)
% 1 given subject = (1|subject)
% categorical - table

% lme_S_F1_I = fitlme(tbl_S_F1_I,'Quiet~SHL~N60~N70');

%% Plots
f1=figure;t1=tiledlayout(1,2);
ax1 = nexttile;
ax2 = nexttile;
f2=figure;t2=tiledlayout(1,2);
ax3 = nexttile;
ax4 = nexttile;
f3=figure;t3=tiledlayout(1,2);
ax5 = nexttile;
ax6 = nexttile;
f4=figure;t4=tiledlayout(1,2);
ax7 = nexttile;
ax8 = nexttile;
f5=figure;t5=tiledlayout(1,2);
ax9 = nexttile;
ax10 = nexttile;
f6=figure;t6=tiledlayout(1,2);
ax11 = nexttile;
ax12 = nexttile;
f7=figure;t7=tiledlayout(1,2);
ax13 = nexttile;
ax14 = nexttile;
f8=figure;t8=tiledlayout(1,2);
ax15 = nexttile;
ax16 = nexttile;
f9=figure;t9=tiledlayout(1,2);
ax17 = nexttile;
ax18 = nexttile;

hold([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10 ax11 ax12 ax13 ax14 ax15 ax16 ax17 ax18],'on')

% Boxplot
% Noise Conditions A_I
% Table_I(and(contains(Table_I.Activity,'Listen')

boxplotGroup(ax1,{F1_S_Q_I(~isnan(F1_S_Q_I)),F1_S_SHL_I(~isnan(F1_S_SHL_I)),F1_S_N60_I(~isnan(F1_S_N60_I)),F1_S_N70_I(~isnan(F1_S_N70_I))},...
'PrimaryLabels', {'Quiet','SHL','N60','N70'}, ...
'SecondaryLabels', {'Speaking'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[QuietColor;SHLColor;N60Color;N70Color],'GroupType','betweenGroups');
boxplotGroup(ax2,{F1_L_Q_I(~isnan(F1_L_Q_I)),F1_L_SHL_I(~isnan(F1_L_SHL_I)),F1_L_N60_I(~isnan(F1_L_N60_I)),F1_L_N70_I(~isnan(F1_L_N70_I))},...
'PrimaryLabels', {'Quiet','SHL','N60','N70'}, ...
'SecondaryLabels', {'Listening'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[QuietColor;SHLColor;N60Color;N70Color],'GroupType','betweenGroups');

boxplotGroup(ax3,{F2_S_Q_I(~isnan(F2_S_Q_I)),F2_S_SHL_I(~isnan(F2_S_SHL_I)),F2_S_N60_I(~isnan(F2_S_N60_I)),F2_S_N70_I(~isnan(F2_S_N70_I))},...
'PrimaryLabels', {'Quiet','SHL','N60','N70'}, ...
'SecondaryLabels', {'Speaking'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[QuietColor;SHLColor;N60Color;N70Color],'GroupType','betweenGroups');
boxplotGroup(ax4,{F2_L_Q_I(~isnan(F2_L_Q_I)),F2_L_SHL_I(~isnan(F2_L_SHL_I)),F2_L_N60_I(~isnan(F2_L_N60_I)),F2_L_N70_I(~isnan(F2_L_N70_I))},...
'PrimaryLabels', {'Quiet','SHL','N60','N70'}, ...
'SecondaryLabels', {'Listening'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[QuietColor;SHLColor;N60Color;N70Color],'GroupType','betweenGroups');

boxplotGroup(ax5,{F3_S_Q_I(~isnan(F3_S_Q_I)),F3_S_SHL_I(~isnan(F3_S_SHL_I)),F3_S_N60_I(~isnan(F3_S_N60_I)),F3_S_N70_I(~isnan(F3_S_N70_I))},...
'PrimaryLabels', {'Quiet','SHL','N60','N70'}, ...
'SecondaryLabels', {'Speaking'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[QuietColor;SHLColor;N60Color;N70Color],'GroupType','betweenGroups');
boxplotGroup(ax6,{F3_L_Q_I(~isnan(F3_L_Q_I)),F3_L_SHL_I(~isnan(F3_L_SHL_I)),F3_L_N60_I(~isnan(F3_L_N60_I)),F3_L_N70_I(~isnan(F3_L_N70_I))},...
'PrimaryLabels', {'Quiet','SHL','N60','N70'}, ...
'SecondaryLabels', {'Listening'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[QuietColor;SHLColor;N60Color;N70Color],'GroupType','betweenGroups');

ax1.YGrid = 'on';
ax1.XTickLabelRotation = 90;
set(ax1,'Color',[SColor,0.04])
ax2.YGrid = 'on';
ax2.XTickLabelRotation = 90;
set(ax2,'Color',[LColor,0.04])
ylim([ax1 ax2],[min([ylim(ax1) ylim(ax2)]) max([ylim(ax1) ylim(ax2)])])

ax3.YGrid = 'on';
ax3.XTickLabelRotation = 90;
set(ax3,'Color',[SColor,0.04])
ax4.YGrid = 'on';
ax4.XTickLabelRotation = 90;
set(ax4,'Color',[LColor,0.04])
ylim([ax3 ax4],[min([ylim(ax3) ylim(ax4)]) max([ylim(ax3) ylim(ax4)])])

ax5.YGrid = 'on';
ax5.XTickLabelRotation = 90;
set(ax5,'Color',[SColor,0.04])
ax6.YGrid = 'on';
ax6.XTickLabelRotation = 90;
set(ax6,'Color',[LColor,0.04])
ylim([ax5 ax6],[min([ylim(ax5) ylim(ax6)]) max([ylim(ax5) ylim(ax6)])])

% Noise Conditions A_II
boxplotGroup(ax7,{F1_S_Q_II(~isnan(F1_S_Q_II)),F1_S_N60_II(~isnan(F1_S_N60_II)),F1_S_N70_II(~isnan(F1_S_N70_II))},...
'PrimaryLabels', {'Quiet','N60','N70'}, ...
'SecondaryLabels', {'Speaking'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[QuietColor;N60Color;N70Color],'GroupType','betweenGroups');
boxplotGroup(ax8,{F1_L_Q_II(~isnan(F1_L_Q_II)),F1_L_N60_II(~isnan(F1_L_N60_II)),F1_L_N70_II(~isnan(F1_L_N70_II))},...
'PrimaryLabels', {'Quiet','N60','N70'}, ...
'SecondaryLabels', {'Listening'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[QuietColor;N60Color;N70Color],'GroupType','betweenGroups');

boxplotGroup(ax9,{F2_S_Q_II(~isnan(F2_S_Q_II)),F2_S_N60_II(~isnan(F2_S_N60_II)),F2_S_N70_II(~isnan(F2_S_N70_II))},...
'PrimaryLabels', {'Quiet','N60','N70'}, ...
'SecondaryLabels', {'Speaking'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[QuietColor;N60Color;N70Color],'GroupType','betweenGroups');
boxplotGroup(ax10,{F2_L_Q_II(~isnan(F2_L_Q_II)),F2_L_N60_II(~isnan(F2_L_N60_II)),F2_L_N70_II(~isnan(F2_L_N70_II))},...
'PrimaryLabels', {'Quiet','N60','N70'}, ...
'SecondaryLabels', {'Listening'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[QuietColor;N60Color;N70Color],'GroupType','betweenGroups');

boxplotGroup(ax11,{F3_S_Q_II(~isnan(F3_S_Q_II)),F3_S_N60_II(~isnan(F3_S_N60_II)),F3_S_N70_II(~isnan(F3_S_N70_II))},...
'PrimaryLabels', {'Quiet','N60','N70'}, ...
'SecondaryLabels', {'Speaking'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[QuietColor;N60Color;N70Color],'GroupType','betweenGroups');
boxplotGroup(ax12,{F3_L_Q_II(~isnan(F3_L_Q_II)),F3_L_N60_II(~isnan(F3_L_N60_II)),F3_L_N70_II(~isnan(F3_L_N70_II))},...
'PrimaryLabels', {'Quiet','N60','N70'}, ...
'SecondaryLabels', {'Listening'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[QuietColor;N60Color;N70Color],'GroupType','betweenGroups');

ax7.YGrid = 'on';
ax7.XTickLabelRotation = 90;
set(ax7,'Color',[SColor,0.04])
ax8.YGrid = 'on';
ax8.XTickLabelRotation = 90;
set(ax8,'Color',[LColor,0.04])
ylim([ax7 ax8],[min([ylim(ax7) ylim(ax8)]) max([ylim(ax7) ylim(ax8)])])

ax9.YGrid = 'on';
ax9.XTickLabelRotation = 90;
set(ax9,'Color',[SColor,0.04])
ax10.YGrid = 'on';
ax10.XTickLabelRotation = 90;
set(ax10,'Color',[LColor,0.04])
ylim([ax9 ax10],[min([ylim(ax9) ylim(ax10)]) max([ylim(ax9) ylim(ax10)])])

ax11.YGrid = 'on';
ax11.XTickLabelRotation = 90;
set(ax11,'Color',[SColor,0.04])
ax12.YGrid = 'on';
ax12.XTickLabelRotation = 90;
set(ax12,'Color',[LColor,0.04])
ylim([ax11 ax12],[min([ylim(ax11) ylim(ax12)]) max([ylim(ax11) ylim(ax12)])])

% Settings A_II
boxplotGroup(ax13,{F1_S_UN_II(~isnan(F1_S_UN_II)),F1_S_AA_II(~isnan(F1_S_AA_II)),F1_S_AB_II(~isnan(F1_S_AB_II))},...
'PrimaryLabels', {'UN','AA','AB'}, ...
'SecondaryLabels', {'Speaking'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[UNColor;AAColor;ABColor],'GroupType','betweenGroups');
boxplotGroup(ax14,{F1_L_UN_II(~isnan(F1_L_UN_II)),F1_L_AA_II(~isnan(F1_L_AA_II)),F1_L_AB_II(~isnan(F1_L_AB_II))},...
'PrimaryLabels', {'UN','AA','AB'}, ...
'SecondaryLabels', {'Listening'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[UNColor;AAColor;ABColor],'GroupType','betweenGroups');

boxplotGroup(ax15,{F2_S_UN_II(~isnan(F2_S_UN_II)),F2_S_AA_II(~isnan(F2_S_AA_II)),F2_S_AB_II(~isnan(F2_S_AB_II))},...
'PrimaryLabels', {'UN','AA','AB'}, ...
'SecondaryLabels', {'Speaking'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[UNColor;AAColor;ABColor],'GroupType','betweenGroups');
boxplotGroup(ax16,{F2_L_UN_II(~isnan(F2_L_UN_II)),F2_L_AA_II(~isnan(F2_L_AA_II)),F2_L_AB_II(~isnan(F2_L_AB_II))},...
'PrimaryLabels', {'UN','AA','AB'}, ...
'SecondaryLabels', {'Listening'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[UNColor;AAColor;ABColor],'GroupType','betweenGroups');

boxplotGroup(ax17,{F3_S_UN_II(~isnan(F3_S_UN_II)),F3_S_AA_II(~isnan(F3_S_AA_II)),F3_S_AB_II(~isnan(F3_S_AB_II))},...
'PrimaryLabels', {'UN','AA','AB'}, ...
'SecondaryLabels', {'Speaking'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[UNColor;AAColor;ABColor],'GroupType','betweenGroups');
boxplotGroup(ax18,{F3_L_UN_II(~isnan(F3_L_UN_II)),F3_L_AA_II(~isnan(F3_L_AA_II)),F3_L_AB_II(~isnan(F3_L_AB_II))},...
'PrimaryLabels', {'UN','AA','AB'}, ...
'SecondaryLabels', {'Listening'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[UNColor;AAColor;ABColor],'GroupType','betweenGroups');

ax13.YGrid = 'on';
ax13.XTickLabelRotation = 90;
set(ax13,'Color',[SColor,0.04])
ax14.YGrid = 'on';
ax14.XTickLabelRotation = 90;
set(ax14,'Color',[LColor,0.04])
ylim([ax13 ax14],[min([ylim(ax13) ylim(ax14)]) max([ylim(ax13) ylim(ax14)])])

ax15.YGrid = 'on';
ax15.XTickLabelRotation = 90;
set(ax15,'Color',[SColor,0.04])
ax16.YGrid = 'on';
ax16.XTickLabelRotation = 90;
set(ax16,'Color',[LColor,0.04])
ylim([ax15 ax16],[min([ylim(ax15) ylim(ax16)]) max([ylim(ax15) ylim(ax16)])])

ax17.YGrid = 'on';
ax17.XTickLabelRotation = 90;
set(ax17,'Color',[SColor,0.04])
ax18.YGrid = 'on';
ax18.XTickLabelRotation = 90;
set(ax18,'Color',[LColor,0.04])
ylim([ax17 ax18],[min([ylim(ax17) ylim(ax18)]) max([ylim(ax17) ylim(ax18)])])

% title(tl11,'Mean Pupil Diameter','FontWeight','bold')
% title(tl12,'Pupil Diameter Slope','FontWeight','bold')
% title(tl13,'Peak Pupil Diameter','FontWeight','bold')

ylabel([ax1 ax5 ax7 ax11 ax13 ax17],'Pupil diameter [mm]')
ylabel([ax3 ax9 ax15],'Slope [mm/s]')