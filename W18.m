%% PBCA-Thesis - Week 18 - AMEND II data - Pupillometry
% Pathing
clear all; clc; close all;
BPath = strsplit(pwd,'PBCA-thesis');
addpath('tools')
addpath([BPath{1} 'Pupil-preprocessing-tools\tools']) % For preprocessing

% Colors
SColor = [53, 155, 67]./255;
LColor = [204, 36, 0]./255;
NHColor = [204, 62, 0]./255;
HIColor = [164, 0, 204]./255;
UNColor = [185, 193, 254]./255;
AAColor = [116, 132, 252]./255;
ABColor = [4, 30, 222]./255;
QuietColor = [204, 152, 0]./255;
SHLColor = [123, 31, 162]./255;
N60Color = [0, 196, 215]./255;
N70Color = [2, 36, 223]./255;

% Important variables
[subDirs_I] = GetSubDirsFirstLevelOnly('data\AMEND_I');
[subDirs_II] = GetSubDirsFirstLevelOnly('data\AMEND_II');
subDirs_II(contains(subDirs_II,{'Pilot 1','Pair01','Pair14'}))=[]; % Only using data from Pair01 to Pair13, others removed
FileNames_I={'P1_N0_B1.mat','P1_N0_B2.mat','P1_N60_B1.mat','P1_N60_B2.mat','P1_Noise60_B1.mat','P1_Noise60_B2.mat','P1_Noise70_B1.mat','P1_Noise70_B2.mat','P2_N0_B1.mat','P2_N0_B2.mat','P2_N60_B1.mat','P2_N60_B2.mat','P2_Noise60_B1.mat','P2_Noise60_B2.mat','P2_Noise70_B1.mat','P2_Noise70_B2.mat'};
FileNames_II={'UNHI_N0.mat','UNHI_N60.mat','UNHI_N70.mat','AAHI_N0.mat','AAHI_N60.mat','AAHI_N70.mat','ABHI_N0.mat','ABHI_N60.mat','ABHI_N70.mat'};
LoadUtt_I=load('data\AMEND_I\utterances1110.mat');
LoadUtt_II=load('data\AMEND_II\Utterances0805.mat');
LoadDelays_I=load('data\AMEND_I\delays1110.mat');
% LoadDelays_II=load('data\AMEND_II\delays1110.mat'); % Should be no delay
Param.Fs = 50; % Sampling frequency of pupil data
Param.RemoveBeforeAndAfter = [35 100]*1e-3; % Samples within the time range before and after NaNs will set NaNs as well.
Param.MinLengthNaNRepair = 5; % Drop values (i.e., change to NaN) before and after NaNs only for contiguous NaNs of at least __ samples. 
LPWinSize = 0.5; % [s]: Window size of hamming-window for low-pass filtering
FilterWidth = round((LPWinSize*Param.Fs)/2); % [samples]: Width of hamming filter used for fixation duration
LPWindow = hamming(round(LPWinSize*Param.Fs));
LPWindow = LPWindow/sum(LPWindow); % Hamming-window
AudFs = 48000; % [Hz], sampling frequency of the audio files
AdapBL = 0.3; % [s], Duration of baseline prior to event
TimeMinWin = 0.5; % [s], Minimum time of a window
TimeInitialMerge = 0.3; % [s], Time threshold for merging windows initially
TimeMerge = 2; % [s], Time threshold for merging windows after rejecting small windows
VisRatio = 0.8; % Percentage of data to visualize in plots
RejectRatio = 0.4; % Rejection threshold based on the ratio of NaNs in data
RejectDelay = 0.5; % [s], Rejection threshold based on delay between timestamps and n-samples
TimeStartW = 0.5; % [s], time before Utt/Lis starts
TimeEndW = 0; % [s], time after Utt/Lis starts
x = 1; % idx to store global values
% NCond_II = numel(FileNames_II)/NPs; % N of conditions
% NTPs = numel(subDirs)*numel(FileNames)/NCond; % Total N of TPs
% TPsOrder = zeros(NTPs,NCond_II); % Vector that will contain the indexes of trials/file for each TP

NCols=60*Param.Fs; % Duration (samples) of each window
NRows=100; % Number of windows per trial
NLayers=numel(FileNames_II)*numel(subDirs_II); % Number of trials

GSW = zeros(NLayers,NRows,NCols); GLW = GSW; % Global Speaking/Listening Windows
SW_NH = GSW; LW_NH = GSW; % NH Speaking/Listening Windows
SW_HI = GSW; LW_HI = GSW; % HI Speaking/Listening Windows
SW_N0 = GSW; LW_N0 = GSW; % N0 Speaking/Listening Windows
SW_N60 = GSW; LW_N60 = GSW; % N60 Speaking/Listening Windows
SW_N70 = GSW; LW_N70 = GSW; % N70 Speaking/Listening Windows
SW_NH_N0 = GSW; LW_NH_N0 = GSW; % NH N0 Speaking/Listening Windows
SW_NH_N60 = GSW; LW_NH_N60 = GSW; % NH N60 Speaking/Listening Windows
SW_NH_N70 = GSW; LW_NH_N70 = GSW; % NH N70 Speaking/Listening Windows
SW_HI_N0 = GSW; LW_HI_N0 = GSW; % HI N0 Speaking/Listening Windows
SW_HI_N60 = GSW; LW_HI_N60 = GSW; % HI N60 Speaking/Listening Windows
SW_HI_N70 = GSW; LW_HI_N70 = GSW; % HI N70 Speaking/Listening Windows
SW_HI_UN = GSW; LW_HI_UN = GSW; % HI UN Speaking/Listening Windows
SW_HI_AA = GSW; LW_HI_AA = GSW; % HI AA Speaking/Listening Windows
SW_HI_AB = GSW; LW_HI_AB = GSW; % HI AB Speaking/Listening Windows
SW_HI_UN_N0 = GSW; LW_HI_UN_N0 = GSW; % HI UN N0 Speaking/Listening Windows
SW_HI_UN_N60 = GSW; LW_HI_UN_N60 = GSW; % HI UN N60 Speaking/Listening Windows
SW_HI_UN_N70 = GSW; LW_HI_UN_N70 = GSW; % HI UN N70 Speaking/Listening Windows
SW_HI_AA_N0 = GSW; LW_HI_AA_N0 = GSW; % HI AA N0 Speaking/Listening Windows
SW_HI_AA_N60 = GSW; LW_HI_AA_N60 = GSW; % HI AA N60 Speaking/Listening Windows
SW_HI_AA_N70 = GSW; LW_HI_AA_N70 = GSW; % HI AA N70 Speaking/Listening Windows
SW_HI_AB_N0 = GSW; LW_HI_AB_N0 = GSW; % HI AB N0 Speaking/Listening Windows
SW_HI_AB_N60 = GSW; LW_HI_AB_N60 = GSW; % HI AB N60 Speaking/Listening Windows
SW_HI_AB_N70 = GSW; LW_HI_AB_N70 = GSW; % HI AB N70 Speaking/Listening Windows

GSW_B = GSW; GLW_B = GSW; % Global Speaking/Listening Windows Adaptive-Baseline corrected
SW_NH_B = GSW; LW_NH_B = GSW; % NH Speaking/Listening Windows Adaptive-Baseline corrected
SW_HI_B = GSW; LW_HI_B = GSW; % HI Speaking/Listening Windows Adaptive-Baseline corrected
SW_N0_B = GSW; LW_N0_B = GSW; % N0 Speaking/Listening Windows Adaptive-Baseline corrected
SW_N60_B = GSW; LW_N60_B = GSW; % N60 Speaking/Listening Windows Adaptive-Baseline corrected
SW_N70_B = GSW; LW_N70_B = GSW; % N70 Speaking/Listening Windows Adaptive-Baseline corrected
SW_NH_N0_B = GSW; LW_NH_N0_B = GSW; % NH N0 Speaking/Listening Windows Adaptive-Baseline corrected
SW_NH_N60_B = GSW; LW_NH_N60_B = GSW; % NH N60 Speaking/Listening Windows Adaptive-Baseline corrected
SW_NH_N70_B = GSW; LW_NH_N70_B = GSW; % NH N70 Speaking/Listening Windows Adaptive-Baseline corrected
SW_HI_N0_B = GSW; LW_HI_N0_B = GSW; % HI N0 Speaking/Listening Windows Adaptive-Baseline corrected
SW_HI_N60_B = GSW; LW_HI_N60_B = GSW; % HI N60 Speaking/Listening Windows Adaptive-Baseline corrected
SW_HI_N70_B = GSW; LW_HI_N70_B = GSW; % HI N70 Speaking/Listening Windows Adaptive-Baseline corrected
SW_HI_UN_B = GSW; LW_HI_UN_B = GSW; % HI UN Speaking/Listening Windows Adaptive-Baseline corrected
SW_HI_AA_B = GSW; LW_HI_AA_B = GSW; % HI AA Speaking/Listening Windows Adaptive-Baseline corrected
SW_HI_AB_B = GSW; LW_HI_AB_B = GSW; % HI AB Speaking/Listening Windows Adaptive-Baseline corrected
SW_HI_UN_N0_B = GSW; LW_HI_UN_N0_B = GSW; % HI UN N0 Speaking/Listening Windows Adaptive-Baseline corrected
SW_HI_UN_N60_B = GSW; LW_HI_UN_N60_B = GSW; % HI UN N60 Speaking/Listening Windows Adaptive-Baseline corrected
SW_HI_UN_N70_B = GSW; LW_HI_UN_N70_B = GSW; % HI UN N70 Speaking/Listening Windows Adaptive-Baseline corrected
SW_HI_AA_N0_B = GSW; LW_HI_AA_N0_B = GSW; % HI AA N0 Speaking/Listening Windows Adaptive-Baseline corrected
SW_HI_AA_N60_B = GSW; LW_HI_AA_N60_B = GSW; % HI AA N60 Speaking/Listening Windows Adaptive-Baseline corrected
SW_HI_AA_N70_B = GSW; LW_HI_AA_N70_B = GSW; % HI AA N70 Speaking/Listening Windows Adaptive-Baseline corrected
SW_HI_AB_N0_B = GSW; LW_HI_AB_N0_B = GSW; % HI AB N0 Speaking/Listening Windows Adaptive-Baseline corrected
SW_HI_AB_N60_B = GSW; LW_HI_AB_N60_B = GSW; % HI AB N60 Speaking/Listening Windows Adaptive-Baseline corrected
SW_HI_AB_N70_B = GSW; LW_HI_AB_N70_B = GSW; % HI AB N70 Speaking/Listening Windows Adaptive-Baseline corrected

GSDur = zeros(NLayers,NRows); % Global Speaking Windows Duration
GLDur = zeros(NLayers,NRows); % Global Listening Windows Duration
SDur_N0 = zeros(NLayers,NRows); % N0 Speaking Windows Duration
LDur_N0 = zeros(NLayers,NRows); % N0 Listening Windows Duration
SDur_N60 = zeros(NLayers,NRows); % N60 Speaking Windows Duration
LDur_N60 = zeros(NLayers,NRows); % N60 Listening Windows Duration
SDur_N70 = zeros(NLayers,NRows); % N70 Speaking Windows Duration
LDur_N70 = zeros(NLayers,NRows); % N70 Listening Windows Duration

% Loop for AMEND II
for q=1:numel(subDirs_II)
    PairIn_II = q;
    for ChosenFolder = {'\NH\','\HI\'}
        PairFolder_II=[pwd,'\data\AMEND_II\',cell2mat(subDirs_II(q)),cell2mat(ChosenFolder)]; % Folder naming changed
        PairFiles_II=dir(PairFolder_II); % Folder naming changed
        PairUtt_II=LoadUtt_II.Utterances(PairIn_II,:);
        
        if isempty(PairUtt_II{1})
            disp(['Warning: No Utterance found for folder "',cell2mat(ChosenFolder),cell2mat(subDirs_II(q)),'"'])
            continue
        end
%         PairDelay_II=LoadDelays_II.TobAudDelay(PairIn_II,:);
    
        for i=1:numel(FileNames_II)
            try
                alldata = load([PairFiles_II(1).folder, '\', cell2mat(FileNames_II(i))]);
            catch ME
                disp(['Warning: File ', PairFiles_II(1).folder, '\', cell2mat(FileNames_II(i)), ' not found (no Gaze data).']);
                continue
            end
            alldata_mat = cell2mat(alldata.data);

            % Replace blanks '[]' for 'NaN' in fields diameterLeft and diameterRight
            [alldata_mat(cellfun(@isempty,{alldata_mat.diameterLeft})).diameterLeft] = deal(NaN);
            [alldata_mat(cellfun(@isempty,{alldata_mat.diameterRight})).diameterRight] = deal(NaN);

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
            LDiam = LDiamConv;
            RDiam = RDiamConv;

            % Decide 'better' eye results
            [Min, idx_decision] = min([sum(LMetadata.Isnan) sum(RMetadata.Isnan)]);
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
                disp(['Warning: File ',PairFiles_II(1).folder, '\', cell2mat(FileNames_II(i)), ' was rejected because it contains too many NaNs (',sprintf('%0.2f',100*DiamNaN/length(Diameter)),'%).'])
                continue
            end
            
            % Add nan padding - Diameter
            NaNPadS = 2*TimeStartW*Param.Fs;
            Diameter = [nan*ones(NaNPadS,1);Diameter];
            
            % Retrieve Utterances
            if contains(ChosenFolder,'NH')
                SpeakKey = 'utteranceCH1';
                ListenKey = 'utteranceCH2';
            elseif contains(ChosenFolder,'HI')
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

            if isempty(SpeakRaw) && isempty(ListenRaw)
                disp(['Warning: File ',PairFiles_II(1).folder, '\', cell2mat(FileNames_II(i)),' was rejected for not having associated Utterance windows.']);
                continue
            end
            
            % SAME PROCESSING AS IN W1.m
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

            % Time-locked indexes (based on Start or End of events)
            WSpeakIdx=[Speak(:,2)-TimeStartW*Param.Fs,Speak(:,2),Speak(:,3)+TimeEndW*Param.Fs];
            WListenIdx=[Listen(:,2)-TimeStartW*Param.Fs,Listen(:,2),Listen(:,2)+TimeEndW*Param.Fs];
            
        %         EWSpeakIdx=[Speak(:,3)-TimeStartW*Param.Fs,Speak(:,3),Speak(:,3)+TimeEndW*Param.Fs];
        %         EWListenIdx=[Listen(:,3)-TimeStartW*Param.Fs,Listen(:,3),Listen(:,3)+TimeEndW*Param.Fs];

        %         t_Diam = linspace(0,length(Diameter)./Param.Fs,length(Diameter));

            % Storing Speaking/Listening by conditions
            if contains(cell2mat(FileNames_II(i)),'N0')
                for j=1:size(Speak,1)
                    SW_N0(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    SW_N0_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_N0_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    SDur_N0(i,j) = Speak(j,1);
                end
                for j=1:size(Listen,1)
                    LW_N0(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                    LW_N0_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_N0_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    LDur_N0(i,j) = Listen(j,1);
                end
            elseif contains(cell2mat(FileNames_II(i)),'N60')
                for j=1:size(Speak,1)
                    SW_N60(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    SW_N60_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_N60_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    SDur_N60(i,j) = Speak(j,1);
                end
                for j=1:size(Listen,1)
                    LW_N60(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                    LW_N60_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_N60_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    LDur_N60(i,j) = Listen(j,1);
                end
            elseif contains(cell2mat(FileNames_II(i)),'N70')
                for j=1:size(Speak,1)
                    SW_N70(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    SW_N70_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_N70_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    SDur_N70(i,j) = Speak(j,1);
                end
                for j=1:size(Listen,1)
                    LW_N70(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                    LW_N70_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_N70_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    LDur_N70(i,j) = Listen(j,1);
                end
            end
            % Either NH or HI
            if contains(ChosenFolder,'NH')
                for j=1:size(Speak,1)
                    SW_NH(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_NH)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    SW_NH_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_NH_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                end
                for j=1:size(Listen,1)
                    LW_NH(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_NH)-length(WListenIdx(j,1):Listen(j,3)))'];
                    LW_NH_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_NH_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                end
                if contains(cell2mat(FileNames_II(i)),'N0')
                    for j=1:size(Speak,1)
                        SW_NH_N0(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_NH_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        SW_NH_N0_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_NH_N0_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    end
                    for j=1:size(Listen,1)
                        LW_NH_N0(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_NH_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                        LW_NH_N0_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_NH_N0_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    end
                elseif contains(cell2mat(FileNames_II(i)),'N60')
                    for j=1:size(Speak,1)
                        SW_NH_N60(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_NH_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        SW_NH_N60_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_NH_N60_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    end
                    for j=1:size(Listen,1)
                        LW_NH_N60(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_NH_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                        LW_NH_N60_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_NH_N60_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    end
                elseif contains(cell2mat(FileNames_II(i)),'N70')
                    for j=1:size(Speak,1)
                        SW_NH_N70(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_NH_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        SW_NH_N70_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_NH_N70_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    end
                    for j=1:size(Listen,1)
                        LW_NH_N70(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_NH_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                        LW_NH_N70_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_NH_N70_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    end
                end
            elseif contains(ChosenFolder,'HI')
                for j=1:size(Speak,1)
                    SW_HI(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_HI)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    SW_HI_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_HI_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                end
                for j=1:size(Listen,1)
                    LW_HI(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_HI)-length(WListenIdx(j,1):Listen(j,3)))'];
                    LW_HI_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_HI_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                end
                if contains(cell2mat(FileNames_II(i)),'N0')
                    for j=1:size(Speak,1)
                        SW_HI_N0(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_HI_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        SW_HI_N0_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_HI_N0_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    end
                    for j=1:size(Listen,1)
                        LW_HI_N0(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_HI_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                        LW_HI_N0_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_HI_N0_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    end
                elseif contains(cell2mat(FileNames_II(i)),'N60')
                    for j=1:size(Speak,1)
                        SW_HI_N60(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_HI_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        SW_HI_N60_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_HI_N60_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    end
                    for j=1:size(Listen,1)
                        LW_HI_N60(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_HI_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                        LW_HI_N60_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_HI_N60_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    end
                elseif contains(cell2mat(FileNames_II(i)),'N70')
                    for j=1:size(Speak,1)
                        SW_HI_N70(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_HI_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        SW_HI_N70_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_HI_N70_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    end
                    for j=1:size(Listen,1)
                        LW_HI_N70(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_HI_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                        LW_HI_N70_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_HI_N70_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    end
                end
                if contains(cell2mat(FileNames_II(i)),'UN')
                    for j=1:size(Speak,1)
                        SW_HI_UN(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_HI_UN)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        SW_HI_UN_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_HI_UN_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    end
                    for j=1:size(Listen,1)
                        LW_HI_UN(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_HI_UN)-length(WListenIdx(j,1):Listen(j,3)))'];
                        LW_HI_UN_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_HI_UN_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    end
                    if contains(cell2mat(FileNames_II(i)),'N0')
                        for j=1:size(Speak,1)
                            SW_HI_UN_N0(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_HI_UN_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            SW_HI_UN_N0_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_HI_UN_N0_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            LW_HI_UN_N0(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_HI_UN_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                            LW_HI_UN_N0_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_HI_UN_N0_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    elseif contains(cell2mat(FileNames_II(i)),'N60')
                        for j=1:size(Speak,1)
                            SW_HI_UN_N60(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_HI_UN_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            SW_HI_UN_N60_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_HI_UN_N60_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            LW_HI_UN_N60(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_HI_UN_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                            LW_HI_UN_N60_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_HI_UN_N60_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    elseif contains(cell2mat(FileNames_II(i)),'N70')
                        for j=1:size(Speak,1)
                            SW_HI_UN_N70(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_HI_UN_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            SW_HI_UN_N70_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_HI_UN_N70_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            LW_HI_UN_N70(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_HI_UN_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                            LW_HI_UN_N70_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_HI_UN_N70_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    end
                elseif contains(cell2mat(FileNames_II(i)),'AA')
                    for j=1:size(Speak,1)
                        SW_HI_AA(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_HI_AA)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        SW_HI_AA_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_HI_AA_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    end
                    for j=1:size(Listen,1)
                        LW_HI_AA(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_HI_AA)-length(WListenIdx(j,1):Listen(j,3)))'];
                        LW_HI_AA_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_HI_AA_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    end
                    if contains(cell2mat(FileNames_II(i)),'N0')
                        for j=1:size(Speak,1)
                            SW_HI_AA_N0(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_HI_AA_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            SW_HI_AA_N0_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_HI_AA_N0_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            LW_HI_AA_N0(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_HI_AA_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                            LW_HI_AA_N0_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_HI_AA_N0_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    elseif contains(cell2mat(FileNames_II(i)),'N60')
                        for j=1:size(Speak,1)
                            SW_HI_AA_N60(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_HI_AA_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            SW_HI_AA_N60_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_HI_AA_N60_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            LW_HI_AA_N60(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_HI_AA_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                            LW_HI_AA_N60_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_HI_AA_N60_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    elseif contains(cell2mat(FileNames_II(i)),'N70')
                        for j=1:size(Speak,1)
                            SW_HI_AA_N70(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_HI_AA_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            SW_HI_AA_N70_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_HI_AA_N70_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            LW_HI_AA_N70(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_HI_AA_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                            LW_HI_AA_N70_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_HI_AA_N70_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    end
                elseif contains(cell2mat(FileNames_II(i)),'AB')
                    for j=1:size(Speak,1)
                        SW_HI_AB(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_HI_AB)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        SW_HI_AB_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_HI_AB_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    end
                    for j=1:size(Listen,1)
                        LW_HI_AB(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_HI_AB)-length(WListenIdx(j,1):Listen(j,3)))'];
                        LW_HI_AB_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_HI_AB_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    end
                    if contains(cell2mat(FileNames_II(i)),'N0')
                        for j=1:size(Speak,1)
                            SW_HI_AB_N0(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_HI_AB_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            SW_HI_AB_N0_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_HI_AB_N0_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            LW_HI_AB_N0(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_HI_AB_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                            LW_HI_AB_N0_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_HI_AB_N0_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    elseif contains(cell2mat(FileNames_II(i)),'N60')
                        for j=1:size(Speak,1)
                            SW_HI_AB_N60(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_HI_AB_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            SW_HI_AB_N60_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_HI_AB_N60_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            LW_HI_AB_N60(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_HI_AB_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                            LW_HI_AB_N60_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_HI_AB_N60_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    elseif contains(cell2mat(FileNames_II(i)),'N70')
                        for j=1:size(Speak,1)
                            SW_HI_AB_N70(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_HI_AB_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            SW_HI_AB_N70_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_HI_AB_N70_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            LW_HI_AB_N70(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_HI_AB_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                            LW_HI_AB_N70_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_HI_AB_N70_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    end
                end
            end

            % Storing Global Speaking/Listening windows
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
%             if contains(cell2mat(FileNames_II(i)),'P2')
%                 TPsOrder(2*q,i-NCond) = x;
%             elseif contains(cell2mat(FileNames_II(i)),'P1')
%                 TPsOrder(2*q-1,i) = x;
%             end

            % Increase index of num of files used
            x=x+1;
            
        end
    end
end
%% Out-of-files-loop calculations
% Clean empty rows and layers, set 0's to NaN
GSW(~any(GSW,[2 3]),:,:)=[];GSW(:,~any(GSW,[1 3]),:)=[];GSW(GSW==0)=NaN;
GSW_B(~any(GSW_B,[2 3]),:,:)=[];GSW_B(:,~any(GSW_B,[1 3]),:)=[];GSW_B(GSW_B==0)=NaN;
GLW(~any(GLW,[2 3]),:,:)=[];GLW(:,~any(GLW,[1 3]),:)=[];GLW(GLW==0)=NaN;
GLW_B(~any(GLW_B,[2 3]),:,:)=[];GLW_B(:,~any(GLW_B,[1 3]),:)=[];GLW_B(GLW_B==0)=NaN;

SW_NH(~any(SW_NH,[2 3]),:,:)=[];SW_NH(:,~any(SW_NH,[1 3]),:)=[];SW_NH(SW_NH==0)=NaN;
SW_NH_B(~any(SW_NH_B,[2 3]),:,:)=[];SW_NH_B(:,~any(SW_NH_B,[1 3]),:)=[];SW_NH_B(SW_NH_B==0)=NaN;
LW_NH(~any(LW_NH,[2 3]),:,:)=[];LW_NH(:,~any(LW_NH,[1 3]),:)=[];LW_NH(LW_NH==0)=NaN;
LW_NH_B(~any(LW_NH_B,[2 3]),:,:)=[];LW_NH_B(:,~any(LW_NH_B,[1 3]),:)=[];LW_NH_B(LW_NH_B==0)=NaN;
SW_HI(~any(SW_HI,[2 3]),:,:)=[];SW_HI(:,~any(SW_HI,[1 3]),:)=[];SW_HI(SW_HI==0)=NaN;
SW_HI_B(~any(SW_HI_B,[2 3]),:,:)=[];SW_HI_B(:,~any(SW_HI_B,[1 3]),:)=[];SW_HI_B(SW_HI_B==0)=NaN;
LW_HI(~any(LW_HI,[2 3]),:,:)=[];LW_HI(:,~any(LW_HI,[1 3]),:)=[];LW_HI(LW_HI==0)=NaN;
LW_HI_B(~any(LW_HI_B,[2 3]),:,:)=[];LW_HI_B(:,~any(LW_HI_B,[1 3]),:)=[];LW_HI_B(LW_HI_B==0)=NaN;

SW_N0(~any(SW_N0,[2 3]),:,:)=[];SW_N0(:,~any(SW_N0,[1 3]),:)=[];SW_N0(SW_N0==0)=NaN;
SW_N0_B(~any(SW_N0_B,[2 3]),:,:)=[];SW_N0_B(:,~any(SW_N0_B,[1 3]),:)=[];SW_N0_B(SW_N0_B==0)=NaN;
LW_N0(~any(LW_N0,[2 3]),:,:)=[];LW_N0(:,~any(LW_N0,[1 3]),:)=[];LW_N0(LW_N0==0)=NaN;
LW_N0_B(~any(LW_N0_B,[2 3]),:,:)=[];LW_N0_B(:,~any(LW_N0_B,[1 3]),:)=[];LW_N0_B(LW_N0_B==0)=NaN;
SW_N60(~any(SW_N60,[2 3]),:,:)=[];SW_N60(:,~any(SW_N60,[1 3]),:)=[];SW_N60(SW_N60==0)=NaN;
SW_N60_B(~any(SW_N60_B,[2 3]),:,:)=[];SW_N60_B(:,~any(SW_N60_B,[1 3]),:)=[];SW_N60_B(SW_N60_B==0)=NaN;
LW_N60(~any(LW_N60,[2 3]),:,:)=[];LW_N60(:,~any(LW_N60,[1 3]),:)=[];LW_N60(LW_N60==0)=NaN;
LW_N60_B(~any(LW_N60_B,[2 3]),:,:)=[];LW_N60_B(:,~any(LW_N60_B,[1 3]),:)=[];LW_N60_B(LW_N60_B==0)=NaN;
SW_N70(~any(SW_N70,[2 3]),:,:)=[];SW_N70(:,~any(SW_N70,[1 3]),:)=[];SW_N70(SW_N70==0)=NaN;
SW_N70_B(~any(SW_N70_B,[2 3]),:,:)=[];SW_N70_B(:,~any(SW_N70_B,[1 3]),:)=[];SW_N70_B(SW_N70_B==0)=NaN;
LW_N70(~any(LW_N70,[2 3]),:,:)=[];LW_N70(:,~any(LW_N70,[1 3]),:)=[];LW_N70(LW_N70==0)=NaN;
LW_N70_B(~any(LW_N70_B,[2 3]),:,:)=[];LW_N70_B(:,~any(LW_N70_B,[1 3]),:)=[];LW_N70_B(LW_N70_B==0)=NaN;

SW_NH_N0(~any(SW_NH_N0,[2 3]),:,:)=[];SW_NH_N0(:,~any(SW_NH_N0,[1 3]),:)=[];SW_NH_N0(SW_NH_N0==0)=NaN;
SW_NH_N0_B(~any(SW_NH_N0_B,[2 3]),:,:)=[];SW_NH_N0_B(:,~any(SW_NH_N0_B,[1 3]),:)=[];SW_NH_N0_B(SW_NH_N0_B==0)=NaN;
LW_NH_N0(~any(LW_NH_N0,[2 3]),:,:)=[];LW_NH_N0(:,~any(LW_NH_N0,[1 3]),:)=[];LW_NH_N0(LW_NH_N0==0)=NaN;
LW_NH_N0_B(~any(LW_NH_N0_B,[2 3]),:,:)=[];LW_NH_N0_B(:,~any(LW_NH_N0_B,[1 3]),:)=[];LW_NH_N0_B(LW_NH_N0_B==0)=NaN;
SW_NH_N60(~any(SW_NH_N60,[2 3]),:,:)=[];SW_NH_N60(:,~any(SW_NH_N60,[1 3]),:)=[];SW_NH_N60(SW_NH_N60==0)=NaN;
SW_NH_N60_B(~any(SW_NH_N60_B,[2 3]),:,:)=[];SW_NH_N60_B(:,~any(SW_NH_N60_B,[1 3]),:)=[];SW_NH_N60_B(SW_NH_N60_B==0)=NaN;
LW_NH_N60(~any(LW_NH_N60,[2 3]),:,:)=[];LW_NH_N60(:,~any(LW_NH_N60,[1 3]),:)=[];LW_NH_N60(LW_NH_N60==0)=NaN;
LW_NH_N60_B(~any(LW_NH_N60_B,[2 3]),:,:)=[];LW_NH_N60_B(:,~any(LW_NH_N60_B,[1 3]),:)=[];LW_NH_N60_B(LW_NH_N60_B==0)=NaN;
SW_NH_N70(~any(SW_NH_N70,[2 3]),:,:)=[];SW_NH_N70(:,~any(SW_NH_N70,[1 3]),:)=[];SW_NH_N70(SW_NH_N70==0)=NaN;
SW_NH_N70_B(~any(SW_NH_N70_B,[2 3]),:,:)=[];SW_NH_N70_B(:,~any(SW_NH_N70_B,[1 3]),:)=[];SW_NH_N70_B(SW_NH_N70_B==0)=NaN;
LW_NH_N70(~any(LW_NH_N70,[2 3]),:,:)=[];LW_NH_N70(:,~any(LW_NH_N70,[1 3]),:)=[];LW_NH_N70(LW_NH_N70==0)=NaN;
LW_NH_N70_B(~any(LW_NH_N70_B,[2 3]),:,:)=[];LW_NH_N70_B(:,~any(LW_NH_N70_B,[1 3]),:)=[];LW_NH_N70_B(LW_NH_N70_B==0)=NaN;

SW_HI_N0(~any(SW_HI_N0,[2 3]),:,:)=[];SW_HI_N0(:,~any(SW_HI_N0,[1 3]),:)=[];SW_HI_N0(SW_HI_N0==0)=NaN;
SW_HI_N0_B(~any(SW_HI_N0_B,[2 3]),:,:)=[];SW_HI_N0_B(:,~any(SW_HI_N0_B,[1 3]),:)=[];SW_HI_N0_B(SW_HI_N0_B==0)=NaN;
LW_HI_N0(~any(LW_HI_N0,[2 3]),:,:)=[];LW_HI_N0(:,~any(LW_HI_N0,[1 3]),:)=[];LW_HI_N0(LW_HI_N0==0)=NaN;
LW_HI_N0_B(~any(LW_HI_N0_B,[2 3]),:,:)=[];LW_HI_N0_B(:,~any(LW_HI_N0_B,[1 3]),:)=[];LW_HI_N0_B(LW_HI_N0_B==0)=NaN;
SW_HI_N60(~any(SW_HI_N60,[2 3]),:,:)=[];SW_HI_N60(:,~any(SW_HI_N60,[1 3]),:)=[];SW_HI_N60(SW_HI_N60==0)=NaN;
SW_HI_N60_B(~any(SW_HI_N60_B,[2 3]),:,:)=[];SW_HI_N60_B(:,~any(SW_HI_N60_B,[1 3]),:)=[];SW_HI_N60_B(SW_HI_N60_B==0)=NaN;
LW_HI_N60(~any(LW_HI_N60,[2 3]),:,:)=[];LW_HI_N60(:,~any(LW_HI_N60,[1 3]),:)=[];LW_HI_N60(LW_HI_N60==0)=NaN;
LW_HI_N60_B(~any(LW_HI_N60_B,[2 3]),:,:)=[];LW_HI_N60_B(:,~any(LW_HI_N60_B,[1 3]),:)=[];LW_HI_N60_B(LW_HI_N60_B==0)=NaN;
SW_HI_N70(~any(SW_HI_N70,[2 3]),:,:)=[];SW_HI_N70(:,~any(SW_HI_N70,[1 3]),:)=[];SW_HI_N70(SW_HI_N70==0)=NaN;
SW_HI_N70_B(~any(SW_HI_N70_B,[2 3]),:,:)=[];SW_HI_N70_B(:,~any(SW_HI_N70_B,[1 3]),:)=[];SW_HI_N70_B(SW_HI_N70_B==0)=NaN;
LW_HI_N70(~any(LW_HI_N70,[2 3]),:,:)=[];LW_HI_N70(:,~any(LW_HI_N70,[1 3]),:)=[];LW_HI_N70(LW_HI_N70==0)=NaN;
LW_HI_N70_B(~any(LW_HI_N70_B,[2 3]),:,:)=[];LW_HI_N70_B(:,~any(LW_HI_N70_B,[1 3]),:)=[];LW_HI_N70_B(LW_HI_N70_B==0)=NaN;

SW_HI_UN(~any(SW_HI_UN,[2 3]),:,:)=[];SW_HI_UN(:,~any(SW_HI_UN,[1 3]),:)=[];SW_HI_UN(SW_HI_UN==0)=NaN;
SW_HI_UN_B(~any(SW_HI_UN_B,[2 3]),:,:)=[];SW_HI_UN_B(:,~any(SW_HI_UN_B,[1 3]),:)=[];SW_HI_UN_B(SW_HI_UN_B==0)=NaN;
LW_HI_UN(~any(LW_HI_UN,[2 3]),:,:)=[];LW_HI_UN(:,~any(LW_HI_UN,[1 3]),:)=[];LW_HI_UN(LW_HI_UN==0)=NaN;
LW_HI_UN_B(~any(LW_HI_UN_B,[2 3]),:,:)=[];LW_HI_UN_B(:,~any(LW_HI_UN_B,[1 3]),:)=[];LW_HI_UN_B(LW_HI_UN_B==0)=NaN;
SW_HI_AA(~any(SW_HI_AA,[2 3]),:,:)=[];SW_HI_AA(:,~any(SW_HI_AA,[1 3]),:)=[];SW_HI_AA(SW_HI_AA==0)=NaN;
SW_HI_AA_B(~any(SW_HI_AA_B,[2 3]),:,:)=[];SW_HI_AA_B(:,~any(SW_HI_AA_B,[1 3]),:)=[];SW_HI_AA_B(SW_HI_AA_B==0)=NaN;
LW_HI_AA(~any(LW_HI_AA,[2 3]),:,:)=[];LW_HI_AA(:,~any(LW_HI_AA,[1 3]),:)=[];LW_HI_AA(LW_HI_AA==0)=NaN;
LW_HI_AA_B(~any(LW_HI_AA_B,[2 3]),:,:)=[];LW_HI_AA_B(:,~any(LW_HI_AA_B,[1 3]),:)=[];LW_HI_AA_B(LW_HI_AA_B==0)=NaN;
SW_HI_AB(~any(SW_HI_AB,[2 3]),:,:)=[];SW_HI_AB(:,~any(SW_HI_AB,[1 3]),:)=[];SW_HI_AB(SW_HI_AB==0)=NaN;
SW_HI_AB_B(~any(SW_HI_AB_B,[2 3]),:,:)=[];SW_HI_AB_B(:,~any(SW_HI_AB_B,[1 3]),:)=[];SW_HI_AB_B(SW_HI_AB_B==0)=NaN;
LW_HI_AB(~any(LW_HI_AB,[2 3]),:,:)=[];LW_HI_AB(:,~any(LW_HI_AB,[1 3]),:)=[];LW_HI_AB(LW_HI_AB==0)=NaN;
LW_HI_AB_B(~any(LW_HI_AB_B,[2 3]),:,:)=[];LW_HI_AB_B(:,~any(LW_HI_AB_B,[1 3]),:)=[];LW_HI_AB_B(LW_HI_AB_B==0)=NaN;

SW_HI_UN_N0(~any(SW_HI_UN_N0,[2 3]),:,:)=[];SW_HI_UN_N0(:,~any(SW_HI_UN_N0,[1 3]),:)=[];SW_HI_UN_N0(SW_HI_UN_N0==0)=NaN;
SW_HI_UN_N0_B(~any(SW_HI_UN_N0_B,[2 3]),:,:)=[];SW_HI_UN_N0_B(:,~any(SW_HI_UN_N0_B,[1 3]),:)=[];SW_HI_UN_N0_B(SW_HI_UN_N0_B==0)=NaN;
LW_HI_UN_N0(~any(LW_HI_UN_N0,[2 3]),:,:)=[];LW_HI_UN_N0(:,~any(LW_HI_UN_N0,[1 3]),:)=[];LW_HI_UN_N0(LW_HI_UN_N0==0)=NaN;
LW_HI_UN_N0_B(~any(LW_HI_UN_N0_B,[2 3]),:,:)=[];LW_HI_UN_N0_B(:,~any(LW_HI_UN_N0_B,[1 3]),:)=[];LW_HI_UN_N0_B(LW_HI_UN_N0_B==0)=NaN;
SW_HI_UN_N60(~any(SW_HI_UN_N60,[2 3]),:,:)=[];SW_HI_UN_N60(:,~any(SW_HI_UN_N60,[1 3]),:)=[];SW_HI_UN_N60(SW_HI_UN_N60==0)=NaN;
SW_HI_UN_N60_B(~any(SW_HI_UN_N60_B,[2 3]),:,:)=[];SW_HI_UN_N60_B(:,~any(SW_HI_UN_N60_B,[1 3]),:)=[];SW_HI_UN_N60_B(SW_HI_UN_N60_B==0)=NaN;
LW_HI_UN_N60(~any(LW_HI_UN_N60,[2 3]),:,:)=[];LW_HI_UN_N60(:,~any(LW_HI_UN_N60,[1 3]),:)=[];LW_HI_UN_N60(LW_HI_UN_N60==0)=NaN;
LW_HI_UN_N60_B(~any(LW_HI_UN_N60_B,[2 3]),:,:)=[];LW_HI_UN_N60_B(:,~any(LW_HI_UN_N60_B,[1 3]),:)=[];LW_HI_UN_N60_B(LW_HI_UN_N60_B==0)=NaN;
SW_HI_UN_N70(~any(SW_HI_UN_N70,[2 3]),:,:)=[];SW_HI_UN_N70(:,~any(SW_HI_UN_N70,[1 3]),:)=[];SW_HI_UN_N70(SW_HI_UN_N70==0)=NaN;
SW_HI_UN_N70_B(~any(SW_HI_UN_N70_B,[2 3]),:,:)=[];SW_HI_UN_N70_B(:,~any(SW_HI_UN_N70_B,[1 3]),:)=[];SW_HI_UN_N70_B(SW_HI_UN_N70_B==0)=NaN;
LW_HI_UN_N70(~any(LW_HI_UN_N70,[2 3]),:,:)=[];LW_HI_UN_N70(:,~any(LW_HI_UN_N70,[1 3]),:)=[];LW_HI_UN_N70(LW_HI_UN_N70==0)=NaN;
LW_HI_UN_N70_B(~any(LW_HI_UN_N70_B,[2 3]),:,:)=[];LW_HI_UN_N70_B(:,~any(LW_HI_UN_N70_B,[1 3]),:)=[];LW_HI_UN_N70_B(LW_HI_UN_N70_B==0)=NaN;
SW_HI_AA_N0(~any(SW_HI_AA_N0,[2 3]),:,:)=[];SW_HI_AA_N0(:,~any(SW_HI_AA_N0,[1 3]),:)=[];SW_HI_AA_N0(SW_HI_AA_N0==0)=NaN;
SW_HI_AA_N0_B(~any(SW_HI_AA_N0_B,[2 3]),:,:)=[];SW_HI_AA_N0_B(:,~any(SW_HI_AA_N0_B,[1 3]),:)=[];SW_HI_AA_N0_B(SW_HI_AA_N0_B==0)=NaN;
LW_HI_AA_N0(~any(LW_HI_AA_N0,[2 3]),:,:)=[];LW_HI_AA_N0(:,~any(LW_HI_AA_N0,[1 3]),:)=[];LW_HI_AA_N0(LW_HI_AA_N0==0)=NaN;
LW_HI_AA_N0_B(~any(LW_HI_AA_N0_B,[2 3]),:,:)=[];LW_HI_AA_N0_B(:,~any(LW_HI_AA_N0_B,[1 3]),:)=[];LW_HI_AA_N0_B(LW_HI_AA_N0_B==0)=NaN;
SW_HI_AA_N60(~any(SW_HI_AA_N60,[2 3]),:,:)=[];SW_HI_AA_N60(:,~any(SW_HI_AA_N60,[1 3]),:)=[];SW_HI_AA_N60(SW_HI_AA_N60==0)=NaN;
SW_HI_AA_N60_B(~any(SW_HI_AA_N60_B,[2 3]),:,:)=[];SW_HI_AA_N60_B(:,~any(SW_HI_AA_N60_B,[1 3]),:)=[];SW_HI_AA_N60_B(SW_HI_AA_N60_B==0)=NaN;
LW_HI_AA_N60(~any(LW_HI_AA_N60,[2 3]),:,:)=[];LW_HI_AA_N60(:,~any(LW_HI_AA_N60,[1 3]),:)=[];LW_HI_AA_N60(LW_HI_AA_N60==0)=NaN;
LW_HI_AA_N60_B(~any(LW_HI_AA_N60_B,[2 3]),:,:)=[];LW_HI_AA_N60_B(:,~any(LW_HI_AA_N60_B,[1 3]),:)=[];LW_HI_AA_N60_B(LW_HI_AA_N60_B==0)=NaN;
SW_HI_AA_N70(~any(SW_HI_AA_N70,[2 3]),:,:)=[];SW_HI_AA_N70(:,~any(SW_HI_AA_N70,[1 3]),:)=[];SW_HI_AA_N70(SW_HI_AA_N70==0)=NaN;
SW_HI_AA_N70_B(~any(SW_HI_AA_N70_B,[2 3]),:,:)=[];SW_HI_AA_N70_B(:,~any(SW_HI_AA_N70_B,[1 3]),:)=[];SW_HI_AA_N70_B(SW_HI_AA_N70_B==0)=NaN;
LW_HI_AA_N70(~any(LW_HI_AA_N70,[2 3]),:,:)=[];LW_HI_AA_N70(:,~any(LW_HI_AA_N70,[1 3]),:)=[];LW_HI_AA_N70(LW_HI_AA_N70==0)=NaN;
LW_HI_AA_N70_B(~any(LW_HI_AA_N70_B,[2 3]),:,:)=[];LW_HI_AA_N70_B(:,~any(LW_HI_AA_N70_B,[1 3]),:)=[];LW_HI_AA_N70_B(LW_HI_AA_N70_B==0)=NaN;
SW_HI_AB_N0(~any(SW_HI_AB_N0,[2 3]),:,:)=[];SW_HI_AB_N0(:,~any(SW_HI_AB_N0,[1 3]),:)=[];SW_HI_AB_N0(SW_HI_AB_N0==0)=NaN;
SW_HI_AB_N0_B(~any(SW_HI_AB_N0_B,[2 3]),:,:)=[];SW_HI_AB_N0_B(:,~any(SW_HI_AB_N0_B,[1 3]),:)=[];SW_HI_AB_N0_B(SW_HI_AB_N0_B==0)=NaN;
LW_HI_AB_N0(~any(LW_HI_AB_N0,[2 3]),:,:)=[];LW_HI_AB_N0(:,~any(LW_HI_AB_N0,[1 3]),:)=[];LW_HI_AB_N0(LW_HI_AB_N0==0)=NaN;
LW_HI_AB_N0_B(~any(LW_HI_AB_N0_B,[2 3]),:,:)=[];LW_HI_AB_N0_B(:,~any(LW_HI_AB_N0_B,[1 3]),:)=[];LW_HI_AB_N0_B(LW_HI_AB_N0_B==0)=NaN;
SW_HI_AB_N60(~any(SW_HI_AB_N60,[2 3]),:,:)=[];SW_HI_AB_N60(:,~any(SW_HI_AB_N60,[1 3]),:)=[];SW_HI_AB_N60(SW_HI_AB_N60==0)=NaN;
SW_HI_AB_N60_B(~any(SW_HI_AB_N60_B,[2 3]),:,:)=[];SW_HI_AB_N60_B(:,~any(SW_HI_AB_N60_B,[1 3]),:)=[];SW_HI_AB_N60_B(SW_HI_AB_N60_B==0)=NaN;
LW_HI_AB_N60(~any(LW_HI_AB_N60,[2 3]),:,:)=[];LW_HI_AB_N60(:,~any(LW_HI_AB_N60,[1 3]),:)=[];LW_HI_AB_N60(LW_HI_AB_N60==0)=NaN;
LW_HI_AB_N60_B(~any(LW_HI_AB_N60_B,[2 3]),:,:)=[];LW_HI_AB_N60_B(:,~any(LW_HI_AB_N60_B,[1 3]),:)=[];LW_HI_AB_N60_B(LW_HI_AB_N60_B==0)=NaN;
SW_HI_AB_N70(~any(SW_HI_AB_N70,[2 3]),:,:)=[];SW_HI_AB_N70(:,~any(SW_HI_AB_N70,[1 3]),:)=[];SW_HI_AB_N70(SW_HI_AB_N70==0)=NaN;
SW_HI_AB_N70_B(~any(SW_HI_AB_N70_B,[2 3]),:,:)=[];SW_HI_AB_N70_B(:,~any(SW_HI_AB_N70_B,[1 3]),:)=[];SW_HI_AB_N70_B(SW_HI_AB_N70_B==0)=NaN;
LW_HI_AB_N70(~any(LW_HI_AB_N70,[2 3]),:,:)=[];LW_HI_AB_N70(:,~any(LW_HI_AB_N70,[1 3]),:)=[];LW_HI_AB_N70(LW_HI_AB_N70==0)=NaN;
LW_HI_AB_N70_B(~any(LW_HI_AB_N70_B,[2 3]),:,:)=[];LW_HI_AB_N70_B(:,~any(LW_HI_AB_N70_B,[1 3]),:)=[];LW_HI_AB_N70_B(LW_HI_AB_N70_B==0)=NaN;

% Calculate means omitting NaNs, LP filtering with mean-padding at start
% end of each group of Speaking/Listening windows
GSW_Mean = ndnanfilter(reshape(mean(GSW,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
GSW_B_Mean = ndnanfilter(reshape(mean(GSW_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
GLW_Mean = ndnanfilter(reshape(mean(GLW,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
GLW_B_Mean = ndnanfilter(reshape(mean(GLW_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

SW_NH_Mean = ndnanfilter(reshape(mean(SW_NH,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_NH_B_Mean = ndnanfilter(reshape(mean(SW_NH_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_NH_Mean = ndnanfilter(reshape(mean(LW_NH,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_NH_B_Mean = ndnanfilter(reshape(mean(LW_NH_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_Mean = ndnanfilter(reshape(mean(SW_HI,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_B_Mean = ndnanfilter(reshape(mean(SW_HI_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_Mean = ndnanfilter(reshape(mean(LW_HI,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_B_Mean = ndnanfilter(reshape(mean(LW_HI_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

SW_N0_Mean = ndnanfilter(reshape(mean(SW_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_N0_B_Mean = ndnanfilter(reshape(mean(SW_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_N0_Mean = ndnanfilter(reshape(mean(LW_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_N0_B_Mean = ndnanfilter(reshape(mean(LW_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_N60_Mean = ndnanfilter(reshape(mean(SW_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_N60_B_Mean = ndnanfilter(reshape(mean(SW_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_N60_Mean = ndnanfilter(reshape(mean(LW_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_N60_B_Mean = ndnanfilter(reshape(mean(LW_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_N70_Mean = ndnanfilter(reshape(mean(SW_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_N70_B_Mean = ndnanfilter(reshape(mean(SW_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_N70_Mean = ndnanfilter(reshape(mean(LW_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_N70_B_Mean = ndnanfilter(reshape(mean(LW_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

SW_NH_N0_Mean = ndnanfilter(reshape(mean(SW_NH_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_NH_N0_B_Mean = ndnanfilter(reshape(mean(SW_NH_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_NH_N0_Mean = ndnanfilter(reshape(mean(LW_NH_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_NH_N0_B_Mean = ndnanfilter(reshape(mean(LW_NH_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_NH_N60_Mean = ndnanfilter(reshape(mean(SW_NH_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_NH_N60_B_Mean = ndnanfilter(reshape(mean(SW_NH_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_NH_N60_Mean = ndnanfilter(reshape(mean(LW_NH_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_NH_N60_B_Mean = ndnanfilter(reshape(mean(LW_NH_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_NH_N70_Mean = ndnanfilter(reshape(mean(SW_NH_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_NH_N70_B_Mean = ndnanfilter(reshape(mean(SW_NH_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_NH_N70_Mean = ndnanfilter(reshape(mean(LW_NH_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_NH_N70_B_Mean = ndnanfilter(reshape(mean(LW_NH_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

SW_HI_N0_Mean = ndnanfilter(reshape(mean(SW_HI_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_N0_B_Mean = ndnanfilter(reshape(mean(SW_HI_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_N0_Mean = ndnanfilter(reshape(mean(LW_HI_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_N0_B_Mean = ndnanfilter(reshape(mean(LW_HI_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_N60_Mean = ndnanfilter(reshape(mean(SW_HI_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_N60_B_Mean = ndnanfilter(reshape(mean(SW_HI_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_N60_Mean = ndnanfilter(reshape(mean(LW_HI_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_N60_B_Mean = ndnanfilter(reshape(mean(LW_HI_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_N70_Mean = ndnanfilter(reshape(mean(SW_HI_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_N70_B_Mean = ndnanfilter(reshape(mean(SW_HI_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_N70_Mean = ndnanfilter(reshape(mean(LW_HI_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_N70_B_Mean = ndnanfilter(reshape(mean(LW_HI_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

SW_HI_UN_Mean = ndnanfilter(reshape(mean(SW_HI_UN,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_UN_B_Mean = ndnanfilter(reshape(mean(SW_HI_UN_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_UN_Mean = ndnanfilter(reshape(mean(LW_HI_UN,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_UN_B_Mean = ndnanfilter(reshape(mean(LW_HI_UN_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_AA_Mean = ndnanfilter(reshape(mean(SW_HI_AA,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_AA_B_Mean = ndnanfilter(reshape(mean(SW_HI_AA_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_AA_Mean = ndnanfilter(reshape(mean(LW_HI_AA,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_AA_B_Mean = ndnanfilter(reshape(mean(LW_HI_AA_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_AB_Mean = ndnanfilter(reshape(mean(SW_HI_AB,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_AB_B_Mean = ndnanfilter(reshape(mean(SW_HI_AB_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_AB_Mean = ndnanfilter(reshape(mean(LW_HI_AB,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_AB_B_Mean = ndnanfilter(reshape(mean(LW_HI_AB_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

SW_HI_UN_N0_Mean = ndnanfilter(reshape(mean(SW_HI_UN_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_UN_N0_B_Mean = ndnanfilter(reshape(mean(SW_HI_UN_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_UN_N0_Mean = ndnanfilter(reshape(mean(LW_HI_UN_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_UN_N0_B_Mean = ndnanfilter(reshape(mean(LW_HI_UN_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_UN_N60_Mean = ndnanfilter(reshape(mean(SW_HI_UN_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_UN_N60_B_Mean = ndnanfilter(reshape(mean(SW_HI_UN_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_UN_N60_Mean = ndnanfilter(reshape(mean(LW_HI_UN_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_UN_N60_B_Mean = ndnanfilter(reshape(mean(LW_HI_UN_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_UN_N70_Mean = ndnanfilter(reshape(mean(SW_HI_UN_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_UN_N70_B_Mean = ndnanfilter(reshape(mean(SW_HI_UN_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_UN_N70_Mean = ndnanfilter(reshape(mean(LW_HI_UN_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_UN_N70_B_Mean = ndnanfilter(reshape(mean(LW_HI_UN_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_AA_N0_Mean = ndnanfilter(reshape(mean(SW_HI_AA_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_AA_N0_B_Mean = ndnanfilter(reshape(mean(SW_HI_AA_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_AA_N0_Mean = ndnanfilter(reshape(mean(LW_HI_AA_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_AA_N0_B_Mean = ndnanfilter(reshape(mean(LW_HI_AA_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_AA_N60_Mean = ndnanfilter(reshape(mean(SW_HI_AA_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_AA_N60_B_Mean = ndnanfilter(reshape(mean(SW_HI_AA_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_AA_N60_Mean = ndnanfilter(reshape(mean(LW_HI_AA_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_AA_N60_B_Mean = ndnanfilter(reshape(mean(LW_HI_AA_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_AA_N70_Mean = ndnanfilter(reshape(mean(SW_HI_AA_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_AA_N70_B_Mean = ndnanfilter(reshape(mean(SW_HI_AA_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_AA_N70_Mean = ndnanfilter(reshape(mean(LW_HI_AA_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_AA_N70_B_Mean = ndnanfilter(reshape(mean(LW_HI_AA_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_AB_N0_Mean = ndnanfilter(reshape(mean(SW_HI_AB_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_AB_N0_B_Mean = ndnanfilter(reshape(mean(SW_HI_AB_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_AB_N0_Mean = ndnanfilter(reshape(mean(LW_HI_AB_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_AB_N0_B_Mean = ndnanfilter(reshape(mean(LW_HI_AB_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_AB_N60_Mean = ndnanfilter(reshape(mean(SW_HI_AB_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_AB_N60_B_Mean = ndnanfilter(reshape(mean(SW_HI_AB_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_AB_N60_Mean = ndnanfilter(reshape(mean(LW_HI_AB_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_AB_N60_B_Mean = ndnanfilter(reshape(mean(LW_HI_AB_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_AB_N70_Mean = ndnanfilter(reshape(mean(SW_HI_AB_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_HI_AB_N70_B_Mean = ndnanfilter(reshape(mean(SW_HI_AB_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_AB_N70_Mean = ndnanfilter(reshape(mean(LW_HI_AB_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_HI_AB_N70_B_Mean = ndnanfilter(reshape(mean(LW_HI_AB_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

% Calculate SEM as: std(X)/sqrt(squeeze(sum(~isnan(X),[1 2]))))
GSW_SEM = (squeeze(std(GSW,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(GSW),[1 2]))))';
GSW_B_SEM = (squeeze(std(GSW_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(GSW_B),[1 2]))))';
GLW_SEM = (squeeze(std(GLW,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(GLW),[1 2]))))';
GLW_B_SEM = (squeeze(std(GLW_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(GLW_B),[1 2]))))';

SW_NH_SEM = (squeeze(std(SW_NH,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_NH),[1 2]))))';
SW_NH_B_SEM = (squeeze(std(SW_NH_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_NH_B),[1 2]))))';
LW_NH_SEM = (squeeze(std(LW_NH,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_NH),[1 2]))))';
LW_NH_B_SEM = (squeeze(std(LW_NH_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_NH_B),[1 2]))))';
SW_HI_SEM = (squeeze(std(SW_HI,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI),[1 2]))))';
SW_HI_B_SEM = (squeeze(std(SW_HI_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_B),[1 2]))))';
LW_HI_SEM = (squeeze(std(LW_HI,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI),[1 2]))))';
LW_HI_B_SEM = (squeeze(std(LW_HI_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_B),[1 2]))))';

SW_N0_SEM = (squeeze(std(SW_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_N0),[1 2]))))';
SW_N0_B_SEM = (squeeze(std(SW_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_N0_B),[1 2]))))';
LW_N0_SEM = (squeeze(std(LW_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_N0),[1 2]))))';
LW_N0_B_SEM = (squeeze(std(LW_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_N0_B),[1 2]))))';
SW_N60_SEM = (squeeze(std(SW_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_N60),[1 2]))))';
SW_N60_B_SEM = (squeeze(std(SW_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_N60_B),[1 2]))))';
LW_N60_SEM = (squeeze(std(LW_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_N60),[1 2]))))';
LW_N60_B_SEM = (squeeze(std(LW_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_N60_B),[1 2]))))';
SW_N70_SEM = (squeeze(std(SW_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_N70),[1 2]))))';
SW_N70_B_SEM = (squeeze(std(SW_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_N70_B),[1 2]))))';
LW_N70_SEM = (squeeze(std(LW_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_N70),[1 2]))))';
LW_N70_B_SEM = (squeeze(std(LW_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_N70_B),[1 2]))))';

SW_NH_N0_SEM = (squeeze(std(SW_NH_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_NH_N0),[1 2]))))';
SW_NH_N0_B_SEM = (squeeze(std(SW_NH_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_NH_N0_B),[1 2]))))';
LW_NH_N0_SEM = (squeeze(std(LW_NH_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_NH_N0),[1 2]))))';
LW_NH_N0_B_SEM = (squeeze(std(LW_NH_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_NH_N0_B),[1 2]))))';
SW_NH_N60_SEM = (squeeze(std(SW_NH_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_NH_N60),[1 2]))))';
SW_NH_N60_B_SEM = (squeeze(std(SW_NH_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_NH_N60_B),[1 2]))))';
LW_NH_N60_SEM = (squeeze(std(LW_NH_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_NH_N60),[1 2]))))';
LW_NH_N60_B_SEM = (squeeze(std(LW_NH_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_NH_N60_B),[1 2]))))';
SW_NH_N70_SEM = (squeeze(std(SW_NH_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_NH_N70),[1 2]))))';
SW_NH_N70_B_SEM = (squeeze(std(SW_NH_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_NH_N70_B),[1 2]))))';
LW_NH_N70_SEM = (squeeze(std(LW_NH_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_NH_N70),[1 2]))))';
LW_NH_N70_B_SEM = (squeeze(std(LW_NH_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_NH_N70_B),[1 2]))))';

SW_HI_N0_SEM = (squeeze(std(SW_HI_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_N0),[1 2]))))';
SW_HI_N0_B_SEM = (squeeze(std(SW_HI_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_N0_B),[1 2]))))';
LW_HI_N0_SEM = (squeeze(std(LW_HI_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_N0),[1 2]))))';
LW_HI_N0_B_SEM = (squeeze(std(LW_HI_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_N0_B),[1 2]))))';
SW_HI_N60_SEM = (squeeze(std(SW_HI_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_N60),[1 2]))))';
SW_HI_N60_B_SEM = (squeeze(std(SW_HI_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_N60_B),[1 2]))))';
LW_HI_N60_SEM = (squeeze(std(LW_HI_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_N60),[1 2]))))';
LW_HI_N60_B_SEM = (squeeze(std(LW_HI_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_N60_B),[1 2]))))';
SW_HI_N70_SEM = (squeeze(std(SW_HI_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_N70),[1 2]))))';
SW_HI_N70_B_SEM = (squeeze(std(SW_HI_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_N70_B),[1 2]))))';
LW_HI_N70_SEM = (squeeze(std(LW_HI_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_N70),[1 2]))))';
LW_HI_N70_B_SEM = (squeeze(std(LW_HI_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_N70_B),[1 2]))))';

SW_HI_UN_SEM = (squeeze(std(SW_HI_UN,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_UN),[1 2]))))';
SW_HI_UN_B_SEM = (squeeze(std(SW_HI_UN_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_UN_B),[1 2]))))';
LW_HI_UN_SEM = (squeeze(std(LW_HI_UN,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_UN),[1 2]))))';
LW_HI_UN_B_SEM = (squeeze(std(LW_HI_UN_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_UN_B),[1 2]))))';
SW_HI_AA_SEM = (squeeze(std(SW_HI_AA,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_AA),[1 2]))))';
SW_HI_AA_B_SEM = (squeeze(std(SW_HI_AA_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_AA_B),[1 2]))))';
LW_HI_AA_SEM = (squeeze(std(LW_HI_AA,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_AA),[1 2]))))';
LW_HI_AA_B_SEM = (squeeze(std(LW_HI_AA_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_AA_B),[1 2]))))';
SW_HI_AB_SEM = (squeeze(std(SW_HI_AB,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_AB),[1 2]))))';
SW_HI_AB_B_SEM = (squeeze(std(SW_HI_AB_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_AB_B),[1 2]))))';
LW_HI_AB_SEM = (squeeze(std(LW_HI_AB,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_AB),[1 2]))))';
LW_HI_AB_B_SEM = (squeeze(std(LW_HI_AB_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_AB_B),[1 2]))))';

SW_HI_UN_N0_SEM = (squeeze(std(SW_HI_UN_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_UN_N0),[1 2]))))';
SW_HI_UN_N0_B_SEM = (squeeze(std(SW_HI_UN_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_UN_N0_B),[1 2]))))';
LW_HI_UN_N0_SEM = (squeeze(std(LW_HI_UN_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_UN_N0),[1 2]))))';
LW_HI_UN_N0_B_SEM = (squeeze(std(LW_HI_UN_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_UN_N0_B),[1 2]))))';
SW_HI_UN_N60_SEM = (squeeze(std(SW_HI_UN_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_UN_N60),[1 2]))))';
SW_HI_UN_N60_B_SEM = (squeeze(std(SW_HI_UN_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_UN_N60_B),[1 2]))))';
LW_HI_UN_N60_SEM = (squeeze(std(LW_HI_UN_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_UN_N60),[1 2]))))';
LW_HI_UN_N60_B_SEM = (squeeze(std(LW_HI_UN_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_UN_N60_B),[1 2]))))';
SW_HI_UN_N70_SEM = (squeeze(std(SW_HI_UN_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_UN_N70),[1 2]))))';
SW_HI_UN_N70_B_SEM = (squeeze(std(SW_HI_UN_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_UN_N70_B),[1 2]))))';
LW_HI_UN_N70_SEM = (squeeze(std(LW_HI_UN_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_UN_N70),[1 2]))))';
LW_HI_UN_N70_B_SEM = (squeeze(std(LW_HI_UN_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_UN_N70_B),[1 2]))))';
SW_HI_AA_N0_SEM = (squeeze(std(SW_HI_AA_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_AA_N0),[1 2]))))';
SW_HI_AA_N0_B_SEM = (squeeze(std(SW_HI_AA_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_AA_N0_B),[1 2]))))';
LW_HI_AA_N0_SEM = (squeeze(std(LW_HI_AA_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_AA_N0),[1 2]))))';
LW_HI_AA_N0_B_SEM = (squeeze(std(LW_HI_AA_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_AA_N0_B),[1 2]))))';
SW_HI_AA_N60_SEM = (squeeze(std(SW_HI_AA_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_AA_N60),[1 2]))))';
SW_HI_AA_N60_B_SEM = (squeeze(std(SW_HI_AA_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_AA_N60_B),[1 2]))))';
LW_HI_AA_N60_SEM = (squeeze(std(LW_HI_AA_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_AA_N60),[1 2]))))';
LW_HI_AA_N60_B_SEM = (squeeze(std(LW_HI_AA_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_AA_N60_B),[1 2]))))';
SW_HI_AA_N70_SEM = (squeeze(std(SW_HI_AA_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_AA_N70),[1 2]))))';
SW_HI_AA_N70_B_SEM = (squeeze(std(SW_HI_AA_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_AA_N70_B),[1 2]))))';
LW_HI_AA_N70_SEM = (squeeze(std(LW_HI_AA_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_AA_N70),[1 2]))))';
LW_HI_AA_N70_B_SEM = (squeeze(std(LW_HI_AA_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_AA_N70_B),[1 2]))))';
SW_HI_AB_N0_SEM = (squeeze(std(SW_HI_AB_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_AB_N0),[1 2]))))';
SW_HI_AB_N0_B_SEM = (squeeze(std(SW_HI_AB_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_AB_N0_B),[1 2]))))';
LW_HI_AB_N0_SEM = (squeeze(std(LW_HI_AB_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_AB_N0),[1 2]))))';
LW_HI_AB_N0_B_SEM = (squeeze(std(LW_HI_AB_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_AB_N0_B),[1 2]))))';
SW_HI_AB_N60_SEM = (squeeze(std(SW_HI_AB_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_AB_N60),[1 2]))))';
SW_HI_AB_N60_B_SEM = (squeeze(std(SW_HI_AB_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_AB_N60_B),[1 2]))))';
LW_HI_AB_N60_SEM = (squeeze(std(LW_HI_AB_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_AB_N60),[1 2]))))';
LW_HI_AB_N60_B_SEM = (squeeze(std(LW_HI_AB_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_AB_N60_B),[1 2]))))';
SW_HI_AB_N70_SEM = (squeeze(std(SW_HI_AB_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_AB_N70),[1 2]))))';
SW_HI_AB_N70_B_SEM = (squeeze(std(SW_HI_AB_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(SW_HI_AB_N70_B),[1 2]))))';
LW_HI_AB_N70_SEM = (squeeze(std(LW_HI_AB_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_AB_N70),[1 2]))))';
LW_HI_AB_N70_B_SEM = (squeeze(std(LW_HI_AB_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(LW_HI_AB_N70_B),[1 2]))))';
%% Global Plots
f1 = figure;tiledlayout(2,2);ax1 = nexttile;ax2 = nexttile;ax3 = nexttile;ax4 = nexttile;
f2 = figure;tiledlayout(2,2);ax5 = nexttile;ax6 = nexttile;ax7 = nexttile;ax8 = nexttile;
f3 = figure;tiledlayout(1,2);
f4 = figure;tiledlayout(1,2);

hold([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8],'on')
grid([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8],'on')

% Plot event onset and baseline markers
xline(ax1,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax2,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax3,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax4,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax5,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax6,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax7,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax8,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')

xline(ax5,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
xline(ax6,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
xline(ax7,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
xline(ax8,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')

% plot(linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3)),GSW_Mean,color=SColor,linewidth=2);hold on;
% GSW_Mean(isnan(GSW_Mean))=0;GSW_SEM(isnan(GSW_SEM))=0;
% fill([linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3)), flipud(linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3))')'],[(GSW_Mean+GSW_SEM), flipud((GSW_Mean-GSW_SEM)')'],SColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')

% Define time vectors
t_SW_NH = linspace(-TimeStartW,size(SW_NH,3)/Param.Fs,size(SW_NH,3));
t_SW_NH_N0 = linspace(-TimeStartW,size(SW_NH_N0,3)/Param.Fs,size(SW_NH_N0,3));
t_SW_NH_N60 = linspace(-TimeStartW,size(SW_NH_N60,3)/Param.Fs,size(SW_NH_N60,3));
t_SW_NH_N70 = linspace(-TimeStartW,size(SW_NH_N70,3)/Param.Fs,size(SW_NH_N70,3));
t_LW_NH = linspace(-TimeStartW,size(LW_NH,3)/Param.Fs,size(LW_NH,3));
t_LW_NH_N0 = linspace(-TimeStartW,size(LW_NH_N0,3)/Param.Fs,size(LW_NH_N0,3));
t_LW_NH_N60 = linspace(-TimeStartW,size(LW_NH_N60,3)/Param.Fs,size(LW_NH_N60,3));
t_LW_NH_N70 = linspace(-TimeStartW,size(LW_NH_N70,3)/Param.Fs,size(LW_NH_N70,3));
t_SW_HI = linspace(-TimeStartW,size(SW_HI,3)/Param.Fs,size(SW_HI,3));
t_SW_HI_N0 = linspace(-TimeStartW,size(SW_HI_N0,3)/Param.Fs,size(SW_HI_N0,3));
t_SW_HI_N60 = linspace(-TimeStartW,size(SW_HI_N60,3)/Param.Fs,size(SW_HI_N60,3));
t_SW_HI_N70 = linspace(-TimeStartW,size(SW_HI_N70,3)/Param.Fs,size(SW_HI_N70,3));
t_LW_HI = linspace(-TimeStartW,size(LW_HI,3)/Param.Fs,size(LW_HI,3));
t_LW_HI_N0 = linspace(-TimeStartW,size(LW_HI_N0,3)/Param.Fs,size(LW_HI_N0,3));
t_LW_HI_N60 = linspace(-TimeStartW,size(LW_HI_N60,3)/Param.Fs,size(LW_HI_N60,3));
t_LW_HI_N70 = linspace(-TimeStartW,size(LW_HI_N70,3)/Param.Fs,size(LW_HI_N70,3));

t_SW_NH_B = linspace(-TimeStartW,size(SW_NH_B,3)/Param.Fs,size(SW_NH_B,3));
t_SW_NH_N0_B = linspace(-TimeStartW,size(SW_NH_N0_B,3)/Param.Fs,size(SW_NH_N0_B,3));
t_SW_NH_N60_B = linspace(-TimeStartW,size(SW_NH_N60_B,3)/Param.Fs,size(SW_NH_N60_B,3));
t_SW_NH_N70_B = linspace(-TimeStartW,size(SW_NH_N70_B,3)/Param.Fs,size(SW_NH_N70_B,3));
t_LW_NH_B = linspace(-TimeStartW,size(LW_NH_B,3)/Param.Fs,size(LW_NH_B,3));
t_LW_NH_N0_B = linspace(-TimeStartW,size(LW_NH_N0_B,3)/Param.Fs,size(LW_NH_N0_B,3));
t_LW_NH_N60_B = linspace(-TimeStartW,size(LW_NH_N60_B,3)/Param.Fs,size(LW_NH_N60_B,3));
t_LW_NH_N70_B = linspace(-TimeStartW,size(LW_NH_N70_B,3)/Param.Fs,size(LW_NH_N70_B,3));
t_SW_HI_B = linspace(-TimeStartW,size(SW_HI_B,3)/Param.Fs,size(SW_HI_B,3));
t_SW_HI_N0_B = linspace(-TimeStartW,size(SW_HI_N0_B,3)/Param.Fs,size(SW_HI_N0_B,3));
t_SW_HI_N60_B = linspace(-TimeStartW,size(SW_HI_N60_B,3)/Param.Fs,size(SW_HI_N60_B,3));
t_SW_HI_N70_B = linspace(-TimeStartW,size(SW_HI_N70_B,3)/Param.Fs,size(SW_HI_N70_B,3));
t_LW_HI_B = linspace(-TimeStartW,size(LW_HI_B,3)/Param.Fs,size(LW_HI_B,3));
t_LW_HI_N0_B = linspace(-TimeStartW,size(LW_HI_N0_B,3)/Param.Fs,size(LW_HI_N0_B,3));
t_LW_HI_N60_B = linspace(-TimeStartW,size(LW_HI_N60_B,3)/Param.Fs,size(LW_HI_N60_B,3));
t_LW_HI_N70_B = linspace(-TimeStartW,size(LW_HI_N70_B,3)/Param.Fs,size(LW_HI_N70_B,3));


plot(ax1,t_SW_NH,SW_NH_Mean,color=SColor,linewidth=2)
plot(ax1,t_SW_NH_N0,SW_NH_N0_Mean,color=QuietColor)
plot(ax1,t_SW_NH_N60,SW_NH_N60_Mean,color=N60Color)
plot(ax1,t_SW_NH_N70,SW_NH_N70_Mean,color=N70Color)
plot(ax2,t_LW_NH,LW_NH_Mean,color=LColor,linewidth=2)
plot(ax2,t_LW_NH_N0,LW_NH_N0_Mean,color=QuietColor)
plot(ax2,t_LW_NH_N60,LW_NH_N60_Mean,color=N60Color)
plot(ax2,t_LW_NH_N70,LW_NH_N70_Mean,color=N70Color)
plot(ax3,t_SW_HI,SW_HI_Mean,color=SColor,linewidth=2)
plot(ax3,t_SW_HI_N0,SW_HI_N0_Mean,color=QuietColor)
plot(ax3,t_SW_HI_N60,SW_HI_N60_Mean,color=N60Color)
plot(ax3,t_SW_HI_N70,SW_HI_N70_Mean,color=N70Color)
plot(ax4,nan,color=SColor,linewidth=2) % plot nans to show color in legend
plot(ax4,t_LW_HI,LW_HI_Mean,color=LColor,linewidth=2)
plot(ax4,t_LW_HI_N0,LW_HI_N0_Mean,color=QuietColor)
plot(ax4,t_LW_HI_N60,LW_HI_N60_Mean,color=N60Color)
plot(ax4,t_LW_HI_N70,LW_HI_N70_Mean,color=N70Color)

plot(ax5,t_SW_NH_B,SW_NH_B_Mean,color=SColor,linewidth=2)
plot(ax5,t_SW_NH_N0_B,SW_NH_N0_B_Mean,color=QuietColor)
plot(ax5,t_SW_NH_N60_B,SW_NH_N60_B_Mean,color=N60Color)
plot(ax5,t_SW_NH_N70_B,SW_NH_N70_B_Mean,color=N70Color)
plot(ax6,t_LW_NH_B,LW_NH_B_Mean,color=LColor,linewidth=2)
plot(ax6,t_LW_NH_N0_B,LW_NH_N0_B_Mean,color=QuietColor)
plot(ax6,t_LW_NH_N60_B,LW_NH_N60_B_Mean,color=N60Color)
plot(ax6,t_LW_NH_N70_B,LW_NH_N70_B_Mean,color=N70Color)
plot(ax7,t_SW_HI_B,SW_HI_B_Mean,color=SColor,linewidth=2)
plot(ax7,t_SW_HI_N0_B,SW_HI_N0_B_Mean,color=QuietColor)
plot(ax7,t_SW_HI_N60_B,SW_HI_N60_B_Mean,color=N60Color)
plot(ax7,t_SW_HI_N70_B,SW_HI_N70_B_Mean,color=N70Color)
plot(ax8,nan,color=SColor,linewidth=2) % plot nans to show color in legend
plot(ax8,t_LW_HI_B,LW_HI_B_Mean,color=LColor,linewidth=2)
plot(ax8,t_LW_HI_N0_B,LW_HI_N0_B_Mean,color=QuietColor)
plot(ax8,t_LW_HI_N60_B,LW_HI_N60_B_Mean,color=N60Color)
plot(ax8,t_LW_HI_N70_B,LW_HI_N70_B_Mean,color=N70Color)

% plot(ax2,nan,color=SColor,linewidth=2) % plot nans to show color in legend
% plot(ax2,linspace(-TimeStartW,size(GLW,3)/Param.Fs,size(GLW,3)),GLW_Mean,color=LColor,linewidth=2)
% plot(ax2,t_LW_NH,LW_NH_Mean,color=NHColor)
% plot(ax2,t_LW_HI,LW_HI_Mean,color=HIColor)
% plot(ax2,linspace(-TimeStartW,size(LW_N0,3)/Param.Fs,size(LW_N0,3)),LW_N0_Mean,color=QuietColor)
% plot(ax2,linspace(-TimeStartW,size(LW_N60,3)/Param.Fs,size(LW_N60,3)),LW_N60_Mean,color=N60Color)
% plot(ax2,linspace(-TimeStartW,size(LW_N70,3)/Param.Fs,size(LW_N70,3)),LW_N70_Mean,color=N70Color)
% plot(ax3,linspace(-TimeStartW,size(GSW_B,3)/Param.Fs,size(GSW_B,3)),GSW_B_Mean,color=SColor,linewidth=2)
% plot(ax3,t_SW_NH_B,SW_NH_B_Mean,color=NHColor)
% plot(ax3,t_SW_HI_B,SW_HI_B_Mean,color=HIColor)
% plot(ax3,linspace(-TimeStartW,size(SW_N0_B,3)/Param.Fs,size(SW_N0_B,3)),SW_N0_B_Mean,color=QuietColor)
% plot(ax3,linspace(-TimeStartW,size(SW_N60_B,3)/Param.Fs,size(SW_N60_B,3)),SW_N60_B_Mean,color=N60Color)
% plot(ax3,linspace(-TimeStartW,size(SW_N70_B,3)/Param.Fs,size(SW_N70_B,3)),SW_N70_B_Mean,color=N70Color)
% plot(ax4,nan,color=SColor,linewidth=2) % plot nans to show color in legend
% plot(ax4,linspace(-TimeStartW,size(GLW_B,3)/Param.Fs,size(GLW_B,3)),GLW_B_Mean,color=LColor,linewidth=2)
% plot(ax4,t_LW_NH_B,LW_NH_B_Mean,color=NHColor)
% plot(ax4,t_LW_HI_B,LW_HI_B_Mean,color=HIColor)
% plot(ax4,linspace(-TimeStartW,size(LW_N0_B,3)/Param.Fs,size(LW_N0_B,3)),LW_N0_B_Mean,color=QuietColor)
% plot(ax4,linspace(-TimeStartW,size(LW_N60_B,3)/Param.Fs,size(LW_N60_B,3)),LW_N60_B_Mean,color=N60Color)
% plot(ax4,linspace(-TimeStartW,size(LW_N70_B,3)/Param.Fs,size(LW_N70_B,3)),LW_N70_B_Mean,color=N70Color)
% 
% plot(ax5,linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3)),GSW_Mean,color=SColor,linewidth=2)
% plot(ax5,t_SW_NH,SW_NH_Mean,color=NHColor)
% plot(ax5,t_SW_HI,SW_HI_Mean,color=HIColor)
% plot(ax5,linspace(-TimeStartW,size(SW_HI_UN,3)/Param.Fs,size(SW_HI_UN,3)),SW_HI_UN_Mean,color=UNColor)
% plot(ax5,linspace(-TimeStartW,size(SW_HI_AA,3)/Param.Fs,size(SW_HI_AA,3)),SW_HI_AA_Mean,color=AAColor)
% plot(ax5,linspace(-TimeStartW,size(SW_HI_AB,3)/Param.Fs,size(SW_HI_AB,3)),SW_HI_AB_Mean,color=ABColor)
% plot(ax6,nan,color=SColor,linewidth=2) % plot nans to show color in legend
% plot(ax6,linspace(-TimeStartW,size(GLW,3)/Param.Fs,size(GLW,3)),GLW_Mean,color=LColor,linewidth=2)
% plot(ax6,t_LW_NH,LW_NH_Mean,color=NHColor)
% plot(ax6,t_LW_HI,LW_HI_Mean,color=HIColor)
% plot(ax6,linspace(-TimeStartW,size(LW_HI_UN,3)/Param.Fs,size(LW_HI_UN,3)),LW_HI_UN_Mean,color=UNColor)
% plot(ax6,linspace(-TimeStartW,size(LW_HI_AA,3)/Param.Fs,size(LW_HI_AA,3)),LW_HI_AA_Mean,color=AAColor)
% plot(ax6,linspace(-TimeStartW,size(LW_HI_AB,3)/Param.Fs,size(LW_HI_AB,3)),LW_HI_AB_Mean,color=ABColor)
% plot(ax7,linspace(-TimeStartW,size(GSW_B,3)/Param.Fs,size(GSW_B,3)),GSW_B_Mean,color=SColor,linewidth=2)
% plot(ax7,t_SW_NH_B,SW_NH_B_Mean,color=NHColor)
% plot(ax7,t_SW_HI_B,SW_HI_B_Mean,color=HIColor)
% plot(ax7,linspace(-TimeStartW,size(SW_HI_UN_B,3)/Param.Fs,size(SW_HI_UN_B,3)),SW_HI_UN_B_Mean,color=UNColor)
% plot(ax7,linspace(-TimeStartW,size(SW_HI_AA_B,3)/Param.Fs,size(SW_HI_AA_B,3)),SW_HI_AA_B_Mean,color=AAColor)
% plot(ax7,linspace(-TimeStartW,size(SW_HI_AB_B,3)/Param.Fs,size(SW_HI_AB_B,3)),SW_HI_AB_B_Mean,color=ABColor)
% plot(ax8,nan,color=SColor,linewidth=2) % plot nans to show color in legend
% plot(ax8,linspace(-TimeStartW,size(GLW_B,3)/Param.Fs,size(GLW_B,3)),GLW_B_Mean,color=LColor,linewidth=2)
% plot(ax8,t_LW_NH_B,LW_NH_B_Mean,color=NHColor)
% plot(ax8,t_LW_HI_B,LW_HI_B_Mean,color=HIColor)
% plot(ax8,linspace(-TimeStartW,size(LW_HI_UN_B,3)/Param.Fs,size(LW_HI_UN_B,3)),LW_HI_UN_B_Mean,color=UNColor)
% plot(ax8,linspace(-TimeStartW,size(LW_HI_AA_B,3)/Param.Fs,size(LW_HI_AA_B,3)),LW_HI_AA_B_Mean,color=AAColor)
% plot(ax8,linspace(-TimeStartW,size(LW_HI_AB_B,3)/Param.Fs,size(LW_HI_AB_B,3)),LW_HI_AB_B_Mean,color=ABColor)

% plot(ax1,t_SW_NH_N0,SW_NH_N0_Mean,color=QuietColor,linestyle=':')
% plot(ax1,t_SW_NH_N60,SW_NH_N60_Mean,color=N60Color,linestyle=':')
% plot(ax1,t_SW_NH_N70,SW_NH_N70_Mean,color=N70Color,linestyle=':')
% plot(ax1,linspace(-TimeStartW,size(SW_HI_UN_N0,3)/Param.Fs,size(SW_HI_UN_N0,3)),SW_HI_UN_N0_Mean,color=QuietColor,linestyle='--')
% plot(ax1,linspace(-TimeStartW,size(SW_HI_UN_N60,3)/Param.Fs,size(SW_HI_UN_N60,3)),SW_HI_UN_N60_Mean,color=N60Color,linestyle='--')
% plot(ax1,linspace(-TimeStartW,size(SW_HI_UN_N70,3)/Param.Fs,size(SW_HI_UN_N70,3)),SW_HI_UN_N70_Mean,color=N70Color,linestyle='--')
% plot(ax1,linspace(-TimeStartW,size(SW_HI_AA_N0,3)/Param.Fs,size(SW_HI_AA_N0,3)),SW_HI_AA_N0_Mean,color=QuietColor,linestyle='-.')
% plot(ax1,linspace(-TimeStartW,size(SW_HI_AA_N60,3)/Param.Fs,size(SW_HI_AA_N60,3)),SW_HI_AA_N60_Mean,color=N60Color,linestyle='-.')
% plot(ax1,linspace(-TimeStartW,size(SW_HI_AA_N70,3)/Param.Fs,size(SW_HI_AA_N70,3)),SW_HI_AA_N70_Mean,color=N70Color,linestyle='-.')
% plot(ax1,linspace(-TimeStartW,size(SW_HI_AB_N0,3)/Param.Fs,size(SW_HI_AB_N0,3)),SW_HI_AB_N0_Mean,color=QuietColor,linestyle='-.')
% plot(ax1,linspace(-TimeStartW,size(SW_HI_AB_N60,3)/Param.Fs,size(SW_HI_AB_N60,3)),SW_HI_AB_N60_Mean,color=N60Color,linestyle='-.')
% plot(ax1,linspace(-TimeStartW,size(SW_HI_AB_N70,3)/Param.Fs,size(SW_HI_AB_N70,3)),SW_HI_AB_N70_Mean,color=N70Color,linestyle='-.')
% plot(ax2,t_LW_NH_N0,LW_NH_N0_Mean,color=QuietColor,linestyle=':')
% plot(ax2,t_LW_NH_N60,LW_NH_N60_Mean,color=N60Color,linestyle=':')
% plot(ax2,t_LW_NH_N70,LW_NH_N70_Mean,color=N70Color,linestyle=':')
% plot(ax2,linspace(-TimeStartW,size(LW_HI_UN_N0,3)/Param.Fs,size(LW_HI_UN_N0,3)),LW_HI_UN_N0_Mean,color=QuietColor,linestyle='--')
% plot(ax2,linspace(-TimeStartW,size(LW_HI_UN_N60,3)/Param.Fs,size(LW_HI_UN_N60,3)),LW_HI_UN_N60_Mean,color=N60Color,linestyle='--')
% plot(ax2,linspace(-TimeStartW,size(LW_HI_UN_N70,3)/Param.Fs,size(LW_HI_UN_N70,3)),LW_HI_UN_N70_Mean,color=N70Color,linestyle='--')
% plot(ax2,linspace(-TimeStartW,size(LW_HI_AA_N0,3)/Param.Fs,size(LW_HI_AA_N0,3)),LW_HI_AA_N0_Mean,color=QuietColor,linestyle='-.')
% plot(ax2,linspace(-TimeStartW,size(LW_HI_AA_N60,3)/Param.Fs,size(LW_HI_AA_N60,3)),LW_HI_AA_N60_Mean,color=N60Color,linestyle='-.')
% plot(ax2,linspace(-TimeStartW,size(LW_HI_AA_N70,3)/Param.Fs,size(LW_HI_AA_N70,3)),LW_HI_AA_N70_Mean,color=N70Color,linestyle='-.')
% plot(ax2,linspace(-TimeStartW,size(LW_HI_AB_N0,3)/Param.Fs,size(LW_HI_AB_N0,3)),LW_HI_AB_N0_Mean,color=QuietColor,linestyle='-.')
% plot(ax2,linspace(-TimeStartW,size(LW_HI_AB_N60,3)/Param.Fs,size(LW_HI_AB_N60,3)),LW_HI_AB_N60_Mean,color=N60Color,linestyle='-.')
% plot(ax2,linspace(-TimeStartW,size(LW_HI_AB_N70,3)/Param.Fs,size(LW_HI_AB_N70,3)),LW_HI_AB_N70_Mean,color=N70Color,linestyle='-.')
% 
% plot(ax3,t_SW_NH_N0_B,SW_NH_N0_B_Mean,color=QuietColor,linestyle=':')
% plot(ax3,t_SW_NH_N60_B,SW_NH_N60_B_Mean,color=N60Color,linestyle=':')
% plot(ax3,t_SW_NH_N70_B,SW_NH_N70_B_Mean,color=N70Color,linestyle=':')
% plot(ax3,linspace(-TimeStartW,size(SW_HI_UN_N0_B,3)/Param.Fs,size(SW_HI_UN_N0_B,3)),SW_HI_UN_N0_B_Mean,color=QuietColor,linestyle='--')
% plot(ax3,linspace(-TimeStartW,size(SW_HI_UN_N60_B,3)/Param.Fs,size(SW_HI_UN_N60_B,3)),SW_HI_UN_N60_B_Mean,color=N60Color,linestyle='--')
% plot(ax3,linspace(-TimeStartW,size(SW_HI_UN_N70_B,3)/Param.Fs,size(SW_HI_UN_N70_B,3)),SW_HI_UN_N70_B_Mean,color=N70Color,linestyle='--')
% plot(ax3,linspace(-TimeStartW,size(SW_HI_AA_N0_B,3)/Param.Fs,size(SW_HI_AA_N0_B,3)),SW_HI_AA_N0_B_Mean,color=QuietColor,linestyle='-.')
% plot(ax3,linspace(-TimeStartW,size(SW_HI_AA_N60_B,3)/Param.Fs,size(SW_HI_AA_N60_B,3)),SW_HI_AA_N60_B_Mean,color=N60Color,linestyle='-.')
% plot(ax3,linspace(-TimeStartW,size(SW_HI_AA_N70_B,3)/Param.Fs,size(SW_HI_AA_N70_B,3)),SW_HI_AA_N70_B_Mean,color=N70Color,linestyle='-.')
% plot(ax3,linspace(-TimeStartW,size(SW_HI_AB_N0_B,3)/Param.Fs,size(SW_HI_AB_N0_B,3)),SW_HI_AB_N0_B_Mean,color=QuietColor,linestyle='-.')
% plot(ax3,linspace(-TimeStartW,size(SW_HI_AB_N60_B,3)/Param.Fs,size(SW_HI_AB_N60_B,3)),SW_HI_AB_N60_B_Mean,color=N60Color,linestyle='-.')
% plot(ax3,linspace(-TimeStartW,size(SW_HI_AB_N70_B,3)/Param.Fs,size(SW_HI_AB_N70_B,3)),SW_HI_AB_N70_B_Mean,color=N70Color,linestyle='-.')
% plot(ax4,t_LW_NH_N0_B,LW_NH_N0_B_Mean,color=QuietColor,linestyle=':')
% plot(ax4,t_LW_NH_N60_B,LW_NH_N60_B_Mean,color=N60Color,linestyle=':')
% plot(ax4,t_LW_NH_N70_B,LW_NH_N70_B_Mean,color=N70Color,linestyle=':')
% plot(ax4,linspace(-TimeStartW,size(LW_HI_UN_N0_B,3)/Param.Fs,size(LW_HI_UN_N0_B,3)),LW_HI_UN_N0_B_Mean,color=QuietColor,linestyle='--')
% plot(ax4,linspace(-TimeStartW,size(LW_HI_UN_N60_B,3)/Param.Fs,size(LW_HI_UN_N60_B,3)),LW_HI_UN_N60_B_Mean,color=N60Color,linestyle='--')
% plot(ax4,linspace(-TimeStartW,size(LW_HI_UN_N70_B,3)/Param.Fs,size(LW_HI_UN_N70_B,3)),LW_HI_UN_N70_B_Mean,color=N70Color,linestyle='--')
% plot(ax4,linspace(-TimeStartW,size(LW_HI_AA_N0_B,3)/Param.Fs,size(LW_HI_AA_N0_B,3)),LW_HI_AA_N0_B_Mean,color=QuietColor,linestyle='-.')
% plot(ax4,linspace(-TimeStartW,size(LW_HI_AA_N60_B,3)/Param.Fs,size(LW_HI_AA_N60_B,3)),LW_HI_AA_N60_B_Mean,color=N60Color,linestyle='-.')
% plot(ax4,linspace(-TimeStartW,size(LW_HI_AA_N70_B,3)/Param.Fs,size(LW_HI_AA_N70_B,3)),LW_HI_AA_N70_B_Mean,color=N70Color,linestyle='-.')
% plot(ax4,linspace(-TimeStartW,size(LW_HI_AB_N0_B,3)/Param.Fs,size(LW_HI_AB_N0_B,3)),LW_HI_AB_N0_B_Mean,color=QuietColor,linestyle='-.')
% plot(ax4,linspace(-TimeStartW,size(LW_HI_AB_N60_B,3)/Param.Fs,size(LW_HI_AB_N60_B,3)),LW_HI_AB_N60_B_Mean,color=N60Color,linestyle='-.')
% plot(ax4,linspace(-TimeStartW,size(LW_HI_AB_N70_B,3)/Param.Fs,size(LW_HI_AB_N70_B,3)),LW_HI_AB_N70_B_Mean,color=N70Color,linestyle='-.')

% Mean and SEM: Set NaNs to 0 in order to use "fill" correctly
GSW_Mean(isnan(GSW_Mean))=0;GSW_SEM(isnan(GSW_SEM))=0;
GLW_Mean(isnan(GLW_Mean))=0;GLW_SEM(isnan(GLW_SEM))=0;
SW_NH_Mean(isnan(SW_NH_Mean))=0;SW_NH_SEM(isnan(SW_NH_SEM))=0;
LW_NH_Mean(isnan(LW_NH_Mean))=0;LW_NH_SEM(isnan(LW_NH_SEM))=0;
SW_N0_Mean(isnan(SW_N0_Mean))=0;SW_N0_SEM(isnan(SW_N0_SEM))=0;
SW_N60_Mean(isnan(SW_N60_Mean))=0;SW_N60_SEM(isnan(SW_N60_SEM))=0;
SW_N70_Mean(isnan(SW_N70_Mean))=0;SW_N70_SEM(isnan(SW_N70_SEM))=0;
LW_N0_Mean(isnan(LW_N0_Mean))=0;LW_N0_SEM(isnan(LW_N0_SEM))=0;
LW_N60_Mean(isnan(LW_N60_Mean))=0;LW_N60_SEM(isnan(LW_N60_SEM))=0;
LW_N70_Mean(isnan(LW_N70_Mean))=0;LW_N70_SEM(isnan(LW_N70_SEM))=0;
SW_NH_N0_Mean(isnan(SW_NH_N0_Mean))=0;SW_NH_N0_SEM(isnan(SW_NH_N0_SEM))=0;
SW_NH_N60_Mean(isnan(SW_NH_N60_Mean))=0;SW_NH_N60_SEM(isnan(SW_NH_N60_SEM))=0;
SW_NH_N70_Mean(isnan(SW_NH_N70_Mean))=0;SW_NH_N70_SEM(isnan(SW_NH_N70_SEM))=0;
LW_NH_N0_Mean(isnan(LW_NH_N0_Mean))=0;LW_NH_N0_SEM(isnan(LW_NH_N0_SEM))=0;
LW_NH_N60_Mean(isnan(LW_NH_N60_Mean))=0;LW_NH_N60_SEM(isnan(LW_NH_N60_SEM))=0;
LW_NH_N70_Mean(isnan(LW_NH_N70_Mean))=0;LW_NH_N70_SEM(isnan(LW_NH_N70_SEM))=0;
SW_HI_Mean(isnan(SW_HI_Mean))=0;SW_HI_SEM(isnan(SW_HI_SEM))=0;
LW_HI_Mean(isnan(LW_HI_Mean))=0;LW_HI_SEM(isnan(LW_HI_SEM))=0;
SW_HI_N0_Mean(isnan(SW_HI_N0_Mean))=0;SW_HI_N0_SEM(isnan(SW_HI_N0_SEM))=0;
SW_HI_N60_Mean(isnan(SW_HI_N60_Mean))=0;SW_HI_N60_SEM(isnan(SW_HI_N60_SEM))=0;
SW_HI_N70_Mean(isnan(SW_HI_N70_Mean))=0;SW_HI_N70_SEM(isnan(SW_HI_N70_SEM))=0;
LW_HI_N0_Mean(isnan(LW_HI_N0_Mean))=0;LW_HI_N0_SEM(isnan(LW_HI_N0_SEM))=0;
LW_HI_N60_Mean(isnan(LW_HI_N60_Mean))=0;LW_HI_N60_SEM(isnan(LW_HI_N60_SEM))=0;
LW_HI_N70_Mean(isnan(LW_HI_N70_Mean))=0;LW_HI_N70_SEM(isnan(LW_HI_N70_SEM))=0;
SW_HI_UN_Mean(isnan(SW_HI_UN_Mean))=0;SW_HI_UN_SEM(isnan(SW_HI_UN_SEM))=0;
LW_HI_UN_Mean(isnan(LW_HI_UN_Mean))=0;LW_HI_UN_SEM(isnan(LW_HI_UN_SEM))=0;
SW_HI_AA_Mean(isnan(SW_HI_AA_Mean))=0;SW_HI_AA_SEM(isnan(SW_HI_AA_SEM))=0;
LW_HI_AA_Mean(isnan(LW_HI_AA_Mean))=0;LW_HI_AA_SEM(isnan(LW_HI_AA_SEM))=0;
SW_HI_AB_Mean(isnan(SW_HI_AB_Mean))=0;SW_HI_AB_SEM(isnan(SW_HI_AB_SEM))=0;
LW_HI_AB_Mean(isnan(LW_HI_AB_Mean))=0;LW_HI_AB_SEM(isnan(LW_HI_AB_SEM))=0;
SW_HI_UN_N0_Mean(isnan(SW_HI_UN_N0_Mean))=0;SW_HI_UN_N0_SEM(isnan(SW_HI_UN_N0_SEM))=0;
SW_HI_UN_N60_Mean(isnan(SW_HI_UN_N60_Mean))=0;SW_HI_UN_N60_SEM(isnan(SW_HI_UN_N60_SEM))=0;
SW_HI_UN_N70_Mean(isnan(SW_HI_UN_N70_Mean))=0;SW_HI_UN_N70_SEM(isnan(SW_HI_UN_N70_SEM))=0;
LW_HI_UN_N0_Mean(isnan(LW_HI_UN_N0_Mean))=0;LW_HI_UN_N0_SEM(isnan(LW_HI_UN_N0_SEM))=0;
LW_HI_UN_N60_Mean(isnan(LW_HI_UN_N60_Mean))=0;LW_HI_UN_N60_SEM(isnan(LW_HI_UN_N60_SEM))=0;
LW_HI_UN_N70_Mean(isnan(LW_HI_UN_N70_Mean))=0;LW_HI_UN_N70_SEM(isnan(LW_HI_UN_N70_SEM))=0;
SW_HI_AA_N0_Mean(isnan(SW_HI_AA_N0_Mean))=0;SW_HI_AA_N0_SEM(isnan(SW_HI_AA_N0_SEM))=0;
SW_HI_AA_N60_Mean(isnan(SW_HI_AA_N60_Mean))=0;SW_HI_AA_N60_SEM(isnan(SW_HI_AA_N60_SEM))=0;
SW_HI_AA_N70_Mean(isnan(SW_HI_AA_N70_Mean))=0;SW_HI_AA_N70_SEM(isnan(SW_HI_AA_N70_SEM))=0;
LW_HI_AA_N0_Mean(isnan(LW_HI_AA_N0_Mean))=0;LW_HI_AA_N0_SEM(isnan(LW_HI_AA_N0_SEM))=0;
LW_HI_AA_N60_Mean(isnan(LW_HI_AA_N60_Mean))=0;LW_HI_AA_N60_SEM(isnan(LW_HI_AA_N60_SEM))=0;
LW_HI_AA_N70_Mean(isnan(LW_HI_AA_N70_Mean))=0;LW_HI_AA_N70_SEM(isnan(LW_HI_AA_N70_SEM))=0;
SW_HI_AB_N0_Mean(isnan(SW_HI_AB_N0_Mean))=0;SW_HI_AB_N0_SEM(isnan(SW_HI_AB_N0_SEM))=0;
SW_HI_AB_N60_Mean(isnan(SW_HI_AB_N60_Mean))=0;SW_HI_AB_N60_SEM(isnan(SW_HI_AB_N60_SEM))=0;
SW_HI_AB_N70_Mean(isnan(SW_HI_AB_N70_Mean))=0;SW_HI_AB_N70_SEM(isnan(SW_HI_AB_N70_SEM))=0;
LW_HI_AB_N0_Mean(isnan(LW_HI_AB_N0_Mean))=0;LW_HI_AB_N0_SEM(isnan(LW_HI_AB_N0_SEM))=0;
LW_HI_AB_N60_Mean(isnan(LW_HI_AB_N60_Mean))=0;LW_HI_AB_N60_SEM(isnan(LW_HI_AB_N60_SEM))=0;
LW_HI_AB_N70_Mean(isnan(LW_HI_AB_N70_Mean))=0;LW_HI_AB_N70_SEM(isnan(LW_HI_AB_N70_SEM))=0;

GSW_B_Mean(isnan(GSW_B_Mean))=0;GSW_B_SEM(isnan(GSW_B_SEM))=0;
GLW_B_Mean(isnan(GLW_B_Mean))=0;GLW_B_SEM(isnan(GLW_B_SEM))=0;
SW_NH_B_Mean(isnan(SW_NH_B_Mean))=0;SW_NH_B_SEM(isnan(SW_NH_B_SEM))=0;
LW_NH_B_Mean(isnan(LW_NH_B_Mean))=0;LW_NH_B_SEM(isnan(LW_NH_B_SEM))=0;
SW_N0_B_Mean(isnan(SW_N0_B_Mean))=0;SW_N0_B_SEM(isnan(SW_N0_B_SEM))=0;
SW_N60_B_Mean(isnan(SW_N60_B_Mean))=0;SW_N60_B_SEM(isnan(SW_N60_B_SEM))=0;
SW_N70_B_Mean(isnan(SW_N70_B_Mean))=0;SW_N70_B_SEM(isnan(SW_N70_B_SEM))=0;
LW_N0_B_Mean(isnan(LW_N0_B_Mean))=0;LW_N0_B_SEM(isnan(LW_N0_B_SEM))=0;
LW_N60_B_Mean(isnan(LW_N60_B_Mean))=0;LW_N60_B_SEM(isnan(LW_N60_B_SEM))=0;
LW_N70_B_Mean(isnan(LW_N70_B_Mean))=0;LW_N70_B_SEM(isnan(LW_N70_B_SEM))=0;
SW_NH_N0_B_Mean(isnan(SW_NH_N0_B_Mean))=0;SW_NH_N0_B_SEM(isnan(SW_NH_N0_B_SEM))=0;
SW_NH_N60_B_Mean(isnan(SW_NH_N60_B_Mean))=0;SW_NH_N60_B_SEM(isnan(SW_NH_N60_B_SEM))=0;
SW_NH_N70_B_Mean(isnan(SW_NH_N70_B_Mean))=0;SW_NH_N70_B_SEM(isnan(SW_NH_N70_B_SEM))=0;
LW_NH_N0_B_Mean(isnan(LW_NH_N0_B_Mean))=0;LW_NH_N0_B_SEM(isnan(LW_NH_N0_B_SEM))=0;
LW_NH_N60_B_Mean(isnan(LW_NH_N60_B_Mean))=0;LW_NH_N60_B_SEM(isnan(LW_NH_N60_B_SEM))=0;
LW_NH_N70_B_Mean(isnan(LW_NH_N70_B_Mean))=0;LW_NH_N70_B_SEM(isnan(LW_NH_N70_B_SEM))=0;
SW_HI_B_Mean(isnan(SW_HI_B_Mean))=0;SW_HI_B_SEM(isnan(SW_HI_B_SEM))=0;
LW_HI_B_Mean(isnan(LW_HI_B_Mean))=0;LW_HI_B_SEM(isnan(LW_HI_B_SEM))=0;
SW_HI_N0_B_Mean(isnan(SW_HI_N0_B_Mean))=0;SW_HI_N0_B_SEM(isnan(SW_HI_N0_B_SEM))=0;
SW_HI_N60_B_Mean(isnan(SW_HI_N60_B_Mean))=0;SW_HI_N60_B_SEM(isnan(SW_HI_N60_B_SEM))=0;
SW_HI_N70_B_Mean(isnan(SW_HI_N70_B_Mean))=0;SW_HI_N70_B_SEM(isnan(SW_HI_N70_B_SEM))=0;
LW_HI_N0_B_Mean(isnan(LW_HI_N0_B_Mean))=0;LW_HI_N0_B_SEM(isnan(LW_HI_N0_B_SEM))=0;
LW_HI_N60_B_Mean(isnan(LW_HI_N60_B_Mean))=0;LW_HI_N60_B_SEM(isnan(LW_HI_N60_B_SEM))=0;
LW_HI_N70_B_Mean(isnan(LW_HI_N70_B_Mean))=0;LW_HI_N70_B_SEM(isnan(LW_HI_N70_B_SEM))=0;
SW_HI_UN_B_Mean(isnan(SW_HI_UN_B_Mean))=0;SW_HI_UN_B_SEM(isnan(SW_HI_UN_B_SEM))=0;
LW_HI_UN_B_Mean(isnan(LW_HI_UN_B_Mean))=0;LW_HI_UN_B_SEM(isnan(LW_HI_UN_B_SEM))=0;
SW_HI_AA_B_Mean(isnan(SW_HI_AA_B_Mean))=0;SW_HI_AA_B_SEM(isnan(SW_HI_AA_B_SEM))=0;
LW_HI_AA_B_Mean(isnan(LW_HI_AA_B_Mean))=0;LW_HI_AA_B_SEM(isnan(LW_HI_AA_B_SEM))=0;
SW_HI_AB_B_Mean(isnan(SW_HI_AB_B_Mean))=0;SW_HI_AB_B_SEM(isnan(SW_HI_AB_B_SEM))=0;
LW_HI_AB_B_Mean(isnan(LW_HI_AB_B_Mean))=0;LW_HI_AB_B_SEM(isnan(LW_HI_AB_B_SEM))=0;
SW_HI_UN_N0_B_Mean(isnan(SW_HI_UN_N0_B_Mean))=0;SW_HI_UN_N0_B_SEM(isnan(SW_HI_UN_N0_B_SEM))=0;
SW_HI_UN_N60_B_Mean(isnan(SW_HI_UN_N60_B_Mean))=0;SW_HI_UN_N60_B_SEM(isnan(SW_HI_UN_N60_B_SEM))=0;
SW_HI_UN_N70_B_Mean(isnan(SW_HI_UN_N70_B_Mean))=0;SW_HI_UN_N70_B_SEM(isnan(SW_HI_UN_N70_B_SEM))=0;
LW_HI_UN_N0_B_Mean(isnan(LW_HI_UN_N0_B_Mean))=0;LW_HI_UN_N0_B_SEM(isnan(LW_HI_UN_N0_B_SEM))=0;
LW_HI_UN_N60_B_Mean(isnan(LW_HI_UN_N60_B_Mean))=0;LW_HI_UN_N60_B_SEM(isnan(LW_HI_UN_N60_B_SEM))=0;
LW_HI_UN_N70_B_Mean(isnan(LW_HI_UN_N70_B_Mean))=0;LW_HI_UN_N70_B_SEM(isnan(LW_HI_UN_N70_B_SEM))=0;
SW_HI_AA_N0_B_Mean(isnan(SW_HI_AA_N0_B_Mean))=0;SW_HI_AA_N0_B_SEM(isnan(SW_HI_AA_N0_B_SEM))=0;
SW_HI_AA_N60_B_Mean(isnan(SW_HI_AA_N60_B_Mean))=0;SW_HI_AA_N60_B_SEM(isnan(SW_HI_AA_N60_B_SEM))=0;
SW_HI_AA_N70_B_Mean(isnan(SW_HI_AA_N70_B_Mean))=0;SW_HI_AA_N70_B_SEM(isnan(SW_HI_AA_N70_B_SEM))=0;
LW_HI_AA_N0_B_Mean(isnan(LW_HI_AA_N0_B_Mean))=0;LW_HI_AA_N0_B_SEM(isnan(LW_HI_AA_N0_B_SEM))=0;
LW_HI_AA_N60_B_Mean(isnan(LW_HI_AA_N60_B_Mean))=0;LW_HI_AA_N60_B_SEM(isnan(LW_HI_AA_N60_B_SEM))=0;
LW_HI_AA_N70_B_Mean(isnan(LW_HI_AA_N70_B_Mean))=0;LW_HI_AA_N70_B_SEM(isnan(LW_HI_AA_N70_B_SEM))=0;
SW_HI_AB_N0_B_Mean(isnan(SW_HI_AB_N0_B_Mean))=0;SW_HI_AB_N0_B_SEM(isnan(SW_HI_AB_N0_B_SEM))=0;
SW_HI_AB_N60_B_Mean(isnan(SW_HI_AB_N60_B_Mean))=0;SW_HI_AB_N60_B_SEM(isnan(SW_HI_AB_N60_B_SEM))=0;
SW_HI_AB_N70_B_Mean(isnan(SW_HI_AB_N70_B_Mean))=0;SW_HI_AB_N70_B_SEM(isnan(SW_HI_AB_N70_B_SEM))=0;
LW_HI_AB_N0_B_Mean(isnan(LW_HI_AB_N0_B_Mean))=0;LW_HI_AB_N0_B_SEM(isnan(LW_HI_AB_N0_B_SEM))=0;
LW_HI_AB_N60_B_Mean(isnan(LW_HI_AB_N60_B_Mean))=0;LW_HI_AB_N60_B_SEM(isnan(LW_HI_AB_N60_B_SEM))=0;
LW_HI_AB_N70_B_Mean(isnan(LW_HI_AB_N70_B_Mean))=0;LW_HI_AB_N70_B_SEM(isnan(LW_HI_AB_N70_B_SEM))=0;

fill(ax1,[t_SW_NH, flipud(t_SW_NH')'],[(SW_NH_Mean+SW_NH_SEM), flipud((SW_NH_Mean-SW_NH_SEM)')'],SColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[t_SW_NH_N0, flipud(t_SW_NH_N0')'],[(SW_NH_N0_Mean+SW_NH_N0_SEM), flipud((SW_NH_N0_Mean-SW_NH_N0_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[t_SW_NH_N60, flipud(t_SW_NH_N60')'],[(SW_NH_N60_Mean+SW_NH_N60_SEM), flipud((SW_NH_N60_Mean-SW_NH_N60_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[t_SW_NH_N70, flipud(t_SW_NH_N70')'],[(SW_NH_N70_Mean+SW_NH_N70_SEM), flipud((SW_NH_N70_Mean-SW_NH_N70_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[t_LW_NH, flipud(t_LW_NH')'],[(LW_NH_Mean+LW_NH_SEM), flipud((LW_NH_Mean-LW_NH_SEM)')'],LColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[t_LW_NH_N0, flipud(t_LW_NH_N0')'],[(LW_NH_N0_Mean+LW_NH_N0_SEM), flipud((LW_NH_N0_Mean-LW_NH_N0_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[t_LW_NH_N60, flipud(t_LW_NH_N60')'],[(LW_NH_N60_Mean+LW_NH_N60_SEM), flipud((LW_NH_N60_Mean-LW_NH_N60_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[t_LW_NH_N70, flipud(t_LW_NH_N70')'],[(LW_NH_N70_Mean+LW_NH_N70_SEM), flipud((LW_NH_N70_Mean-LW_NH_N70_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[t_SW_HI, flipud(t_SW_HI')'],[(SW_HI_Mean+SW_HI_SEM), flipud((SW_HI_Mean-SW_HI_SEM)')'],SColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[t_SW_HI_N0, flipud(t_SW_HI_N0')'],[(SW_HI_N0_Mean+SW_HI_N0_SEM), flipud((SW_HI_N0_Mean-SW_HI_N0_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[t_SW_HI_N60, flipud(t_SW_HI_N60')'],[(SW_HI_N60_Mean+SW_HI_N60_SEM), flipud((SW_HI_N60_Mean-SW_HI_N60_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[t_SW_HI_N70, flipud(t_SW_HI_N70')'],[(SW_HI_N70_Mean+SW_HI_N70_SEM), flipud((SW_HI_N70_Mean-SW_HI_N70_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[t_LW_HI, flipud(t_LW_HI')'],[(LW_HI_Mean+LW_HI_SEM), flipud((LW_HI_Mean-LW_HI_SEM)')'],LColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[t_LW_HI_N0, flipud(t_LW_HI_N0')'],[(LW_HI_N0_Mean+LW_HI_N0_SEM), flipud((LW_HI_N0_Mean-LW_HI_N0_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[t_LW_HI_N60, flipud(t_LW_HI_N60')'],[(LW_HI_N60_Mean+LW_HI_N60_SEM), flipud((LW_HI_N60_Mean-LW_HI_N60_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[t_LW_HI_N70, flipud(t_LW_HI_N70')'],[(LW_HI_N70_Mean+LW_HI_N70_SEM), flipud((LW_HI_N70_Mean-LW_HI_N70_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')

fill(ax5,[t_SW_NH_B, flipud(t_SW_NH_B')'],[(SW_NH_B_Mean+SW_NH_B_SEM), flipud((SW_NH_B_Mean-SW_NH_B_SEM)')'],SColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax5,[t_SW_NH_N0_B, flipud(t_SW_NH_N0_B')'],[(SW_NH_N0_B_Mean+SW_NH_N0_B_SEM), flipud((SW_NH_N0_B_Mean-SW_NH_N0_B_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax5,[t_SW_NH_N60_B, flipud(t_SW_NH_N60_B')'],[(SW_NH_N60_B_Mean+SW_NH_N60_B_SEM), flipud((SW_NH_N60_B_Mean-SW_NH_N60_B_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax5,[t_SW_NH_N70_B, flipud(t_SW_NH_N70_B')'],[(SW_NH_N70_B_Mean+SW_NH_N70_B_SEM), flipud((SW_NH_N70_B_Mean-SW_NH_N70_B_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax6,[t_LW_NH_B, flipud(t_LW_NH_B')'],[(LW_NH_B_Mean+LW_NH_B_SEM), flipud((LW_NH_B_Mean-LW_NH_B_SEM)')'],LColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax6,[t_LW_NH_N0_B, flipud(t_LW_NH_N0_B')'],[(LW_NH_N0_B_Mean+LW_NH_N0_B_SEM), flipud((LW_NH_N0_B_Mean-LW_NH_N0_B_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax6,[t_LW_NH_N60_B, flipud(t_LW_NH_N60_B')'],[(LW_NH_N60_B_Mean+LW_NH_N60_B_SEM), flipud((LW_NH_N60_B_Mean-LW_NH_N60_B_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax6,[t_LW_NH_N70_B, flipud(t_LW_NH_N70_B')'],[(LW_NH_N70_B_Mean+LW_NH_N70_B_SEM), flipud((LW_NH_N70_B_Mean-LW_NH_N70_B_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax7,[t_SW_HI_B, flipud(t_SW_HI_B')'],[(SW_HI_B_Mean+SW_HI_B_SEM), flipud((SW_HI_B_Mean-SW_HI_B_SEM)')'],SColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax7,[t_SW_HI_N0_B, flipud(t_SW_HI_N0_B')'],[(SW_HI_N0_B_Mean+SW_HI_N0_B_SEM), flipud((SW_HI_N0_B_Mean-SW_HI_N0_B_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax7,[t_SW_HI_N60_B, flipud(t_SW_HI_N60_B')'],[(SW_HI_N60_B_Mean+SW_HI_N60_B_SEM), flipud((SW_HI_N60_B_Mean-SW_HI_N60_B_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax7,[t_SW_HI_N70_B, flipud(t_SW_HI_N70_B')'],[(SW_HI_N70_B_Mean+SW_HI_N70_B_SEM), flipud((SW_HI_N70_B_Mean-SW_HI_N70_B_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax8,[t_LW_HI_B, flipud(t_LW_HI_B')'],[(LW_HI_B_Mean+LW_HI_B_SEM), flipud((LW_HI_B_Mean-LW_HI_B_SEM)')'],LColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax8,[t_LW_HI_N0_B, flipud(t_LW_HI_N0_B')'],[(LW_HI_N0_B_Mean+LW_HI_N0_B_SEM), flipud((LW_HI_N0_B_Mean-LW_HI_N0_B_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax8,[t_LW_HI_N60_B, flipud(t_LW_HI_N60_B')'],[(LW_HI_N60_B_Mean+LW_HI_N60_B_SEM), flipud((LW_HI_N60_B_Mean-LW_HI_N60_B_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
fill(ax8,[t_LW_HI_N70_B, flipud(t_LW_HI_N70_B')'],[(LW_HI_N70_B_Mean+LW_HI_N70_B_SEM), flipud((LW_HI_N70_B_Mean-LW_HI_N70_B_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')


% fill(ax1,[linspace(-TimeStartW,size(SW_HI_UN_N0,3)/Param.Fs,size(SW_HI_UN_N0,3)), flipud(linspace(-TimeStartW,size(SW_HI_UN_N0,3)/Param.Fs,size(SW_HI_UN_N0,3))')'],[(SW_HI_UN_N0_Mean+SW_HI_UN_N0_SEM), flipud((SW_HI_UN_N0_Mean-SW_HI_UN_N0_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax1,[linspace(-TimeStartW,size(SW_HI_UN_N60,3)/Param.Fs,size(SW_HI_UN_N60,3)), flipud(linspace(-TimeStartW,size(SW_HI_UN_N60,3)/Param.Fs,size(SW_HI_UN_N60,3))')'],[(SW_HI_UN_N60_Mean+SW_HI_UN_N60_SEM), flipud((SW_HI_UN_N60_Mean-SW_HI_UN_N60_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax1,[linspace(-TimeStartW,size(SW_HI_UN_N70,3)/Param.Fs,size(SW_HI_UN_N70,3)), flipud(linspace(-TimeStartW,size(SW_HI_UN_N70,3)/Param.Fs,size(SW_HI_UN_N70,3))')'],[(SW_HI_UN_N70_Mean+SW_HI_UN_N70_SEM), flipud((SW_HI_UN_N70_Mean-SW_HI_UN_N70_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax1,[linspace(-TimeStartW,size(SW_HI_AA_N0,3)/Param.Fs,size(SW_HI_AA_N0,3)), flipud(linspace(-TimeStartW,size(SW_HI_AA_N0,3)/Param.Fs,size(SW_HI_AA_N0,3))')'],[(SW_HI_AA_N0_Mean+SW_HI_AA_N0_SEM), flipud((SW_HI_AA_N0_Mean-SW_HI_AA_N0_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax1,[linspace(-TimeStartW,size(SW_HI_AA_N60,3)/Param.Fs,size(SW_HI_AA_N60,3)), flipud(linspace(-TimeStartW,size(SW_HI_AA_N60,3)/Param.Fs,size(SW_HI_AA_N60,3))')'],[(SW_HI_AA_N60_Mean+SW_HI_AA_N60_SEM), flipud((SW_HI_AA_N60_Mean-SW_HI_AA_N60_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax1,[linspace(-TimeStartW,size(SW_HI_AA_N70,3)/Param.Fs,size(SW_HI_AA_N70,3)), flipud(linspace(-TimeStartW,size(SW_HI_AA_N70,3)/Param.Fs,size(SW_HI_AA_N70,3))')'],[(SW_HI_AA_N70_Mean+SW_HI_AA_N70_SEM), flipud((SW_HI_AA_N70_Mean-SW_HI_AA_N70_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax1,[linspace(-TimeStartW,size(SW_HI_AB_N0,3)/Param.Fs,size(SW_HI_AB_N0,3)), flipud(linspace(-TimeStartW,size(SW_HI_AB_N0,3)/Param.Fs,size(SW_HI_AB_N0,3))')'],[(SW_HI_AB_N0_Mean+SW_HI_AB_N0_SEM), flipud((SW_HI_AB_N0_Mean-SW_HI_AB_N0_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax1,[linspace(-TimeStartW,size(SW_HI_AB_N60,3)/Param.Fs,size(SW_HI_AB_N60,3)), flipud(linspace(-TimeStartW,size(SW_HI_AB_N60,3)/Param.Fs,size(SW_HI_AB_N60,3))')'],[(SW_HI_AB_N60_Mean+SW_HI_AB_N60_SEM), flipud((SW_HI_AB_N60_Mean-SW_HI_AB_N60_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax1,[linspace(-TimeStartW,size(SW_HI_AB_N70,3)/Param.Fs,size(SW_HI_AB_N70,3)), flipud(linspace(-TimeStartW,size(SW_HI_AB_N70,3)/Param.Fs,size(SW_HI_AB_N70,3))')'],[(SW_HI_AB_N70_Mean+SW_HI_AB_N70_SEM), flipud((SW_HI_AB_N70_Mean-SW_HI_AB_N70_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax2,[linspace(-TimeStartW,size(GLW,3)/Param.Fs,size(GLW,3)), flipud(linspace(-TimeStartW,size(GLW,3)/Param.Fs,size(GLW,3))')'],[(GLW_Mean+GLW_SEM), flipud((GLW_Mean-GLW_SEM)')'],SColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax2,[linspace(-TimeStartW,size(LW_N0,3)/Param.Fs,size(LW_N0,3)), flipud(linspace(-TimeStartW,size(LW_N0,3)/Param.Fs,size(LW_N0,3))')'],[(LW_N0_Mean+LW_N0_SEM), flipud((LW_N0_Mean-LW_N0_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax2,[linspace(-TimeStartW,size(LW_N60,3)/Param.Fs,size(LW_N60,3)), flipud(linspace(-TimeStartW,size(LW_N60,3)/Param.Fs,size(LW_N60,3))')'],[(LW_N60_Mean+LW_N60_SEM), flipud((LW_N60_Mean-LW_N60_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax2,[linspace(-TimeStartW,size(LW_N70,3)/Param.Fs,size(LW_N70,3)), flipud(linspace(-TimeStartW,size(LW_N70,3)/Param.Fs,size(LW_N70,3))')'],[(LW_N70_Mean+LW_N70_SEM), flipud((LW_N70_Mean-LW_N70_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax2,[t_LW_NH_N0, flipud(t_LW_NH_N0')'],[(LW_NH_N0_Mean+LW_NH_N0_SEM), flipud((LW_NH_N0_Mean-LW_NH_N0_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax2,[t_LW_NH_N60, flipud(t_LW_NH_N60')'],[(LW_NH_N60_Mean+LW_NH_N60_SEM), flipud((LW_NH_N60_Mean-LW_NH_N60_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax2,[t_LW_NH_N70, flipud(t_LW_NH_N70')'],[(LW_NH_N70_Mean+LW_NH_N70_SEM), flipud((LW_NH_N70_Mean-LW_NH_N70_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax2,[linspace(-TimeStartW,size(LW_HI_UN_N0,3)/Param.Fs,size(LW_HI_UN_N0,3)), flipud(linspace(-TimeStartW,size(LW_HI_UN_N0,3)/Param.Fs,size(LW_HI_UN_N0,3))')'],[(LW_HI_UN_N0_Mean+LW_HI_UN_N0_SEM), flipud((LW_HI_UN_N0_Mean-LW_HI_UN_N0_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax2,[linspace(-TimeStartW,size(LW_HI_UN_N60,3)/Param.Fs,size(LW_HI_UN_N60,3)), flipud(linspace(-TimeStartW,size(LW_HI_UN_N60,3)/Param.Fs,size(LW_HI_UN_N60,3))')'],[(LW_HI_UN_N60_Mean+LW_HI_UN_N60_SEM), flipud((LW_HI_UN_N60_Mean-LW_HI_UN_N60_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax2,[linspace(-TimeStartW,size(LW_HI_UN_N70,3)/Param.Fs,size(LW_HI_UN_N70,3)), flipud(linspace(-TimeStartW,size(LW_HI_UN_N70,3)/Param.Fs,size(LW_HI_UN_N70,3))')'],[(LW_HI_UN_N70_Mean+LW_HI_UN_N70_SEM), flipud((LW_HI_UN_N70_Mean-LW_HI_UN_N70_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax2,[linspace(-TimeStartW,size(LW_HI_AA_N0,3)/Param.Fs,size(LW_HI_AA_N0,3)), flipud(linspace(-TimeStartW,size(LW_HI_AA_N0,3)/Param.Fs,size(LW_HI_AA_N0,3))')'],[(LW_HI_AA_N0_Mean+LW_HI_AA_N0_SEM), flipud((LW_HI_AA_N0_Mean-LW_HI_AA_N0_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax2,[linspace(-TimeStartW,size(LW_HI_AA_N60,3)/Param.Fs,size(LW_HI_AA_N60,3)), flipud(linspace(-TimeStartW,size(LW_HI_AA_N60,3)/Param.Fs,size(LW_HI_AA_N60,3))')'],[(LW_HI_AA_N60_Mean+LW_HI_AA_N60_SEM), flipud((LW_HI_AA_N60_Mean-LW_HI_AA_N60_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax2,[linspace(-TimeStartW,size(LW_HI_AA_N70,3)/Param.Fs,size(LW_HI_AA_N70,3)), flipud(linspace(-TimeStartW,size(LW_HI_AA_N70,3)/Param.Fs,size(LW_HI_AA_N70,3))')'],[(LW_HI_AA_N70_Mean+LW_HI_AA_N70_SEM), flipud((LW_HI_AA_N70_Mean-LW_HI_AA_N70_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax2,[linspace(-TimeStartW,size(LW_HI_AB_N0,3)/Param.Fs,size(LW_HI_AB_N0,3)), flipud(linspace(-TimeStartW,size(LW_HI_AB_N0,3)/Param.Fs,size(LW_HI_AB_N0,3))')'],[(LW_HI_AB_N0_Mean+LW_HI_AB_N0_SEM), flipud((LW_HI_AB_N0_Mean-LW_HI_AB_N0_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax2,[linspace(-TimeStartW,size(LW_HI_AB_N60,3)/Param.Fs,size(LW_HI_AB_N60,3)), flipud(linspace(-TimeStartW,size(LW_HI_AB_N60,3)/Param.Fs,size(LW_HI_AB_N60,3))')'],[(LW_HI_AB_N60_Mean+LW_HI_AB_N60_SEM), flipud((LW_HI_AB_N60_Mean-LW_HI_AB_N60_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax2,[linspace(-TimeStartW,size(LW_HI_AB_N70,3)/Param.Fs,size(LW_HI_AB_N70,3)), flipud(linspace(-TimeStartW,size(LW_HI_AB_N70,3)/Param.Fs,size(LW_HI_AB_N70,3))')'],[(LW_HI_AB_N70_Mean+LW_HI_AB_N70_SEM), flipud((LW_HI_AB_N70_Mean-LW_HI_AB_N70_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% 
% fill(ax3,[linspace(-TimeStartW,size(GSW_B,3)/Param.Fs,size(GSW_B,3)), flipud(linspace(-TimeStartW,size(GSW_B,3)/Param.Fs,size(GSW_B,3))')'],[(GSW_B_Mean+GSW_B_SEM), flipud((GSW_B_Mean-GSW_B_SEM)')'],SColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax3,[linspace(-TimeStartW,size(SW_N0_B,3)/Param.Fs,size(SW_N0_B,3)), flipud(linspace(-TimeStartW,size(SW_N0_B,3)/Param.Fs,size(SW_N0_B,3))')'],[(SW_N0_B_Mean+SW_N0_B_SEM), flipud((SW_N0_B_Mean-SW_N0_B_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax3,[linspace(-TimeStartW,size(SW_N60_B,3)/Param.Fs,size(SW_N60_B,3)), flipud(linspace(-TimeStartW,size(SW_N60_B,3)/Param.Fs,size(SW_N60_B,3))')'],[(SW_N60_B_Mean+SW_N60_B_SEM), flipud((SW_N60_B_Mean-SW_N60_B_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax3,[linspace(-TimeStartW,size(SW_N70_B,3)/Param.Fs,size(SW_N70_B,3)), flipud(linspace(-TimeStartW,size(SW_N70_B,3)/Param.Fs,size(SW_N70_B,3))')'],[(SW_N70_B_Mean+SW_N70_B_SEM), flipud((SW_N70_B_Mean-SW_N70_B_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax3,[t_SW_NH_N0_B, flipud(t_SW_NH_N0_B')'],[(SW_NH_N0_B_Mean+SW_NH_N0_B_SEM), flipud((SW_NH_N0_B_Mean-SW_NH_N0_B_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax3,[t_SW_NH_N60_B, flipud(t_SW_NH_N60_B')'],[(SW_NH_N60_B_Mean+SW_NH_N60_B_SEM), flipud((SW_NH_N60_B_Mean-SW_NH_N60_B_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax3,[t_SW_NH_N70_B, flipud(t_SW_NH_N70_B')'],[(SW_NH_N70_B_Mean+SW_NH_N70_B_SEM), flipud((SW_NH_N70_B_Mean-SW_NH_N70_B_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax3,[linspace(-TimeStartW,size(SW_HI_UN_N0_B,3)/Param.Fs,size(SW_HI_UN_N0_B,3)), flipud(linspace(-TimeStartW,size(SW_HI_UN_N0_B,3)/Param.Fs,size(SW_HI_UN_N0_B,3))')'],[(SW_HI_UN_N0_B_Mean+SW_HI_UN_N0_B_SEM), flipud((SW_HI_UN_N0_B_Mean-SW_HI_UN_N0_B_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax3,[linspace(-TimeStartW,size(SW_HI_UN_N60_B,3)/Param.Fs,size(SW_HI_UN_N60_B,3)), flipud(linspace(-TimeStartW,size(SW_HI_UN_N60_B,3)/Param.Fs,size(SW_HI_UN_N60_B,3))')'],[(SW_HI_UN_N60_B_Mean+SW_HI_UN_N60_B_SEM), flipud((SW_HI_UN_N60_B_Mean-SW_HI_UN_N60_B_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax3,[linspace(-TimeStartW,size(SW_HI_UN_N70_B,3)/Param.Fs,size(SW_HI_UN_N70_B,3)), flipud(linspace(-TimeStartW,size(SW_HI_UN_N70_B,3)/Param.Fs,size(SW_HI_UN_N70_B,3))')'],[(SW_HI_UN_N70_B_Mean+SW_HI_UN_N70_B_SEM), flipud((SW_HI_UN_N70_B_Mean-SW_HI_UN_N70_B_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax3,[linspace(-TimeStartW,size(SW_HI_AA_N0_B,3)/Param.Fs,size(SW_HI_AA_N0_B,3)), flipud(linspace(-TimeStartW,size(SW_HI_AA_N0_B,3)/Param.Fs,size(SW_HI_AA_N0_B,3))')'],[(SW_HI_AA_N0_B_Mean+SW_HI_AA_N0_B_SEM), flipud((SW_HI_AA_N0_B_Mean-SW_HI_AA_N0_B_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax3,[linspace(-TimeStartW,size(SW_HI_AA_N60_B,3)/Param.Fs,size(SW_HI_AA_N60_B,3)), flipud(linspace(-TimeStartW,size(SW_HI_AA_N60_B,3)/Param.Fs,size(SW_HI_AA_N60_B,3))')'],[(SW_HI_AA_N60_B_Mean+SW_HI_AA_N60_B_SEM), flipud((SW_HI_AA_N60_B_Mean-SW_HI_AA_N60_B_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax3,[linspace(-TimeStartW,size(SW_HI_AA_N70_B,3)/Param.Fs,size(SW_HI_AA_N70_B,3)), flipud(linspace(-TimeStartW,size(SW_HI_AA_N70_B,3)/Param.Fs,size(SW_HI_AA_N70_B,3))')'],[(SW_HI_AA_N70_B_Mean+SW_HI_AA_N70_B_SEM), flipud((SW_HI_AA_N70_B_Mean-SW_HI_AA_N70_B_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax3,[linspace(-TimeStartW,size(SW_HI_AB_N0_B,3)/Param.Fs,size(SW_HI_AB_N0_B,3)), flipud(linspace(-TimeStartW,size(SW_HI_AB_N0_B,3)/Param.Fs,size(SW_HI_AB_N0_B,3))')'],[(SW_HI_AB_N0_B_Mean+SW_HI_AB_N0_B_SEM), flipud((SW_HI_AB_N0_B_Mean-SW_HI_AB_N0_B_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax3,[linspace(-TimeStartW,size(SW_HI_AB_N60_B,3)/Param.Fs,size(SW_HI_AB_N60_B,3)), flipud(linspace(-TimeStartW,size(SW_HI_AB_N60_B,3)/Param.Fs,size(SW_HI_AB_N60_B,3))')'],[(SW_HI_AB_N60_B_Mean+SW_HI_AB_N60_B_SEM), flipud((SW_HI_AB_N60_B_Mean-SW_HI_AB_N60_B_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax3,[linspace(-TimeStartW,size(SW_HI_AB_N70_B,3)/Param.Fs,size(SW_HI_AB_N70_B,3)), flipud(linspace(-TimeStartW,size(SW_HI_AB_N70_B,3)/Param.Fs,size(SW_HI_AB_N70_B,3))')'],[(SW_HI_AB_N70_B_Mean+SW_HI_AB_N70_B_SEM), flipud((SW_HI_AB_N70_B_Mean-SW_HI_AB_N70_B_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax4,[linspace(-TimeStartW,size(GLW_B,3)/Param.Fs,size(GLW_B,3)), flipud(linspace(-TimeStartW,size(GLW_B,3)/Param.Fs,size(GLW_B,3))')'],[(GLW_B_Mean+GLW_B_SEM), flipud((GLW_B_Mean-GLW_B_SEM)')'],SColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax4,[linspace(-TimeStartW,size(LW_N0_B,3)/Param.Fs,size(LW_N0_B,3)), flipud(linspace(-TimeStartW,size(LW_N0_B,3)/Param.Fs,size(LW_N0_B,3))')'],[(LW_N0_B_Mean+LW_N0_B_SEM), flipud((LW_N0_B_Mean-LW_N0_B_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax4,[linspace(-TimeStartW,size(LW_N60_B,3)/Param.Fs,size(LW_N60_B,3)), flipud(linspace(-TimeStartW,size(LW_N60_B,3)/Param.Fs,size(LW_N60_B,3))')'],[(LW_N60_B_Mean+LW_N60_B_SEM), flipud((LW_N60_B_Mean-LW_N60_B_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax4,[linspace(-TimeStartW,size(LW_N70_B,3)/Param.Fs,size(LW_N70_B,3)), flipud(linspace(-TimeStartW,size(LW_N70_B,3)/Param.Fs,size(LW_N70_B,3))')'],[(LW_N70_B_Mean+LW_N70_B_SEM), flipud((LW_N70_B_Mean-LW_N70_B_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax4,[t_LW_NH_N0_B, flipud(t_LW_NH_N0_B')'],[(LW_NH_N0_B_Mean+LW_NH_N0_B_SEM), flipud((LW_NH_N0_B_Mean-LW_NH_N0_B_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax4,[t_LW_NH_N60_B, flipud(t_LW_NH_N60_B')'],[(LW_NH_N60_B_Mean+LW_NH_N60_B_SEM), flipud((LW_NH_N60_B_Mean-LW_NH_N60_B_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax4,[t_LW_NH_N70_B, flipud(t_LW_NH_N70_B')'],[(LW_NH_N70_B_Mean+LW_NH_N70_B_SEM), flipud((LW_NH_N70_B_Mean-LW_NH_N70_B_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax4,[linspace(-TimeStartW,size(LW_HI_UN_N0_B,3)/Param.Fs,size(LW_HI_UN_N0_B,3)), flipud(linspace(-TimeStartW,size(LW_HI_UN_N0_B,3)/Param.Fs,size(LW_HI_UN_N0_B,3))')'],[(LW_HI_UN_N0_B_Mean+LW_HI_UN_N0_B_SEM), flipud((LW_HI_UN_N0_B_Mean-LW_HI_UN_N0_B_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax4,[linspace(-TimeStartW,size(LW_HI_UN_N60_B,3)/Param.Fs,size(LW_HI_UN_N60_B,3)), flipud(linspace(-TimeStartW,size(LW_HI_UN_N60_B,3)/Param.Fs,size(LW_HI_UN_N60_B,3))')'],[(LW_HI_UN_N60_B_Mean+LW_HI_UN_N60_B_SEM), flipud((LW_HI_UN_N60_B_Mean-LW_HI_UN_N60_B_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax4,[linspace(-TimeStartW,size(LW_HI_UN_N70_B,3)/Param.Fs,size(LW_HI_UN_N70_B,3)), flipud(linspace(-TimeStartW,size(LW_HI_UN_N70_B,3)/Param.Fs,size(LW_HI_UN_N70_B,3))')'],[(LW_HI_UN_N70_B_Mean+LW_HI_UN_N70_B_SEM), flipud((LW_HI_UN_N70_B_Mean-LW_HI_UN_N70_B_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax4,[linspace(-TimeStartW,size(LW_HI_AA_N0_B,3)/Param.Fs,size(LW_HI_AA_N0_B,3)), flipud(linspace(-TimeStartW,size(LW_HI_AA_N0_B,3)/Param.Fs,size(LW_HI_AA_N0_B,3))')'],[(LW_HI_AA_N0_B_Mean+LW_HI_AA_N0_B_SEM), flipud((LW_HI_AA_N0_B_Mean-LW_HI_AA_N0_B_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax4,[linspace(-TimeStartW,size(LW_HI_AA_N60_B,3)/Param.Fs,size(LW_HI_AA_N60_B,3)), flipud(linspace(-TimeStartW,size(LW_HI_AA_N60_B,3)/Param.Fs,size(LW_HI_AA_N60_B,3))')'],[(LW_HI_AA_N60_B_Mean+LW_HI_AA_N60_B_SEM), flipud((LW_HI_AA_N60_B_Mean-LW_HI_AA_N60_B_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax4,[linspace(-TimeStartW,size(LW_HI_AA_N70_B,3)/Param.Fs,size(LW_HI_AA_N70_B,3)), flipud(linspace(-TimeStartW,size(LW_HI_AA_N70_B,3)/Param.Fs,size(LW_HI_AA_N70_B,3))')'],[(LW_HI_AA_N70_B_Mean+LW_HI_AA_N70_B_SEM), flipud((LW_HI_AA_N70_B_Mean-LW_HI_AA_N70_B_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax4,[linspace(-TimeStartW,size(LW_HI_AB_N0_B,3)/Param.Fs,size(LW_HI_AB_N0_B,3)), flipud(linspace(-TimeStartW,size(LW_HI_AB_N0_B,3)/Param.Fs,size(LW_HI_AB_N0_B,3))')'],[(LW_HI_AB_N0_B_Mean+LW_HI_AB_N0_B_SEM), flipud((LW_HI_AB_N0_B_Mean-LW_HI_AB_N0_B_SEM)')'],QuietColor,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax4,[linspace(-TimeStartW,size(LW_HI_AB_N60_B,3)/Param.Fs,size(LW_HI_AB_N60_B,3)), flipud(linspace(-TimeStartW,size(LW_HI_AB_N60_B,3)/Param.Fs,size(LW_HI_AB_N60_B,3))')'],[(LW_HI_AB_N60_B_Mean+LW_HI_AB_N60_B_SEM), flipud((LW_HI_AB_N60_B_Mean-LW_HI_AB_N60_B_SEM)')'],N60Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')
% fill(ax4,[linspace(-TimeStartW,size(LW_HI_AB_N70_B,3)/Param.Fs,size(LW_HI_AB_N70_B,3)), flipud(linspace(-TimeStartW,size(LW_HI_AB_N70_B,3)/Param.Fs,size(LW_HI_AB_N70_B,3))')'],[(LW_HI_AB_N70_B_Mean+LW_HI_AB_N70_B_SEM), flipud((LW_HI_AB_N70_B_Mean-LW_HI_AB_N70_B_SEM)')'],N70Color,'FaceAlpha',.2,'Edgecolor','none','handlevisibility' ,'off')

xlabel([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8],'Time [s]')
ylabel([ax1 ax2 ax5 ax6],'Pupil diameter [mm]')
ylabel([ax3 ax4 ax7 ax8],'Pupil baseline difference [mm]')
xlim([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8],[-TimeStartW 3])
xlim(ax1,[-TimeStartW, min(t_SW_NH(cumsum(squeeze(sum(~isnan(SW_NH),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(SW_NH),[1 2])))))])
xlim(ax2,[-TimeStartW, min(t_LW_NH(cumsum(squeeze(sum(~isnan(LW_NH),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(LW_NH),[1 2])))))])
xlim(ax3,[-TimeStartW, min(t_SW_HI(cumsum(squeeze(sum(~isnan(SW_HI),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(SW_HI),[1 2])))))])
xlim(ax4,[-TimeStartW, min(t_SW_NH(cumsum(squeeze(sum(~isnan(LW_HI),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(LW_HI),[1 2])))))])
xlim(ax5,[-TimeStartW, min(t_SW_NH_B(cumsum(squeeze(sum(~isnan(SW_NH_B),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(SW_NH_B),[1 2])))))])
xlim(ax6,[-TimeStartW, min(t_LW_NH_B(cumsum(squeeze(sum(~isnan(LW_NH_B),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(LW_NH_B),[1 2])))))])
xlim(ax7,[-TimeStartW, min(t_SW_HI_B(cumsum(squeeze(sum(~isnan(SW_HI_B),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(SW_HI_B),[1 2])))))])
xlim(ax8,[-TimeStartW, min(t_SW_NH_B(cumsum(squeeze(sum(~isnan(LW_HI_B),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(LW_HI_B),[1 2])))))])
ylmax5_8=cell2mat(ylim([ax5 ax6 ax7 ax8]));
[max5_8,idx5_8] = max(diff(ylmax5_8,1,2));
% ylim([ax5 ax6 ax7 ax8],ylmax5_8(idx5_8,:)) % Setting ylim to the max out of all
lgd4=legend(ax4,'Speaking','Listening','N0','N60','N70','Location','southeastoutside');
lgd4.Title.String = 'Types of windows:';
lgd8=legend(ax8,'Speaking','Listening','N0','N60','N70','Location','southeastoutside');
lgd8.Title.String = 'Types of windows:';
title([ax1 ax5],'NH - Speaking windows')
title([ax2 ax6],'NH - Listening windows')
title([ax3 ax7],'HI - Speaking windows')
title([ax4 ax8],'HI - Listening windows')
sgtitle(f1,'Global response separated by noise conditions')
sgtitle(f2,'Global adaptive baselined response separated by noise conditions')
sgtitle(f3,'')
sgtitle(f4,'')

%% NOTES
sum(sum(sum(~isnan(GSW),3)>0,2)); % Used to know how many windows (non-nan rows) are inside each variable (Layers,Rows,Columns)