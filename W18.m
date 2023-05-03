%% PBCA-Thesis - Week 18 - AMEND II data - Pupillometry
% Pathing
clear all; clc; close all;
BPath = strsplit(pwd,'PBCA-thesis');
addpath('tools')
addpath([BPath{1} 'Pupil-preprocessing-tools\tools']) % For preprocessing

% Colors
SColor = [53, 155, 67]./255;
LColor = [204, 36, 0]./255;
QuietColor = [204, 152, 0]./255;
SHLColor = [123, 31, 162]./255;
N60Color = [0, 196, 215]./255;
N70Color = [2, 36, 223]./255;

% Important variables
[subDirs_I] = GetSubDirsFirstLevelOnly('data\AMEND_I');
[subDirs_II] = GetSubDirsFirstLevelOnly('data\AMEND_II');
subDirs_II(contains(subDirs_II,{'Pilot 1','Pair14'}))=[]; % Only using data from Pair01 to Pair13, others removed
FileNames_I={'P1_N0_B1.mat','P1_N0_B2.mat','P1_N60_B1.mat','P1_N60_B2.mat','P1_Noise60_B1.mat','P1_Noise60_B2.mat','P1_Noise70_B1.mat','P1_Noise70_B2.mat','P2_N0_B1.mat','P2_N0_B2.mat','P2_N60_B1.mat','P2_N60_B2.mat','P2_Noise60_B1.mat','P2_Noise60_B2.mat','P2_Noise70_B1.mat','P2_Noise70_B2.mat'};
FileNames_II={'UNHI_N0.mat','UNHI_N60.mat','UNHI_N70.mat','AAHI_N0.mat','AAHI_N60.mat','AAHI_N70.mat','ABHI_N0.mat','ABHI_N60.mat','ABHI_N70.mat'};
LoadUtt_I=load('data\AMEND_I\utterances1110.mat');
LoadUtt_II=load('data\AMEND_II\utterances.mat');
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

figure;tiledlayout(1,2);ax1 = nexttile;ax2 = nexttile;
figure;tiledlayout(1,2);ax3 = nexttile;ax4 = nexttile;

hold([ax1 ax2 ax3 ax4],'on')
grid([ax1 ax2 ax3 ax4],'on')

% Loop for AMEND II
for q=1:numel(subDirs_II)
    PairIn_II = q;
    for ChosenFolder = {'\NH\','\HI\'}
        PairFolder_II=[pwd,'\data\AMEND_II\',cell2mat(subDirs_II(q)),cell2mat(ChosenFolder)]; % Folder naming changed
        PairFiles_II=dir(PairFolder_II); % Folder naming changed
        try
            PairUtt_II=LoadUtt_II.Utterances(PairIn_II,:);
        catch ME
            disp(['Warning: No Utterance found for folder "',cell2mat(subDirs_II(q)),'"'])
            continue
        end
%         PairDelay_II=LoadDelays_II.TobAudDelay(PairIn_II,:);
    
        for i=1:numel(FileNames_II)
            if isempty(PairUtt_II{q})
                disp(['Warning: File ', PairFiles_II(1).folder, '\', cell2mat(FileNames_II(i)), ' was skipped, utterance not found.'])
            end
            
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

            % Extract delay [s] (Duration_timeStamps - Duration_nsamples)
            EyeAudDelay=alldata_mat(end).timeStamp-alldata_mat(1).timeStamp-length(alldata_mat)/Param.Fs;

            % Skip file if the difference in duration from the number of
            % samples and the duration given by timestamps is bigger than 0.5 s
            if EyeAudDelay > RejectDelay
                disp(['Warning: File ',PairFiles_II(1).folder, '\', cell2mat(FileNames_II(i)), ' was rejected, too much delay (',sprintf('%0.2f',EyeAudDelay),'s).']);
                continue
            end

            LDiamRaw = [alldata_mat.diameterLeft];
            RDiamRaw = [alldata_mat.diameterRight];

            % Preprocessing - Setting outliers as NaNs (remove artifacts)
            LThreshOut = [mean(LDiamRaw,'omitnan')-std(LDiamRaw,'omitnan'),mean(LDiamRaw,'omitnan')+std(LDiamRaw,'omitnan')];
            RThreshOut = [mean(RDiamRaw,'omitnan')-std(RDiamRaw,'omitnan'),mean(RDiamRaw,'omitnan')+std(RDiamRaw,'omitnan')];
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
            
            try
                SpeakRaw = PairUtt_II{1,SpeCond}.(SpeakKey);
                ListenRaw = PairUtt_II{1,SpeCond}.(ListenKey);
                binResUtt = PairUtt_II{1,SpeCond}.binRes;
            catch ME % if isempty(SpeakRaw) && isempty(ListenRaw)
                disp(['Warning: File ',PairFiles_II(1).folder, '\', cell2mat(FileNames_II(i)),' was rejected for not having associated Utterance windows.']);
                continue
            end

            % SAME PROCESSING AS IN W1.m
            % Downsample (rounding) Utt from 250 Hz (1/binRes) to 50 Hz, no
            % shift in time nor delay added this time
            SpeakRaw(:,2:3)=round((SpeakRaw(:,2:3)*binResUtt)*Param.Fs);
            ListenRaw(:,2:3)=round((ListenRaw(:,2:3)*binResUtt)*Param.Fs);

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
                if contains(cell2mat(FileNames_II(i)),'UN')
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

% Calculate SEM as: std(X)/sqrt(numel(X(~isnan(X)))
GSW_SEM = (reshape(std(GSW,0,[1 2],'omitnan'),[],1)/sqrt(numel(GSW(~isnan(GSW)))))';
GSW_B_SEM = (reshape(std(GSW_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(GSW_B(~isnan(GSW_B)))))';
GLW_SEM = (reshape(std(GLW,0,[1 2],'omitnan'),[],1)/sqrt(numel(GLW(~isnan(GLW)))))';
GLW_B_SEM = (reshape(std(GLW_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(GLW_B(~isnan(GLW_B)))))';

SW_NH_SEM = (reshape(std(SW_NH,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_NH(~isnan(SW_NH)))))';
SW_NH_B_SEM = (reshape(std(SW_NH_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_NH_B(~isnan(SW_NH_B)))))';
LW_NH_SEM = (reshape(std(LW_NH,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_NH(~isnan(LW_NH)))))';
LW_NH_B_SEM = (reshape(std(LW_NH_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_NH_B(~isnan(LW_NH_B)))))';
SW_HI_SEM = (reshape(std(SW_HI,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_HI(~isnan(SW_HI)))))';
SW_HI_B_SEM = (reshape(std(SW_HI_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_HI_B(~isnan(SW_HI_B)))))';
LW_HI_SEM = (reshape(std(LW_HI,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_HI(~isnan(LW_HI)))))';
LW_HI_B_SEM = (reshape(std(LW_HI_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_HI_B(~isnan(LW_HI_B)))))';

SW_N0_SEM = (reshape(std(SW_N0,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_N0(~isnan(SW_N0)))))';
SW_N0_B_SEM = (reshape(std(SW_N0_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_N0_B(~isnan(SW_N0_B)))))';
LW_N0_SEM = (reshape(std(LW_N0,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_N0(~isnan(LW_N0)))))';
LW_N0_B_SEM = (reshape(std(LW_N0_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_N0_B(~isnan(LW_N0_B)))))';
SW_N60_SEM = (reshape(std(SW_N60,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_N60(~isnan(SW_N60)))))';
SW_N60_B_SEM = (reshape(std(SW_N60_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_N60_B(~isnan(SW_N60_B)))))';
LW_N60_SEM = (reshape(std(LW_N60,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_N60(~isnan(LW_N60)))))';
LW_N60_B_SEM = (reshape(std(LW_N60_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_N60_B(~isnan(LW_N60_B)))))';
SW_N70_SEM = (reshape(std(SW_N70,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_N70(~isnan(SW_N70)))))';
SW_N70_B_SEM = (reshape(std(SW_N70_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_N70_B(~isnan(SW_N70_B)))))';
LW_N70_SEM = (reshape(std(LW_N70,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_N70(~isnan(LW_N70)))))';
LW_N70_B_SEM = (reshape(std(LW_N70_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_N70_B(~isnan(LW_N70_B)))))';

SW_NH_N0_SEM = (reshape(std(SW_NH_N0,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_NH_N0(~isnan(SW_NH_N0)))))';
SW_NH_N0_B_SEM = (reshape(std(SW_NH_N0_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_NH_N0_B(~isnan(SW_NH_N0_B)))))';
LW_NH_N0_SEM = (reshape(std(LW_NH_N0,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_NH_N0(~isnan(LW_NH_N0)))))';
LW_NH_N0_B_SEM = (reshape(std(LW_NH_N0_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_NH_N0_B(~isnan(LW_NH_N0_B)))))';
SW_NH_N60_SEM = (reshape(std(SW_NH_N60,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_NH_N60(~isnan(SW_NH_N60)))))';
SW_NH_N60_B_SEM = (reshape(std(SW_NH_N60_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_NH_N60_B(~isnan(SW_NH_N60_B)))))';
LW_NH_N60_SEM = (reshape(std(LW_NH_N60,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_NH_N60(~isnan(LW_NH_N60)))))';
LW_NH_N60_B_SEM = (reshape(std(LW_NH_N60_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_NH_N60_B(~isnan(LW_NH_N60_B)))))';
SW_NH_N70_SEM = (reshape(std(SW_NH_N70,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_NH_N70(~isnan(SW_NH_N70)))))';
SW_NH_N70_B_SEM = (reshape(std(SW_NH_N70_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_NH_N70_B(~isnan(SW_NH_N70_B)))))';
LW_NH_N70_SEM = (reshape(std(LW_NH_N70,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_NH_N70(~isnan(LW_NH_N70)))))';
LW_NH_N70_B_SEM = (reshape(std(LW_NH_N70_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_NH_N70_B(~isnan(LW_NH_N70_B)))))';

SW_HI_UN_N0_SEM = (reshape(std(SW_HI_UN_N0,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_HI_UN_N0(~isnan(SW_HI_UN_N0)))))';
SW_HI_UN_N0_B_SEM = (reshape(std(SW_HI_UN_N0_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_HI_UN_N0_B(~isnan(SW_HI_UN_N0_B)))))';
LW_HI_UN_N0_SEM = (reshape(std(LW_HI_UN_N0,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_HI_UN_N0(~isnan(LW_HI_UN_N0)))))';
LW_HI_UN_N0_B_SEM = (reshape(std(LW_HI_UN_N0_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_HI_UN_N0_B(~isnan(LW_HI_UN_N0_B)))))';
SW_HI_UN_N60_SEM = (reshape(std(SW_HI_UN_N60,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_HI_UN_N60(~isnan(SW_HI_UN_N60)))))';
SW_HI_UN_N60_B_SEM = (reshape(std(SW_HI_UN_N60_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_HI_UN_N60_B(~isnan(SW_HI_UN_N60_B)))))';
LW_HI_UN_N60_SEM = (reshape(std(LW_HI_UN_N60,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_HI_UN_N60(~isnan(LW_HI_UN_N60)))))';
LW_HI_UN_N60_B_SEM = (reshape(std(LW_HI_UN_N60_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_HI_UN_N60_B(~isnan(LW_HI_UN_N60_B)))))';
SW_HI_UN_N70_SEM = (reshape(std(SW_HI_UN_N70,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_HI_UN_N70(~isnan(SW_HI_UN_N70)))))';
SW_HI_UN_N70_B_SEM = (reshape(std(SW_HI_UN_N70_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_HI_UN_N70_B(~isnan(SW_HI_UN_N70_B)))))';
LW_HI_UN_N70_SEM = (reshape(std(LW_HI_UN_N70,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_HI_UN_N70(~isnan(LW_HI_UN_N70)))))';
LW_HI_UN_N70_B_SEM = (reshape(std(LW_HI_UN_N70_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_HI_UN_N70_B(~isnan(LW_HI_UN_N70_B)))))';
SW_HI_AA_N0_SEM = (reshape(std(SW_HI_AA_N0,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_HI_AA_N0(~isnan(SW_HI_AA_N0)))))';
SW_HI_AA_N0_B_SEM = (reshape(std(SW_HI_AA_N0_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_HI_AA_N0_B(~isnan(SW_HI_AA_N0_B)))))';
LW_HI_AA_N0_SEM = (reshape(std(LW_HI_AA_N0,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_HI_AA_N0(~isnan(LW_HI_AA_N0)))))';
LW_HI_AA_N0_B_SEM = (reshape(std(LW_HI_AA_N0_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_HI_AA_N0_B(~isnan(LW_HI_AA_N0_B)))))';
SW_HI_AA_N60_SEM = (reshape(std(SW_HI_AA_N60,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_HI_AA_N60(~isnan(SW_HI_AA_N60)))))';
SW_HI_AA_N60_B_SEM = (reshape(std(SW_HI_AA_N60_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_HI_AA_N60_B(~isnan(SW_HI_AA_N60_B)))))';
LW_HI_AA_N60_SEM = (reshape(std(LW_HI_AA_N60,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_HI_AA_N60(~isnan(LW_HI_AA_N60)))))';
LW_HI_AA_N60_B_SEM = (reshape(std(LW_HI_AA_N60_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_HI_AA_N60_B(~isnan(LW_HI_AA_N60_B)))))';
SW_HI_AA_N70_SEM = (reshape(std(SW_HI_AA_N70,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_HI_AA_N70(~isnan(SW_HI_AA_N70)))))';
SW_HI_AA_N70_B_SEM = (reshape(std(SW_HI_AA_N70_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_HI_AA_N70_B(~isnan(SW_HI_AA_N70_B)))))';
LW_HI_AA_N70_SEM = (reshape(std(LW_HI_AA_N70,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_HI_AA_N70(~isnan(LW_HI_AA_N70)))))';
LW_HI_AA_N70_B_SEM = (reshape(std(LW_HI_AA_N70_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_HI_AA_N70_B(~isnan(LW_HI_AA_N70_B)))))';
SW_HI_AB_N0_SEM = (reshape(std(SW_HI_AB_N0,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_HI_AB_N0(~isnan(SW_HI_AB_N0)))))';
SW_HI_AB_N0_B_SEM = (reshape(std(SW_HI_AB_N0_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_HI_AB_N0_B(~isnan(SW_HI_AB_N0_B)))))';
LW_HI_AB_N0_SEM = (reshape(std(LW_HI_AB_N0,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_HI_AB_N0(~isnan(LW_HI_AB_N0)))))';
LW_HI_AB_N0_B_SEM = (reshape(std(LW_HI_AB_N0_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_HI_AB_N0_B(~isnan(LW_HI_AB_N0_B)))))';
SW_HI_AB_N60_SEM = (reshape(std(SW_HI_AB_N60,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_HI_AB_N60(~isnan(SW_HI_AB_N60)))))';
SW_HI_AB_N60_B_SEM = (reshape(std(SW_HI_AB_N60_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_HI_AB_N60_B(~isnan(SW_HI_AB_N60_B)))))';
LW_HI_AB_N60_SEM = (reshape(std(LW_HI_AB_N60,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_HI_AB_N60(~isnan(LW_HI_AB_N60)))))';
LW_HI_AB_N60_B_SEM = (reshape(std(LW_HI_AB_N60_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_HI_AB_N60_B(~isnan(LW_HI_AB_N60_B)))))';
SW_HI_AB_N70_SEM = (reshape(std(SW_HI_AB_N70,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_HI_AB_N70(~isnan(SW_HI_AB_N70)))))';
SW_HI_AB_N70_B_SEM = (reshape(std(SW_HI_AB_N70_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_HI_AB_N70_B(~isnan(SW_HI_AB_N70_B)))))';
LW_HI_AB_N70_SEM = (reshape(std(LW_HI_AB_N70,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_HI_AB_N70(~isnan(LW_HI_AB_N70)))))';
LW_HI_AB_N70_B_SEM = (reshape(std(LW_HI_AB_N70_B,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_HI_AB_N70_B(~isnan(LW_HI_AB_N70_B)))))';

%% Global Plots
% Plot event onset and baseline markers
xline(ax1,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax2,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax3,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax4,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')

% xline(ax1,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
% xline(ax2,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
xline(ax3,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
xline(ax4,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')

plot(ax1,linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3)),GSW_Mean,color=SColor,linewidth=2)
plot(ax1,linspace(-TimeStartW,size(SW_N0,3)/Param.Fs,size(SW_N0,3)),SW_N0_Mean,color=QuietColor)
plot(ax1,linspace(-TimeStartW,size(SW_N60,3)/Param.Fs,size(SW_N60,3)),SW_N60_Mean,color=N60Color)
plot(ax1,linspace(-TimeStartW,size(SW_N70,3)/Param.Fs,size(SW_N70,3)),SW_N70_Mean,color=N70Color)
plot(ax1,linspace(-TimeStartW,size(SW_NH_N0,3)/Param.Fs,size(SW_NH_N0,3)),SW_NH_N0_Mean,color=QuietColor,linestyle=':')
plot(ax1,linspace(-TimeStartW,size(SW_NH_N60,3)/Param.Fs,size(SW_NH_N60,3)),SW_NH_N60_Mean,color=N60Color,linestyle=':')
plot(ax1,linspace(-TimeStartW,size(SW_NH_N70,3)/Param.Fs,size(SW_NH_N70,3)),SW_NH_N70_Mean,color=N70Color,linestyle=':')
plot(ax1,linspace(-TimeStartW,size(SW_HI_UN_N0,3)/Param.Fs,size(SW_HI_UN_N0,3)),SW_HI_UN_N0_Mean,color=QuietColor,linestyle='--')
plot(ax1,linspace(-TimeStartW,size(SW_HI_UN_N60,3)/Param.Fs,size(SW_HI_UN_N60,3)),SW_HI_UN_N60_Mean,color=N60Color,linestyle='--')
plot(ax1,linspace(-TimeStartW,size(SW_HI_UN_N70,3)/Param.Fs,size(SW_HI_UN_N70,3)),SW_HI_UN_N70_Mean,color=N70Color,linestyle='--')
plot(ax1,linspace(-TimeStartW,size(SW_HI_AA_N0,3)/Param.Fs,size(SW_HI_AA_N0,3)),SW_HI_AA_N0_Mean,color=QuietColor,linestyle='-.')
plot(ax1,linspace(-TimeStartW,size(SW_HI_AA_N60,3)/Param.Fs,size(SW_HI_AA_N60,3)),SW_HI_AA_N60_Mean,color=N60Color,linestyle='-.')
plot(ax1,linspace(-TimeStartW,size(SW_HI_AA_N70,3)/Param.Fs,size(SW_HI_AA_N70,3)),SW_HI_AA_N70_Mean,color=N70Color,linestyle='-.')
plot(ax1,linspace(-TimeStartW,size(SW_HI_AB_N0,3)/Param.Fs,size(SW_HI_AB_N0,3)),SW_HI_AB_N0_Mean,color=QuietColor,linestyle='-.')
plot(ax1,linspace(-TimeStartW,size(SW_HI_AB_N60,3)/Param.Fs,size(SW_HI_AB_N60,3)),SW_HI_AB_N60_Mean,color=N60Color,linestyle='-.')
plot(ax1,linspace(-TimeStartW,size(SW_HI_AB_N70,3)/Param.Fs,size(SW_HI_AB_N70,3)),SW_HI_AB_N70_Mean,color=N70Color,linestyle='-.')
plot(ax2,nan,color=SColor,linewidth=2) % plot nans to show color in legend
plot(ax2,linspace(-TimeStartW,size(GLW,3)/Param.Fs,size(GLW,3)),GLW_Mean,color=LColor,linewidth=2)
plot(ax2,linspace(-TimeStartW,size(LW_N0,3)/Param.Fs,size(LW_N0,3)),LW_N0_Mean,color=QuietColor)
plot(ax2,linspace(-TimeStartW,size(LW_N60,3)/Param.Fs,size(LW_N60,3)),LW_N60_Mean,color=N60Color)
plot(ax2,linspace(-TimeStartW,size(LW_N70,3)/Param.Fs,size(LW_N70,3)),LW_N70_Mean,color=N70Color)
plot(ax2,linspace(-TimeStartW,size(LW_NH_N0,3)/Param.Fs,size(LW_NH_N0,3)),LW_NH_N0_Mean,color=QuietColor,linestyle=':')
plot(ax2,linspace(-TimeStartW,size(LW_NH_N60,3)/Param.Fs,size(LW_NH_N60,3)),LW_NH_N60_Mean,color=N60Color,linestyle=':')
plot(ax2,linspace(-TimeStartW,size(LW_NH_N70,3)/Param.Fs,size(LW_NH_N70,3)),LW_NH_N70_Mean,color=N70Color,linestyle=':')
plot(ax2,linspace(-TimeStartW,size(LW_HI_UN_N0,3)/Param.Fs,size(LW_HI_UN_N0,3)),LW_HI_UN_N0_Mean,color=QuietColor,linestyle='--')
plot(ax2,linspace(-TimeStartW,size(LW_HI_UN_N60,3)/Param.Fs,size(LW_HI_UN_N60,3)),LW_HI_UN_N60_Mean,color=N60Color,linestyle='--')
plot(ax2,linspace(-TimeStartW,size(LW_HI_UN_N70,3)/Param.Fs,size(LW_HI_UN_N70,3)),LW_HI_UN_N70_Mean,color=N70Color,linestyle='--')
plot(ax2,linspace(-TimeStartW,size(LW_HI_AA_N0,3)/Param.Fs,size(LW_HI_AA_N0,3)),LW_HI_AA_N0_Mean,color=QuietColor,linestyle='-.')
plot(ax2,linspace(-TimeStartW,size(LW_HI_AA_N60,3)/Param.Fs,size(LW_HI_AA_N60,3)),LW_HI_AA_N60_Mean,color=N60Color,linestyle='-.')
plot(ax2,linspace(-TimeStartW,size(LW_HI_AA_N70,3)/Param.Fs,size(LW_HI_AA_N70,3)),LW_HI_AA_N70_Mean,color=N70Color,linestyle='-.')
plot(ax2,linspace(-TimeStartW,size(LW_HI_AB_N0,3)/Param.Fs,size(LW_HI_AB_N0,3)),LW_HI_AB_N0_Mean,color=QuietColor,linestyle='-.')
plot(ax2,linspace(-TimeStartW,size(LW_HI_AB_N60,3)/Param.Fs,size(LW_HI_AB_N60,3)),LW_HI_AB_N60_Mean,color=N60Color,linestyle='-.')
plot(ax2,linspace(-TimeStartW,size(LW_HI_AB_N70,3)/Param.Fs,size(LW_HI_AB_N70,3)),LW_HI_AB_N70_Mean,color=N70Color,linestyle='-.')

plot(ax3,linspace(-TimeStartW,size(GSW_B,3)/Param.Fs,size(GSW_B,3)),GSW_B_Mean,color=SColor,linewidth=2)
plot(ax3,linspace(-TimeStartW,size(SW_N0_B,3)/Param.Fs,size(SW_N0_B,3)),SW_N0_B_Mean,color=QuietColor)
plot(ax3,linspace(-TimeStartW,size(SW_N60_B,3)/Param.Fs,size(SW_N60_B,3)),SW_N60_B_Mean,color=N60Color)
plot(ax3,linspace(-TimeStartW,size(SW_N70_B,3)/Param.Fs,size(SW_N70_B,3)),SW_N70_B_Mean,color=N70Color)
plot(ax3,linspace(-TimeStartW,size(SW_NH_N0_B,3)/Param.Fs,size(SW_NH_N0_B,3)),SW_NH_N0_B_Mean,color=QuietColor,linestyle=':')
plot(ax3,linspace(-TimeStartW,size(SW_NH_N60_B,3)/Param.Fs,size(SW_NH_N60_B,3)),SW_NH_N60_B_Mean,color=N60Color,linestyle=':')
plot(ax3,linspace(-TimeStartW,size(SW_NH_N70_B,3)/Param.Fs,size(SW_NH_N70_B,3)),SW_NH_N70_B_Mean,color=N70Color,linestyle=':')
plot(ax3,linspace(-TimeStartW,size(SW_HI_UN_N0_B,3)/Param.Fs,size(SW_HI_UN_N0_B,3)),SW_HI_UN_N0_B_Mean,color=QuietColor,linestyle='--')
plot(ax3,linspace(-TimeStartW,size(SW_HI_UN_N60_B,3)/Param.Fs,size(SW_HI_UN_N60_B,3)),SW_HI_UN_N60_B_Mean,color=N60Color,linestyle='--')
plot(ax3,linspace(-TimeStartW,size(SW_HI_UN_N70_B,3)/Param.Fs,size(SW_HI_UN_N70_B,3)),SW_HI_UN_N70_B_Mean,color=N70Color,linestyle='--')
plot(ax3,linspace(-TimeStartW,size(SW_HI_AA_N0_B,3)/Param.Fs,size(SW_HI_AA_N0_B,3)),SW_HI_AA_N0_B_Mean,color=QuietColor,linestyle='-.')
plot(ax3,linspace(-TimeStartW,size(SW_HI_AA_N60_B,3)/Param.Fs,size(SW_HI_AA_N60_B,3)),SW_HI_AA_N60_B_Mean,color=N60Color,linestyle='-.')
plot(ax3,linspace(-TimeStartW,size(SW_HI_AA_N70_B,3)/Param.Fs,size(SW_HI_AA_N70_B,3)),SW_HI_AA_N70_B_Mean,color=N70Color,linestyle='-.')
plot(ax3,linspace(-TimeStartW,size(SW_HI_AB_N0_B,3)/Param.Fs,size(SW_HI_AB_N0_B,3)),SW_HI_AB_N0_B_Mean,color=QuietColor,linestyle='-.')
plot(ax3,linspace(-TimeStartW,size(SW_HI_AB_N60_B,3)/Param.Fs,size(SW_HI_AB_N60_B,3)),SW_HI_AB_N60_B_Mean,color=N60Color,linestyle='-.')
plot(ax3,linspace(-TimeStartW,size(SW_HI_AB_N70_B,3)/Param.Fs,size(SW_HI_AB_N70_B,3)),SW_HI_AB_N70_B_Mean,color=N70Color,linestyle='-.')
plot(ax4,nan,color=SColor,linewidth=2) % plot nans to show color in legend
plot(ax4,linspace(-TimeStartW,size(GLW_B,3)/Param.Fs,size(GLW_B,3)),GLW_B_Mean,color=LColor,linewidth=2)
plot(ax4,linspace(-TimeStartW,size(LW_N0_B,3)/Param.Fs,size(LW_N0_B,3)),LW_N0_B_Mean,color=QuietColor)
plot(ax4,linspace(-TimeStartW,size(LW_N60_B,3)/Param.Fs,size(LW_N60_B,3)),LW_N60_B_Mean,color=N60Color)
plot(ax4,linspace(-TimeStartW,size(LW_N70_B,3)/Param.Fs,size(LW_N70_B,3)),LW_N70_B_Mean,color=N70Color)
plot(ax4,linspace(-TimeStartW,size(LW_NH_N0_B,3)/Param.Fs,size(LW_NH_N0_B,3)),LW_NH_N0_B_Mean,color=QuietColor,linestyle=':')
plot(ax4,linspace(-TimeStartW,size(LW_NH_N60_B,3)/Param.Fs,size(LW_NH_N60_B,3)),LW_NH_N60_B_Mean,color=N60Color,linestyle=':')
plot(ax4,linspace(-TimeStartW,size(LW_NH_N70_B,3)/Param.Fs,size(LW_NH_N70_B,3)),LW_NH_N70_B_Mean,color=N70Color,linestyle=':')
plot(ax4,linspace(-TimeStartW,size(LW_HI_UN_N0_B,3)/Param.Fs,size(LW_HI_UN_N0_B,3)),LW_HI_UN_N0_B_Mean,color=QuietColor,linestyle='--')
plot(ax4,linspace(-TimeStartW,size(LW_HI_UN_N60_B,3)/Param.Fs,size(LW_HI_UN_N60_B,3)),LW_HI_UN_N60_B_Mean,color=N60Color,linestyle='--')
plot(ax4,linspace(-TimeStartW,size(LW_HI_UN_N70_B,3)/Param.Fs,size(LW_HI_UN_N70_B,3)),LW_HI_UN_N70_B_Mean,color=N70Color,linestyle='--')
plot(ax4,linspace(-TimeStartW,size(LW_HI_AA_N0_B,3)/Param.Fs,size(LW_HI_AA_N0_B,3)),LW_HI_AA_N0_B_Mean,color=QuietColor,linestyle='-.')
plot(ax4,linspace(-TimeStartW,size(LW_HI_AA_N60_B,3)/Param.Fs,size(LW_HI_AA_N60_B,3)),LW_HI_AA_N60_B_Mean,color=N60Color,linestyle='-.')
plot(ax4,linspace(-TimeStartW,size(LW_HI_AA_N70_B,3)/Param.Fs,size(LW_HI_AA_N70_B,3)),LW_HI_AA_N70_B_Mean,color=N70Color,linestyle='-.')
plot(ax4,linspace(-TimeStartW,size(LW_HI_AB_N0_B,3)/Param.Fs,size(LW_HI_AB_N0_B,3)),LW_HI_AB_N0_B_Mean,color=QuietColor,linestyle='-.')
plot(ax4,linspace(-TimeStartW,size(LW_HI_AB_N60_B,3)/Param.Fs,size(LW_HI_AB_N60_B,3)),LW_HI_AB_N60_B_Mean,color=N60Color,linestyle='-.')
plot(ax4,linspace(-TimeStartW,size(LW_HI_AB_N70_B,3)/Param.Fs,size(LW_HI_AB_N70_B,3)),LW_HI_AB_N70_B_Mean,color=N70Color,linestyle='-.')

% Mean and SEM: Set NaNs to 0 in order to use "fill" correctly
GSW_Mean(isnan(GSW_Mean))=0;GSW_SEM(isnan(GSW_SEM))=0;
GLW_Mean(isnan(GLW_Mean))=0;GLW_SEM(isnan(GLW_SEM))=0;
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

fill(ax1,[linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3)), flipud(linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3))')'],[(GSW_Mean+GSW_SEM), flipud((GSW_Mean-GSW_SEM)')'],SColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_N0,3)/Param.Fs,size(SW_N0,3)), flipud(linspace(-TimeStartW,size(SW_N0,3)/Param.Fs,size(SW_N0,3))')'],[(SW_N0_Mean+SW_N0_SEM), flipud((SW_N0_Mean-SW_N0_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_N60,3)/Param.Fs,size(SW_N60,3)), flipud(linspace(-TimeStartW,size(SW_N60,3)/Param.Fs,size(SW_N60,3))')'],[(SW_N60_Mean+SW_N60_SEM), flipud((SW_N60_Mean-SW_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_N70,3)/Param.Fs,size(SW_N70,3)), flipud(linspace(-TimeStartW,size(SW_N70,3)/Param.Fs,size(SW_N70,3))')'],[(SW_N70_Mean+SW_N70_SEM), flipud((SW_N70_Mean-SW_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_NH_N0,3)/Param.Fs,size(SW_NH_N0,3)), flipud(linspace(-TimeStartW,size(SW_NH_N0,3)/Param.Fs,size(SW_NH_N0,3))')'],[(SW_NH_N0_Mean+SW_NH_N0_SEM), flipud((SW_NH_N0_Mean-SW_NH_N0_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_NH_N60,3)/Param.Fs,size(SW_NH_N60,3)), flipud(linspace(-TimeStartW,size(SW_NH_N60,3)/Param.Fs,size(SW_NH_N60,3))')'],[(SW_NH_N60_Mean+SW_NH_N60_SEM), flipud((SW_NH_N60_Mean-SW_NH_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_NH_N70,3)/Param.Fs,size(SW_NH_N70,3)), flipud(linspace(-TimeStartW,size(SW_NH_N70,3)/Param.Fs,size(SW_NH_N70,3))')'],[(SW_NH_N70_Mean+SW_NH_N70_SEM), flipud((SW_NH_N70_Mean-SW_NH_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_HI_UN_N0,3)/Param.Fs,size(SW_HI_UN_N0,3)), flipud(linspace(-TimeStartW,size(SW_HI_UN_N0,3)/Param.Fs,size(SW_HI_UN_N0,3))')'],[(SW_HI_UN_N0_Mean+SW_HI_UN_N0_SEM), flipud((SW_HI_UN_N0_Mean-SW_HI_UN_N0_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_HI_UN_N60,3)/Param.Fs,size(SW_HI_UN_N60,3)), flipud(linspace(-TimeStartW,size(SW_HI_UN_N60,3)/Param.Fs,size(SW_HI_UN_N60,3))')'],[(SW_HI_UN_N60_Mean+SW_HI_UN_N60_SEM), flipud((SW_HI_UN_N60_Mean-SW_HI_UN_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_HI_UN_N70,3)/Param.Fs,size(SW_HI_UN_N70,3)), flipud(linspace(-TimeStartW,size(SW_HI_UN_N70,3)/Param.Fs,size(SW_HI_UN_N70,3))')'],[(SW_HI_UN_N70_Mean+SW_HI_UN_N70_SEM), flipud((SW_HI_UN_N70_Mean-SW_HI_UN_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_HI_AA_N0,3)/Param.Fs,size(SW_HI_AA_N0,3)), flipud(linspace(-TimeStartW,size(SW_HI_AA_N0,3)/Param.Fs,size(SW_HI_AA_N0,3))')'],[(SW_HI_AA_N0_Mean+SW_HI_AA_N0_SEM), flipud((SW_HI_AA_N0_Mean-SW_HI_AA_N0_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_HI_AA_N60,3)/Param.Fs,size(SW_HI_AA_N60,3)), flipud(linspace(-TimeStartW,size(SW_HI_AA_N60,3)/Param.Fs,size(SW_HI_AA_N60,3))')'],[(SW_HI_AA_N60_Mean+SW_HI_AA_N60_SEM), flipud((SW_HI_AA_N60_Mean-SW_HI_AA_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_HI_AA_N70,3)/Param.Fs,size(SW_HI_AA_N70,3)), flipud(linspace(-TimeStartW,size(SW_HI_AA_N70,3)/Param.Fs,size(SW_HI_AA_N70,3))')'],[(SW_HI_AA_N70_Mean+SW_HI_AA_N70_SEM), flipud((SW_HI_AA_N70_Mean-SW_HI_AA_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_HI_AB_N0,3)/Param.Fs,size(SW_HI_AB_N0,3)), flipud(linspace(-TimeStartW,size(SW_HI_AB_N0,3)/Param.Fs,size(SW_HI_AB_N0,3))')'],[(SW_HI_AB_N0_Mean+SW_HI_AB_N0_SEM), flipud((SW_HI_AB_N0_Mean-SW_HI_AB_N0_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_HI_AB_N60,3)/Param.Fs,size(SW_HI_AB_N60,3)), flipud(linspace(-TimeStartW,size(SW_HI_AB_N60,3)/Param.Fs,size(SW_HI_AB_N60,3))')'],[(SW_HI_AB_N60_Mean+SW_HI_AB_N60_SEM), flipud((SW_HI_AB_N60_Mean-SW_HI_AB_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_HI_AB_N70,3)/Param.Fs,size(SW_HI_AB_N70,3)), flipud(linspace(-TimeStartW,size(SW_HI_AB_N70,3)/Param.Fs,size(SW_HI_AB_N70,3))')'],[(SW_HI_AB_N70_Mean+SW_HI_AB_N70_SEM), flipud((SW_HI_AB_N70_Mean-SW_HI_AB_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(GLW,3)/Param.Fs,size(GLW,3)), flipud(linspace(-TimeStartW,size(GLW,3)/Param.Fs,size(GLW,3))')'],[(GLW_Mean+GLW_SEM), flipud((GLW_Mean-GLW_SEM)')'],SColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_N0,3)/Param.Fs,size(LW_N0,3)), flipud(linspace(-TimeStartW,size(LW_N0,3)/Param.Fs,size(LW_N0,3))')'],[(LW_N0_Mean+LW_N0_SEM), flipud((LW_N0_Mean-LW_N0_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_N60,3)/Param.Fs,size(LW_N60,3)), flipud(linspace(-TimeStartW,size(LW_N60,3)/Param.Fs,size(LW_N60,3))')'],[(LW_N60_Mean+LW_N60_SEM), flipud((LW_N60_Mean-LW_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_N70,3)/Param.Fs,size(LW_N70,3)), flipud(linspace(-TimeStartW,size(LW_N70,3)/Param.Fs,size(LW_N70,3))')'],[(LW_N70_Mean+LW_N70_SEM), flipud((LW_N70_Mean-LW_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_NH_N0,3)/Param.Fs,size(LW_NH_N0,3)), flipud(linspace(-TimeStartW,size(LW_NH_N0,3)/Param.Fs,size(LW_NH_N0,3))')'],[(LW_NH_N0_Mean+LW_NH_N0_SEM), flipud((LW_NH_N0_Mean-LW_NH_N0_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_NH_N60,3)/Param.Fs,size(LW_NH_N60,3)), flipud(linspace(-TimeStartW,size(LW_NH_N60,3)/Param.Fs,size(LW_NH_N60,3))')'],[(LW_NH_N60_Mean+LW_NH_N60_SEM), flipud((LW_NH_N60_Mean-LW_NH_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_NH_N70,3)/Param.Fs,size(LW_NH_N70,3)), flipud(linspace(-TimeStartW,size(LW_NH_N70,3)/Param.Fs,size(LW_NH_N70,3))')'],[(LW_NH_N70_Mean+LW_NH_N70_SEM), flipud((LW_NH_N70_Mean-LW_NH_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_HI_UN_N0,3)/Param.Fs,size(LW_HI_UN_N0,3)), flipud(linspace(-TimeStartW,size(LW_HI_UN_N0,3)/Param.Fs,size(LW_HI_UN_N0,3))')'],[(LW_HI_UN_N0_Mean+LW_HI_UN_N0_SEM), flipud((LW_HI_UN_N0_Mean-LW_HI_UN_N0_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_HI_UN_N60,3)/Param.Fs,size(LW_HI_UN_N60,3)), flipud(linspace(-TimeStartW,size(LW_HI_UN_N60,3)/Param.Fs,size(LW_HI_UN_N60,3))')'],[(LW_HI_UN_N60_Mean+LW_HI_UN_N60_SEM), flipud((LW_HI_UN_N60_Mean-LW_HI_UN_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_HI_UN_N70,3)/Param.Fs,size(LW_HI_UN_N70,3)), flipud(linspace(-TimeStartW,size(LW_HI_UN_N70,3)/Param.Fs,size(LW_HI_UN_N70,3))')'],[(LW_HI_UN_N70_Mean+LW_HI_UN_N70_SEM), flipud((LW_HI_UN_N70_Mean-LW_HI_UN_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_HI_AA_N0,3)/Param.Fs,size(LW_HI_AA_N0,3)), flipud(linspace(-TimeStartW,size(LW_HI_AA_N0,3)/Param.Fs,size(LW_HI_AA_N0,3))')'],[(LW_HI_AA_N0_Mean+LW_HI_AA_N0_SEM), flipud((LW_HI_AA_N0_Mean-LW_HI_AA_N0_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_HI_AA_N60,3)/Param.Fs,size(LW_HI_AA_N60,3)), flipud(linspace(-TimeStartW,size(LW_HI_AA_N60,3)/Param.Fs,size(LW_HI_AA_N60,3))')'],[(LW_HI_AA_N60_Mean+LW_HI_AA_N60_SEM), flipud((LW_HI_AA_N60_Mean-LW_HI_AA_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_HI_AA_N70,3)/Param.Fs,size(LW_HI_AA_N70,3)), flipud(linspace(-TimeStartW,size(LW_HI_AA_N70,3)/Param.Fs,size(LW_HI_AA_N70,3))')'],[(LW_HI_AA_N70_Mean+LW_HI_AA_N70_SEM), flipud((LW_HI_AA_N70_Mean-LW_HI_AA_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_HI_AB_N0,3)/Param.Fs,size(LW_HI_AB_N0,3)), flipud(linspace(-TimeStartW,size(LW_HI_AB_N0,3)/Param.Fs,size(LW_HI_AB_N0,3))')'],[(LW_HI_AB_N0_Mean+LW_HI_AB_N0_SEM), flipud((LW_HI_AB_N0_Mean-LW_HI_AB_N0_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_HI_AB_N60,3)/Param.Fs,size(LW_HI_AB_N60,3)), flipud(linspace(-TimeStartW,size(LW_HI_AB_N60,3)/Param.Fs,size(LW_HI_AB_N60,3))')'],[(LW_HI_AB_N60_Mean+LW_HI_AB_N60_SEM), flipud((LW_HI_AB_N60_Mean-LW_HI_AB_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_HI_AB_N70,3)/Param.Fs,size(LW_HI_AB_N70,3)), flipud(linspace(-TimeStartW,size(LW_HI_AB_N70,3)/Param.Fs,size(LW_HI_AB_N70,3))')'],[(LW_HI_AB_N70_Mean+LW_HI_AB_N70_SEM), flipud((LW_HI_AB_N70_Mean-LW_HI_AB_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')

fill(ax3,[linspace(-TimeStartW,size(GSW_B,3)/Param.Fs,size(GSW_B,3)), flipud(linspace(-TimeStartW,size(GSW_B,3)/Param.Fs,size(GSW_B,3))')'],[(GSW_B_Mean+GSW_B_SEM), flipud((GSW_B_Mean-GSW_B_SEM)')'],SColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SW_N0_B,3)/Param.Fs,size(SW_N0_B,3)), flipud(linspace(-TimeStartW,size(SW_N0_B,3)/Param.Fs,size(SW_N0_B,3))')'],[(SW_N0_B_Mean+SW_N0_B_SEM), flipud((SW_N0_B_Mean-SW_N0_B_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SW_N60_B,3)/Param.Fs,size(SW_N60_B,3)), flipud(linspace(-TimeStartW,size(SW_N60_B,3)/Param.Fs,size(SW_N60_B,3))')'],[(SW_N60_B_Mean+SW_N60_B_SEM), flipud((SW_N60_B_Mean-SW_N60_B_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SW_N70_B,3)/Param.Fs,size(SW_N70_B,3)), flipud(linspace(-TimeStartW,size(SW_N70_B,3)/Param.Fs,size(SW_N70_B,3))')'],[(SW_N70_B_Mean+SW_N70_B_SEM), flipud((SW_N70_B_Mean-SW_N70_B_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SW_NH_N0_B,3)/Param.Fs,size(SW_NH_N0_B,3)), flipud(linspace(-TimeStartW,size(SW_NH_N0_B,3)/Param.Fs,size(SW_NH_N0_B,3))')'],[(SW_NH_N0_B_Mean+SW_NH_N0_B_SEM), flipud((SW_NH_N0_B_Mean-SW_NH_N0_B_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SW_NH_N60_B,3)/Param.Fs,size(SW_NH_N60_B,3)), flipud(linspace(-TimeStartW,size(SW_NH_N60_B,3)/Param.Fs,size(SW_NH_N60_B,3))')'],[(SW_NH_N60_B_Mean+SW_NH_N60_B_SEM), flipud((SW_NH_N60_B_Mean-SW_NH_N60_B_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SW_NH_N70_B,3)/Param.Fs,size(SW_NH_N70_B,3)), flipud(linspace(-TimeStartW,size(SW_NH_N70_B,3)/Param.Fs,size(SW_NH_N70_B,3))')'],[(SW_NH_N70_B_Mean+SW_NH_N70_B_SEM), flipud((SW_NH_N70_B_Mean-SW_NH_N70_B_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SW_HI_UN_N0_B,3)/Param.Fs,size(SW_HI_UN_N0_B,3)), flipud(linspace(-TimeStartW,size(SW_HI_UN_N0_B,3)/Param.Fs,size(SW_HI_UN_N0_B,3))')'],[(SW_HI_UN_N0_B_Mean+SW_HI_UN_N0_B_SEM), flipud((SW_HI_UN_N0_B_Mean-SW_HI_UN_N0_B_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SW_HI_UN_N60_B,3)/Param.Fs,size(SW_HI_UN_N60_B,3)), flipud(linspace(-TimeStartW,size(SW_HI_UN_N60_B,3)/Param.Fs,size(SW_HI_UN_N60_B,3))')'],[(SW_HI_UN_N60_B_Mean+SW_HI_UN_N60_B_SEM), flipud((SW_HI_UN_N60_B_Mean-SW_HI_UN_N60_B_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SW_HI_UN_N70_B,3)/Param.Fs,size(SW_HI_UN_N70_B,3)), flipud(linspace(-TimeStartW,size(SW_HI_UN_N70_B,3)/Param.Fs,size(SW_HI_UN_N70_B,3))')'],[(SW_HI_UN_N70_B_Mean+SW_HI_UN_N70_B_SEM), flipud((SW_HI_UN_N70_B_Mean-SW_HI_UN_N70_B_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SW_HI_AA_N0_B,3)/Param.Fs,size(SW_HI_AA_N0_B,3)), flipud(linspace(-TimeStartW,size(SW_HI_AA_N0_B,3)/Param.Fs,size(SW_HI_AA_N0_B,3))')'],[(SW_HI_AA_N0_B_Mean+SW_HI_AA_N0_B_SEM), flipud((SW_HI_AA_N0_B_Mean-SW_HI_AA_N0_B_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SW_HI_AA_N60_B,3)/Param.Fs,size(SW_HI_AA_N60_B,3)), flipud(linspace(-TimeStartW,size(SW_HI_AA_N60_B,3)/Param.Fs,size(SW_HI_AA_N60_B,3))')'],[(SW_HI_AA_N60_B_Mean+SW_HI_AA_N60_B_SEM), flipud((SW_HI_AA_N60_B_Mean-SW_HI_AA_N60_B_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SW_HI_AA_N70_B,3)/Param.Fs,size(SW_HI_AA_N70_B,3)), flipud(linspace(-TimeStartW,size(SW_HI_AA_N70_B,3)/Param.Fs,size(SW_HI_AA_N70_B,3))')'],[(SW_HI_AA_N70_B_Mean+SW_HI_AA_N70_B_SEM), flipud((SW_HI_AA_N70_B_Mean-SW_HI_AA_N70_B_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SW_HI_AB_N0_B,3)/Param.Fs,size(SW_HI_AB_N0_B,3)), flipud(linspace(-TimeStartW,size(SW_HI_AB_N0_B,3)/Param.Fs,size(SW_HI_AB_N0_B,3))')'],[(SW_HI_AB_N0_B_Mean+SW_HI_AB_N0_B_SEM), flipud((SW_HI_AB_N0_B_Mean-SW_HI_AB_N0_B_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SW_HI_AB_N60_B,3)/Param.Fs,size(SW_HI_AB_N60_B,3)), flipud(linspace(-TimeStartW,size(SW_HI_AB_N60_B,3)/Param.Fs,size(SW_HI_AB_N60_B,3))')'],[(SW_HI_AB_N60_B_Mean+SW_HI_AB_N60_B_SEM), flipud((SW_HI_AB_N60_B_Mean-SW_HI_AB_N60_B_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SW_HI_AB_N70_B,3)/Param.Fs,size(SW_HI_AB_N70_B,3)), flipud(linspace(-TimeStartW,size(SW_HI_AB_N70_B,3)/Param.Fs,size(SW_HI_AB_N70_B,3))')'],[(SW_HI_AB_N70_B_Mean+SW_HI_AB_N70_B_SEM), flipud((SW_HI_AB_N70_B_Mean-SW_HI_AB_N70_B_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(GLW_B,3)/Param.Fs,size(GLW_B,3)), flipud(linspace(-TimeStartW,size(GLW_B,3)/Param.Fs,size(GLW_B,3))')'],[(GLW_B_Mean+GLW_B_SEM), flipud((GLW_B_Mean-GLW_B_SEM)')'],SColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LW_N0_B,3)/Param.Fs,size(LW_N0_B,3)), flipud(linspace(-TimeStartW,size(LW_N0_B,3)/Param.Fs,size(LW_N0_B,3))')'],[(LW_N0_B_Mean+LW_N0_B_SEM), flipud((LW_N0_B_Mean-LW_N0_B_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LW_N60_B,3)/Param.Fs,size(LW_N60_B,3)), flipud(linspace(-TimeStartW,size(LW_N60_B,3)/Param.Fs,size(LW_N60_B,3))')'],[(LW_N60_B_Mean+LW_N60_B_SEM), flipud((LW_N60_B_Mean-LW_N60_B_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LW_N70_B,3)/Param.Fs,size(LW_N70_B,3)), flipud(linspace(-TimeStartW,size(LW_N70_B,3)/Param.Fs,size(LW_N70_B,3))')'],[(LW_N70_B_Mean+LW_N70_B_SEM), flipud((LW_N70_B_Mean-LW_N70_B_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LW_NH_N0_B,3)/Param.Fs,size(LW_NH_N0_B,3)), flipud(linspace(-TimeStartW,size(LW_NH_N0_B,3)/Param.Fs,size(LW_NH_N0_B,3))')'],[(LW_NH_N0_B_Mean+LW_NH_N0_B_SEM), flipud((LW_NH_N0_B_Mean-LW_NH_N0_B_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LW_NH_N60_B,3)/Param.Fs,size(LW_NH_N60_B,3)), flipud(linspace(-TimeStartW,size(LW_NH_N60_B,3)/Param.Fs,size(LW_NH_N60_B,3))')'],[(LW_NH_N60_B_Mean+LW_NH_N60_B_SEM), flipud((LW_NH_N60_B_Mean-LW_NH_N60_B_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LW_NH_N70_B,3)/Param.Fs,size(LW_NH_N70_B,3)), flipud(linspace(-TimeStartW,size(LW_NH_N70_B,3)/Param.Fs,size(LW_NH_N70_B,3))')'],[(LW_NH_N70_B_Mean+LW_NH_N70_B_SEM), flipud((LW_NH_N70_B_Mean-LW_NH_N70_B_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LW_HI_UN_N0_B,3)/Param.Fs,size(LW_HI_UN_N0_B,3)), flipud(linspace(-TimeStartW,size(LW_HI_UN_N0_B,3)/Param.Fs,size(LW_HI_UN_N0_B,3))')'],[(LW_HI_UN_N0_B_Mean+LW_HI_UN_N0_B_SEM), flipud((LW_HI_UN_N0_B_Mean-LW_HI_UN_N0_B_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LW_HI_UN_N60_B,3)/Param.Fs,size(LW_HI_UN_N60_B,3)), flipud(linspace(-TimeStartW,size(LW_HI_UN_N60_B,3)/Param.Fs,size(LW_HI_UN_N60_B,3))')'],[(LW_HI_UN_N60_B_Mean+LW_HI_UN_N60_B_SEM), flipud((LW_HI_UN_N60_B_Mean-LW_HI_UN_N60_B_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LW_HI_UN_N70_B,3)/Param.Fs,size(LW_HI_UN_N70_B,3)), flipud(linspace(-TimeStartW,size(LW_HI_UN_N70_B,3)/Param.Fs,size(LW_HI_UN_N70_B,3))')'],[(LW_HI_UN_N70_B_Mean+LW_HI_UN_N70_B_SEM), flipud((LW_HI_UN_N70_B_Mean-LW_HI_UN_N70_B_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LW_HI_AA_N0_B,3)/Param.Fs,size(LW_HI_AA_N0_B,3)), flipud(linspace(-TimeStartW,size(LW_HI_AA_N0_B,3)/Param.Fs,size(LW_HI_AA_N0_B,3))')'],[(LW_HI_AA_N0_B_Mean+LW_HI_AA_N0_B_SEM), flipud((LW_HI_AA_N0_B_Mean-LW_HI_AA_N0_B_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LW_HI_AA_N60_B,3)/Param.Fs,size(LW_HI_AA_N60_B,3)), flipud(linspace(-TimeStartW,size(LW_HI_AA_N60_B,3)/Param.Fs,size(LW_HI_AA_N60_B,3))')'],[(LW_HI_AA_N60_B_Mean+LW_HI_AA_N60_B_SEM), flipud((LW_HI_AA_N60_B_Mean-LW_HI_AA_N60_B_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LW_HI_AA_N70_B,3)/Param.Fs,size(LW_HI_AA_N70_B,3)), flipud(linspace(-TimeStartW,size(LW_HI_AA_N70_B,3)/Param.Fs,size(LW_HI_AA_N70_B,3))')'],[(LW_HI_AA_N70_B_Mean+LW_HI_AA_N70_B_SEM), flipud((LW_HI_AA_N70_B_Mean-LW_HI_AA_N70_B_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LW_HI_AB_N0_B,3)/Param.Fs,size(LW_HI_AB_N0_B,3)), flipud(linspace(-TimeStartW,size(LW_HI_AB_N0_B,3)/Param.Fs,size(LW_HI_AB_N0_B,3))')'],[(LW_HI_AB_N0_B_Mean+LW_HI_AB_N0_B_SEM), flipud((LW_HI_AB_N0_B_Mean-LW_HI_AB_N0_B_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LW_HI_AB_N60_B,3)/Param.Fs,size(LW_HI_AB_N60_B,3)), flipud(linspace(-TimeStartW,size(LW_HI_AB_N60_B,3)/Param.Fs,size(LW_HI_AB_N60_B,3))')'],[(LW_HI_AB_N60_B_Mean+LW_HI_AB_N60_B_SEM), flipud((LW_HI_AB_N60_B_Mean-LW_HI_AB_N60_B_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LW_HI_AB_N70_B,3)/Param.Fs,size(LW_HI_AB_N70_B,3)), flipud(linspace(-TimeStartW,size(LW_HI_AB_N70_B,3)/Param.Fs,size(LW_HI_AB_N70_B,3))')'],[(LW_HI_AB_N70_B_Mean+LW_HI_AB_N70_B_SEM), flipud((LW_HI_AB_N70_B_Mean-LW_HI_AB_N70_B_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')

xlabel([ax1 ax2 ax3 ax4],'Time [s]')
ylabel([ax1 ax2],'Pupil diameter [mm]')
ylabel([ax3 ax4],'Pupil baseline difference [mm]')
xlim([ax1 ax2 ax3 ax4],[-TimeStartW 3])
lgd2=legend(ax2,'Speaking','Listening','Speaking','Listening','N0','N60','N70','NH_N0','NH_N60','NH_N70','HI_UN_N0','HI_UN_N60','HI_UN_N70','HI_AA_N0','HI_AA_N60','HI_AA_N70','HI_AB_N0','HI_AB_N60','HI_AB_N70','Location','southeastoutside');
lgd2.Title.String = 'Types of windows:';
lgd4=legend(ax4,'Speaking','Listening','Speaking','Listening','N0','N60','N70','NH_N0','NH_N60','NH_N70','HI_UN_N0','HI_UN_N60','HI_UN_N70','HI_AA_N0','HI_AA_N60','HI_AA_N70','HI_AB_N0','HI_AB_N60','HI_AB_N70','Location','southeastoutside');
lgd4.Title.String = 'Types of windows:';
title(ax1,'Global Speaking-evoked')
title(ax2,'Global Listening-evoked')
title(ax3,'Global adaptive baselined Speaking-evoked')
title(ax4,'Global adaptive baselined Listening-evoked')