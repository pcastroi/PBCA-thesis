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

% Variables
NCols=60*Param.Fs; % Duration (samples) of each window
NRows=100; % Number of windows per trial
NLayers=numel(FileNames)*numel(subDirs); % Number of trials

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

GSW_B = zeros(NLayers,NRows,NCols); GLW_B = GSW; % Global Speaking/Listening Windows Adaptive-Baseline corrected
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
SDur_N60 = zeros(NLayers,NRows); % SHL Speaking Windows Duration
LDur_N60 = zeros(NLayers,NRows); % SHL Listening Windows Duration
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

            % Time-locked indexes (based on Start or End of events)
            WSpeakIdx=[Speak(:,2)-TimeStartW*Param.Fs,Speak(:,2),Speak(:,3)+TimeEndW*Param.Fs];
            WListenIdx=[Listen(:,2)-TimeStartW*Param.Fs,Listen(:,2),Listen(:,2)+TimeEndW*Param.Fs];
        %         EWSpeakIdx=[Speak(:,3)-TimeStartW*Param.Fs,Speak(:,3),Speak(:,3)+TimeEndW*Param.Fs];
        %         EWListenIdx=[Listen(:,3)-TimeStartW*Param.Fs,Listen(:,3),Listen(:,3)+TimeEndW*Param.Fs];

        %         t_Diam = linspace(0,length(Diameter)./Param.Fs,length(Diameter));

            % Storing Speaking/Listening by conditions
            if contains(cell2mat(FileNames_II(i)),'UN')
                if contains(cell2mat(FileNames_II(i)),'N0')
                    for j=1:size(Speak,1)
                        SW_N0(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        SW_B_N0(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_B_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        SDur_N0(i,j) = Speak(j,1);
                    end
                    for j=1:size(Listen,1)
                        LW_N0(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                        LW_B_N0(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_B_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                        LDur_N0(i,j) = Listen(j,1);
                    end
                elseif contains(cell2mat(FileNames_II(i)),'N60')
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
                elseif contains(cell2mat(FileNames_II(i)),'N70')
                    for j=1:size(Speak,1)
                        SW_N70(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(SW_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        SW_B_N70(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SW_B_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        SDur_N70(i,j) = Speak(j,1);
                    end
                    for j=1:size(Listen,1)
                        LW_N70(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(LW_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                        LW_B_N70(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LW_B_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                        LDur_N70(i,j) = Listen(j,1);
                    end
                end
            elseif contains(cell2mat(FileNames_II(i)),'AA')
                
            elseif contains(cell2mat(FileNames_II(i)),'AB')
                    
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
            if contains(cell2mat(FileNames_II(i)),'P2')
                TPsOrder(2*q,i-NCond) = x;
            elseif contains(cell2mat(FileNames_II(i)),'P1')
                TPsOrder(2*q-1,i) = x;
            end

            % Increase index of num of files used
            x=x+1;
            
            % Plots
            figure
            plot(linspace(0,size(Diameter,1)/Param.Fs,size(Diameter,1)),Diameter,color='black');
            grid on
            yl=ylim();
            ylim(yl);
            % Plot rectangles (Utterance and listening time windows)
%             arrayfun(@(i)rectangle('Position', [startStopS(i,1),yl(1),widthS(i),range(yl)], ...
%             'EdgeColor', 'none', 'FaceColor', [0 1 0 .2]), 1:size(startStopS,1))
%             arrayfun(@(i)rectangle('Position', [startStopL(i,1),yl(1),widthL(i),range(yl)], ...
%             'EdgeColor', 'none', 'FaceColor', [1 0 1 .2]), 1:size(startStopL,1))
%             eline1=line(NaN,NaN,'LineWidth',2,'Color',[0 1 0 .2]);
%             eline2=line(NaN,NaN,'LineWidth',2,'Color',[1 0 1 .2]);
%             legend(['Baselined diameter (', eyeChosen,' Eye)'],'Speaking windows','Listening windows')
            sgtitle(strrep(strrep([PairFiles_II(1).folder,'\',PairFiles_II(i).name],'_','-'),'\','\\'))
            xlabel('Time [s]')
            ylabel('Pupil diameter [mm]');

        end
    end
end