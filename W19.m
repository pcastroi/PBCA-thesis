%% PBCA-Thesis - Week 19 - AMEND II data - CHECK FILE BY FILE
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
Param.RemoveBeforeAndAfter = [0 0]*1e-3; % Samples within the time range before and after NaNs will set NaNs as well.
Param.MinLengthNaNRepair = 0; % Drop values (i.e., change to NaN) before and after NaNs only for contiguous NaNs of at least __ samples. 
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
ChosenFolder = {'\NH\'};
% for q=1:numel(subDirs_II)
q=4;
i=9;
    PairIn_II = q;
%     for ChosenFolder = {'\NH\','\HI\'}
        PairFolder_II=[pwd,'\data\AMEND_II\',cell2mat(subDirs_II(q)),cell2mat(ChosenFolder)]; % Folder naming changed
        PairFiles_II=dir(PairFolder_II); % Folder naming changed
        PairUtt_II=LoadUtt_II.Utterances(PairIn_II,:);
        
        if isempty(PairUtt_II{1})
            disp(['Warning: No Utterance found for folder "',cell2mat(ChosenFolder),cell2mat(subDirs_II(q)),'"'])
%             continue
        end
%         PairDelay_II=LoadDelays_II.TobAudDelay(PairIn_II,:);
    
%         for i=1:numel(FileNames_II)
            try
%                 alldata = load([PairFiles_II(1).folder, '\', cell2mat(FileNames_II(i))]);
                alldata = load(['C:\git\PBCA-thesis\data\AMEND_II\Pair05\NH', '\', 'ABHI_N70.mat']);             
            catch ME
                disp(['Warning: File ', PairFiles_II(1).folder, '\', cell2mat(FileNames_II(i)), ' not found (no Gaze data).']);
%                 continue
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
%                 continue
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
            
            LDiamPre=LDiamRaw;LDiamPre(~LvalOut)=NaN;
            RDiamPre=RDiamRaw;RDiamPre(~RvalOut)=NaN;

            % Processing - Interpolating NaNs
            [LDiam,LMetadata] = preprocpupil(LDiamPre,Param);
            [RDiam,RMetadata] = preprocpupil(RDiamPre,Param);

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
            [Min, idx_decision] = min([sum(LMetadata.Isnan) sum(RMetadata.Isnan)]);
            if idx_decision == 1
                Diameter = LDiam;
                DiameterRaw = LDiamRaw';
                DiameterPre = LDiamPre';
                eyeChosen = 'Left';
                DiamNaN = sum(LMetadata.Isnan);
            elseif idx_decision == 2
                Diameter = RDiam;
                DiameterRaw = RDiamRaw';
                DiameterPre = RDiamPre';
                eyeChosen = 'Right';
                DiamNaN = sum(RMetadata.Isnan);
            end

            if DiamNaN/length(Diameter) >= RejectRatio
                disp(['Warning: File ',PairFiles_II(1).folder, '\', cell2mat(FileNames_II(i)), ' was rejected because it contains too many NaNs (',sprintf('%0.2f',100*DiamNaN/length(Diameter)),'%).'])
%                 continue
            end
            
            % Add nan padding - Diameter
            NaNPadS = 2*TimeStartW*Param.Fs;
            Diameter = [nan*ones(NaNPadS,1);Diameter];
            DiameterRaw = [nan*ones(NaNPadS,1);DiameterRaw];
            DiameterPre = [nan*ones(NaNPadS,1);DiameterPre];
            
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
%                 continue
            end
            
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

            % Time-locked indexes (based on Start or End of events)
            WSpeakIdx=[Speak(:,2)-TimeStartW*Param.Fs,Speak(:,2),Speak(:,3)+TimeEndW*Param.Fs];
            WListenIdx=[Listen(:,2)-TimeStartW*Param.Fs,Listen(:,2),Listen(:,2)+TimeEndW*Param.Fs];
            
        %         EWSpeakIdx=[Speak(:,3)-TimeStartW*Param.Fs,Speak(:,3),Speak(:,3)+TimeEndW*Param.Fs];
        %         EWListenIdx=[Listen(:,3)-TimeStartW*Param.Fs,Listen(:,3),Listen(:,3)+TimeEndW*Param.Fs];

            s_Diam = linspace(1,length(Diameter),length(Diameter));

            % Increase index of num of files used
            x=x+1;
            
            figure('Position',[10 10 1500 400],'Name',strrep(strrep([PairFiles_II(1).folder, '\', cell2mat(FileNames_II(i))],'_','-'),'\','\\'))
            hold on
            grid on
            plot(s_Diam,Diameter,color='black')
            scatter(s_Diam,DiameterPre,'k*')
            scatter(s_Diam,DiameterRaw,'r.')
            startStopS = s_Diam(Speak(:,2:3));widthS = startStopS(:,2)-startStopS(:,1);
            startStopL = s_Diam(Listen(:,2:3));widthL = startStopL(:,2)-startStopL(:,1);
            ylim('tight')
            xlim('tight')
            yl=ylim();
            arrayfun(@(i)rectangle('Position', [startStopS(i,1),yl(1),widthS(i),range(yl)],'EdgeColor', 'none', 'FaceColor', [0 1 0 .2]), 1:size(startStopS,1))
            arrayfun(@(i)rectangle('Position', [startStopL(i,1),yl(1),widthL(i),range(yl)],'EdgeColor', 'none', 'FaceColor', [1 0 1 .2]), 1:size(startStopL,1))
            eline1=line(NaN,NaN,'LineWidth',4,'Color',[0 1 0 .2]);
            eline2=line(NaN,NaN,'LineWidth',4,'Color',[1 0 1 .2]);
            legend(['Diameter (', eyeChosen,' Eye)'],['Diameter Raw (', eyeChosen,' Eye)'],'Speaking windows','Listening windows','Location','southeastoutside')
            title(['Delay = ', sprintf('%0.2f',EyeAudDelay), 's | NaNs = ', sprintf('%0.2f',100*DiamNaN/length(Diameter)), '%'])
            
            ABC_cases=[11080,6605,3025];
            figure;tiledlayout(1,3);ax1 = nexttile;ax2 = nexttile;ax3 = nexttile;
            hold([ax1 ax2 ax3],'on')
            grid([ax1 ax2 ax3],'on')
            w=50;
            DiameterRaw(3015)=NaN;DiameterRaw(3024)=NaN;DiameterRaw(3025)=NaN;DiameterRaw(3026)=NaN;
            t1 = linspace(0,1,length(ABC_cases(1)-w:ABC_cases(1)+w));
            scatter(ax1,t1,DiameterRaw(ABC_cases(1)-w:ABC_cases(1)+w),100,'k.')
            xd1 = find(DiameterRaw(ABC_cases(1)-w:ABC_cases(1)+w,:)>3.0995);
            scatter(ax1,t1(xd1),DiameterRaw(xd1+ABC_cases(1)-w-1),100,'r.')
            
            t2 = linspace(0,1,length(ABC_cases(2)-w:ABC_cases(2)+w));
            plot(ax2,t2,Diameter(ABC_cases(2)-w:ABC_cases(2)+w),color='black',Linewidth=1);
            scatter(ax2,linspace(0,1,length(ABC_cases(2)-w:ABC_cases(2)+w)),DiameterRaw(ABC_cases(2)-w:ABC_cases(2)+w),100,'k.')
            xd2 = find(DiameterRaw(ABC_cases(2)-w:ABC_cases(2)+w,:)>3.8);
            scatter(ax2,t2(xd2),DiameterRaw(xd2+ABC_cases(2)-w-1),100,'r.')
            
            t3 = linspace(0,1,length(ABC_cases(3)-w:ABC_cases(3)+w));
            scatter(ax3,t3,DiameterRaw(ABC_cases(3)-w:ABC_cases(3)+w),100,'k.')
            xd3 = find(DiameterRaw(ABC_cases(3)-w:ABC_cases(3)+w,:)>3.4);
            scatter(ax3,t3(xd3),DiameterRaw(xd3+ABC_cases(3)-w-1),100,'r.')
            
            title(ax1,'(1)')
            title(ax2,'(2)')
            title(ax3,'(3)')
            xlabel([ax1 ax2 ax3],'Time [s]')
            ylabel(ax1,'Pupil diameter [mm]');
            
%         end
%     end
% end