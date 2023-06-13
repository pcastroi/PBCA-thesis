%% PBCA-Thesis - Week 9, 10 & 11 - Fixation duration? Not using diameter anymore
% Pathing
clear all; clc; close all;
BPath = strsplit(pwd,'PBCA-thesis');
addpath('tools')
addpath([BPath{1} '\Pupil-preprocessing-tools\tools']) % For preprocessing

[subDirs] = GetSubDirsFirstLevelOnly('data\AMEND_I');
FileNames={'P1_Quiet_B1.mat','P1_Quiet_B2.mat','P1_SHL_B1.mat','P1_SHL_B2.mat','P1_Noise60_B1.mat','P1_Noise60_B2.mat','P1_Noise70_B1.mat','P1_Noise70_B2.mat','P2_Quiet_B1.mat','P2_Quiet_B2.mat','P2_SHL_B1.mat','P2_SHL_B2.mat','P2_Noise60_B1.mat','P2_Noise60_B2.mat','P2_Noise70_B1.mat','P2_Noise70_B2.mat'};
LoadUtt=load('data\AMEND_I\utterances1110.mat');
LoadDelays=load('data\AMEND_I\delays1110.mat');
LoadTPsOrder=load('data\AMEND_I\TPsOrder_I.mat');

% Parameters for processing
Param.Fs = 50; % Sampling frequency of pupil data
Param.Fs = 50; % Sampling frequency of pupil data
Param.RemoveBeforeAndAfter = [0 0]*1e-3; % Samples within the time range before and after NaNs will set NaNs as well.
Param.MinLengthNaNRepair = 0; % Drop values (i.e., change to NaN) before and after NaNs only for contiguous NaNs of at least __ samples. 
MaxVisualAngle = 1.5; % [degrees], critical visual angle for Fixation definition
LPWinSize = 0.5; % [s]: Window size of hamming-window for low-pass filtering
FilterWidth = round((LPWinSize*Param.Fs)/2); % [samples]: Width of hamming filter used for Fixation duration
LPWindow = hamming(round(LPWinSize*Param.Fs));
LPWindow = LPWindow/sum(LPWindow); % Hamming-window

NCols=60*Param.Fs; % Duration (samples) of each window
NRows=100; % Number of windows per trial
NLayers=numel(FileNames)*numel(subDirs); % Number of trials
NPs = 2; % N of TPs per Pair
NCond = numel(FileNames)/NPs; % N of conditions
NTPs = numel(subDirs)*numel(FileNames)/NCond; % Total N of TPs

TimeStartW = 0.5; % [s], time before Utt/Lis starts
AdapBL = 0.3; % [s], Baseline period
TimeEndW = 0; % [s], time after Utt/Lis starts
TimeStart = 20; % [s], time at which simultaneous recording started
TimeMinWin = 0.5; % [s], Minimum time of a window
TimeInitialMerge = 0.3; % [s], Time threshold for merging windows initially
TimeMerge = 2; % [s], Time threshold for merging windows after rejecting small windows
RejectRatio = 0.4; % Rejection threshold based on the ratio of NaNs in data
RejectDelay = 0.5; % [s], Rejection threshold based on delay between timestamps and n-samples
x=1; % idx to store global values

% Colors
SColor = [53, 155, 67]./255;
LColor = [204, 36, 0]./255;
QuietColor = [204, 152, 0]./255;
SHLColor = [123, 31, 162]./255;
N60Color = [0, 196, 215]./255;
N70Color = [2, 36, 223]./255;

% Variables
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

GSDur = zeros(200,400);
GLDur = zeros(200,400);
SDur_Quiet = zeros(200,400);
LDur_Quiet = zeros(200,400);
SDur_SHL = zeros(200,400);
LDur_SHL = zeros(200,400);
SDur_N60 = zeros(200,400);
LDur_N60 = zeros(200,400);
SDur_N70 = zeros(200,400);
LDur_N70 = zeros(200,400);

for q=1:numel(subDirs)
    PairIn = q;
    PairFiles=dir(['data\AMEND_I\Main',sprintf('%d',PairIn),'\*.mat']);
    PairUtt=LoadUtt.Utterances(PairIn,:);
    PairDelay=LoadDelays.TobAudDelay(PairIn,:);
    
    for i=1:numel(FileNames)
        if contains(cell2mat(FileNames(i)),'P2')
            if LoadTPsOrder.TPsOrder(2*q,i-NCond) ~= x
                disp(['Warning: File ',PairFiles(1).folder, '\', cell2mat(FileNames(i)),' was rejected in AMEND I (P2) analysis.']);
                continue
            end
        elseif contains(cell2mat(FileNames(i)),'P1')
            if LoadTPsOrder.TPsOrder(2*q-1,i) ~= x
                disp(['Warning: File ',PairFiles(1).folder, '\', cell2mat(FileNames(i)),' was rejected in AMEND I (P1) analysis.']);
                continue
            end
        end
        alldata = load([PairFiles(1).folder, '\', cell2mat(FileNames(i))]);
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

        % GAZE 3D
        % In case gaze3d is transposed from origin [x;y;z] -> transpose to [x,y,z]
        if size([alldata_mat(:,cellfun(@(xd) ~any(isnan(xd)),{alldata_mat.gaze3d})).gaze3d],1) == 3
            for k=1:size(alldata_mat,2)
               alldata_mat(k).gaze3d=alldata_mat(k).gaze3d'; 
            end
        end

        % Replace blanks '[]' for 'NaN' in gaze3d
        [alldata_mat(cellfun(@isempty,{alldata_mat.gaze3d})).gaze3d] = deal(NaN);

        % Replace 'NaN' to '[NaN,NaN,NaN]' in gaze3d
        [alldata_mat(:,cellfun(@(xd) any(isnan(xd)),{alldata_mat.gaze3d})).gaze3d] = deal([NaN,NaN,NaN]);

        gaze3d = vertcat(alldata_mat.gaze3d);

        GazeX = gaze3d(:,1);
        GazeY = gaze3d(:,2);
        GazeZ = gaze3d(:,3);

        % Preprocessing - PUPILS pipeline
        % 1) Calculate angular velocity
        [az, el] = calcAngular(GazeX,GazeY,GazeZ);
        gazeAz_velocity = [0;diff(az)/(1/Param.Fs)];
        gazeEl_velocity = [0;diff(el)/(1/Param.Fs)];
        gaze_vel=sqrt(gazeAz_velocity.^2.*cosd(el).^2+gazeEl_velocity.^2);

        % 2) Select preprocessing options
        options = struct;
        options.fs = Param.Fs;            % sampling frequency (Hz)
        options.blink_rule = 'std';       % Rule for blink detection 'std' / 'vel'
        options.pre_blink_t   = 100;      % region to interpolate before blink (ms)
        options.post_blink_t  = 200;      % region to interpolate after blink (ms)
        options.xy_units = 'mm';          % xy coordinate units 'px' / 'mm' / 'cm'
        options.vel_threshold =  30;      % velocity threshold for saccade detection
        options.min_sacc_duration = 10;   % minimum saccade duration (ms)
        options.interpolate_saccades = 0; % Specify whether saccadic distortions should be interpolated 1-yes 0-noB
        options.pre_sacc_t   = 50;        % Region to interpolate before saccade (ms)
        options.post_sacc_t  = 100;       % Region to interpolate after saccade (ms)
        options.low_pass_fc   = 10;       % Low-pass filter cut-off frequency (Hz)

        % 3) Using PUPILS toolbox 
        %   -IN [samples]: [timestamps, X, Y, Z, Diameter, gaze_vel]
        %   -OUT[samples]: [timestamps, X, Y, Z, Diameter, gaze_vel, blinks, saccades, fixations, interpolated, denoised]
        data = [timestamps*Param.Fs GazeX GazeY DiameterRaw gaze_vel];
        [proc_data, proc_info] = processPupilData(data, options);

        % 4) Fixation duration - from Fixation samples (from PUPILS)
        GazeX_fix = GazeX; GazeY_fix = GazeY; GazeZ_fix = GazeZ;
        GazeX_fix(proc_data(:,8)==0)=NaN; GazeY_fix(proc_data(:,8)==0)=NaN; GazeZ_fix(proc_data(:,8)==0)=NaN;
        Fixation = fix_duration([GazeX_fix,GazeY_fix,GazeZ_fix],MaxVisualAngle,Param.Fs);        

        % 5) Filtering - Fixation duration with hamming window (F = 3)
        Fixation = ndnanfilter(Fixation,'hamming',FilterWidth);
        
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
            RepB = 0;
        elseif contains(cell2mat(FileNames(i)),'B2')
            RepB = 1;
        end

        if contains(cell2mat(FileNames(i)),'Quiet')
            CondKey = RepB + 1;
        elseif contains(cell2mat(FileNames(i)),'SHL')
            CondKey = RepB + 3;
        elseif contains(cell2mat(FileNames(i)),'Noise60')
            CondKey = RepB + 5;
        elseif contains(cell2mat(FileNames(i)),'Noise70')
            CondKey = RepB + 7;
        end

        SpeakRaw = PairUtt{1,CondKey}.(SpeakKey);
        ListenRaw = PairUtt{1,CondKey}.(ListenKey);
        binResUtt = PairUtt{1,CondKey}.binRes;
        SDelayRaw = PairDelay{1,CondKey}.(SDelayKey);
        LDelayRaw = PairDelay{1,CondKey}.(LDelayKey);
        
        if or(SDelayRaw < 0,LDelayRaw < 0)
            SDelayRaw=[0,0];
            LDelayRaw=[0,0];
        end
        
        binResDel = PairDelay{1,CondKey}.binRes;
        
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
        WSpeakIdx=[Speak(:,2)-TimeStartW*Param.Fs,Speak(:,2),Speak(:,2)+TimeEndW*Param.Fs];
        WListenIdx=[Listen(:,2)-TimeStartW*Param.Fs,Listen(:,2),Listen(:,2)+TimeEndW*Param.Fs];
%         EWSpeakIdx=[Speak(:,3)-TimeStartW*Param.Fs,Speak(:,3),Speak(:,3)+TimeEndW*Param.Fs];
%         EWListenIdx=[Listen(:,3)-TimeStartW*Param.Fs,Listen(:,3),Listen(:,3)+TimeEndW*Param.Fs];
        
        % Storing Speaking/Listening by conditions
        if contains(cell2mat(FileNames(i)),'Quiet')
            for j=1:size(Speak,1)
                SW_Quiet(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3)-1);NaN*ones(1,length(SW_Quiet)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SDur_Quiet(i,j) = Speak(j,1);
            end
            for j=1:size(Listen,1)
                LW_Quiet(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3)-1);NaN*ones(1,length(LW_Quiet)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LDur_Quiet(i,j) = Listen(j,1);
            end
        elseif contains(cell2mat(FileNames(i)),'SHL')
            for j=1:size(Speak,1)
                SW_SHL(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3)-1);NaN*ones(1,length(SW_SHL)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SDur_SHL(i,j) = Speak(j,1);
            end
            for j=1:size(Listen,1)
                LW_SHL(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3)-1);NaN*ones(1,length(LW_SHL)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LDur_SHL(i,j) = Listen(j,1);
            end
        elseif contains(cell2mat(FileNames(i)),'Noise60')
            for j=1:size(Speak,1)
                SW_N60(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3)-1);NaN*ones(1,length(SW_N60)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SDur_N60(i,j) = Speak(j,1);
            end
            for j=1:size(Listen,1)
                LW_N60(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3)-1);NaN*ones(1,length(LW_N60)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LDur_N60(i,j) = Listen(j,1);
            end
        elseif contains(cell2mat(FileNames(i)),'Noise70')
            for j=1:size(Speak,1)
                SW_N70(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3)-1);NaN*ones(1,length(SW_N70)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SDur_N70(i,j) = Speak(j,1);
            end
            for j=1:size(Listen,1)
                LW_N70(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3)-1);NaN*ones(1,length(LW_N70)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LDur_N70(x,j) = Listen(j,1);
            end
        end
        
        SW = zeros(size(Speak,1),ceil(Param.Fs*(TimeStartW+max(Speak(:,1)))));
        for j=1:size(Speak,1)
            % Add nan-padding when necessary
            SW(j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3)-1);NaN*ones(1,length(SW)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
            GSW(x,j,1:length(SW(j,:))) = SW(j,:);
            GSDur(x,j) = Speak(j,1);
        end
        
        LW = zeros(size(Listen,1),ceil(Param.Fs*(TimeStartW+max(Listen(:,1)))));
        for j=1:size(Listen,1)
            % Add nan-padding when necessary
            LW(j,:)=[Fixation(WListenIdx(j,1):Listen(j,3)-1);NaN*ones(1,length(LW)-length(WListenIdx(j,1):Listen(j,3)-1))'];
            GLW(x,j,1:length(LW(j,:))) = LW(j,:);
            GLDur(x,j) = Listen(j,1);
        end
        
        % Increase index of num of files used
        x=x+1;
    end
end

%% Avg across TPs
% Clean empty rows and layers, set 0's to NaN
GSW(~any(GSW,[2 3]),:,:)=[];GSW(:,~any(GSW,[1 3]),:)=[];GSW(GSW==0)=NaN;
GLW(~any(GLW,[2 3]),:,:)=[];GLW(:,~any(GLW,[1 3]),:)=[];GLW(GLW==0)=NaN;

SW_Quiet(~any(SW_Quiet,[2 3]),:,:)=[];SW_Quiet(:,~any(SW_Quiet,[1 3]),:)=[];SW_Quiet(SW_Quiet==0)=NaN;
LW_Quiet(~any(LW_Quiet,[2 3]),:,:)=[];LW_Quiet(:,~any(LW_Quiet,[1 3]),:)=[];LW_Quiet(LW_Quiet==0)=NaN;
SW_SHL(~any(SW_SHL,[2 3]),:,:)=[];SW_SHL(:,~any(SW_SHL,[1 3]),:)=[];SW_SHL(SW_SHL==0)=NaN;
LW_SHL(~any(LW_SHL,[2 3]),:,:)=[];LW_SHL(:,~any(LW_SHL,[1 3]),:)=[];LW_SHL(LW_SHL==0)=NaN;
SW_N60(~any(SW_N60,[2 3]),:,:)=[];SW_N60(:,~any(SW_N60,[1 3]),:)=[];SW_N60(SW_N60==0)=NaN;
LW_N60(~any(LW_N60,[2 3]),:,:)=[];LW_N60(:,~any(LW_N60,[1 3]),:)=[];LW_N60(LW_N60==0)=NaN;
SW_N70(~any(SW_N70,[2 3]),:,:)=[];SW_N70(:,~any(SW_N70,[1 3]),:)=[];SW_N70(SW_N70==0)=NaN;
LW_N70(~any(LW_N70,[2 3]),:,:)=[];LW_N70(:,~any(LW_N70,[1 3]),:)=[];LW_N70(LW_N70==0)=NaN;

TP_GSW = nan*ones(NTPs,NCols); TP_GLW = TP_GSW; % Global TP Speaking/Listening Windows
TP_SW_Quiet = TP_GSW; TP_LW_Quiet = TP_GSW; % Quiet TP Speaking/Listening Windows
TP_SW_SHL = TP_GSW; TP_LW_SHL = TP_GSW; % SHL TP Speaking/Listening Windows
TP_SW_N60 = TP_GSW; TP_LW_N60 = TP_GSW; % N60 TP Speaking/Listening Windows
TP_SW_N70 = TP_GSW; TP_LW_N70 = TP_GSW; % N70 TP Speaking/Listening Windows
for i=1:NTPs    
    if ~isempty(nonzeros(LoadTPsOrder.TPsOrder(i,:)))
        TP_GSW(i,:) = (reshape(mean(GSW(min(nonzeros(LoadTPsOrder.TPsOrder(i,:))):max(nonzeros(LoadTPsOrder.TPsOrder(i,:))),:,:),[1 2],'omitnan'),[],1))';
        TP_GLW(i,:) = (reshape(mean(GLW(min(nonzeros(LoadTPsOrder.TPsOrder(i,:))):max(nonzeros(LoadTPsOrder.TPsOrder(i,:))),:,:),[1 2],'omitnan'),[],1))';
        
        % Quiet condition
        if ~isempty(nonzeros(LoadTPsOrder.TPsOrder(i,1:2)))
            TPCondIdx = nonzeros(LoadTPsOrder.TPsOrder(i,1:2));
            if length(TPCondIdx)==1
                TPCondIdx = TPCondIdx*ones(2,1);
            end
            
            TP_SW_Quiet(i,:) = (reshape(mean(GSW(TPCondIdx(1):TPCondIdx(2),:,:),[1 2],'omitnan'),[],1))';
            TP_LW_Quiet(i,:) = (reshape(mean(GLW(TPCondIdx(1):TPCondIdx(2),:,:),[1 2],'omitnan'),[],1))';
        end
        
        % SHL condition
        if ~isempty(nonzeros(LoadTPsOrder.TPsOrder(i,3:4)))
            TPCondIdx = nonzeros(LoadTPsOrder.TPsOrder(i,3:4));
            if length(TPCondIdx)==1
                TPCondIdx = TPCondIdx*ones(2,1);
            end
            
            TP_SW_SHL(i,:) = (reshape(mean(GSW(TPCondIdx(1):TPCondIdx(2),:,:),[1 2],'omitnan'),[],1))';
            TP_LW_SHL(i,:) = (reshape(mean(GLW(TPCondIdx(1):TPCondIdx(2),:,:),[1 2],'omitnan'),[],1))';
        end
        
        % N60 condition
        if ~isempty(nonzeros(LoadTPsOrder.TPsOrder(i,5:6)))
            TPCondIdx = nonzeros(LoadTPsOrder.TPsOrder(i,5:6));
            if length(TPCondIdx)==1
                TPCondIdx = TPCondIdx*ones(2,1);
            end
            
            TP_SW_N60(i,:) = (reshape(mean(GSW(TPCondIdx(1):TPCondIdx(2),:,:),[1 2],'omitnan'),[],1))';
            TP_LW_N60(i,:) = (reshape(mean(GLW(TPCondIdx(1):TPCondIdx(2),:,:),[1 2],'omitnan'),[],1))';
        end
        
        % N70 condition
        if ~isempty(nonzeros(LoadTPsOrder.TPsOrder(i,7:8)))
            TPCondIdx = nonzeros(LoadTPsOrder.TPsOrder(i,7:8));
            if length(TPCondIdx)==1
                TPCondIdx = TPCondIdx*ones(2,1);
            end
            
            TP_SW_N70(i,:) = (reshape(mean(GSW(TPCondIdx(1):TPCondIdx(2),:,:),[1 2],'omitnan'),[],1))';
            TP_LW_N70(i,:) = (reshape(mean(GLW(TPCondIdx(1):TPCondIdx(2),:,:),[1 2],'omitnan'),[],1))';
        end
    end
end



TP_GSW_Mean = ndnanfilter(mean(TP_GSW,'omitnan'),'hamming',FilterWidth);
TP_GLW_Mean = ndnanfilter(mean(TP_GLW,'omitnan'),'hamming',FilterWidth);
TP_GSW_SEM = (std(TP_GSW,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(GSW),[1 2])))');
TP_GLW_SEM = (std(TP_GLW,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(GLW),[1 2])))');

TP_SW_Quiet_Mean = ndnanfilter(mean(TP_SW_Quiet,'omitnan'),'hamming',FilterWidth);
TP_LW_Quiet_Mean = ndnanfilter(mean(TP_LW_Quiet,'omitnan'),'hamming',FilterWidth);
TP_SW_Quiet_SEM = (std(TP_SW_Quiet,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(SW_Quiet),[1 2])))');
TP_LW_Quiet_SEM = (std(TP_LW_Quiet,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(LW_Quiet),[1 2])))');

TP_SW_SHL_Mean = ndnanfilter(mean(TP_SW_SHL,'omitnan'),'hamming',FilterWidth);
TP_LW_SHL_Mean = ndnanfilter(mean(TP_LW_SHL,'omitnan'),'hamming',FilterWidth);
TP_SW_SHL_SEM = (std(TP_SW_SHL,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(SW_SHL),[1 2])))');
TP_LW_SHL_SEM = (std(TP_LW_SHL,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(LW_SHL),[1 2])))');

TP_SW_N60_Mean = ndnanfilter(mean(TP_SW_N60,'omitnan'),'hamming',FilterWidth);
TP_LW_N60_Mean = ndnanfilter(mean(TP_LW_N60,'omitnan'),'hamming',FilterWidth);
TP_SW_N60_SEM = (std(TP_SW_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(SW_N60),[1 2])))');
TP_LW_N60_SEM = (std(TP_LW_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(LW_N60),[1 2])))');

TP_SW_N70_Mean = ndnanfilter(mean(TP_SW_N70,'omitnan'),'hamming',FilterWidth);
TP_LW_N70_Mean = ndnanfilter(mean(TP_LW_N70,'omitnan'),'hamming',FilterWidth);
TP_SW_N70_SEM = (std(TP_SW_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(SW_N70),[1 2])))');
TP_LW_N70_SEM = (std(TP_LW_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(LW_N70),[1 2])))');

%% Global Plots
figure;t1=tiledlayout(1,2);ax1 = nexttile;ax2 = nexttile;
figure;t2=tiledlayout(1,2);ax5 = nexttile;ax6 = nexttile;
hold([ax1 ax2 ax5 ax6],'on')
grid([ax1 ax2 ax5 ax6],'on')

% Global Speak/Listening durations
data5 = [mean(nonzeros(GSDur));mean(nonzeros(GLDur))];
b5 = bar(ax5,data5,'grouped','FaceColor','flat');
[ngroups,nbars] = size(data5);
x5 = nan(nbars, ngroups);
y5 = x5;
b5.FaceColor = 'flat';
ColBarW = {SColor,LColor};
for e = 1:nbars
    x5(e,:) = b5(e).XEndPoints;
    y5(e,:) = b5(e).YEndPoints;
end
for t = 1:length(ColBarW)
    b5.CData(t,:) = cell2mat(ColBarW(t));
end
errorbar(ax5,x5',data5,[std(nonzeros(GSDur))/sqrt(numel(nonzeros(GSDur)));std(nonzeros(GLDur))/sqrt(numel(nonzeros(GLDur)))],'k','linestyle','none','handlevisibility' ,'off')
text(ax5,x5,y5,num2str(data5,'%0.2f'),'HorizontalAlignment','center','VerticalAlignment','top','Color',[1 1 1 1],'FontWeight','bold')

bar(ax6,nan,'FaceColor',SColor);bar(ax6,nan,'FaceColor',LColor)
data6 = [mean(nonzeros(SDur_Quiet)),mean(nonzeros(SDur_SHL)),mean(nonzeros(SDur_N60)),mean(nonzeros(SDur_N70));mean(nonzeros(LDur_Quiet)),mean(nonzeros(LDur_SHL)),mean(nonzeros(LDur_N60)),mean(nonzeros(LDur_N70))];
b6 = bar(ax6,data6,'grouped');
[ngroups,nbars] = size(data6);
x6 = nan(nbars, ngroups);
ColBarCond = {QuietColor,SHLColor,N60Color,N70Color};
for e = 1:nbars
    x6(e,:) = b6(e).XEndPoints;
end
for t = 1:length(ColBarCond)
    b6(t).FaceColor = cell2mat(ColBarCond(t));
end
errorbar(ax6,x6',data6,[std(nonzeros(SDur_Quiet))/sqrt(numel(nonzeros(SDur_Quiet))),std(nonzeros(SDur_SHL))/sqrt(numel(nonzeros(SDur_SHL))),std(nonzeros(SDur_N60))/sqrt(numel(nonzeros(SDur_N60))),std(nonzeros(SDur_N70))/sqrt(numel(nonzeros(SDur_N70)));std(nonzeros(LDur_Quiet))/sqrt(numel(nonzeros(LDur_Quiet))),std(nonzeros(LDur_SHL))/sqrt(numel(nonzeros(LDur_SHL))),std(nonzeros(LDur_N60))/sqrt(numel(nonzeros(LDur_N60))),std(nonzeros(LDur_N70))/sqrt(numel(nonzeros(LDur_N70)))],'k','linestyle','none','handlevisibility' ,'off')

xticks([ax5 ax6],1:2)
xticklabels([ax5 ax6],{'Speaking','Listening'})
xlim(ax5,ax6.XLim)
ylim(ax5,ax6.YLim)
lgd6=legend(ax6,'Speaking','Listening','Quiet','SHL','N60','N70','Location','southeastoutside');
lgd6.Title.String = 'Types of Windows:';
title(ax5,'Global average window duration')
title(ax6,'Condition-based average window duration')
ylabel([ax5 ax6],'Time [s]')

xline(ax1,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax2,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')

% plot(ax1,linspace(-TimeStartW,size(TP_GSW,2)/Param.Fs,size(TP_GSW,2)),TP_GSW_Mean,color=SColor,linewidth=2)
plot(ax1,linspace(-TimeStartW,size(TP_SW_Quiet,2)/Param.Fs,size(TP_SW_Quiet,2)),TP_SW_Quiet_Mean,color=QuietColor)
plot(ax1,linspace(-TimeStartW,size(TP_SW_SHL,2)/Param.Fs,size(TP_SW_SHL,2)),TP_SW_SHL_Mean,color=SHLColor)
plot(ax1,linspace(-TimeStartW,size(TP_SW_N60,2)/Param.Fs,size(TP_SW_N60,2)),TP_SW_N60_Mean,color=N60Color)
plot(ax1,linspace(-TimeStartW,size(TP_SW_N70,2)/Param.Fs,size(TP_SW_N70,2)),TP_SW_N70_Mean,color=N70Color)

% plot(ax2,nan,color=SColor,linewidth=2) % plot nans to show color in legend
% plot(ax2,linspace(-TimeStartW,size(TP_GLW,2)/Param.Fs,size(TP_GLW,2)),TP_GLW_Mean,color=LColor,linewidth=2)
plot(ax2,linspace(-TimeStartW,size(TP_LW_Quiet,2)/Param.Fs,size(TP_LW_Quiet,2)),TP_LW_Quiet_Mean,color=QuietColor)
plot(ax2,linspace(-TimeStartW,size(TP_LW_SHL,2)/Param.Fs,size(TP_LW_SHL,2)),TP_LW_SHL_Mean,color=SHLColor)
plot(ax2,linspace(-TimeStartW,size(TP_LW_N60,2)/Param.Fs,size(TP_LW_N60,2)),TP_LW_N60_Mean,color=N60Color)
plot(ax2,linspace(-TimeStartW,size(TP_LW_N70,2)/Param.Fs,size(TP_LW_N70,2)),TP_LW_N70_Mean,color=N70Color)

% plot(ax1,linspace(-TimeStartW,size(TP_SW_SHL,2)/Param.Fs,size(TP_SW_SHL,2)),mean([TP_SW_Quiet_Mean;TP_SW_SHL_Mean;TP_SW_N60_Mean;TP_SW_N70_Mean],1,'omitnan'),'k--')
% plot(ax2,linspace(-TimeStartW,size(TP_LW_SHL,2)/Param.Fs,size(TP_LW_SHL,2)),mean([TP_LW_Quiet_Mean;TP_LW_SHL_Mean;TP_LW_N60_Mean;TP_LW_N70_Mean],1,'omitnan'),'k--')

TP_GSW_Mean(isnan(TP_GSW_Mean))=0;TP_GSW_SEM(isnan(TP_GSW_SEM))=0;
TP_SW_Quiet_Mean(isnan(TP_SW_Quiet_Mean))=0;TP_SW_Quiet_SEM(isnan(TP_SW_Quiet_SEM))=0;
TP_SW_SHL_Mean(isnan(TP_SW_SHL_Mean))=0;TP_SW_SHL_SEM(isnan(TP_SW_SHL_SEM))=0;
TP_SW_N60_Mean(isnan(TP_SW_N60_Mean))=0;TP_SW_N60_SEM(isnan(TP_SW_N60_SEM))=0;
TP_SW_N70_Mean(isnan(TP_SW_N70_Mean))=0;TP_SW_N70_SEM(isnan(TP_SW_N70_SEM))=0;
TP_GLW_Mean(isnan(TP_GLW_Mean))=0;TP_GLW_SEM(isnan(TP_GLW_SEM))=0;
TP_LW_Quiet_Mean(isnan(TP_LW_Quiet_Mean))=0;TP_LW_Quiet_SEM(isnan(TP_LW_Quiet_SEM))=0;
TP_LW_SHL_Mean(isnan(TP_LW_SHL_Mean))=0;TP_LW_SHL_SEM(isnan(TP_LW_SHL_SEM))=0;
TP_LW_N60_Mean(isnan(TP_LW_N60_Mean))=0;TP_LW_N60_SEM(isnan(TP_LW_N60_SEM))=0;
TP_LW_N70_Mean(isnan(TP_LW_N70_Mean))=0;TP_LW_N70_SEM(isnan(TP_LW_N70_SEM))=0;

% fill(ax1,[linspace(-TimeStartW,size(TP_GSW,2)/Param.Fs,size(TP_GSW,2)), flipud(linspace(-TimeStartW,size(TP_GSW,2)/Param.Fs,size(TP_GSW,2))')'],[(TP_GSW_Mean+TP_GSW_SEM), flipud((TP_GSW_Mean-TP_GSW_SEM)')'],SColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(TP_SW_Quiet,2)/Param.Fs,size(TP_SW_Quiet,2)), flipud(linspace(-TimeStartW,size(TP_SW_Quiet,2)/Param.Fs,size(TP_SW_Quiet,2))')'],[(TP_SW_Quiet_Mean+TP_SW_Quiet_SEM), flipud((TP_SW_Quiet_Mean-TP_SW_Quiet_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(TP_SW_SHL,2)/Param.Fs,size(TP_SW_SHL,2)), flipud(linspace(-TimeStartW,size(TP_SW_SHL,2)/Param.Fs,size(TP_SW_SHL,2))')'],[(TP_SW_SHL_Mean+TP_SW_SHL_SEM), flipud((TP_SW_SHL_Mean-TP_SW_SHL_SEM)')'],SHLColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(TP_SW_N60,2)/Param.Fs,size(TP_SW_N60,2)), flipud(linspace(-TimeStartW,size(TP_SW_N60,2)/Param.Fs,size(TP_SW_N60,2))')'],[(TP_SW_N60_Mean+TP_SW_N60_SEM), flipud((TP_SW_N60_Mean-TP_SW_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(TP_SW_N70,2)/Param.Fs,size(TP_SW_N70,2)), flipud(linspace(-TimeStartW,size(TP_SW_N70,2)/Param.Fs,size(TP_SW_N70,2))')'],[(TP_SW_N70_Mean+TP_SW_N70_SEM), flipud((TP_SW_N70_Mean-TP_SW_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')

% fill(ax2,[linspace(-TimeStartW,size(TP_GLW,2)/Param.Fs,size(TP_GLW,2)), flipud(linspace(-TimeStartW,size(TP_GLW,2)/Param.Fs,size(TP_GLW,2))')'],[(TP_GLW_Mean+TP_GLW_SEM), flipud((TP_GLW_Mean-TP_GLW_SEM)')'],LColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(TP_LW_Quiet,2)/Param.Fs,size(TP_LW_Quiet,2)), flipud(linspace(-TimeStartW,size(TP_LW_Quiet,2)/Param.Fs,size(TP_LW_Quiet,2))')'],[(TP_LW_Quiet_Mean+TP_LW_Quiet_SEM), flipud((TP_LW_Quiet_Mean-TP_LW_Quiet_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(TP_LW_SHL,2)/Param.Fs,size(TP_LW_SHL,2)), flipud(linspace(-TimeStartW,size(TP_LW_SHL,2)/Param.Fs,size(TP_LW_SHL,2))')'],[(TP_LW_SHL_Mean+TP_LW_SHL_SEM), flipud((TP_LW_SHL_Mean-TP_LW_SHL_SEM)')'],SHLColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(TP_LW_N60,2)/Param.Fs,size(TP_LW_N60,2)), flipud(linspace(-TimeStartW,size(TP_LW_N60,2)/Param.Fs,size(TP_LW_N60,2))')'],[(TP_LW_N60_Mean+TP_LW_N60_SEM), flipud((TP_LW_N60_Mean-TP_LW_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(TP_LW_N70,2)/Param.Fs,size(TP_LW_N70,2)), flipud(linspace(-TimeStartW,size(TP_LW_N70,2)/Param.Fs,size(TP_LW_N70,2))')'],[(TP_LW_N70_Mean+TP_LW_N70_SEM), flipud((TP_LW_N70_Mean-TP_LW_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')

xlabel([ax1 ax2],'Time [s]')
ylabel(ax1,'Fixation duration [s]')
xlim([ax1 ax2],[-TimeStartW 3])
ylim([ax1 ax2],[min([ylim(ax1),ylim(ax2)]), max([ylim(ax1),ylim(ax2)])])

% lgd2=legend(ax2,'Speaking','Listening','Quiet','SHL','N60','N70','Location','southeastoutside');
% lgd2.Title.String = 'Types of windows:';
% lgd4=legend(ax4,'Speaking','Listening','Quiet','SHL','N60','N70','Location','southeastoutside');
% lgd4.Title.String = 'Types of windows:';
% title(ax1,'Global Speaking-evoked averaged across TPs response')
% title(ax2,'Global Listening-evoked averaged across TPs response')

set(ax1,'Color',[SColor,0.04])
set(ax2,'Color',[LColor,0.04])


