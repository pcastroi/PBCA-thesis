%% PBCA-Thesis - Week 4, Week 5, Week 6 - Groups and Features -- Figures. Fixed or adaptive baseline?
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
WP = [-0.5,3]; % [s], Duration of analysis window
AdapBL = 0.1; % [s], Duration of baseline prior to event
BL = [15,20]; % [s]
BLPeriod = [0,20]; % [s]
BLStartEnd = BLPeriod*Param.Fs + 1; % [samples]
TimeMinWin = 0.5; % [s], Minimum time of a window
TimeInitialMerge = 0.3; % [s], Time threshold for merging windows initially
TimeMerge = 2; % [s], Time threshold for merging windows after rejecting small windows
RejectRatio = 0.4; % Rejection threshold based on the ratio of NaNs in data
RejectDelay = 0.5; % [s], Rejection threshold based on delay between timestamps and n-samples

% Preallocate groups
G1 = zeros(200,BLStartEnd(2)); % Initially Speaking: [0,20] s
G2 = G1; % Initially Listening: [0,20] s
G3 = zeros(200,(WP(2)-WP(1))*Param.Fs+1); % Initially Speaking: [-3,6] s with event starting at 0 s.
G4 = G3; % Initially Listening: [-3,6] s with event starting at 0 s.
G5 = G3; % Initially Speaking: Baselined [15,20] s: [-3,6] s with event starting at 0 s.
G6 = G3; % Initially Listening: Baselined [15,20] s: [-3,6] s with event starting at 0 s.
G7 = zeros(200,200,(WP(2)-WP(1))*Param.Fs+1); % All Speaking windows: [-3,6] s with event starting at 0 s.
G8 = G7; % All Listening windows: [-3,6] s with event starting at 0 s.
G9 = G7; % All Speaking windows: Baselined [15,20] s: [-3,6] s with event starting at 0 s.
G10 = G7; % All Listening windows: Baselined [15,20] s: [-3,6] s with event starting at 0 s.
G11 = G7; % All Speaking windows: Baselined 1 s prior to event: [-3,6] s with event starting at 0 s.
G12 = G7; % All Listening windows: Baselined 1 s prior to event: [-3,6] s with event starting at 0 s.
G13 = G7; % Quiet Speaking windows: [-3,6] s with event starting at 0 s.
G14 = G7; % Quiet Listening windows: [-3,6] s with event starting at 0 s.
G15 = G7; % SHL Speaking windows: [-3,6] s with event starting at 0 s.
G16 = G7; % SHL Listening windows: [-3,6] s with event starting at 0 s.
G17 = G7; % N60 Speaking windows: [-3,6] s with event starting at 0 s.
G18 = G7; % N60 Listening windows: [-3,6] s with event starting at 0 s.
G19 = G7; % N70 Speaking windows: [-3,6] s with event starting at 0 s.
G20 = G7; % N70 Listening windows: [-3,6] s with event starting at 0 s.

% Preallocate features
F1_S_Q=zeros(200,1); % Feature 1 (Quiet): Mean pupil size (Speaking)
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
[subDirs] = GetSubDirsFirstLevelOnly('data');
LoadDelays=load('data\delays1110.mat');

for q=1:numel(subDirs)
    PairIn = q;
    PairFiles=dir(['data\Main',sprintf('%d',PairIn),'\*.mat']);
    PairUtt=load('data\utterances1110.mat');
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
        
        %% Plots
        if x==1
            figure;tiledlayout(1,2);
            ax1 = nexttile;
            ax2 = nexttile;
            figure;tiledlayout(2,2);
            ax3 = nexttile;
            ax4 = nexttile;
            ax5 = nexttile;
            ax6 = nexttile;
            figure;tiledlayout(3,2);
            ax7 = nexttile;
            ax8 = nexttile;
            ax9 = nexttile;
            ax10 = nexttile;
            ax11 = nexttile;
            ax12 = nexttile;
            figure;tiledlayout(1,3);
            ax13 = nexttile;
            ax14 = nexttile;
            ax15 = nexttile;
            figure;tiledlayout(1,2);
            ax16 = nexttile;
            ax17 = nexttile;
            figure;tiledlayout(1,2);
            ax18 = nexttile;
            ax19 = nexttile;
            figure;tiledlayout(1,2);
            ax20 = nexttile;
            ax21 = nexttile;
        end
        x=x+1;
        hold([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10 ax11 ax12 ax13 ax14 ax15 ax16 ax17 ax18 ax19 ax20 ax21],'on')
        
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
        
        if Speak(1,2) <= Listen(1,2) % Initially Speaking (or Speaking+Listening at same time)
            G1(x,:)=Diameter(BLStartEnd(1):BLStartEnd(2));
            G3(x,:)=Diameter(Speak(1,2)+WP(1)*Param.Fs:Speak(1,2)+WP(2)*Param.Fs);
            G5(x,:)=Diameter(Speak(1,2)+WP(1)*Param.Fs:Speak(1,2)+WP(2)*Param.Fs)-mean(Diameter(BL(1)*Param.Fs:BL(2)*Param.Fs));
%             plot(ax1,t_Diam(BLStartEnd(1):BLStartEnd(2)),G1(x,:),color=[0 0 0 0.2],linewidth=0.5)
%             plot(ax3,linspace(WP(1),WP(2),size(G3,2)),G3(x,:),color=[0 0 0 0.2],linewidth=0.5)
%             plot(ax5,linspace(WP(1),WP(2),size(G5,2)),G5(x,:),color=[0 0 0 0.2],linewidth=0.5)
        elseif Speak(1,2) > Listen(1,2) % Initially Listening/Non-speaking
            G2(x,:)=Diameter(BLStartEnd(1):BLStartEnd(2));
            G4(x,:)=Diameter(Listen(1,2)+WP(1)*Param.Fs:Listen(1,2)+WP(2)*Param.Fs);
            G6(x,:)=Diameter(Listen(1,2)+WP(1)*Param.Fs:Listen(1,2)+WP(2)*Param.Fs)-mean(Diameter(BL(1)*Param.Fs:BL(2)*Param.Fs));
%             plot(ax2,t_Diam(BLStartEnd(1):BLStartEnd(2)),G2(x,:),color=[0 0 0 0.2],linewidth=0.5)
%             plot(ax4,linspace(WP(1),WP(2),size(G4,2)),G4(x,:),color=[0 0 0 0.2],linewidth=0.5)
%             plot(ax6,linspace(WP(1),WP(2),size(G6,2)),G6(x,:),color=[0 0 0 0.2],linewidth=0.5)
        end
        
        if Speak ~= Listen(1,:)+1
            for j=1:size(Speak,1)
                if Speak(j,2)+WP(2)*Param.Fs <= length(Diameter) % When reaching 'end of tobii data'
                    G7(x,j,:)=Diameter(Speak(j,2)+WP(1)*Param.Fs:Speak(j,2)+WP(2)*Param.Fs);
                    G9(x,j,:)=Diameter(Speak(j,2)+WP(1)*Param.Fs:Speak(j,2)+WP(2)*Param.Fs)-mean(Diameter(BL(1)*Param.Fs:BL(2)*Param.Fs));
                    G11(x,j,:)=Diameter(Speak(j,2)+WP(1)*Param.Fs:Speak(j,2)+WP(2)*Param.Fs)-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));
%                     plot(ax7,linspace(WP(1),WP(2),size(G7,3)),reshape(G7(x,j,:),[],1),color=[0 0 0 0.01],linewidth=0.5)
                end
            end
        end
        
        if Listen ~= Speak(1,:)+1
            for j=1:size(Listen,1)
                if Listen(j,2)+WP(2)*Param.Fs <= length(Diameter) % When reaching 'end of tobii data'
                    G8(x,j,:)=Diameter(Listen(j,2)+WP(1)*Param.Fs:Listen(j,2)+WP(2)*Param.Fs);
                    G10(x,j,:)=Diameter(Listen(j,2)+WP(1)*Param.Fs:Listen(j,2)+WP(2)*Param.Fs)-mean(Diameter(BL(1)*Param.Fs:BL(2)*Param.Fs));
                    G12(x,j,:)=Diameter(Listen(j,2)+WP(1)*Param.Fs:Listen(j,2)+WP(2)*Param.Fs)-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));
%                     plot(ax8,linspace(WP(1),WP(2),size(G8,3)),reshape(G8(x,j,:),[],1),color=[0 0 0 0.01],linewidth=0.5)
                end
            end
        end
        
        % Extract features for different conditions
        temp=zeros(1000,1);
        temp2=zeros(1000,2);
        temp3=temp;
        if contains(PairFiles(i).name,'Quiet')
            for j=1:size(Speak,1)
                temp(j)=mean(Diameter(Speak(j,2):Speak(j,3)));
                temp2(j,:)=polyfit(t_Diam(Speak(j,2):Speak(j,3)),Diameter(Speak(j,2):Speak(j,3)),1);
                temp3(j)=max(Diameter(Speak(j,2):Speak(j,3)));
                if Speak(j,2)+WP(2)*Param.Fs <= length(Diameter) % When reaching 'end of tobii data'
                    G13(x,j,:)=Diameter(Speak(j,2)+WP(1)*Param.Fs:Speak(j,2)+WP(2)*Param.Fs);
                    G13_FB(x,j,:)=Diameter(Speak(j,2)+WP(1)*Param.Fs:Speak(j,2)+WP(2)*Param.Fs)-mean(Diameter(BL(1)*Param.Fs:BL(2)*Param.Fs));
                    G13_AB(x,j,:)=Diameter(Speak(j,2)+WP(1)*Param.Fs:Speak(j,2)+WP(2)*Param.Fs)-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));
                end
            end
            F1_S_Q(x)=mean(nonzeros(temp));
            F2_S_Q(x)=mean(nonzeros(temp2(:,2)));
            F3_S_Q(x)=mean(nonzeros(temp3));
            for j=1:size(Listen,1)
                temp(j)=mean(Diameter(Listen(j,2):Listen(j,3)));
                temp2(j,:)=polyfit(t_Diam(Listen(j,2):Listen(j,3)),Diameter(Listen(j,2):Listen(j,3)),1);
                temp3(j)=max(Diameter(Listen(j,2):Listen(j,3)));
                if Listen(j,2)+WP(2)*Param.Fs <= length(Diameter) % When reaching 'end of tobii data'
                    G14(x,j,:)=Diameter(Listen(j,2)+WP(1)*Param.Fs:Listen(j,2)+WP(2)*Param.Fs);
                    G14_FB(x,j,:)=Diameter(Listen(j,2)+WP(1)*Param.Fs:Listen(j,2)+WP(2)*Param.Fs)-mean(Diameter(BL(1)*Param.Fs:BL(2)*Param.Fs));
                    G14_AB(x,j,:)=Diameter(Listen(j,2)+WP(1)*Param.Fs:Listen(j,2)+WP(2)*Param.Fs)-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));
                end
            end
            F1_L_Q(x)=mean(nonzeros(temp));
            F2_L_Q(x)=mean(nonzeros(temp2(:,2))); % Check if this makes sense
            F3_L_Q(x)=mean(nonzeros(temp3));
        elseif contains(PairFiles(i).name,'SHL')
            for j=1:size(Speak,1)
                temp(j)=mean(Diameter(Speak(j,2):Speak(j,3)));
                temp2(j,:)=polyfit(t_Diam(Speak(j,2):Speak(j,3)),Diameter(Speak(j,2):Speak(j,3)),1);
                temp3(j)=max(Diameter(Speak(j,2):Speak(j,3)));
                if Speak(j,2)+WP(2)*Param.Fs <= length(Diameter) % When reaching 'end of tobii data'
                    G15(x,j,:)=Diameter(Speak(j,2)+WP(1)*Param.Fs:Speak(j,2)+WP(2)*Param.Fs);
                    G15_FB(x,j,:)=Diameter(Speak(j,2)+WP(1)*Param.Fs:Speak(j,2)+WP(2)*Param.Fs)-mean(Diameter(BL(1)*Param.Fs:BL(2)*Param.Fs));
                    G15_AB(x,j,:)=Diameter(Speak(j,2)+WP(1)*Param.Fs:Speak(j,2)+WP(2)*Param.Fs)-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));
                end
            end
            F1_S_SHL(x)=mean(nonzeros(temp));
            F2_S_SHL(x)=mean(nonzeros(temp2(:,2)));
            F3_S_SHL(x)=mean(nonzeros(temp3));
            for j=1:size(Listen,1)
                temp(j)=mean(Diameter(Listen(j,2):Listen(j,3)));
                temp2(j,:)=polyfit(t_Diam(Listen(j,2):Listen(j,3)),Diameter(Listen(j,2):Listen(j,3)),1);
                temp3(j)=max(Diameter(Listen(j,2):Listen(j,3)));
                if Listen(j,2)+WP(2)*Param.Fs <= length(Diameter) % When reaching 'end of tobii data'
                    G16(x,j,:)=Diameter(Listen(j,2)+WP(1)*Param.Fs:Listen(j,2)+WP(2)*Param.Fs);
                    G16_FB(x,j,:)=Diameter(Listen(j,2)+WP(1)*Param.Fs:Listen(j,2)+WP(2)*Param.Fs)-mean(Diameter(BL(1)*Param.Fs:BL(2)*Param.Fs));
                    G16_AB(x,j,:)=Diameter(Listen(j,2)+WP(1)*Param.Fs:Listen(j,2)+WP(2)*Param.Fs)-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));
                end
            end
            F1_L_SHL(x)=mean(nonzeros(temp));
            F2_L_SHL(x)=mean(nonzeros(temp2(:,2)));
            F3_L_SHL(x)=mean(nonzeros(temp3));
        elseif contains(PairFiles(i).name,'Noise60')
            for j=1:size(Speak,1)
                temp(j)=mean(Diameter(Speak(j,2):Speak(j,3)));
                temp2(j,:)=polyfit(t_Diam(Speak(j,2):Speak(j,3)),Diameter(Speak(j,2):Speak(j,3)),1);
                temp3(j)=max(Diameter(Speak(j,2):Speak(j,3)));
                if Speak(j,2)+WP(2)*Param.Fs <= length(Diameter) % When reaching 'end of tobii data'
                    G17(x,j,:)=Diameter(Speak(j,2)+WP(1)*Param.Fs:Speak(j,2)+WP(2)*Param.Fs);
                    G17_FB(x,j,:)=Diameter(Speak(j,2)+WP(1)*Param.Fs:Speak(j,2)+WP(2)*Param.Fs)-mean(Diameter(BL(1)*Param.Fs:BL(2)*Param.Fs));
                    G17_AB(x,j,:)=Diameter(Speak(j,2)+WP(1)*Param.Fs:Speak(j,2)+WP(2)*Param.Fs)-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));
                end
            end
            F1_S_N60(x)=mean(nonzeros(temp));
            F2_S_N60(x)=mean(nonzeros(temp2(:,2)));
            F3_S_N60(x)=mean(nonzeros(temp3));
            for j=1:size(Listen,1)
                temp(j)=mean(Diameter(Listen(j,2):Listen(j,3)));
                temp2(j,:)=polyfit(t_Diam(Listen(j,2):Listen(j,3)),Diameter(Listen(j,2):Listen(j,3)),1);
                temp3(j)=max(Diameter(Listen(j,2):Listen(j,3)));
                if Listen(j,2)+WP(2)*Param.Fs <= length(Diameter) % When reaching 'end of tobii data'
                    G18(x,j,:)=Diameter(Listen(j,2)+WP(1)*Param.Fs:Listen(j,2)+WP(2)*Param.Fs);
                    G18_FB(x,j,:)=Diameter(Listen(j,2)+WP(1)*Param.Fs:Listen(j,2)+WP(2)*Param.Fs)-mean(Diameter(BL(1)*Param.Fs:BL(2)*Param.Fs));
                    G18_AB(x,j,:)=Diameter(Listen(j,2)+WP(1)*Param.Fs:Listen(j,2)+WP(2)*Param.Fs)-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));
                end
            end
            F1_L_N60(x)=mean(nonzeros(temp));
            F2_L_N60(x)=mean(nonzeros(temp2(:,2)));
            F3_L_N60(x)=mean(nonzeros(temp3));
        elseif contains(PairFiles(i).name,'Noise70')
            for j=1:size(Speak,1)
                temp(j)=mean(Diameter(Speak(j,2):Speak(j,3)));
                temp2(j,:)=polyfit(t_Diam(Speak(j,2):Speak(j,3)),Diameter(Speak(j,2):Speak(j,3)),1);
                temp3(j)=max(Diameter(Speak(j,2):Speak(j,3)));
                if Speak(j,2)+WP(2)*Param.Fs <= length(Diameter) % When reaching 'end of tobii data'
                    G19(x,j,:)=Diameter(Speak(j,2)+WP(1)*Param.Fs:Speak(j,2)+WP(2)*Param.Fs);
                    G19_FB(x,j,:)=Diameter(Speak(j,2)+WP(1)*Param.Fs:Speak(j,2)+WP(2)*Param.Fs)-mean(Diameter(BL(1)*Param.Fs:BL(2)*Param.Fs));
                    G19_AB(x,j,:)=Diameter(Speak(j,2)+WP(1)*Param.Fs:Speak(j,2)+WP(2)*Param.Fs)-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));
                end
            end
            F1_S_N70(x)=mean(nonzeros(temp));
            F2_S_N70(x)=mean(nonzeros(temp2(:,2)));
            F3_S_N70(x)=mean(nonzeros(temp3));
            for j=1:size(Listen,1)
                temp(j)=mean(Diameter(Listen(j,2):Listen(j,3)));
                temp2(j,:)=polyfit(t_Diam(Listen(j,2):Listen(j,3)),Diameter(Listen(j,2):Listen(j,3)),1);
                temp3(j)=max(Diameter(Listen(j,2):Listen(j,3)));
                if Listen(j,2)+WP(2)*Param.Fs <= length(Diameter) % When reaching 'end of tobii data'
                    G20(x,j,:)=Diameter(Listen(j,2)+WP(1)*Param.Fs:Listen(j,2)+WP(2)*Param.Fs);
                    G20_FB(x,j,:)=Diameter(Listen(j,2)+WP(1)*Param.Fs:Listen(j,2)+WP(2)*Param.Fs)-mean(Diameter(BL(1)*Param.Fs:BL(2)*Param.Fs));
                    G20_AB(x,j,:)=Diameter(Listen(j,2)+WP(1)*Param.Fs:Listen(j,2)+WP(2)*Param.Fs)-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));
                end
            end
            F1_L_N70(x)=mean(nonzeros(temp));
            F2_L_N70(x)=mean(nonzeros(temp2(:,2)));
            F3_L_N70(x)=mean(nonzeros(temp3));
        end
        
        % Last file = \Main12\P1_SHL_B2.mat -> default rejection criteria
        if contains([PairFiles(i).folder, '\', PairFiles(i).name],'\Main12\P1_SHL_B2.mat')  
            title(ax1,'G1: BL Period Initially Speaking')
            title(ax2,'G2: BL Period Initially Listening')
            title(ax3,'G3: Speaking-Evoked')
            title(ax4,'G4: Listening-Evoked')
            title(ax5,'G5: Baselined Speaking-Evoked')
            title(ax6,'G6: Baselined Listening-Evoked')
            title(ax7,'G7: Global Speaking-evoked response')
            title(ax8,'G8: Global Listening-evoked response')
            title(ax9,['G9: Baselined [',num2str(BL(1)),',',num2str(BL(2)),'] global Speaking-evoked response'])
            title(ax10,['G10: Baselined [',num2str(BL(1)),',',num2str(BL(2)),'] global Listening-evoked response'])
            title(ax11,'G11: Baselined (adaptive) global Speaking-evoked response')
            title(ax12,'G12: Baselined (adaptive) global Listening-evoked response')
            title(ax13,'Feature 1: Mean Pupil Size')
            title(ax14,'Feature 2: Mean Slope')
            title(ax15,'Feature 3: Mean Peak Pupil Size')
            title(ax16,'Global Speaking-evoked response')
            title(ax17,'Global Listening-evoked response')
            title(ax18,'Global fixed [15,20] s baselined Speaking-evoked response')
            title(ax19,'Global fixed [15,20] s baselined Listening-evoked response')
            title(ax20,'Global adaptive baselined Speaking-evoked response')
            title(ax21,'Global adaptive baselined Listening-evoked response')
            xlabel([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10 ax11 ax12 ax16 ax17 ax18 ax19 ax20 ax21],'Time [s]')
            xlabel([ax13 ax14 ax15],'Conditions')
            ylabel([ax1 ax2 ax3 ax4 ax7 ax8 ax13 ax15 ax16 ax17],'Pupil diameter [mm]')
            ylabel([ax5 ax6 ax9 ax10 ax11 ax12 ax18 ax19 ax20 ax21],'Pupil baseline difference [mm]')
            ylabel(ax14,'Slope [mm/s]')
            xlim([ax1 ax2],[0, 20])
            xlim([ax3 ax4 ax5 ax6 ax16 ax17 ax18 ax19 ax20 ax21],[WP(1), WP(2)])
            grid([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10 ax11 ax12 ax13 ax14 ax15 ax16 ax17 ax18 ax19 ax20 ax21],'on')
            xline(ax3,0,'--')
            xline(ax4,0,'--')
            xline(ax5,0,'--')
            xline(ax6,0,'--')
            xline(ax7,0,'--')
            xline(ax8,0,'--')
            xline(ax9,0,'--')
            xline(ax10,0,'--')
            xline(ax11,0,'--')
            xline(ax12,0,'--')
            xline(ax16,0,'--','handlevisibility' ,'off')
            xline(ax17,0,'--','handlevisibility' ,'off')
            xline(ax18,0,'--','handlevisibility' ,'off')
            xline(ax19,0,'--','handlevisibility' ,'off')
            xline(ax20,0,'--','handlevisibility' ,'off')
            xline(ax21,0,'--','handlevisibility' ,'off')
            G1(~any(G1,2),:)=[];
            G2(~any(G2,2),:)=[];
            G3(~any(G3,2),:)=[];
            G4(~any(G4,2),:)=[];
            G5(~any(G5,2),:)=[];
            G6(~any(G6,2),:)=[];
            G7(~any(G7,[2 3]),:,:)=[];G7(:,~any(G7,[1 3]),:)=[];
            G8(~any(G8,[2 3]),:,:)=[];G8(:,~any(G8,[1 3]),:)=[];
            G9(~any(G9,[2 3]),:,:)=[];G9(:,~any(G9,[1 3]),:)=[];
            G10(~any(G10,[2 3]),:,:)=[];G10(:,~any(G10,[1 3]),:)=[];
            G11(~any(G11,[2 3]),:,:)=[];G11(:,~any(G11,[1 3]),:)=[];
            G12(~any(G12,[2 3]),:,:)=[];G12(:,~any(G12,[1 3]),:)=[];
            G13(~any(G13,[2 3]),:,:)=[];G13(:,~any(G13,[1 3]),:)=[];
            G14(~any(G14,[2 3]),:,:)=[];G14(:,~any(G14,[1 3]),:)=[];
            G15(~any(G15,[2 3]),:,:)=[];G15(:,~any(G15,[1 3]),:)=[];
            G16(~any(G16,[2 3]),:,:)=[];G16(:,~any(G16,[1 3]),:)=[];
            G17(~any(G17,[2 3]),:,:)=[];G17(:,~any(G17,[1 3]),:)=[];
            G18(~any(G18,[2 3]),:,:)=[];G18(:,~any(G18,[1 3]),:)=[];
            G19(~any(G19,[2 3]),:,:)=[];G19(:,~any(G19,[1 3]),:)=[];
            G20(~any(G20,[2 3]),:,:)=[];G20(:,~any(G20,[1 3]),:)=[];
            G13_FB(~any(G13_FB,[2 3]),:,:)=[];G13_FB(:,~any(G13_FB,[1 3]),:)=[];
            G14_FB(~any(G14_FB,[2 3]),:,:)=[];G14_FB(:,~any(G14_FB,[1 3]),:)=[];
            G15_FB(~any(G15_FB,[2 3]),:,:)=[];G15_FB(:,~any(G15_FB,[1 3]),:)=[];
            G16_FB(~any(G16_FB,[2 3]),:,:)=[];G16_FB(:,~any(G16_FB,[1 3]),:)=[];
            G17_FB(~any(G17_FB,[2 3]),:,:)=[];G17_FB(:,~any(G17_FB,[1 3]),:)=[];
            G18_FB(~any(G18_FB,[2 3]),:,:)=[];G18_FB(:,~any(G18_FB,[1 3]),:)=[];
            G19_FB(~any(G19_FB,[2 3]),:,:)=[];G19_FB(:,~any(G19_FB,[1 3]),:)=[];
            G20_FB(~any(G20_FB,[2 3]),:,:)=[];G20_FB(:,~any(G20_FB,[1 3]),:)=[];
            G13_AB(~any(G13_AB,[2 3]),:,:)=[];G13_AB(:,~any(G13_AB,[1 3]),:)=[];
            G14_AB(~any(G14_AB,[2 3]),:,:)=[];G14_AB(:,~any(G14_AB,[1 3]),:)=[];
            G15_AB(~any(G15_AB,[2 3]),:,:)=[];G15_AB(:,~any(G15_AB,[1 3]),:)=[];
            G16_AB(~any(G16_AB,[2 3]),:,:)=[];G16_AB(:,~any(G16_AB,[1 3]),:)=[];
            G17_AB(~any(G17_AB,[2 3]),:,:)=[];G17_AB(:,~any(G17_AB,[1 3]),:)=[];
            G18_AB(~any(G18_AB,[2 3]),:,:)=[];G18_AB(:,~any(G18_AB,[1 3]),:)=[];
            G19_AB(~any(G19_AB,[2 3]),:,:)=[];G19_AB(:,~any(G19_AB,[1 3]),:)=[];
            G20_AB(~any(G20_AB,[2 3]),:,:)=[];G20_AB(:,~any(G20_AB,[1 3]),:)=[];
            G7(G7==0)=NaN;
            G8(G8==0)=NaN;
            G9(G9==0)=NaN;
            G10(G10==0)=NaN;
            G11(G11==0)=NaN;
            G12(G12==0)=NaN;
            G13(G13==0)=NaN;
            G14(G14==0)=NaN;
            G15(G15==0)=NaN;
            G16(G16==0)=NaN;
            G17(G17==0)=NaN;
            G18(G18==0)=NaN;
            G19(G19==0)=NaN;
            G20(G20==0)=NaN;
            G13_FB(G13_FB==0)=NaN;
            G14_FB(G14_FB==0)=NaN;
            G15_FB(G15_FB==0)=NaN;
            G16_FB(G16_FB==0)=NaN;
            G17_FB(G17_FB==0)=NaN;
            G18_FB(G18_FB==0)=NaN;
            G19_FB(G19_FB==0)=NaN;
            G20_FB(G20_FB==0)=NaN;
            G13_AB(G13_AB==0)=NaN;
            G14_AB(G14_AB==0)=NaN;
            G15_AB(G15_AB==0)=NaN;
            G16_AB(G16_AB==0)=NaN;
            G17_AB(G17_AB==0)=NaN;
            G18_AB(G18_AB==0)=NaN;
            G19_AB(G19_AB==0)=NaN;
            G20_AB(G20_AB==0)=NaN;
            G1_Mean = mean(G1,1);
            G2_Mean = mean(G2,1);
            G3_Mean = mean(G3,1);
            G4_Mean = mean(G4,1);
            G5_Mean = mean(G5,1);
            G6_Mean = mean(G6,1);
            G7_Mean = reshape(mean(mean(G7,'omitnan'),'omitnan'),[],1)';
            G8_Mean = reshape(mean(mean(G8,'omitnan'),'omitnan'),[],1)';
            G9_Mean = reshape(mean(mean(G9,'omitnan'),'omitnan'),[],1)';
            G10_Mean = reshape(mean(mean(G10,'omitnan'),'omitnan'),[],1)';
            G11_Mean = reshape(mean(mean(G11,'omitnan'),'omitnan'),[],1)';
            G12_Mean = reshape(mean(mean(G12,'omitnan'),'omitnan'),[],1)';
            G13_Mean = reshape(mean(mean(G13,'omitnan'),'omitnan'),[],1)';
            G14_Mean = reshape(mean(mean(G14,'omitnan'),'omitnan'),[],1)';
            G15_Mean = reshape(mean(mean(G15,'omitnan'),'omitnan'),[],1)';
            G16_Mean = reshape(mean(mean(G16,'omitnan'),'omitnan'),[],1)';
            G17_Mean = reshape(mean(mean(G17,'omitnan'),'omitnan'),[],1)';
            G18_Mean = reshape(mean(mean(G18,'omitnan'),'omitnan'),[],1)';
            G19_Mean = reshape(mean(mean(G19,'omitnan'),'omitnan'),[],1)';
            G20_Mean = reshape(mean(mean(G20,'omitnan'),'omitnan'),[],1)';
            G13FB_Mean = reshape(mean(mean(G13_FB,'omitnan'),'omitnan'),[],1)';
            G14FB_Mean = reshape(mean(mean(G14_FB,'omitnan'),'omitnan'),[],1)';
            G15FB_Mean = reshape(mean(mean(G15_FB,'omitnan'),'omitnan'),[],1)';
            G16FB_Mean = reshape(mean(mean(G16_FB,'omitnan'),'omitnan'),[],1)';
            G17FB_Mean = reshape(mean(mean(G17_FB,'omitnan'),'omitnan'),[],1)';
            G18FB_Mean = reshape(mean(mean(G18_FB,'omitnan'),'omitnan'),[],1)';
            G19FB_Mean = reshape(mean(mean(G19_FB,'omitnan'),'omitnan'),[],1)';
            G20FB_Mean = reshape(mean(mean(G20_FB,'omitnan'),'omitnan'),[],1)';
            G13AB_Mean = reshape(mean(mean(G13_AB,'omitnan'),'omitnan'),[],1)';
            G14AB_Mean = reshape(mean(mean(G14_AB,'omitnan'),'omitnan'),[],1)';
            G15AB_Mean = reshape(mean(mean(G15_AB,'omitnan'),'omitnan'),[],1)';
            G16AB_Mean = reshape(mean(mean(G16_AB,'omitnan'),'omitnan'),[],1)';
            G17AB_Mean = reshape(mean(mean(G17_AB,'omitnan'),'omitnan'),[],1)';
            G18AB_Mean = reshape(mean(mean(G18_AB,'omitnan'),'omitnan'),[],1)';
            G19AB_Mean = reshape(mean(mean(G19_AB,'omitnan'),'omitnan'),[],1)';
            G20AB_Mean = reshape(mean(mean(G20_AB,'omitnan'),'omitnan'),[],1)';
            G1_SEM = std(G1,[],1)/sqrt(size(G1,1));
            G2_SEM = std(G2,[],1)/sqrt(size(G2,1));
            G3_SEM = std(G3,[],1)/sqrt(size(G3,1));
            G4_SEM = std(G4,[],1)/sqrt(size(G4,1));
            G5_SEM = std(G5,[],1)/sqrt(size(G5,1));
            G6_SEM = std(G6,[],1)/sqrt(size(G6,1));
            G7_SEM = (reshape(std(std(G7,'omitnan'),'omitnan'),[],1)/sqrt(numel(G7(~isnan(G7)))))';
            G8_SEM = (reshape(std(std(G8,'omitnan'),'omitnan'),[],1)/sqrt(numel(G8(~isnan(G8)))))';
            G9_SEM = (reshape(std(std(G9,'omitnan'),'omitnan'),[],1)/sqrt(numel(G9(~isnan(G9)))))';
            G10_SEM = (reshape(std(std(G10,'omitnan'),'omitnan'),[],1)/sqrt(numel(G10(~isnan(G10)))))';
            G11_SEM = (reshape(std(std(G11,'omitnan'),'omitnan'),[],1)/sqrt(numel(G11(~isnan(G11)))))';
            G12_SEM = (reshape(std(std(G12,'omitnan'),'omitnan'),[],1)/sqrt(numel(G12(~isnan(G12)))))';
            G13_SEM = (reshape(std(std(G13,'omitnan'),'omitnan'),[],1)/sqrt(numel(G13(~isnan(G13)))))';
            G14_SEM = (reshape(std(std(G14,'omitnan'),'omitnan'),[],1)/sqrt(numel(G14(~isnan(G14)))))';
            G15_SEM = (reshape(std(std(G15,'omitnan'),'omitnan'),[],1)/sqrt(numel(G15(~isnan(G15)))))';
            G16_SEM = (reshape(std(std(G16,'omitnan'),'omitnan'),[],1)/sqrt(numel(G16(~isnan(G16)))))';
            G17_SEM = (reshape(std(std(G17,'omitnan'),'omitnan'),[],1)/sqrt(numel(G17(~isnan(G17)))))';
            G18_SEM = (reshape(std(std(G18,'omitnan'),'omitnan'),[],1)/sqrt(numel(G18(~isnan(G18)))))';
            G19_SEM = (reshape(std(std(G19,'omitnan'),'omitnan'),[],1)/sqrt(numel(G19(~isnan(G19)))))';
            G20_SEM = (reshape(std(std(G20,'omitnan'),'omitnan'),[],1)/sqrt(numel(G20(~isnan(G20)))))';
            G13FB_SEM = (reshape(std(std(G13_FB,'omitnan'),'omitnan'),[],1)/sqrt(numel(G13_FB(~isnan(G13_FB)))))';
            G14FB_SEM = (reshape(std(std(G14_FB,'omitnan'),'omitnan'),[],1)/sqrt(numel(G14_FB(~isnan(G14_FB)))))';
            G15FB_SEM = (reshape(std(std(G15_FB,'omitnan'),'omitnan'),[],1)/sqrt(numel(G15_FB(~isnan(G15_FB)))))';
            G16FB_SEM = (reshape(std(std(G16_FB,'omitnan'),'omitnan'),[],1)/sqrt(numel(G16_FB(~isnan(G16_FB)))))';
            G17FB_SEM = (reshape(std(std(G17_FB,'omitnan'),'omitnan'),[],1)/sqrt(numel(G17_FB(~isnan(G17_FB)))))';
            G18FB_SEM = (reshape(std(std(G18_FB,'omitnan'),'omitnan'),[],1)/sqrt(numel(G18_FB(~isnan(G18_FB)))))';
            G19FB_SEM = (reshape(std(std(G19_FB,'omitnan'),'omitnan'),[],1)/sqrt(numel(G19_FB(~isnan(G19_FB)))))';
            G20FB_SEM = (reshape(std(std(G20_FB,'omitnan'),'omitnan'),[],1)/sqrt(numel(G20_FB(~isnan(G20_FB)))))';
            G13AB_SEM = (reshape(std(std(G13_AB,'omitnan'),'omitnan'),[],1)/sqrt(numel(G13_AB(~isnan(G13_AB)))))';
            G14AB_SEM = (reshape(std(std(G14_AB,'omitnan'),'omitnan'),[],1)/sqrt(numel(G14_AB(~isnan(G14_AB)))))';
            G15AB_SEM = (reshape(std(std(G15_AB,'omitnan'),'omitnan'),[],1)/sqrt(numel(G15_AB(~isnan(G15_AB)))))';
            G16AB_SEM = (reshape(std(std(G16_AB,'omitnan'),'omitnan'),[],1)/sqrt(numel(G16_AB(~isnan(G16_AB)))))';
            G17AB_SEM = (reshape(std(std(G17_AB,'omitnan'),'omitnan'),[],1)/sqrt(numel(G17_AB(~isnan(G17_AB)))))';
            G18AB_SEM = (reshape(std(std(G18_AB,'omitnan'),'omitnan'),[],1)/sqrt(numel(G18_AB(~isnan(G18_AB)))))';
            G19AB_SEM = (reshape(std(std(G19_AB,'omitnan'),'omitnan'),[],1)/sqrt(numel(G19_AB(~isnan(G19_AB)))))';
            G20AB_SEM = (reshape(std(std(G20_AB,'omitnan'),'omitnan'),[],1)/sqrt(numel(G20_AB(~isnan(G20_AB)))))';
            plot(ax1,t_Diam(1:size(G1_Mean,2)),G1_Mean,"Color",'r',"LineWidth",2)
            plot(ax2,t_Diam(1:size(G2_Mean,2)),G2_Mean,"Color",'r',"LineWidth",2)
            plot(ax3,linspace(WP(1),WP(2),size(G3,2)),G3_Mean,"Color",'r',"LineWidth",2)
            plot(ax4,linspace(WP(1),WP(2),size(G4,2)),G4_Mean,"Color",'r',"LineWidth",2)
            plot(ax5,linspace(WP(1),WP(2),size(G5,2)),G5_Mean,"Color",'r',"LineWidth",2)
            plot(ax6,linspace(WP(1),WP(2),size(G6,2)),G6_Mean,"Color",'r',"LineWidth",2)
            plot(ax7,linspace(WP(1),WP(2),size(G7,3)),G7_Mean,"Color",'r',"LineWidth",2)
            plot(ax8,linspace(WP(1),WP(2),size(G8,3)),G8_Mean,"Color",'r',"LineWidth",2)
            plot(ax9,linspace(WP(1),WP(2),size(G9,3)),G9_Mean,"Color",'r',"LineWidth",2)
            plot(ax10,linspace(WP(1),WP(2),size(G10,3)),G10_Mean,"Color",'r',"LineWidth",2)
            plot(ax11,linspace(WP(1),WP(2),size(G11,3)),G11_Mean,"Color",'r',"LineWidth",2)
            plot(ax12,linspace(WP(1),WP(2),size(G12,3)),G12_Mean,"Color",'r',"LineWidth",2)
            plot(ax16,linspace(WP(1),WP(2),size(G13,3)),G13_Mean,"Color",[204 37 41]./255,"LineWidth",2)
            plot(ax16,linspace(WP(1),WP(2),size(G15,3)),G15_Mean,"Color",[62 150 81]./255,"LineWidth",2)
            plot(ax16,linspace(WP(1),WP(2),size(G17,3)),G17_Mean,"Color",[146 36 40]./255,"LineWidth",2)
            plot(ax16,linspace(WP(1),WP(2),size(G19,3)),G19_Mean,"Color",[107 76 154]./255,"LineWidth",2)
            plot(ax17,linspace(WP(1),WP(2),size(G14,3)),G14_Mean,"Color",[204 37 41]./255,"LineWidth",2)
            plot(ax17,linspace(WP(1),WP(2),size(G16,3)),G16_Mean,"Color",[62 150 81]./255,"LineWidth",2)
            plot(ax17,linspace(WP(1),WP(2),size(G18,3)),G18_Mean,"Color",[146 36 40]./255,"LineWidth",2)
            plot(ax17,linspace(WP(1),WP(2),size(G20,3)),G20_Mean,"Color",[107 76 154]./255,"LineWidth",2)
            plot(ax18,linspace(WP(1),WP(2),size(G13_FB,3)),G13FB_Mean,"Color",[204 37 41]./255,"LineWidth",2)
            plot(ax18,linspace(WP(1),WP(2),size(G15_FB,3)),G15FB_Mean,"Color",[62 150 81]./255,"LineWidth",2)
            plot(ax18,linspace(WP(1),WP(2),size(G17_FB,3)),G17FB_Mean,"Color",[146 36 40]./255,"LineWidth",2)
            plot(ax18,linspace(WP(1),WP(2),size(G19_FB,3)),G19FB_Mean,"Color",[107 76 154]./255,"LineWidth",2)
            plot(ax19,linspace(WP(1),WP(2),size(G14_FB,3)),G14FB_Mean,"Color",[204 37 41]./255,"LineWidth",2)
            plot(ax19,linspace(WP(1),WP(2),size(G16_FB,3)),G16FB_Mean,"Color",[62 150 81]./255,"LineWidth",2)
            plot(ax19,linspace(WP(1),WP(2),size(G18_FB,3)),G18FB_Mean,"Color",[146 36 40]./255,"LineWidth",2)
            plot(ax19,linspace(WP(1),WP(2),size(G20_FB,3)),G20FB_Mean,"Color",[107 76 154]./255,"LineWidth",2)
            plot(ax20,linspace(WP(1),WP(2),size(G13_AB,3)),G13AB_Mean,"Color",[204 37 41]./255,"LineWidth",2)
            plot(ax20,linspace(WP(1),WP(2),size(G15_AB,3)),G15AB_Mean,"Color",[62 150 81]./255,"LineWidth",2)
            plot(ax20,linspace(WP(1),WP(2),size(G17_AB,3)),G17AB_Mean,"Color",[146 36 40]./255,"LineWidth",2)
            plot(ax20,linspace(WP(1),WP(2),size(G19_AB,3)),G19AB_Mean,"Color",[107 76 154]./255,"LineWidth",2)
            plot(ax21,linspace(WP(1),WP(2),size(G14_AB,3)),G14AB_Mean,"Color",[204 37 41]./255,"LineWidth",2)
            plot(ax21,linspace(WP(1),WP(2),size(G16_AB,3)),G16AB_Mean,"Color",[62 150 81]./255,"LineWidth",2)
            plot(ax21,linspace(WP(1),WP(2),size(G18_AB,3)),G18AB_Mean,"Color",[146 36 40]./255,"LineWidth",2)
            plot(ax21,linspace(WP(1),WP(2),size(G20_AB,3)),G20AB_Mean,"Color",[107 76 154]./255,"LineWidth",2)
            fill(ax1,[t_Diam(1:size(G1_Mean,2)), flipud(t_Diam(1:size(G1_Mean,2))')'],[(G1_Mean+G1_SEM), flipud((G1_Mean-G1_SEM)')'],[154 0 0]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax2,[t_Diam(1:size(G2_Mean,2)), flipud(t_Diam(1:size(G2_Mean,2))')'],[(G2_Mean+G2_SEM), flipud((G2_Mean-G2_SEM)')'],[154 0 0]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax3,[linspace(WP(1),WP(2),size(G3,2)), flipud(linspace(WP(1),WP(2),size(G3,2))')'],[(G3_Mean+G3_SEM), flipud((G3_Mean-G3_SEM)')'],[154 0 0]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax4,[linspace(WP(1),WP(2),size(G4,2)), flipud(linspace(WP(1),WP(2),size(G4,2))')'],[(G4_Mean+G4_SEM), flipud((G4_Mean-G4_SEM)')'],[154 0 0]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax5,[linspace(WP(1),WP(2),size(G5,2)), flipud(linspace(WP(1),WP(2),size(G5,2))')'],[(G5_Mean+G5_SEM), flipud((G5_Mean-G5_SEM)')'],[154 0 0]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax6,[linspace(WP(1),WP(2),size(G6,2)), flipud(linspace(WP(1),WP(2),size(G6,2))')'],[(G6_Mean+G6_SEM), flipud((G6_Mean-G6_SEM)')'],[154 0 0]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax7,[linspace(WP(1),WP(2),size(G7,3)), flipud(linspace(WP(1),WP(2),size(G7,3))')'],[(G7_Mean+G7_SEM), flipud((G7_Mean-G7_SEM)')'],[154 0 0]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax8,[linspace(WP(1),WP(2),size(G8,3)), flipud(linspace(WP(1),WP(2),size(G8,3))')'],[(G8_Mean+G8_SEM), flipud((G8_Mean-G8_SEM)')'],[154 0 0]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax9,[linspace(WP(1),WP(2),size(G9,3)), flipud(linspace(WP(1),WP(2),size(G9,3))')'],[(G9_Mean+G9_SEM), flipud((G9_Mean-G9_SEM)')'],[154 0 0]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax10,[linspace(WP(1),WP(2),size(G10,3)), flipud(linspace(WP(1),WP(2),size(G10,3))')'],[(G10_Mean+G10_SEM), flipud((G10_Mean-G10_SEM)')'],[154 0 0]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax11,[linspace(WP(1),WP(2),size(G11,3)), flipud(linspace(WP(1),WP(2),size(G11,3))')'],[(G11_Mean+G11_SEM), flipud((G11_Mean-G11_SEM)')'],[154 0 0]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax12,[linspace(WP(1),WP(2),size(G12,3)), flipud(linspace(WP(1),WP(2),size(G12,3))')'],[(G12_Mean+G12_SEM), flipud((G12_Mean-G12_SEM)')'],[154 0 0]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax16,[linspace(WP(1),WP(2),size(G13,3)), flipud(linspace(WP(1),WP(2),size(G13,3))')'],[(G13_Mean+G13_SEM), flipud((G13_Mean-G13_SEM)')'],[211 94 96]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax16,[linspace(WP(1),WP(2),size(G15,3)), flipud(linspace(WP(1),WP(2),size(G15,3))')'],[(G15_Mean+G15_SEM), flipud((G15_Mean-G15_SEM)')'],[132 186 91]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax16,[linspace(WP(1),WP(2),size(G17,3)), flipud(linspace(WP(1),WP(2),size(G17,3))')'],[(G17_Mean+G17_SEM), flipud((G17_Mean-G17_SEM)')'],[171 104 87]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax16,[linspace(WP(1),WP(2),size(G19,3)), flipud(linspace(WP(1),WP(2),size(G19,3))')'],[(G19_Mean+G19_SEM), flipud((G19_Mean-G19_SEM)')'],[144 103 167]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax17,[linspace(WP(1),WP(2),size(G14,3)), flipud(linspace(WP(1),WP(2),size(G14,3))')'],[(G14_Mean+G14_SEM), flipud((G14_Mean-G14_SEM)')'],[211 94 96]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax17,[linspace(WP(1),WP(2),size(G16,3)), flipud(linspace(WP(1),WP(2),size(G16,3))')'],[(G16_Mean+G16_SEM), flipud((G16_Mean-G16_SEM)')'],[132 186 91]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax17,[linspace(WP(1),WP(2),size(G18,3)), flipud(linspace(WP(1),WP(2),size(G18,3))')'],[(G18_Mean+G18_SEM), flipud((G18_Mean-G18_SEM)')'],[171 104 87]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax17,[linspace(WP(1),WP(2),size(G20,3)), flipud(linspace(WP(1),WP(2),size(G20,3))')'],[(G20_Mean+G20_SEM), flipud((G20_Mean-G20_SEM)')'],[144 103 167]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            lgd17=legend(ax17,'Quiet','SHL','N60','N70','Location','southeastoutside');
            lgd17.Title.String = 'Conditions:';
            fill(ax18,[linspace(WP(1),WP(2),size(G13_FB,3)), flipud(linspace(WP(1),WP(2),size(G13_FB,3))')'],[(G13FB_Mean+G13FB_SEM), flipud((G13FB_Mean-G13FB_SEM)')'],[211 94 96]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax18,[linspace(WP(1),WP(2),size(G15_FB,3)), flipud(linspace(WP(1),WP(2),size(G15_FB,3))')'],[(G15FB_Mean+G15FB_SEM), flipud((G15FB_Mean-G15FB_SEM)')'],[132 186 91]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax18,[linspace(WP(1),WP(2),size(G17_FB,3)), flipud(linspace(WP(1),WP(2),size(G17_FB,3))')'],[(G17FB_Mean+G17FB_SEM), flipud((G17FB_Mean-G17FB_SEM)')'],[171 104 87]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax18,[linspace(WP(1),WP(2),size(G19_FB,3)), flipud(linspace(WP(1),WP(2),size(G19_FB,3))')'],[(G19FB_Mean+G19FB_SEM), flipud((G19FB_Mean-G19FB_SEM)')'],[144 103 167]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax19,[linspace(WP(1),WP(2),size(G14_FB,3)), flipud(linspace(WP(1),WP(2),size(G14_FB,3))')'],[(G14FB_Mean+G14FB_SEM), flipud((G14FB_Mean-G14FB_SEM)')'],[211 94 96]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax19,[linspace(WP(1),WP(2),size(G16_FB,3)), flipud(linspace(WP(1),WP(2),size(G16_FB,3))')'],[(G16FB_Mean+G16FB_SEM), flipud((G16FB_Mean-G16FB_SEM)')'],[132 186 91]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax19,[linspace(WP(1),WP(2),size(G18_FB,3)), flipud(linspace(WP(1),WP(2),size(G18_FB,3))')'],[(G18FB_Mean+G18FB_SEM), flipud((G18FB_Mean-G18FB_SEM)')'],[171 104 87]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax19,[linspace(WP(1),WP(2),size(G20_FB,3)), flipud(linspace(WP(1),WP(2),size(G20_FB,3))')'],[(G20FB_Mean+G20FB_SEM), flipud((G20FB_Mean-G20FB_SEM)')'],[144 103 167]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            lgd19=legend(ax19,'Quiet','SHL','N60','N70','Location','southeastoutside');
            lgd19.Title.String = 'Conditions:';
            fill(ax20,[linspace(WP(1),WP(2),size(G13_AB,3)), flipud(linspace(WP(1),WP(2),size(G13_AB,3))')'],[(G13AB_Mean+G13AB_SEM), flipud((G13AB_Mean-G13AB_SEM)')'],[211 94 96]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax20,[linspace(WP(1),WP(2),size(G15_AB,3)), flipud(linspace(WP(1),WP(2),size(G15_AB,3))')'],[(G15AB_Mean+G15AB_SEM), flipud((G15AB_Mean-G15AB_SEM)')'],[132 186 91]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax20,[linspace(WP(1),WP(2),size(G17_AB,3)), flipud(linspace(WP(1),WP(2),size(G17_AB,3))')'],[(G17AB_Mean+G17AB_SEM), flipud((G17AB_Mean-G17AB_SEM)')'],[171 104 87]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax20,[linspace(WP(1),WP(2),size(G19_AB,3)), flipud(linspace(WP(1),WP(2),size(G19_AB,3))')'],[(G19AB_Mean+G19AB_SEM), flipud((G19AB_Mean-G19AB_SEM)')'],[144 103 167]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax21,[linspace(WP(1),WP(2),size(G14_AB,3)), flipud(linspace(WP(1),WP(2),size(G14_AB,3))')'],[(G14AB_Mean+G14AB_SEM), flipud((G14AB_Mean-G14AB_SEM)')'],[211 94 96]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax21,[linspace(WP(1),WP(2),size(G16_AB,3)), flipud(linspace(WP(1),WP(2),size(G16_AB,3))')'],[(G16AB_Mean+G16AB_SEM), flipud((G16AB_Mean-G16AB_SEM)')'],[132 186 91]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax21,[linspace(WP(1),WP(2),size(G18_AB,3)), flipud(linspace(WP(1),WP(2),size(G18_AB,3))')'],[(G18AB_Mean+G18AB_SEM), flipud((G18AB_Mean-G18AB_SEM)')'],[171 104 87]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax21,[linspace(WP(1),WP(2),size(G20_AB,3)), flipud(linspace(WP(1),WP(2),size(G20_AB,3))')'],[(G20AB_Mean+G20AB_SEM), flipud((G20AB_Mean-G20AB_SEM)')'],[144 103 167]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            lgd21=legend(ax21,'Quiet','SHL','N60','N70','Location','southeastoutside');
            lgd21.Title.String = 'Conditions:';
            % Plot features
            data13=[mean(nonzeros(F1_S_Q)),mean(nonzeros(F1_L_Q));mean(nonzeros(F1_S_SHL)),mean(nonzeros(F1_L_SHL));mean(nonzeros(F1_S_N60)),mean(nonzeros(F1_L_N60));mean(nonzeros(F1_S_N70)),mean(nonzeros(F1_L_N70))];
            b13=bar(ax13,data13,'grouped');
            [ngroups,nbars] = size(data13);
            x13 = nan(nbars, ngroups);
            for e = 1:nbars
                x13(e,:) = b13(e).XEndPoints;
            end
            errorbar(ax13,x13',data13,[std(nonzeros(F1_S_Q)),std(nonzeros(F1_L_Q));std(nonzeros(F1_S_SHL)),std(nonzeros(F1_L_SHL));std(nonzeros(F1_S_N60)),std(nonzeros(F1_L_N60));std(nonzeros(F1_S_N70)),std(nonzeros(F1_L_N70))],'k','linestyle','none','handlevisibility' ,'off')
            data14=[mean(nonzeros(F2_S_Q)),mean(nonzeros(F2_L_Q));mean(nonzeros(F2_S_SHL)),mean(nonzeros(F2_L_SHL));mean(nonzeros(F2_S_N60)),mean(nonzeros(F2_L_N60));mean(nonzeros(F2_S_N70)),mean(nonzeros(F2_L_N70))];
            b14=bar(ax14,data14,'grouped');
            [ngroups,nbars] = size(data14);
            x14 = nan(nbars, ngroups);
            for e = 1:nbars
                x14(e,:) = b14(e).XEndPoints;
            end
            errorbar(ax14,x14',data14,[std(nonzeros(F2_S_Q)),std(nonzeros(F2_L_Q));std(nonzeros(F2_S_SHL)),std(nonzeros(F2_L_SHL));std(nonzeros(F2_S_N60)),std(nonzeros(F2_L_N60));std(nonzeros(F2_S_N70)),std(nonzeros(F2_L_N70))],'k','linestyle','none','handlevisibility' ,'off')
            data15=[mean(nonzeros(F3_S_Q)),mean(nonzeros(F3_L_Q));mean(nonzeros(F3_S_SHL)),mean(nonzeros(F3_L_SHL));mean(nonzeros(F3_S_N60)),mean(nonzeros(F3_L_N60));mean(nonzeros(F3_S_N70)),mean(nonzeros(F3_L_N70))];
            b15=bar(ax15,data15,'grouped');
            [ngroups,nbars] = size(data15);
            x15 = nan(nbars, ngroups);
            for e = 1:nbars
                x15(e,:) = b15(e).XEndPoints;
            end
            errorbar(ax15,x15',data15,[std(nonzeros(F3_S_Q)),std(nonzeros(F3_L_Q));std(nonzeros(F3_S_SHL)),std(nonzeros(F3_L_SHL));std(nonzeros(F3_S_N60)),std(nonzeros(F3_L_N60));std(nonzeros(F3_S_N70)),std(nonzeros(F3_L_N70))],'k','linestyle','none','handlevisibility' ,'off')
            lgd15=legend(ax15,'Speaking','Non-Speaking','Location','southeastoutside');
            lgd15.Title.String = 'Types of Windows:';
            xticks([ax13 ax14 ax15],1:4)
            xticklabels([ax13 ax14 ax15],{'Quiet','SHL','N60','N70'})
% xlim([ax2 ax1],[9.5, 10.5]);xlim([ax3 ax4],[-0.5, 0.5]);ylim([ax1 ax2 ax3 ax4],[3.4, 3.7])
        end
    end
end
