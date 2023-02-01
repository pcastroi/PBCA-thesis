%% PBCA-Thesis - Week 4 - Fixed or adaptive baseline?
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
AudFs=48000;
WP=[-3,6]; % [s]
BL=[15,20]; % [s]
BLPeriod=[0,20]; % [s]
BLStartEnd=BLPeriod*Param.Fs+1; % [samples]
TimeMinWin = 0.5; % [s], Minimum time of a window
TimeMergeGap = 0.3; % [s], Time threshold for merging windows
RejectRatio = 0.3; % Rejection threshold based on the ratio of NaNs in data
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
            disp(['Warning: No associated Delay data for file ', PairFiles(i).folder, '\', PairFiles(i).name, '.']);
            continue
        end
        
        if or(SDelayRaw < 0,LDelayRaw < 0)
            SDelayRaw=[0,0];
            LDelayRaw=[0,0];
        end
        
        if isempty(SpeakRaw) && isempty(ListenRaw)
            disp(['Warning: No associated Utterance/Listening data for file ', PairFiles(i).folder, '\', PairFiles(i).name, '.']);
            continue
        end
        
        SpeakRaw(:,2:3)=round((SpeakRaw(:,2:3)*binResUtt+BLPeriod(2))*Param.Fs+SDelayRaw(1)/2);
        ListenRaw(:,2:3)=round((ListenRaw(:,2:3)*binResUtt+BLPeriod(2))*Param.Fs+LDelayRaw(1)/2);
        
        % Merge windows if gap <= TimeMergeGap
        [SpeakM] = MergeWin(SpeakRaw, Param.Fs, TimeMergeGap);
        [ListenM] = MergeWin(ListenRaw, Param.Fs, TimeMergeGap);
                
        % Discard windows if duration is < TimeMinWin
        Speak = SpeakM(SpeakM(:,1)>TimeMinWin,:);
        Listen = ListenM(ListenM(:,1)>TimeMinWin,:);
        
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
            
        end
        x=x+1;
        hold([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10 ax11 ax12 ax13 ax14 ax15],'on')
        if Speak(1,2) < Listen(1,2) % Initially Speaking
            G1(x,:)=Diameter(BLStartEnd(1):BLStartEnd(2));
            G3(x,:)=Diameter(Speak(1,2)+WP(1)*Param.Fs:Speak(1,2)+WP(2)*Param.Fs);
            G5(x,:)=Diameter(Speak(1,2)+WP(1)*Param.Fs:Speak(1,2)+WP(2)*Param.Fs)-mean(Diameter(BL(1)*Param.Fs:BL(2)*Param.Fs));
            plot(ax1,t_Diam(BLStartEnd(1):BLStartEnd(2)),G1(x,:),color=[0 0 0 0.2],linewidth=0.5)
            plot(ax3,linspace(WP(1),WP(2),size(G3,2)),G3(x,:),color=[0 0 0 0.2],linewidth=0.5)
            plot(ax5,linspace(WP(1),WP(2),size(G5,2)),G5(x,:),color=[0 0 0 0.2],linewidth=0.5)
        elseif Speak(1,2) > Listen(1,2) % Initially Listening/Non-speaking
            G2(x,:)=Diameter(BLStartEnd(1):BLStartEnd(2));
            G4(x,:)=Diameter(Listen(1,2)+WP(1)*Param.Fs:Listen(1,2)+WP(2)*Param.Fs);
            G6(x,:)=Diameter(Listen(1,2)+WP(1)*Param.Fs:Listen(1,2)+WP(2)*Param.Fs)-mean(Diameter(BL(1)*Param.Fs:BL(2)*Param.Fs));
            plot(ax2,t_Diam(BLStartEnd(1):BLStartEnd(2)),G2(x,:),color=[0 0 0 0.2],linewidth=0.5)
            plot(ax4,linspace(WP(1),WP(2),size(G4,2)),G4(x,:),color=[0 0 0 0.2],linewidth=0.5)
            plot(ax6,linspace(WP(1),WP(2),size(G6,2)),G6(x,:),color=[0 0 0 0.2],linewidth=0.5)
        else
            continue
        end
        
        for j=1:size(Speak,1)
            if Speak(j,2)+WP(2)*Param.Fs <= length(Diameter) % When reaching 'end of tobii data'
                G7(x,j,:)=Diameter(Speak(j,2)+WP(1)*Param.Fs:Speak(j,2)+WP(2)*Param.Fs);
                G9(x,j,:)=Diameter(Speak(j,2)+WP(1)*Param.Fs:Speak(j,2)+WP(2)*Param.Fs)-mean(Diameter(BL(1)*Param.Fs:BL(2)*Param.Fs));
                G11(x,j,:)=Diameter(Speak(j,2)+WP(1)*Param.Fs:Speak(j,2)+WP(2)*Param.Fs)-mean(Diameter(Speak(j,2)-Param.Fs:Speak(j,2)));
%                 plot(ax7,linspace(WP(1),WP(2),size(G7,3)),reshape(G7(x,j,:),[],1),color=[0 0 0 0.01],linewidth=0.5)
            end
        end
        
        for j=1:size(Listen,1)
            if Listen(j,2)+WP(2)*Param.Fs <= length(Diameter) % When reaching 'end of tobii data'
                G8(x,j,:)=Diameter(Listen(j,2)+WP(1)*Param.Fs:Listen(j,2)+WP(2)*Param.Fs);
                G10(x,j,:)=Diameter(Listen(j,2)+WP(1)*Param.Fs:Listen(j,2)+WP(2)*Param.Fs)-mean(Diameter(BL(1)*Param.Fs:BL(2)*Param.Fs));
                G12(x,j,:)=Diameter(Listen(j,2)+WP(1)*Param.Fs:Listen(j,2)+WP(2)*Param.Fs)-mean(Diameter(Listen(j,2)-Param.Fs:Listen(j,2)));
%                 plot(ax8,linspace(WP(1),WP(2),size(G8,3)),reshape(G8(x,j,:),[],1),color=[0 0 0 0.01],linewidth=0.5)
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
            end
            F1_S_Q(x)=mean(nonzeros(temp));
            F2_S_Q(x)=mean(nonzeros(temp2(:,2)));
            F3_S_Q(x)=mean(nonzeros(temp3));
            for j=1:size(Listen,1)
                temp(j)=mean(Diameter(Listen(j,2):Listen(j,3)));
                temp2(j,:)=polyfit(t_Diam(Listen(j,2):Listen(j,3)),Diameter(Listen(j,2):Listen(j,3)),1);
                temp3(j)=max(Diameter(Listen(j,2):Listen(j,3)));
            end
            F1_L_Q(x)=mean(nonzeros(temp));
            F2_L_Q(x)=mean(nonzeros(temp2(:,2))); % Check if this makes sense
            F3_L_Q(x)=mean(nonzeros(temp3));
        elseif contains(PairFiles(i).name,'SHL')
            for j=1:size(Speak,1)
                temp(j)=mean(Diameter(Speak(j,2):Speak(j,3)));
                temp2(j,:)=polyfit(t_Diam(Speak(j,2):Speak(j,3)),Diameter(Speak(j,2):Speak(j,3)),1);
                temp3(j)=max(Diameter(Speak(j,2):Speak(j,3)));
            end
            F1_S_SHL(x)=mean(nonzeros(temp));
            F2_S_SHL(x)=mean(nonzeros(temp2(:,2)));
            F3_S_SHL(x)=mean(nonzeros(temp3));
            for j=1:size(Listen,1)
                temp(j)=mean(Diameter(Listen(j,2):Listen(j,3)));
                temp2(j,:)=polyfit(t_Diam(Listen(j,2):Listen(j,3)),Diameter(Listen(j,2):Listen(j,3)),1);
                temp3(j)=max(Diameter(Listen(j,2):Listen(j,3)));
            end
            F1_L_SHL(x)=mean(nonzeros(temp));
            F2_L_SHL(x)=mean(nonzeros(temp2(:,2)));
            F3_L_SHL(x)=mean(nonzeros(temp3));
        elseif contains(PairFiles(i).name,'Noise60')
            for j=1:size(Speak,1)
                temp(j)=mean(Diameter(Speak(j,2):Speak(j,3)));
                temp2(j,:)=polyfit(t_Diam(Speak(j,2):Speak(j,3)),Diameter(Speak(j,2):Speak(j,3)),1);
                temp3(j)=max(Diameter(Speak(j,2):Speak(j,3)));
            end
            F1_S_N60(x)=mean(nonzeros(temp));
            F2_S_N60(x)=mean(nonzeros(temp2(:,2)));
            F3_S_N60(x)=mean(nonzeros(temp3));
            for j=1:size(Listen,1)
                temp(j)=mean(Diameter(Listen(j,2):Listen(j,3)));
                temp2(j,:)=polyfit(t_Diam(Listen(j,2):Listen(j,3)),Diameter(Listen(j,2):Listen(j,3)),1);
                temp3(j)=max(Diameter(Listen(j,2):Listen(j,3)));
            end
            F1_L_N60(x)=mean(nonzeros(temp));
            F2_L_N60(x)=mean(nonzeros(temp2(:,2)));
            F3_L_N60(x)=mean(nonzeros(temp3));
        elseif contains(PairFiles(i).name,'Noise70')
            for j=1:size(Speak,1)
                temp(j)=mean(Diameter(Speak(j,2):Speak(j,3)));
                temp2(j,:)=polyfit(t_Diam(Speak(j,2):Speak(j,3)),Diameter(Speak(j,2):Speak(j,3)),1);
                temp3(j)=max(Diameter(Speak(j,2):Speak(j,3)));
            end
            F1_S_N70(x)=mean(nonzeros(temp));
            F2_S_N70(x)=mean(nonzeros(temp2(:,2)));
            F3_S_N70(x)=mean(nonzeros(temp3));
            for j=1:size(Listen,1)
                temp(j)=mean(Diameter(Listen(j,2):Listen(j,3)));
                temp2(j,:)=polyfit(t_Diam(Listen(j,2):Listen(j,3)),Diameter(Listen(j,2):Listen(j,3)),1);
                temp3(j)=max(Diameter(Listen(j,2):Listen(j,3)));
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
            xlabel([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10 ax11 ax12],'Time [s]')
            ylabel([ax1 ax2 ax3 ax4 ax7 ax8],'Pupil diameter [mm]')
            ylabel([ax5 ax6 ax9 ax10 ax11 ax12],'Pupil baseline difference [mm]')
            xlim([ax1 ax2],[0, 20])
            xlim([ax3 ax4 ax5 ax6],[-3, 6])
            grid([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10 ax11 ax12],'on')
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
            G7(G7==0)=NaN;
            G8(G8==0)=NaN;
            G9(G9==0)=NaN;
            G10(G10==0)=NaN;
            G11(G11==0)=NaN;
            G12(G12==0)=NaN;
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
            % Plot features
            bar(ax13,[mean(nonzeros(F1_S_Q)),mean(nonzeros(F1_L_Q));mean(nonzeros(F1_S_SHL)),mean(nonzeros(F1_L_SHL));mean(nonzeros(F1_S_N60)),mean(nonzeros(F1_L_N60));mean(nonzeros(F1_S_N70)),mean(nonzeros(F1_L_N70))])
%             errorbar(ax13,[std(nonzeros(F1_S_Q)),std(nonzeros(F1_L_Q));std(nonzeros(F1_S_SHL)),std(nonzeros(F1_L_SHL));std(nonzeros(F1_S_N60)),std(nonzeros(F1_L_N60));std(nonzeros(F1_S_N70)),std(nonzeros(F1_L_N70))])
            bar(ax14,[mean(nonzeros(F2_S_Q)),mean(nonzeros(F2_L_Q));mean(nonzeros(F2_S_SHL)),mean(nonzeros(F2_L_SHL));mean(nonzeros(F2_S_N60)),mean(nonzeros(F2_L_N60));mean(nonzeros(F2_S_N70)),mean(nonzeros(F2_L_N70))])
            
            bar(ax15,[mean(nonzeros(F3_S_Q)),mean(nonzeros(F3_L_Q));mean(nonzeros(F3_S_SHL)),mean(nonzeros(F3_L_SHL));mean(nonzeros(F3_S_N60)),mean(nonzeros(F3_L_N60));mean(nonzeros(F3_S_N70)),mean(nonzeros(F3_L_N70))])
            
% xlim([ax2 ax1],[9.5, 10.5]);xlim([ax3 ax4],[-0.5, 0.5]);ylim([ax1 ax2 ax3 ax4],[3.4, 3.7])
        end
    end
end
