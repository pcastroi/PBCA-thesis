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
BLPeriod=[0,20]; % [s]
BLStartEnd=BLPeriod*Param.Fs+1; % [samples]
G1 = zeros(200,BLStartEnd(2));
G2 = G1;
G3 = G1;
G4 = G1;
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
        
        
        if isempty(SpeakRaw) && isempty(ListenRaw)
            disp(['Warning: No associated Utterance/Listening data for file ', PairFiles(i).folder, '\', PairFiles(i).name, '.']);
            continue
        end
        
        SpeakRaw(:,2:3)=round((SpeakRaw(:,2:3)*binResUtt+BLPeriod(2))*Param.Fs+SDelayRaw(1)/2);
        ListenRaw(:,2:3)=round((ListenRaw(:,2:3)*binResUtt+BLPeriod(2))*Param.Fs+LDelayRaw(1)/2);
        
        t_Diam = linspace(0,length(Diameter)./Param.Fs,length(Diameter));
        
        % Plots
        if x==1
            tiledlayout(2,2);
            ax1 = nexttile;
            ax2 = nexttile;
            ax3 = nexttile;
            ax4 = nexttile;
        end
        x=x+1;
        hold([ax1 ax2 ax3 ax4],'on')
        if SpeakRaw(1,2) > ListenRaw(1,2)
            plot(ax1,t_Diam(BLStartEnd(1):BLStartEnd(2)),Diameter(BLStartEnd(1):BLStartEnd(2)),color=[0 0 0 0.2],linewidth=0.5)
            plot(ax3,linspace(-10,10,size(G3,2)),Diameter(SpeakRaw(1,2)-BLPeriod(2)/2*Param.Fs:SpeakRaw(1,2)+BLPeriod(2)/2*Param.Fs),color=[0 0 0 0.2],linewidth=0.5)
            G1(x,:)=Diameter(BLStartEnd(1):BLStartEnd(2));
            G3(x,:)=Diameter(SpeakRaw(1,2)-BLPeriod(2)/2*Param.Fs:SpeakRaw(1,2)+BLPeriod(2)/2*Param.Fs);
        elseif SpeakRaw(1,2) < ListenRaw(1,2)
            plot(ax2,t_Diam(BLStartEnd(1):BLStartEnd(2)),Diameter(BLStartEnd(1):BLStartEnd(2)),color=[0 0 0 0.2],linewidth=0.5)
            plot(ax4,linspace(-10,10,size(G4,2)),Diameter(ListenRaw(1,2)-BLPeriod(2)/2*Param.Fs:ListenRaw(1,2)+BLPeriod(2)/2*Param.Fs),color=[0 0 0 0.2],linewidth=0.5)
            G2(x,:)=Diameter(BLStartEnd(1):BLStartEnd(2));
            G4(x,:)=Diameter(ListenRaw(1,2)-BLPeriod(2)/2*Param.Fs:ListenRaw(1,2)+BLPeriod(2)/2*Param.Fs);
        else
            plot(ax1,t_Diam(BLStartEnd(1):BLStartEnd(2)),Diameter(BLStartEnd(1):BLStartEnd(2)),color=[0 0 0 0.2],linewidth=0.5)
            plot(ax2,t_Diam(BLStartEnd(1):BLStartEnd(2)),Diameter(BLStartEnd(1):BLStartEnd(2)),color=[0 0 0 0.2],linewidth=0.5)
            plot(ax3,linspace(-10,10,size(G3,2)),Diameter(SpeakRaw(1,2)-BLPeriod(2)/2*Param.Fs:SpeakRaw(1,2)+BLPeriod(2)/2*Param.Fs),color=[0 0 0 0.2],linewidth=0.5)
            plot(ax4,linspace(-10,10,size(G4,2)),Diameter(ListenRaw(1,2)-BLPeriod(2)/2*Param.Fs:ListenRaw(1,2)+BLPeriod(2)/2*Param.Fs),color=[0 0 0 0.2],linewidth=0.5)
            G1(x,:)=Diameter(BLStartEnd(1):BLStartEnd(2));
            G2(x,:)=Diameter(BLStartEnd(1):BLStartEnd(2));
            G3(x,:)=Diameter(SpeakRaw(1,2)-BLPeriod(2)/2*Param.Fs:SpeakRaw(1,2)+BLPeriod(2)/2*Param.Fs);
            G4(x,:)=Diameter(ListenRaw(1,2)-BLPeriod(2)/2*Param.Fs:ListenRaw(1,2)+BLPeriod(2)/2*Param.Fs);
        end
        title(ax1,'G1: BL Period Initially Speaking')
        title(ax2,'G2: BL Period Initially Listening')
        title(ax3,'G3: Speaking-Evoked')
        title(ax4,'G4: Listening-Evoked')
        xlabel([ax1 ax2 ax3 ax4],'Time [s]')
        ylabel([ax1 ax2 ax3 ax4],'Pupil diameter [mm]')
        xlim([ax1 ax2],[0, 20])
        xlim([ax3 ax4],[-10, 10])
        grid([ax1 ax2 ax3 ax4],'on')
        
        % Last file
        if contains([PairFiles(i).folder, '\', PairFiles(i).name],'\Main12\P2_SHL_B2.mat')
            xline(ax3,0,'--')
            xline(ax4,0,'--')
            G1(~any(G1,2),:)=[];
            G2(~any(G2,2),:)=[];
            G3(~any(G3,2),:)=[];
            G4(~any(G4,2),:)=[];
            G1_Mean = mean(G1,1);
            G2_Mean = mean(G2,1);
            G3_Mean = mean(G3,1);
            G4_Mean = mean(G4,1);
            G1_SEM = std(G1,[],1)/sqrt(size(G1,1));
            G2_SEM = std(G2,[],1)/sqrt(size(G2,1));
            G3_SEM = std(G3,[],1)/sqrt(size(G3,1));
            G4_SEM = std(G4,[],1)/sqrt(size(G4,1));
            plot(ax1,t_Diam(1:size(G1_Mean,2)),G1_Mean,"Color",'r',"LineWidth",3)
            plot(ax2,t_Diam(1:size(G2_Mean,2)),G2_Mean,"Color",'r',"LineWidth",3)
            plot(ax3,linspace(-10,10,size(G3,2)),G3_Mean,"Color",'r',"LineWidth",3)
            plot(ax4,linspace(-10,10,size(G4,2)),G4_Mean,"Color",'r',"LineWidth",3)
            fill(ax1,[t_Diam(1:size(G1_Mean,2)), flipud(t_Diam(1:size(G1_Mean,2))')'],[(G1_Mean+G1_SEM), flipud((G1_Mean-G1_SEM)')'],[154 0 0]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax2,[t_Diam(1:size(G2_Mean,2)), flipud(t_Diam(1:size(G2_Mean,2))')'],[(G2_Mean+G2_SEM), flipud((G2_Mean-G2_SEM)')'],[154 0 0]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax3,[linspace(-10,10,size(G3,2)), flipud(linspace(-10,10,size(G3,2))')'],[(G3_Mean+G3_SEM), flipud((G3_Mean-G3_SEM)')'],[154 0 0]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
            fill(ax4,[linspace(-10,10,size(G4,2)), flipud(linspace(-10,10,size(G4,2))')'],[(G4_Mean+G4_SEM), flipud((G4_Mean-G4_SEM)')'],[154 0 0]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
% xlim([ax2 ax1],[9.5, 10.5]);xlim([ax3 ax4],[-0.5, 0.5]);ylim([ax1 ax2 ax3 ax4],[3.4, 3.7])
        end
        
    end
end
