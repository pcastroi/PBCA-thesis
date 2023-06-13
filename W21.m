%% PBCA-Thesis - Week 21 - Pupil Features for AMEND I - Copy of W4.m
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

[subDirs] = GetSubDirsFirstLevelOnly('data\AMEND_I');
FileNames={'P1_Quiet_B1.mat','P1_Quiet_B2.mat','P1_SHL_B1.mat','P1_SHL_B2.mat','P1_Noise60_B1.mat','P1_Noise60_B2.mat','P1_Noise70_B1.mat','P1_Noise70_B2.mat','P2_Quiet_B1.mat','P2_Quiet_B2.mat','P2_SHL_B1.mat','P2_SHL_B2.mat','P2_Noise60_B1.mat','P2_Noise60_B2.mat','P2_Noise70_B1.mat','P2_Noise70_B2.mat'};
LoadDelays=load('data\AMEND_I\delays1110.mat');
LoadTPsOrder=load('data\AMEND_I\TPsOrder_I.mat');
LoadUtt=load('data\AMEND_I\utterances1110.mat');

NPs = 2; % N of TPs per Pair
NCond = numel(FileNames)/NPs; % N of conditions
NTPs = numel(subDirs)*numel(FileNames)/NCond; % Total N of TPs
NCols=200; % N of windows
NRows=numel(FileNames)*numel(subDirs); % Number of trials

% Preallocate features - Phase 1, Phase 2
F1_S_Q=nan*ones(NRows,NCols); % Feature 1 (Quiet): Mean pupil size (Speaking)
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
            SpeB = 0;
        elseif contains(cell2mat(FileNames(i)),'B2')
            SpeB = 1;
        end

        if contains(cell2mat(FileNames(i)),'Quiet')
            SpeCond = SpeB + 1;
        elseif contains(cell2mat(FileNames(i)),'SHL')
            SpeCond = SpeB + 3;
        elseif contains(cell2mat(FileNames(i)),'Noise60')
            SpeCond = SpeB + 5;
        elseif contains(cell2mat(FileNames(i)),'Noise70')
            SpeCond = SpeB + 7;
        end

        SpeakRaw = PairUtt{1,SpeCond}.(SpeakKey);
        ListenRaw = PairUtt{1,SpeCond}.(ListenKey);
        binResUtt = PairUtt{1,SpeCond}.binRes;
        SDelayRaw = PairDelay{1,SpeCond}.(SDelayKey);
        LDelayRaw = PairDelay{1,SpeCond}.(LDelayKey);

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
        
        Speak = Speak(Speak(:,1)>TimeMinWin,:);
        Listen = Listen(Listen(:,1)>TimeMinWin,:);
        
        % Time-locked indexes (cannot be bigger than diameter itself)
%         SIdxEnd=Speak(:,3);
%         LIdxEnd=Listen(:,3);
        SIdxEnd=Speak(:,2)+TimeEndW*Param.Fs;
        SIdxEnd(SIdxEnd>length(Diameter)) = length(Diameter);
        LIdxEnd=Listen(:,2)+TimeEndW*Param.Fs;
        LIdxEnd(LIdxEnd>length(Diameter)) = length(Diameter);
        
        t_Diam = linspace(0,length(Diameter)./Param.Fs,length(Diameter));
        
        % Extract features for different conditions
        temp=zeros(NCols,1);
        temp2=zeros(NCols,2);
        temp3=temp;
        if contains(cell2mat(FileNames(i)),'Quiet')
            for j=1:size(Speak,1)
                temp(j,:)=mean(Diameter(Speak(j,2):SIdxEnd(j)));
                temp2(j,:)=polyfit(t_Diam(Speak(j,2):SIdxEnd(j)),Diameter(Speak(j,2):SIdxEnd(j)),1);
                temp3(j,:)=max(Diameter(Speak(j,2):SIdxEnd(j)));
            end
            F1_S_Q(x,1:size(temp(temp(:,1)>0,1),1))=temp(temp(:,1)>0,1);
            F2_S_Q(x,1:size(temp2(temp2(:,2)~=0,2),1))=temp2(temp2(:,2)~=0,2);
            F3_S_Q(x,1:size(temp3(temp3(:,1)>0,1),1))=temp3(temp3(:,1)>0,1);
            for j=1:size(Listen,1)
                temp(j,:)=mean(Diameter(Listen(j,2):LIdxEnd(j)));
                temp2(j,:)=polyfit(t_Diam(Listen(j,2):LIdxEnd(j)),Diameter(Listen(j,2):LIdxEnd(j)),1);
                temp3(j,:)=max(Diameter(Listen(j,2):LIdxEnd(j)));
            end
            F1_L_Q(x,1:size(temp(temp(:,1)>0,1),1))=temp(temp(:,1)>0,1);
            F2_L_Q(x,1:size(temp2(temp2(:,2)~=0,2),1))=temp2(temp2(:,2)~=0,2);
            F3_L_Q(x,1:size(temp3(temp3(:,1)>0,1),1))=temp3(temp3(:,1)>0,1);
        elseif contains(cell2mat(FileNames(i)),'SHL')
            for j=1:size(Speak,1)
                temp(j,:)=mean(Diameter(Speak(j,2):SIdxEnd(j)));
                temp2(j,:)=polyfit(t_Diam(Speak(j,2):SIdxEnd(j)),Diameter(Speak(j,2):SIdxEnd(j)),1);
                temp3(j,:)=max(Diameter(Speak(j,2):SIdxEnd(j)));
            end
            F1_S_SHL(x,1:size(temp(temp(:,1)>0,1),1))=temp(temp(:,1)>0,1);
            F2_S_SHL(x,1:size(temp2(temp2(:,2)~=0,2),1))=temp2(temp2(:,2)~=0,2);
            F3_S_SHL(x,1:size(temp3(temp3(:,1)>0,1),1))=temp3(temp3(:,1)>0,1);
            for j=1:size(Listen,1)
                temp(j,:)=mean(Diameter(Listen(j,2):LIdxEnd(j)));
                temp2(j,:)=polyfit(t_Diam(Listen(j,2):LIdxEnd(j)),Diameter(Listen(j,2):LIdxEnd(j)),1);
                temp3(j,:)=max(Diameter(Listen(j,2):LIdxEnd(j)));
            end
            F1_L_SHL(x,1:size(temp(temp(:,1)>0,1),1))=temp(temp(:,1)>0,1);
            F2_L_SHL(x,1:size(temp2(temp2(:,2)~=0,2),1))=temp2(temp2(:,2)~=0,2);
            F3_L_SHL(x,1:size(temp3(temp3(:,1)>0,1),1))=temp3(temp3(:,1)>0,1);
        elseif contains(cell2mat(FileNames(i)),'Noise60')
            for j=1:size(Speak,1)
                temp(j,:)=mean(Diameter(Speak(j,2):SIdxEnd(j)));
                temp2(j,:)=polyfit(t_Diam(Speak(j,2):SIdxEnd(j)),Diameter(Speak(j,2):SIdxEnd(j)),1);
                temp3(j,:)=max(Diameter(Speak(j,2):SIdxEnd(j)));
            end
            F1_S_N60(x,1:size(temp(temp(:,1)>0,1),1))=temp(temp(:,1)>0,1);
            F2_S_N60(x,1:size(temp2(temp2(:,2)~=0,2),1))=temp2(temp2(:,2)~=0,2);
            F3_S_N60(x,1:size(temp3(temp3(:,1)>0,1),1))=temp3(temp3(:,1)>0,1);
            for j=1:size(Listen,1)
                temp(j,:)=mean(Diameter(Listen(j,2):LIdxEnd(j)));
                temp2(j,:)=polyfit(t_Diam(Listen(j,2):LIdxEnd(j)),Diameter(Listen(j,2):LIdxEnd(j)),1);
                temp3(j,:)=max(Diameter(Listen(j,2):LIdxEnd(j)));
            end
            F1_L_N60(x,1:size(temp(temp(:,1)>0,1),1))=temp(temp(:,1)>0,1);
            F2_L_N60(x,1:size(temp2(temp2(:,2)~=0,2),1))=temp2(temp2(:,2)~=0,2);
            F3_L_N60(x,1:size(temp3(temp3(:,1)>0,1),1))=temp3(temp3(:,1)>0,1);
            
        elseif contains(cell2mat(FileNames(i)),'Noise70')
            for j=1:size(Speak,1)
                temp(j,:)=mean(Diameter(Speak(j,2):SIdxEnd(j)));
                temp2(j,:)=polyfit(t_Diam(Speak(j,2):SIdxEnd(j)),Diameter(Speak(j,2):SIdxEnd(j)),1);
                temp3(j,:)=max(Diameter(Speak(j,2):SIdxEnd(j)));
            end
            F1_S_N70(x,1:size(temp(temp(:,1)>0,1),1))=temp(temp(:,1)>0,1);
            F2_S_N70(x,1:size(temp2(temp2(:,2)~=0,2),1))=temp2(temp2(:,2)~=0,2);
            F3_S_N70(x,1:size(temp3(temp3(:,1)>0,1),1))=temp3(temp3(:,1)>0,1);
            for j=1:size(Listen,1)
                temp(j,:)=mean(Diameter(Listen(j,2):LIdxEnd(j)));
                temp2(j,:)=polyfit(t_Diam(Listen(j,2):LIdxEnd(j)),Diameter(Listen(j,2):LIdxEnd(j)),1);
                temp3(j,:)=max(Diameter(Listen(j,2):LIdxEnd(j)));
            end
            F1_L_N70(x,1:size(temp(temp(:,1)>0,1),1))=temp(temp(:,1)>0,1);
            F2_L_N70(x,1:size(temp2(temp2(:,2)~=0,2),1))=temp2(temp2(:,2)~=0,2);
            F3_L_N70(x,1:size(temp3(temp3(:,1)>0,1),1))=temp3(temp3(:,1)>0,1);
        end
        x=x+1;
    end
end
%% Plots
f1=figure;tl1=tiledlayout(1,2);
ax11 = nexttile;
ax12 = nexttile;
f2=figure;tl2=tiledlayout(1,2);
ax21 = nexttile;
ax22 = nexttile;
f3=figure;tl3=tiledlayout(1,2);
ax31 = nexttile;
ax32 = nexttile;
hold([ax11 ax12 ax21 ax22 ax31 ax32],'on')

% Boxplot
h11=boxplotGroup(ax11,{F1_S_Q(logical(~isnan(F1_S_Q(:,:,:)))),F1_S_SHL(logical(~isnan(F1_S_SHL(:,:,:)))),F1_S_N60(logical(~isnan(F1_S_N60(:,:,:)))),F1_S_N70(logical(~isnan(F1_S_N70(:,:,:))))},...
'PrimaryLabels', {'Quiet','SHL','N60','N70'}, ...
'SecondaryLabels', {'Speaking'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[QuietColor;SHLColor;N60Color;N70Color],'GroupType','betweenGroups');
h12=boxplotGroup(ax12,{F1_L_Q(logical(~isnan(F1_L_Q(:,:,:)))),F1_L_SHL(logical(~isnan(F1_L_SHL(:,:,:)))),F1_L_N60(logical(~isnan(F1_L_N60(:,:,:)))),F1_L_N70(logical(~isnan(F1_L_N70(:,:,:))))},...
'PrimaryLabels', {'Quiet','SHL','N60','N70'}, ...
'SecondaryLabels', {'Listening'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[QuietColor;SHLColor;N60Color;N70Color],'GroupType','betweenGroups');

h21=boxplotGroup(ax21,{F2_S_Q(logical(~isnan(F2_S_Q(:,:,:)))),F2_S_SHL(logical(~isnan(F2_S_SHL(:,:,:)))),F2_S_N60(logical(~isnan(F2_S_N60(:,:,:)))),F2_S_N70(logical(~isnan(F2_S_N70(:,:,:))))},...
'PrimaryLabels', {'Quiet','SHL','N60','N70'}, ...
'SecondaryLabels', {'Speaking'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[QuietColor;SHLColor;N60Color;N70Color],'GroupType','betweenGroups');
h22=boxplotGroup(ax22,{F2_L_Q(logical(~isnan(F2_L_Q(:,:,:)))),F2_L_SHL(logical(~isnan(F2_L_SHL(:,:,:)))),F2_L_N60(logical(~isnan(F2_L_N60(:,:,:)))),F2_L_N70(logical(~isnan(F2_L_N70(:,:,:))))},...
'PrimaryLabels', {'Quiet','SHL','N60','N70'}, ...
'SecondaryLabels', {'Listening'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[QuietColor;SHLColor;N60Color;N70Color],'GroupType','betweenGroups');

h31=boxplotGroup(ax31,{F3_S_Q(logical(~isnan(F3_S_Q(:,:,:)))),F3_S_SHL(logical(~isnan(F3_S_SHL(:,:,:)))),F3_S_N60(logical(~isnan(F3_S_N60(:,:,:)))),F3_S_N70(logical(~isnan(F3_S_N70(:,:,:))))},...
'PrimaryLabels', {'Quiet','SHL','N60','N70'}, ...
'SecondaryLabels', {'Speaking'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[QuietColor;SHLColor;N60Color;N70Color],'GroupType','betweenGroups');
h32=boxplotGroup(ax32,{F3_L_Q(logical(~isnan(F3_L_Q(:,:,:)))),F3_L_SHL(logical(~isnan(F3_L_SHL(:,:,:)))),F3_L_N60(logical(~isnan(F3_L_N60(:,:,:)))),F3_L_N70(logical(~isnan(F3_L_N70(:,:,:))))},...
'PrimaryLabels', {'Quiet','SHL','N60','N70'}, ...
'SecondaryLabels', {'Listening'}, ...
'interGroupSpace',2,'groupLabelType','Vertical', ...
'PlotStyle','Compact','BoxStyle','filled',...
'Colors',[QuietColor;SHLColor;N60Color;N70Color],'GroupType','betweenGroups');

ax11.YGrid = 'on';
ax11.XTickLabelRotation = 90;
set(ax11,'Color',[SColor,0.1])
ax12.YGrid = 'on';
ax12.XTickLabelRotation = 90;
set(ax12,'Color',[LColor,0.1])
ylim([ax11 ax12],[min([ylim(ax11) ylim(ax12)]) max([ylim(ax11) ylim(ax12)])])

ax21.YGrid = 'on';
ax21.XTickLabelRotation = 90;
set(ax21,'Color',[SColor,0.1])
ax22.YGrid = 'on';
ax22.XTickLabelRotation = 90;
set(ax22,'Color',[LColor,0.1])
ylim([ax21 ax22],[min([ylim(ax21) ylim(ax22)]) max([ylim(ax21) ylim(ax22)])])

ax31.YGrid = 'on';
ax31.XTickLabelRotation = 90;
set(ax31,'Color',[SColor,0.1])
ax32.YGrid = 'on';
ax32.XTickLabelRotation = 90;
set(ax32,'Color',[LColor,0.1])
ylim([ax31 ax32],[min([ylim(ax31) ylim(ax32)]) max([ylim(ax31) ylim(ax32)])])

title(tl1,'Mean Pupil Diameter','FontWeight','bold')
title(tl2,'Pupil Diameter Slope','FontWeight','bold')
title(tl3,'Peak Pupil Diameter','FontWeight','bold')

ylabel([ax11 ax31],'Pupil diameter [mm]')
ylabel([ax21],'Slope [mm/s]')