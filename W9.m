%% PBCA-Thesis - Week 9, 10 & 11 - Fixation duration? Not using diameter anymore
% Pathing
clear all; clc; close all;
BPath = strsplit(pwd,'PBCA-thesis');
addpath('tools')
addpath([BPath{1} 'Pupil-preprocessing-tools\tools']) % For preprocessing

[subDirs] = GetSubDirsFirstLevelOnly('data\AMEND_I');
FileNames={'P1_Quiet_B1.mat','P1_Quiet_B2.mat','P1_SHL_B1.mat','P1_SHL_B2.mat','P1_Noise60_B1.mat','P1_Noise60_B2.mat','P1_Noise70_B1.mat','P1_Noise70_B2.mat','P2_Quiet_B1.mat','P2_Quiet_B2.mat','P2_SHL_B1.mat','P2_SHL_B2.mat','P2_Noise60_B1.mat','P2_Noise60_B2.mat','P2_Noise70_B1.mat','P2_Noise70_B2.mat'};
LoadUtt=load('data\AMEND_I\utterances1110.mat');
LoadDelays=load('data\AMEND_I\delays1110.mat');

% Parameters for processing
Param.Fs = 50; % Sampling frequency of pupil data
Param.Preblink = 0.1; % [s], set to NaN, time before blink
Param.Postblink = 0.2; % [s], set to NaN, time after blink
Param.BlinkThresh = 3; % [samples], threshold of samples in between artifacts or blinks
MaxVisualAngle = 0.5; % [degrees], critical visual angle for fixation definition
LPWinSize = 0.5; % [s]: Window size of hamming-window for low-pass filtering
FilterWidth = round((LPWinSize*Param.Fs)/2); % [samples]: Width of hamming filter used for fixation duration
LPWindow = hamming(round(LPWinSize*Param.Fs));
LPWindow = LPWindow/sum(LPWindow); % Hamming-window

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
GSW = zeros(200,200,4000); % Global Speaking Windows
GLW = zeros(200,200,4000); % Global Listening Windows
SW_Quiet = zeros(200,200,4000); % Quiet Speaking Windows
LW_Quiet = zeros(200,200,4000); % Quiet Listening Windows
SW_SHL = zeros(200,200,4000); % SHL Speaking Windows
LW_SHL = zeros(200,200,4000); % SHL Listening Windows
SW_N60 = zeros(200,200,4000); % N60 Speaking Windows
LW_N60 = zeros(200,200,4000); % N60 Listening Windows
SW_N70 = zeros(200,200,4000); % N70 Speaking Windows
LW_N70 = zeros(200,200,4000); % N70 Listening Windows
GSWB = zeros(200,200,4000); % Global Speaking Windows Baselined
GLWB = zeros(200,200,4000); % Global Listening Windows Baselined
SWB_Quiet = zeros(200,200,4000); % Quiet Speaking Windows Baselined
LWB_Quiet = zeros(200,200,4000); % Quiet Listening Windows Baselined
SWB_SHL = zeros(200,200,4000); % SHL Speaking Windows Baselined
LWB_SHL = zeros(200,200,4000); % SHL Listening Windows Baselined
SWB_N60 = zeros(200,200,4000); % N60 Speaking Windows Baselined
LWB_N60 = zeros(200,200,4000); % N60 Listening Windows Baselined
SWB_N70 = zeros(200,200,4000); % N70 Speaking Windows Baselined
LWB_N70 = zeros(200,200,4000); % N70 Listening Windows Baselined

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
        try
            alldata = load([PairFiles(1).folder, '\', cell2mat(FileNames(i))]);
        catch ME
            disp(['Warning: File ', PairFiles(1).folder, '\', cell2mat(FileNames(i)), ' not found (no Gaze data).']);
            continue
        end
        alldata_mat = cell2mat(alldata.data);
        
        % In case gaze3d is transposed from origin [x;y;z] -> transpose to [x,y,z]
        if size([alldata_mat(:,cellfun(@(xd) ~any(isnan(xd)),{alldata_mat.gaze3d})).gaze3d],1) == 3
            for k=1:size(alldata_mat,2)
               alldata_mat(k).gaze3d=alldata_mat(k).gaze3d'; 
            end
        end
        
        % Replace blanks '[]' for 'NaN' in gaze3d
        [alldata_mat(cellfun(@isempty,{alldata_mat.gaze3d})).gaze3d] = deal(NaN);
        
        % Replace 'NaN' to '[NaN,NaN,NaN]'
        [alldata_mat(:,cellfun(@(xd) any(isnan(xd)),{alldata_mat.gaze3d})).gaze3d] = deal([NaN,NaN,NaN]);
        
        % Extract delay [s] (Duration_timeStamps - Duration_nsamples)
        EyeAudDelay=alldata_mat(end).timeStamp-alldata_mat(1).timeStamp-length(alldata_mat)/Param.Fs;
        
        % Skip file if the difference in duration from the number of
        % samples and the duration given by timestamps is bigger than 0.5 s
        if EyeAudDelay > RejectDelay
            disp(['Warning: File ',PairFiles(1).folder, '\', cell2mat(FileNames(i)), ' was rejected, too much delay (',sprintf('%0.2f',EyeAudDelay),'s).']);
            continue
        end
        
        gaze3d = vertcat(alldata_mat.gaze3d);
        
        GazeXRaw = gaze3d(:,1);
        GazeYRaw = gaze3d(:,2);
        GazeZRaw = gaze3d(:,3);
        
        GazeX = GazeXRaw;
        GazeY = GazeYRaw;
        GazeZ = GazeZRaw;
        
        % NaN's appear at the same time at [X,Y,Z], we only look at X.
        XNan = find(isnan(GazeXRaw));
        
        % Reject a file if gaze3d has too many NaNs
        if length(XNan)/length(gaze3d) >= RejectRatio
            disp(['Warning: File ',PairFiles(1).folder, '\', cell2mat(FileNames(i)), ' was rejected because gaze3d contains too many NaNs (',sprintf('%0.2f',100*length(XNan)/length(gaze3d)),'%).'])
            continue
        end
        
%%%%%%%%%%%%%%%%%%% NOT CURRENTLY USED IN THE SCRIPT!!! %%%%%%%%%%%%%%%%%%%
%         % Preprocessing - 100 ms before and 200 ms after a blink -> NaN
%         if ~isempty(XNan)
%             % Pre-blink
%             Pre = [XNan(1);XNan(find(diff(XNan,1) > 1) + 1)];
%             for h=1:length(Pre)
%                 if Pre(h) > Param.Preblink*Param.Fs % Check start delimiter
%                    for k=1:Param.Preblink*Param.Fs
%                        GazeX(Pre(h)-k)=NaN;
%                        GazeY(Pre(h)-k)=NaN;
%                        GazeZ(Pre(h)-k)=NaN;
%                    end
%                 end
%             end
% 
%             % Post-blink
%             Post = XNan(diff(XNan,1) > Param.BlinkThresh);
%             for h=1:length(Post)
%                 if Post(h) < length(GazeX) + Param.Postblink*Param.Fs % Check end delimiter
%                    for k=1:Param.Postblink*Param.Fs
%                        GazeX(Post(h)+k)=NaN;
%                        GazeY(Post(h)+k)=NaN;
%                        GazeZ(Post(h)+k)=NaN;
%                    end
%                 end
%             end
%         end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Preprocessing - Setting outliers as NaNs (remove artifacts)
        XThresh = [mean(GazeX,'omitnan')-std(GazeX,'omitnan'),mean(GazeX,'omitnan')+std(GazeX,'omitnan')];
        YThresh = [mean(GazeY,'omitnan')-std(GazeY,'omitnan'),mean(GazeY,'omitnan')+std(GazeY,'omitnan')];
        ZThresh = [mean(GazeZ,'omitnan')-std(GazeZ,'omitnan'),mean(GazeZ,'omitnan')+std(GazeZ,'omitnan')];
        F_NOutl = 1; % number of outliers per file

        for s=1:length(GazeX)
            if GazeX(s) < XThresh(1) || GazeX(s) > XThresh(2)
                GazeX(s)=NaN;
                F_NOutl = F_NOutl + 1;
            elseif GazeY(s) < YThresh(1) || GazeY(s) > YThresh(2)
                GazeY(s)=NaN;
                F_NOutl = F_NOutl + 1;
            elseif GazeZ(s) < ZThresh(1) || GazeZ(s) > ZThresh(2)
                GazeZ(s)=NaN;
                F_NOutl = F_NOutl + 1;
            end
        end
        
        % Reject a file if its fixation contains too many outliers
        if F_NOutl/numel(gaze3d) >= RejectRatio
            disp(['Warning: File ',PairFiles(1).folder, '\', cell2mat(FileNames(i)), ' was rejected because its fixation duration contains too many outliers (',sprintf('%0.2f',(100*F_NOutl)/numel(gaze3d)),'%).'])
            continue
        end

        % Fixation duration
        fixation = fix_duration([GazeX,GazeY,GazeZ],MaxVisualAngle,Param.Fs);
        
%         figure;plot(linspace(0,length(fixation)/Param.Fs,length(fixation)),fixation,'LineWidth',1.5,'Color','k','DisplayName','Original fixation');hold on;for k=1:10;plot(linspace(0,length(fixation)/Param.Fs,length(fixation)),ndnanfilter(fixation,'hamming',k),'DisplayName',['Hamming window: F=',num2str(k)]);end;legend;grid on;yline(FThresh(1),'--');yline(FThresh(2),'--');xlabel('Time [s]');ylabel('Fixation duration [s]');title(['File: ',PairFiles(1).folder, '\', cell2mat(FileNames(i))])
        
        % Filtering - Fixation duration with hamming window (F = 3)
        fixation = ndnanfilter(fixation,'hamming',FilterWidth);
        
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
        try
            SDelayRaw = PairDelay{1,CondKey}.(SDelayKey);
            LDelayRaw = PairDelay{1,CondKey}.(LDelayKey);
        catch ME
            disp(['Warning: File ', PairFiles(1).folder, '\', cell2mat(FileNames(i)), ' is missing the associated Delay data.']);
            continue
        end
        
        if or(SDelayRaw < 0,LDelayRaw < 0)
            SDelayRaw=[0,0];
            LDelayRaw=[0,0];
        end
        
        binResDel = PairDelay{1,CondKey}.binRes;

        if isempty(SpeakRaw) && isempty(ListenRaw)
            disp(['Warning: File ', PairFiles(1).folder, '\', cell2mat(FileNames(i)), ' is missing the associated Utterance/Listening data.']);
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
        SpeakD = SpeakMI(SpeakMI(:,1)>TimeMinWin,:);
        ListenD = ListenMI(ListenMI(:,1)>TimeMinWin,:);
        
        % Merge again if duration between windows <= TimeMerge (2 s)
        SpeakM = merge_windows(SpeakD, Param.Fs, TimeMerge);
        ListenM = merge_windows(ListenD, Param.Fs, TimeMerge);
        
        % Discard windows if duration is < 2*TimeMinWin (1 s)
        Speak = SpeakM(SpeakM(:,1)>2*TimeMinWin,:);
        Listen = ListenM(ListenM(:,1)>2*TimeMinWin,:);
        
        % Time-locked indexes (based on Start or End of events)
        WSpeakIdx=[Speak(:,2)-TimeStartW*Param.Fs,Speak(:,2),Speak(:,2)+TimeEndW*Param.Fs];
        WListenIdx=[Listen(:,2)-TimeStartW*Param.Fs,Listen(:,2),Listen(:,2)+TimeEndW*Param.Fs];
%         EWSpeakIdx=[Speak(:,3)-TimeStartW*Param.Fs,Speak(:,3),Speak(:,3)+TimeEndW*Param.Fs];
%         EWListenIdx=[Listen(:,3)-TimeStartW*Param.Fs,Listen(:,3),Listen(:,3)+TimeEndW*Param.Fs];
        
        % Storing Speaking/Listening by conditions
        if contains(cell2mat(FileNames(i)),'Quiet')
            for j=1:size(Speak,1)
                SW_Quiet(x,j,:)=[fixation(WSpeakIdx(j,1):Speak(j,3)-1);NaN*ones(1,length(SW_Quiet)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SWB_Quiet(x,j,:)=[fixation(WSpeakIdx(j,1):Speak(j,3)-1)-mean(fixation(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SWB_Quiet)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SDur_Quiet(i,j) = Speak(j,1);
            end
            for j=1:size(Listen,1)
                LW_Quiet(x,j,:)=[fixation(WListenIdx(j,1):Listen(j,3)-1);NaN*ones(1,length(LW_Quiet)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LWB_Quiet(x,j,:)=[fixation(WListenIdx(j,1):Listen(j,3)-1)-mean(fixation(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LWB_Quiet)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LDur_Quiet(i,j) = Listen(j,1);
            end
        elseif contains(cell2mat(FileNames(i)),'SHL')
            for j=1:size(Speak,1)
                SW_SHL(x,j,:)=[fixation(WSpeakIdx(j,1):Speak(j,3)-1);NaN*ones(1,length(SW_SHL)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SWB_SHL(x,j,:)=[fixation(WSpeakIdx(j,1):Speak(j,3)-1)-mean(fixation(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SWB_SHL)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SDur_SHL(i,j) = Speak(j,1);
            end
            for j=1:size(Listen,1)
                LW_SHL(x,j,:)=[fixation(WListenIdx(j,1):Listen(j,3)-1);NaN*ones(1,length(LW_SHL)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LWB_SHL(x,j,:)=[fixation(WListenIdx(j,1):Listen(j,3)-1)-mean(fixation(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LWB_SHL)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LDur_SHL(i,j) = Listen(j,1);
            end
        elseif contains(cell2mat(FileNames(i)),'Noise60')
            for j=1:size(Speak,1)
                SW_N60(x,j,:)=[fixation(WSpeakIdx(j,1):Speak(j,3)-1);NaN*ones(1,length(SW_N60)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SWB_N60(x,j,:)=[fixation(WSpeakIdx(j,1):Speak(j,3)-1)-mean(fixation(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SWB_N60)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SDur_N60(i,j) = Speak(j,1);
            end
            for j=1:size(Listen,1)
                LW_N60(x,j,:)=[fixation(WListenIdx(j,1):Listen(j,3)-1);NaN*ones(1,length(LW_N60)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LWB_N60(x,j,:)=[fixation(WListenIdx(j,1):Listen(j,3)-1)-mean(fixation(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LWB_N60)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LDur_N60(i,j) = Listen(j,1);
            end
        elseif contains(cell2mat(FileNames(i)),'Noise70')
            for j=1:size(Speak,1)
                SW_N70(x,j,:)=[fixation(WSpeakIdx(j,1):Speak(j,3)-1);NaN*ones(1,length(SW_N70)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SWB_N70(x,j,:)=[fixation(WSpeakIdx(j,1):Speak(j,3)-1)-mean(fixation(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SWB_N70)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
                SDur_N70(i,j) = Speak(j,1);
            end
            for j=1:size(Listen,1)
                LW_N70(x,j,:)=[fixation(WListenIdx(j,1):Listen(j,3)-1);NaN*ones(1,length(LW_N70)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LWB_N70(x,j,:)=[fixation(WListenIdx(j,1):Listen(j,3)-1)-mean(fixation(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LWB_N70)-length(WListenIdx(j,1):Listen(j,3)-1))'];
                LDur_N70(x,j) = Listen(j,1);
            end
        end
        
        SW = zeros(size(Speak,1),ceil(Param.Fs*(TimeStartW+max(Speak(:,1)))));
        SWB = SW;
        for j=1:size(Speak,1)
            % Add nan-padding when necessary
            SW(j,:)=[fixation(WSpeakIdx(j,1):Speak(j,3)-1);NaN*ones(1,length(SW)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
            SWB(j,:)=[fixation(WSpeakIdx(j,1):Speak(j,3)-1)-mean(fixation(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(SWB)-length(WSpeakIdx(j,1):Speak(j,3)-1))'];
            GSW(x,j,1:length(SW(j,:))) = SW(j,:);
            GSWB(x,j,1:length(SWB(j,:))) = SWB(j,:);
            GSDur(x,j) = Speak(j,1);
        end
        
        LW = zeros(size(Listen,1),ceil(Param.Fs*(TimeStartW+max(Listen(:,1)))));
        LWB = LW;
        for j=1:size(Listen,1)
            % Add nan-padding when necessary
            LW(j,:)=[fixation(WListenIdx(j,1):Listen(j,3)-1);NaN*ones(1,length(LW)-length(WListenIdx(j,1):Listen(j,3)-1))'];
            LWB(j,:)=[fixation(WListenIdx(j,1):Listen(j,3)-1)-mean(fixation(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(LWB)-length(WListenIdx(j,1):Listen(j,3)-1))'];
            GLW(x,j,1:length(LW(j,:))) = LW(j,:);
            GLWB(x,j,1:length(LWB(j,:))) = LWB(j,:);
            GLDur(x,j) = Listen(j,1);
        end
        
        % Increase index of num of files used
        x=x+1;
    end
end

%% Global Plots

figure;tiledlayout(1,2);ax1 = nexttile;ax2 = nexttile;
figure;tiledlayout(1,2);ax3 = nexttile;ax4 = nexttile;
figure;tiledlayout(1,2);ax5 = nexttile;ax6 = nexttile;
hold([ax1 ax2 ax3 ax4 ax5 ax6],'on')
grid([ax1 ax2 ax3 ax4 ax5 ax6],'on')

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


% Global fixation duration responses
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
GSWB(~any(GSWB,[2 3]),:,:)=[];GSWB(:,~any(GSWB,[1 3]),:)=[];GSWB(GSWB==0)=NaN;
GLWB(~any(GLWB,[2 3]),:,:)=[];GLWB(:,~any(GLWB,[1 3]),:)=[];GLWB(GLWB==0)=NaN;
SWB_Quiet(~any(SWB_Quiet,[2 3]),:,:)=[];SWB_Quiet(:,~any(SWB_Quiet,[1 3]),:)=[];SWB_Quiet(SWB_Quiet==0)=NaN;
LWB_Quiet(~any(LWB_Quiet,[2 3]),:,:)=[];LWB_Quiet(:,~any(LWB_Quiet,[1 3]),:)=[];LWB_Quiet(LWB_Quiet==0)=NaN;
SWB_SHL(~any(SWB_SHL,[2 3]),:,:)=[];SWB_SHL(:,~any(SWB_SHL,[1 3]),:)=[];SWB_SHL(SWB_SHL==0)=NaN;
LWB_SHL(~any(LWB_SHL,[2 3]),:,:)=[];LWB_SHL(:,~any(LWB_SHL,[1 3]),:)=[];LWB_SHL(LWB_SHL==0)=NaN;
SWB_N60(~any(SWB_N60,[2 3]),:,:)=[];SWB_N60(:,~any(SWB_N60,[1 3]),:)=[];SWB_N60(SWB_N60==0)=NaN;
LWB_N60(~any(LWB_N60,[2 3]),:,:)=[];LWB_N60(:,~any(LWB_N60,[1 3]),:)=[];LWB_N60(LWB_N60==0)=NaN;
SWB_N70(~any(SWB_N70,[2 3]),:,:)=[];SWB_N70(:,~any(SWB_N70,[1 3]),:)=[];SWB_N70(SWB_N70==0)=NaN;
LWB_N70(~any(LWB_N70,[2 3]),:,:)=[];LWB_N70(:,~any(LWB_N70,[1 3]),:)=[];LWB_N70(LWB_N70==0)=NaN;

% Derive mean and LP filter (hamming - 0.5 s) - omitting NaNs
GSW_Mean = ndnanfilter(reshape(mean(GSW,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
GLW_Mean = ndnanfilter(reshape(mean(GLW,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_Quiet_Mean = ndnanfilter(reshape(mean(SW_Quiet,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_Quiet_Mean = ndnanfilter(reshape(mean(LW_Quiet,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_SHL_Mean = ndnanfilter(reshape(mean(SW_SHL,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_SHL_Mean = ndnanfilter(reshape(mean(LW_SHL,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_N60_Mean = ndnanfilter(reshape(mean(SW_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_N60_Mean = ndnanfilter(reshape(mean(LW_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SW_N70_Mean = ndnanfilter(reshape(mean(SW_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LW_N70_Mean = ndnanfilter(reshape(mean(LW_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
GSWB_Mean = ndnanfilter(reshape(mean(GSWB,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
GLWB_Mean = ndnanfilter(reshape(mean(GLWB,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SWB_Quiet_Mean = ndnanfilter(reshape(mean(SWB_Quiet,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LWB_Quiet_Mean = ndnanfilter(reshape(mean(LWB_Quiet,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SWB_SHL_Mean = ndnanfilter(reshape(mean(SWB_SHL,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LWB_SHL_Mean = ndnanfilter(reshape(mean(LWB_SHL,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SWB_N60_Mean = ndnanfilter(reshape(mean(SWB_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LWB_N60_Mean = ndnanfilter(reshape(mean(LWB_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
SWB_N70_Mean = ndnanfilter(reshape(mean(SWB_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
LWB_N70_Mean = ndnanfilter(reshape(mean(LWB_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

GSW_SEM = (reshape(std(GSW,0,[1 2],'omitnan'),[],1)/sqrt(numel(GSW(~isnan(GSW)))))';
GLW_SEM = (reshape(std(GLW,0,[1 2],'omitnan'),[],1)/sqrt(numel(GLW(~isnan(GLW)))))';
SW_Quiet_SEM = (reshape(std(SW_Quiet,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_Quiet(~isnan(SW_Quiet)))))';
LW_Quiet_SEM = (reshape(std(LW_Quiet,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_Quiet(~isnan(LW_Quiet)))))';
SW_SHL_SEM = (reshape(std(SW_SHL,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_SHL(~isnan(SW_SHL)))))';
LW_SHL_SEM = (reshape(std(LW_SHL,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_SHL(~isnan(LW_SHL)))))';
SW_N60_SEM = (reshape(std(SW_N60,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_N60(~isnan(SW_N60)))))';
LW_N60_SEM = (reshape(std(LW_N60,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_N60(~isnan(LW_N60)))))';
SW_N70_SEM = (reshape(std(SW_N70,0,[1 2],'omitnan'),[],1)/sqrt(numel(SW_N70(~isnan(SW_N70)))))';
LW_N70_SEM = (reshape(std(LW_N70,0,[1 2],'omitnan'),[],1)/sqrt(numel(LW_N70(~isnan(LW_N70)))))';
GSWB_SEM = (reshape(std(GSWB,0,[1 2],'omitnan'),[],1)/sqrt(numel(GSWB(~isnan(GSWB)))))';
GLWB_SEM = (reshape(std(GLWB,0,[1 2],'omitnan'),[],1)/sqrt(numel(GLWB(~isnan(GLWB)))))';
SWB_Quiet_SEM = (reshape(std(SWB_Quiet,0,[1 2],'omitnan'),[],1)/sqrt(numel(SWB_Quiet(~isnan(SWB_Quiet)))))';
LWB_Quiet_SEM = (reshape(std(LWB_Quiet,0,[1 2],'omitnan'),[],1)/sqrt(numel(LWB_Quiet(~isnan(LWB_Quiet)))))';
SWB_SHL_SEM = (reshape(std(SWB_SHL,0,[1 2],'omitnan'),[],1)/sqrt(numel(SWB_SHL(~isnan(SWB_SHL)))))';
LWB_SHL_SEM = (reshape(std(LWB_SHL,0,[1 2],'omitnan'),[],1)/sqrt(numel(LWB_SHL(~isnan(LWB_SHL)))))';
SWB_N60_SEM = (reshape(std(SWB_N60,0,[1 2],'omitnan'),[],1)/sqrt(numel(SWB_N60(~isnan(SWB_N60)))))';
LWB_N60_SEM = (reshape(std(LWB_N60,0,[1 2],'omitnan'),[],1)/sqrt(numel(LWB_N60(~isnan(LWB_N60)))))';
SWB_N70_SEM = (reshape(std(SWB_N70,0,[1 2],'omitnan'),[],1)/sqrt(numel(SWB_N70(~isnan(SWB_N70)))))';
LWB_N70_SEM = (reshape(std(LWB_N70,0,[1 2],'omitnan'),[],1)/sqrt(numel(LWB_N70(~isnan(LWB_N70)))))';

xline(ax1,0,'--','handlevisibility','off')
xline(ax2,0,'--','handlevisibility','off')
xline(ax3,0,'--','handlevisibility','off')
xline(ax4,0,'--','handlevisibility','off')

plot(ax1,linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3)),GSW_Mean,color=SColor,linewidth=2)
plot(ax1,linspace(-TimeStartW,size(SW_Quiet,3)/Param.Fs,size(SW_Quiet,3)),SW_Quiet_Mean,color=QuietColor)
plot(ax1,linspace(-TimeStartW,size(SW_SHL,3)/Param.Fs,size(SW_SHL,3)),SW_SHL_Mean,color=SHLColor)
plot(ax1,linspace(-TimeStartW,size(SW_N60,3)/Param.Fs,size(SW_N60,3)),SW_N60_Mean,color=N60Color)
plot(ax1,linspace(-TimeStartW,size(SW_N70,3)/Param.Fs,size(SW_N70,3)),SW_N70_Mean,color=N70Color)

plot(ax2,nan,color=SColor,linewidth=2) % plot nans to show color in legend
plot(ax2,linspace(-TimeStartW,size(GLW,3)/Param.Fs,size(GLW,3)),GLW_Mean,color=LColor,linewidth=2)
plot(ax2,linspace(-TimeStartW,size(LW_Quiet,3)/Param.Fs,size(LW_Quiet,3)),LW_Quiet_Mean,color=QuietColor)
plot(ax2,linspace(-TimeStartW,size(LW_SHL,3)/Param.Fs,size(LW_SHL,3)),LW_SHL_Mean,color=SHLColor)
plot(ax2,linspace(-TimeStartW,size(LW_N60,3)/Param.Fs,size(LW_N60,3)),LW_N60_Mean,color=N60Color)
plot(ax2,linspace(-TimeStartW,size(LW_N70,3)/Param.Fs,size(LW_N70,3)),LW_N70_Mean,color=N70Color)

% plot(ax1,linspace(-TimeStartW,size(SW_SHL,3)/Param.Fs,size(SW_SHL,3)),mean([SW_Quiet_Mean;SW_SHL_Mean;SW_N60_Mean;SW_N70_Mean],1,'omitnan'),'k--')
% plot(ax2,linspace(-TimeStartW,size(LW_SHL,3)/Param.Fs,size(LW_SHL,3)),mean([LW_Quiet_Mean;LW_SHL_Mean;LW_N60_Mean;LW_N70_Mean],1,'omitnan'),'k--')

plot(ax3,linspace(-TimeStartW,size(GSWB,3)/Param.Fs,size(GSWB,3)),GSWB_Mean,color=SColor,linewidth=2)
plot(ax3,linspace(-TimeStartW,size(SWB_Quiet,3)/Param.Fs,size(SWB_Quiet,3)),SWB_Quiet_Mean,color=QuietColor)
plot(ax3,linspace(-TimeStartW,size(SWB_SHL,3)/Param.Fs,size(SWB_SHL,3)),SWB_SHL_Mean,color=SHLColor)
plot(ax3,linspace(-TimeStartW,size(SWB_N60,3)/Param.Fs,size(SWB_N60,3)),SWB_N60_Mean,color=N60Color)
plot(ax3,linspace(-TimeStartW,size(SWB_N70,3)/Param.Fs,size(SWB_N70,3)),SWB_N70_Mean,color=N70Color)

plot(ax4,nan,color=SColor,linewidth=2) % plot nans to show color in legend
plot(ax4,linspace(-TimeStartW,size(GLWB,3)/Param.Fs,size(GLWB,3)),GLWB_Mean,color=LColor,linewidth=2)
plot(ax4,linspace(-TimeStartW,size(LWB_Quiet,3)/Param.Fs,size(LWB_Quiet,3)),LWB_Quiet_Mean,color=QuietColor)
plot(ax4,linspace(-TimeStartW,size(LWB_SHL,3)/Param.Fs,size(LWB_SHL,3)),LWB_SHL_Mean,color=SHLColor)
plot(ax4,linspace(-TimeStartW,size(LWB_N60,3)/Param.Fs,size(LWB_N60,3)),LWB_N60_Mean,color=N60Color)
plot(ax4,linspace(-TimeStartW,size(LWB_N70,3)/Param.Fs,size(LWB_N70,3)),LWB_N70_Mean,color=N70Color)

GSW_Mean(isnan(GSW_Mean))=0;GSW_SEM(isnan(GSW_SEM))=0;
SW_Quiet_Mean(isnan(SW_Quiet_Mean))=0;SW_Quiet_SEM(isnan(SW_Quiet_SEM))=0;
SW_SHL_Mean(isnan(SW_SHL_Mean))=0;SW_SHL_SEM(isnan(SW_SHL_SEM))=0;
SW_N60_Mean(isnan(SW_N60_Mean))=0;SW_N60_SEM(isnan(SW_N60_SEM))=0;
SW_N70_Mean(isnan(SW_N70_Mean))=0;SW_N70_SEM(isnan(SW_N70_SEM))=0;

GLW_Mean(isnan(GLW_Mean))=0;GLW_SEM(isnan(GLW_SEM))=0;
LW_Quiet_Mean(isnan(LW_Quiet_Mean))=0;LW_Quiet_SEM(isnan(LW_Quiet_SEM))=0;
LW_SHL_Mean(isnan(LW_SHL_Mean))=0;LW_SHL_SEM(isnan(LW_SHL_SEM))=0;
LW_N60_Mean(isnan(LW_N60_Mean))=0;LW_N60_SEM(isnan(LW_N60_SEM))=0;
LW_N70_Mean(isnan(LW_N70_Mean))=0;LW_N70_SEM(isnan(LW_N70_SEM))=0;

GSWB_Mean(isnan(GSWB_Mean))=0;GSWB_SEM(isnan(GSWB_SEM))=0;
SWB_Quiet_Mean(isnan(SWB_Quiet_Mean))=0;SWB_Quiet_SEM(isnan(SWB_Quiet_SEM))=0;
SWB_SHL_Mean(isnan(SWB_SHL_Mean))=0;SWB_SHL_SEM(isnan(SWB_SHL_SEM))=0;
SWB_N60_Mean(isnan(SWB_N60_Mean))=0;SWB_N60_SEM(isnan(SWB_N60_SEM))=0;
SWB_N70_Mean(isnan(SWB_N70_Mean))=0;SWB_N70_SEM(isnan(SWB_N70_SEM))=0;

GLWB_Mean(isnan(GLWB_Mean))=0;GLWB_SEM(isnan(GLWB_SEM))=0;
LWB_Quiet_Mean(isnan(LWB_Quiet_Mean))=0;LWB_Quiet_SEM(isnan(LWB_Quiet_SEM))=0;
LWB_SHL_Mean(isnan(LWB_SHL_Mean))=0;LWB_SHL_SEM(isnan(LWB_SHL_SEM))=0;
LWB_N60_Mean(isnan(LWB_N60_Mean))=0;LWB_N60_SEM(isnan(LWB_N60_SEM))=0;
LWB_N70_Mean(isnan(LWB_N70_Mean))=0;LWB_N70_SEM(isnan(LWB_N70_SEM))=0;

fill(ax1,[linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3)), flipud(linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3))')'],[(GSW_Mean+GSW_SEM), flipud((GSW_Mean-GSW_SEM)')'],SColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_Quiet,3)/Param.Fs,size(SW_Quiet,3)), flipud(linspace(-TimeStartW,size(SW_Quiet,3)/Param.Fs,size(SW_Quiet,3))')'],[(SW_Quiet_Mean+SW_Quiet_SEM), flipud((SW_Quiet_Mean-SW_Quiet_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_SHL,3)/Param.Fs,size(SW_SHL,3)), flipud(linspace(-TimeStartW,size(SW_SHL,3)/Param.Fs,size(SW_SHL,3))')'],[(SW_SHL_Mean+SW_SHL_SEM), flipud((SW_SHL_Mean-SW_SHL_SEM)')'],SHLColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_N60,3)/Param.Fs,size(SW_N60,3)), flipud(linspace(-TimeStartW,size(SW_N60,3)/Param.Fs,size(SW_N60,3))')'],[(SW_N60_Mean+SW_N60_SEM), flipud((SW_N60_Mean-SW_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[linspace(-TimeStartW,size(SW_N70,3)/Param.Fs,size(SW_N70,3)), flipud(linspace(-TimeStartW,size(SW_N70,3)/Param.Fs,size(SW_N70,3))')'],[(SW_N70_Mean+SW_N70_SEM), flipud((SW_N70_Mean-SW_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')

fill(ax2,[linspace(-TimeStartW,size(GLW,3)/Param.Fs,size(GLW,3)), flipud(linspace(-TimeStartW,size(GLW,3)/Param.Fs,size(GLW,3))')'],[(GLW_Mean+GLW_SEM), flipud((GLW_Mean-GLW_SEM)')'],LColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_Quiet,3)/Param.Fs,size(LW_Quiet,3)), flipud(linspace(-TimeStartW,size(LW_Quiet,3)/Param.Fs,size(LW_Quiet,3))')'],[(LW_Quiet_Mean+LW_Quiet_SEM), flipud((LW_Quiet_Mean-LW_Quiet_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_SHL,3)/Param.Fs,size(LW_SHL,3)), flipud(linspace(-TimeStartW,size(LW_SHL,3)/Param.Fs,size(LW_SHL,3))')'],[(LW_SHL_Mean+LW_SHL_SEM), flipud((LW_SHL_Mean-LW_SHL_SEM)')'],SHLColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_N60,3)/Param.Fs,size(LW_N60,3)), flipud(linspace(-TimeStartW,size(LW_N60,3)/Param.Fs,size(LW_N60,3))')'],[(LW_N60_Mean+LW_N60_SEM), flipud((LW_N60_Mean-LW_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[linspace(-TimeStartW,size(LW_N70,3)/Param.Fs,size(LW_N70,3)), flipud(linspace(-TimeStartW,size(LW_N70,3)/Param.Fs,size(LW_N70,3))')'],[(LW_N70_Mean+LW_N70_SEM), flipud((LW_N70_Mean-LW_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')

fill(ax3,[linspace(-TimeStartW,size(GSWB,3)/Param.Fs,size(GSWB,3)), flipud(linspace(-TimeStartW,size(GSWB,3)/Param.Fs,size(GSWB,3))')'],[(GSWB_Mean+GSWB_SEM), flipud((GSWB_Mean-GSWB_SEM)')'],SColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SWB_Quiet,3)/Param.Fs,size(SW_Quiet,3)), flipud(linspace(-TimeStartW,size(SWB_Quiet,3)/Param.Fs,size(SWB_Quiet,3))')'],[(SWB_Quiet_Mean+SWB_Quiet_SEM), flipud((SWB_Quiet_Mean-SWB_Quiet_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SWB_SHL,3)/Param.Fs,size(SW_SHL,3)), flipud(linspace(-TimeStartW,size(SWB_SHL,3)/Param.Fs,size(SWB_SHL,3))')'],[(SWB_SHL_Mean+SWB_SHL_SEM), flipud((SWB_SHL_Mean-SWB_SHL_SEM)')'],SHLColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SWB_N60,3)/Param.Fs,size(SW_N60,3)), flipud(linspace(-TimeStartW,size(SWB_N60,3)/Param.Fs,size(SWB_N60,3))')'],[(SWB_N60_Mean+SWB_N60_SEM), flipud((SWB_N60_Mean-SWB_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[linspace(-TimeStartW,size(SWB_N70,3)/Param.Fs,size(SW_N70,3)), flipud(linspace(-TimeStartW,size(SWB_N70,3)/Param.Fs,size(SWB_N70,3))')'],[(SWB_N70_Mean+SWB_N70_SEM), flipud((SWB_N70_Mean-SWB_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')

fill(ax4,[linspace(-TimeStartW,size(GLWB,3)/Param.Fs,size(GLWB,3)), flipud(linspace(-TimeStartW,size(GLWB,3)/Param.Fs,size(GLWB,3))')'],[(GLWB_Mean+GLWB_SEM), flipud((GLWB_Mean-GLWB_SEM)')'],LColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LWB_Quiet,3)/Param.Fs,size(LWB_Quiet,3)), flipud(linspace(-TimeStartW,size(LWB_Quiet,3)/Param.Fs,size(LWB_Quiet,3))')'],[(LWB_Quiet_Mean+LWB_Quiet_SEM), flipud((LWB_Quiet_Mean-LWB_Quiet_SEM)')'],QuietColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LWB_SHL,3)/Param.Fs,size(LWB_SHL,3)), flipud(linspace(-TimeStartW,size(LWB_SHL,3)/Param.Fs,size(LWB_SHL,3))')'],[(LWB_SHL_Mean+LWB_SHL_SEM), flipud((LWB_SHL_Mean-LWB_SHL_SEM)')'],SHLColor,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LWB_N60,3)/Param.Fs,size(LWB_N60,3)), flipud(linspace(-TimeStartW,size(LWB_N60,3)/Param.Fs,size(LWB_N60,3))')'],[(LWB_N60_Mean+LWB_N60_SEM), flipud((LWB_N60_Mean-LWB_N60_SEM)')'],N60Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[linspace(-TimeStartW,size(LWB_N70,3)/Param.Fs,size(LWB_N70,3)), flipud(linspace(-TimeStartW,size(LWB_N70,3)/Param.Fs,size(LWB_N70,3))')'],[(LWB_N70_Mean+LWB_N70_SEM), flipud((LWB_N70_Mean-LWB_N70_SEM)')'],N70Color,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')

xlabel([ax1 ax2 ax3 ax4],'Time [s]')
ylabel([ax1 ax2],'Fixation duration [s]')
ylabel([ax3 ax4],'Fixation duration difference from baseline [s]')
xlim([ax1 ax2 ax3 ax4],[-TimeStartW 3])
lgd2=legend(ax2,'Speaking','Listening','Quiet','SHL','N60','N70','Location','southeastoutside');
lgd2.Title.String = 'Types of windows:';
lgd4=legend(ax4,'Speaking','Listening','Quiet','SHL','N60','N70','Location','southeastoutside');
lgd4.Title.String = 'Types of windows:';
title(ax1,'Global Speaking-evoked response')
title(ax2,'Global Listening-evoked response')
title(ax3,'Global adaptive baselined Speaking-evoked response')
title(ax4,'Global adaptive baselined Listening-evoked response')