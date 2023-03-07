%% PBCA-Thesis - Week 9 - Fixation duration? Not using diameter anymore
% Pathing
clear all; clc; close all;
BPath = strsplit(pwd,'PBCA-thesis');
addpath('tools')
addpath([BPath{1} 'Pupil-preprocessing-tools\tools']) % For preprocessing

[subDirs] = GetSubDirsFirstLevelOnly('data');
FileNames={'P1_Quiet_B1.mat','P1_Quiet_B2.mat','P1_SHL_B1.mat','P1_SHL_B2.mat','P1_Noise60_B1.mat','P1_Noise60_B2.mat','P1_Noise70_B1.mat','P1_Noise70_B2.mat','P2_Quiet_B1.mat','P2_Quiet_B2.mat','P2_SHL_B1.mat','P2_SHL_B2.mat','P2_Noise60_B1.mat','P2_Noise60_B2.mat','P2_Noise70_B1.mat','P2_Noise70_B2.mat'};
LoadUtt=load('data\utterances1110.mat');
LoadDelays=load('data\delays1110.mat');

% Parameters for processing
Param.Fs = 50; % Sampling frequency of pupil data
Param.Preblink = 0.1; % [s], set to NaN, time before blink
Param.Postblink = 0.2; % [s], set to NaN, time after blink
Param.BlinkThresh = 3; % [samples], threshold of samples in between artifacts or blinks
MaxVisualAngle = 0.5; % [degrees], critical visual angle for fixation definition

TimeStartW = 0.5; % [s], time before Utt/Lis starts
TimeEndW = 0; % [s], time after Utt/Lis starts
TimeStart = 20; % [s], time at which simultaneous recording started
TimeMinWin = 0.5; % [s], Minimum time of a window
TimeInitialMerge = 0.3; % [s], Time threshold for merging windows initially
TimeMerge = 2; % [s], Time threshold for merging windows after rejecting small windows
RejectRatio = 0.4; % Rejection threshold based on the ratio of NaNs in data
RejectDelay = 0.5; % [s], Rejection threshold based on delay between timestamps and n-samples

for q=1:numel(subDirs)
    PairIn = q;
    PairFiles=dir(['data\Main',sprintf('%d',PairIn),'\*.mat']);
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
        
        gaze3d = vertcat(alldata_mat.gaze3d);
        
        GazeXRaw = gaze3d(:,1);
        GazeYRaw = gaze3d(:,2);
        GazeZRaw = gaze3d(:,3);
        
        GazeX = GazeXRaw;
        GazeY = GazeYRaw;
        GazeZ = GazeZRaw;
        
        % NaN's appear at the same time at [X,Y,Z], we only look at X.
        XNan = find(isnan(GazeXRaw));
        
        % Reject a file if the chosen eye data has too many NaNs
        if length(XNan)/length(gaze3d) >= RejectRatio
            disp(['Warning: File ',PairFiles(1).folder, '\', cell2mat(FileNames(i)), ' was rejected because it contains too many NaNs (',sprintf('%0.2f',100*length(XNan)/length(gaze3d)),'%).'])
            continue
        end
        
        % Extract delay [s] (Duration_timeStamps - Duration_nsamples)
        EyeAudDelay=alldata_mat(end).timeStamp-alldata_mat(1).timeStamp-length(alldata_mat)/Param.Fs;
        
        % Skip file if the difference in duration from the number of
        % samples and the duration given by timestamps is bigger than 0.5 s
        if EyeAudDelay > RejectDelay
            disp(['Warning: File ',PairFiles(1).folder, '\', cell2mat(FileNames(i)), ' was rejected, too much delay (',sprintf('%0.2f',EyeAudDelay),'s).']);
            continue
        end
        
        % Preprocessing - 100 ms before and 200 ms after a blink -> NaN
        if ~isempty(XNan)
            
            % Pre-blink
            Pre = [XNan(1);XNan(find(diff(XNan,1) > 1) + 1)];
            for h=1:length(Pre)
                if Pre(h) > Param.Preblink*Param.Fs % Check start delimiter
                   for k=1:Param.Preblink*Param.Fs
                       GazeX(Pre(h)-k)=NaN;
                       GazeY(Pre(h)-k)=NaN;
                       GazeZ(Pre(h)-k)=NaN;
                   end
                end
            end

            % Post-blink
            Post = XNan(diff(XNan,1) > Param.BlinkThresh);
            for h=1:length(Post)
                if Post(h) < length(GazeXRaw) + Param.Postblink*Param.Fs % Check end delimiter
                   for k=1:Param.Postblink*Param.Fs
                       GazeX(Post(h)+k)=NaN;
                       GazeY(Post(h)+k)=NaN;
                       GazeZ(Post(h)+k)=NaN;
                   end
                end
            end
        end
        % Fixation duration
        fixationRaw = fix_duration([GazeXRaw,GazeYRaw,GazeZRaw],MaxVisualAngle,Param.Fs);
        fixation = fix_duration([GazeX,GazeY,GazeZ],MaxVisualAngle,Param.Fs);
        
        % Time-locked indexes (based on Start or End of events)
        SWSpeakIdx=[Speak(:,2)-TimeStartW*Param.Fs,Speak(:,2),Speak(:,2)+TimeEndW*Param.Fs];
        SWListenIdx=[Listen(:,2)-TimeStartW*Param.Fs,Listen(:,2),Listen(:,2)+TimeEndW*Param.Fs];
        EWSpeakIdx=[Speak(:,3)-TimeStartW*Param.Fs,Speak(:,3),Speak(:,3)+TimeEndW*Param.Fs];
        EWListenIdx=[Listen(:,3)-TimeStartW*Param.Fs,Listen(:,3),Listen(:,3)+TimeEndW*Param.Fs];
        
        % Plots
        figure;tiledlayout(1,2);ax1 = nexttile;ax2 = nexttile;
        hold([ax1 ax2],'on')
        SW = zeros(size(Speak,1),ceil(Param.Fs*(TimeStartW+max(Speak(:,1)))));
        for j=1:size(Speak,1)
            if SWSpeakIdx(j,3)-1 <= length(fixationRaw)
                plot(ax1,linspace(-TimeStartW,Speak(j,1),length(SWSpeakIdx(j,1):Speak(j,3)-1)),fixationRaw(SWSpeakIdx(j,1):Speak(j,3)-1),color=[0 0 0 0.01],linewidth=0.5)
                % Add nan-padding when necessary
                SW(j,:)=[fixationRaw(SWSpeakIdx(j,1):Speak(j,3)-1);NaN*ones(1,length(SW)-length(SWSpeakIdx(j,1):Speak(j,3)-1))'];
            end
        end
        SW_Mean = mean(SW,1);
        SW_SEM = std(SW,[],1)/sqrt(size(SW,1));
        plot(ax1,linspace(-TimeStartW,max(Speak(:,1)),length(SW)),mean(SW,1,'omitnan'),color=[1 0 0 0.8],LineWidth=1.5)
        fill(ax1,[linspace(-TimeStartW,max(Speak(:,1)),length(SW)), flipud(linspace(-TimeStartW,max(Speak(:,1)),length(SW)))'],[(SW_Mean+SW_SEM), flipud((SW_Mean-SW_SEM)')'],[154 0 0]./255,'FaceAlpha',.5,'Edgecolor','none','handlevisibility' ,'off')
        xline(ax1,0,"--")
        
        LW = zeros(size(Listen,1),ceil(Param.Fs*(TimeStartW+max(Listen(:,1)))));
        for j=1:size(Listen,1)
            if SWListenIdx(j,3)-1 <= length(fixationRaw)
                plot(ax2,linspace(-TimeStartW,Listen(j,1),length(SWListenIdx(j,1):Listen(j,3)-1)),fixationRaw(SWListenIdx(j,1):Listen(j,3)-1),color=[0 0 0 0.01],linewidth=0.5)
                % Add nan-padding when necessary
                LW(j,:)=[fixationRaw(SWListenIdx(j,1):Listen(j,3)-1);NaN*ones(1,length(LW)-length(SWListenIdx(j,1):Listen(j,3)-1))'];
            end
        end
        plot(ax2,linspace(-TimeStartW,max(Listen(:,1)),length(LW)),mean(LW,1,'omitnan'),color=[1 0 0 0.8],LineWidth=1.5)
        xline(ax2,0,"--")
        
%         xticks(ax1,[-TimeStartW:TimeStartW:max(Speak(:,1))])
%         xticks(ax2,[-TimeStartW:TimeStartW:max(Listen(:,1))])
        title(ax1,'Fixation duration during Speaking windows')
        title(ax2,'Baselined listening-evoked pupil response')
        xlabel([ax1 ax2],'Time [s]')
        ylabel([ax1 ax2],'Fixation duration [s]');
        
        
%         figure;plot(linspace(0,length(fixation)./Param.Fs,length(fixation)),fixation,'LineWidth',2);hold on;plot(linspace(0,length(fixationRaw)./Param.Fs,length(fixationRaw)),fixationRaw);title([PairFiles(1).folder, '\', cell2mat(FileNames(i))],'interpreter','none')
     
%         backimag = uint8(zeros(600,600,3));         % backimag is the background image we have; let's put black
%         [dimY, dimX, ~] = size(backimag);
%         X = round(rescale(GazeX,1,dimX));           % GazeX from [-1,1] -> convert to [1,600]
%         Y = round(rescale(GazeY,1,dimY));           % GazeY from [-1,1] -> convert to [1,600]
%         W = ones(size(GazeX));                      % we could have weights on each gaze (eg to normalise for different number of trials per subject)
%         heatmap = mat2cell(hot,256,[1 1 1]);        % we choose a heatmap colouring, eg "hot", and convert it to "cell"
%         sigma = 20;                                  % the variance parameter for the gaussian kernel
%         % Create "mask"
%         origmask = ones(dimX, dimY)*0.1;
%         for k = 1:size(GazeX,1)
%             origmask(Y(k), X(k)) = origmask(Y(k), X(k)) + W(k);
%         end
%         % Filter using a gaussian kernel
%         mask = imgaussfilt(origmask, sigma);
%         % Normalise total mass of heatmap
%         mask = rescale(mask);
%         % Colour the background image with the heatmap
%         newImage = backimag;
%         for rgbInd = 1:3
%             thisHeat = heatmap{rgbInd}( floor(mask*255) + 1 );
%             newImage(:,:,rgbInd) = (newImage(:,:,rgbInd) + uint8(thisHeat*255));
%         end
%         figure; imshow(newImage); set(gca, 'ydir', 'normal')
%         figure; imshow(origmask); set(gca, 'ydir', 'normal')
%         figure; plot(GazeX,GazeY); xlim([-1 1]); ylim([-1 1])
        
%         [alldata_mat(:,cellfun(@(xd) any(isnan(xd)),{alldata_mat.gaze3d})).gaze3d] % all non-nans in gaze3d
               
    end
end
