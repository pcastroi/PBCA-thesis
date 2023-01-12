%% PBCA-Thesis - Week 1&2 - Processing existing pupil data
% Pathing
clear all; clc; close all;
BPath = strsplit(pwd,'PBCA-thesis');
addpath('tools')
addpath([BPath{1} 'Pupil-preprocessing-tools\tools']) % For preprocessing

% Ask for pair number (pair of test subjects)
% while true
%     PairIn = str2double(input('Enter participant pair (from 1 to 12): ','s'));
%     if fix(PairIn) == PairIn && PairIn >= 1 && PairIn <= 12
%         break;
%     end
%     disp('Number must be an integer between 1 and 12.');
% end

[subDirs] = GetSubDirsFirstLevelOnly('data');
LoadUtt=load('data\utterances1110.mat');

% Parameters for processing
Param.Fs = 50; % Sampling frequency of pupil data
Param.RemoveBeforeAndAfter = [35 100]*1e-3; % Samples within the time range before and after NaNs will set NaNs as well.
Param.MinLengthNaNRepair = 5; % Drop values (i.e., change to NaN) before and after NaNs only for contiguous NaNs of at least __ samples. 
LPWinSize = 1; % [s]: Window size of hamming-window for low-pass filtering
LPWindow = hamming(round(LPWinSize*Param.Fs));
LPWindow = LPWindow/sum(LPWindow); % Hamming-window
TimeStartW=1; % [s], time before Utt/Lis starts
TimeEndW=2; % [s], time after Utt/Lis starts
TimeStart=20; % [s], time at which simultaneous recording started
TimeBL=[10,15]; % [s], time chosen for baseline
RejectRatio = 0.3; % Limit of rejection (ratio of NaNs in data)

% Run through all subfolders once
for q=1:numel(subDirs)
    PairIn = q;
    
    % Files and Utterances: different conditions
    PairFiles=dir(['data\Main',sprintf('%d',PairIn),'\*.mat']);
    PairUtt=LoadUtt.Utterances(PairIn,:);
    
    % Parameters
    AvgLPupilSize = zeros(1,numel(PairFiles)); % [mm], avg Listening pupil size (per file)
    AvgLPupilSlope = zeros(1,numel(PairFiles)); % [mm/s], avg Listening baselined pupil slope (per file)
    AvgSPupilSize = zeros(1,numel(PairFiles)); % [mm], avg Speaking pupil size (per file)
    AvgSPupilSlope = zeros(1,numel(PairFiles)); % [mm/s], avg Speaking baselined pupil slope (per file)
    AvgTimeSpe = zeros(1,numel(PairFiles)); % [s], avg time of Utterance (per file)
    AvgTimeLis = zeros(1,numel(PairFiles)); % [s], avg time of Listening (per file)
    
    % Run once per file
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
        if alldata_mat(end).timeStamp-alldata_mat(1).timeStamp-length(alldata_mat)/Param.Fs > 0.5
            disp(['Warning: File ',PairFiles(i).folder, '\', PairFiles(i).name, ' was rejected, too much delay (',sprintf('%0.5f',EyeAudDelay),'s total).']);
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

        % Low-Pass Filtering ('valid' to avoid peak/dip at start/finish)
        LDiam = conv(LDiam,LPWindow,'valid'); 
        RDiam = conv(RDiam,LPWindow,'valid');

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
        
        % Reject a file if the chosen eye data has too many NaNs
        if DiamNaN/length(Diameter) >= RejectRatio
            disp(['Warning: File ',PairFiles(i).folder, '\', PairFiles(i).name, ' was rejected because it contains too many NaNs (',num2str(100*DiamNaN/length(Diameter)),'%).'])
            continue
        end
        % Baseline the Diameter (from start to finish)
        BLDiam = mean(Diameter(TimeBL(1)*Param.Fs:TimeBL(2)*Param.Fs))-Diameter;

        % Retrieve Utterances
        if contains(PairFiles(i).name,'P1')
            SpeakKey = 'utteranceCH1';
            ListenKey = 'utteranceCH2';
        elseif contains(PairFiles(i).name,'P2')
            SpeakKey = 'utteranceCH2';
            ListenKey = 'utteranceCH1';
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

        Speak = PairUtt{1,SpeCond}.(SpeakKey);
        Listen = PairUtt{1,SpeCond}.(ListenKey);
        binRes = PairUtt{1,SpeCond}.binRes;

        % Downsample (rounding) Utt from 250 Hz (1/binRes) to 50 Hz, shift
        % in time to account for the time at which the audio recording
        % started (from 0 to 20 s only eye data) and add the relative
        % delay from eye tracker and audio (changes from file to file)
        Speak(:,2:3)=round((Speak(:,2:3)*binRes+TimeStart+EyeAudDelay)*Param.Fs);
        Listen(:,2:3)=round((Listen(:,2:3)*binRes+TimeStart+EyeAudDelay)*Param.Fs);

        % Time-locked indexes (based on Start or End of events)
    %     SWSpeakIdx=[Speak(:,2)-TimeStartW*Param.Fs,Speak(:,2),Speak(:,2)+TimeEndW*Param.Fs];
    %     SWListenIdx=[Listen(:,2)-TimeStartW*Param.Fs,Listen(:,2),Listen(:,2)+TimeEndW*Param.Fs];
    %     EWSpeakIdx=[Speak(:,3)-TimeStartW*Param.Fs,Speak(:,3),Speak(:,3)+TimeEndW*Param.Fs];
    %     EWListenIdx=[Listen(:,3)-TimeStartW*Param.Fs,Listen(:,3),Listen(:,3)+TimeEndW*Param.Fs];

        % Time vectors and idx for plotting
        t_Diam = linspace(1,length(BLDiam)./Param.Fs,length(BLDiam));
        startStopU = t_Diam(Speak(:,2:3)); 
        widthU = startStopU(:,2)-startStopU(:,1);
        startStopL = t_Diam(Listen(:,2:3)); 
        widthL = startStopL(:,2)-startStopL(:,1);

        % Features (of each window): mean slope, mean pupil size
        % Raw
        % Speaking
        if ~isempty(Speak)
            for j=1:size(Speak,1)
                AvgSPupilSize(i)=AvgSPupilSize(i)+mean(Diameter(Speak(j,2):Speak(j,3)));
                AvgSPupilSlope(i)=AvgSPupilSlope(i)+...
                    mean(diff(Diameter(Speak(j,2):Speak(j,3)))./diff(t_Diam(Speak(j,2):Speak(j,3)))');
            end
            AvgSPupilSize(i)=AvgSPupilSize(i)./j;
            AvgSPupilSlope(i)=AvgSPupilSlope(i)./j;
        else
            AvgSPupilSize(i)=NaN;
            AvgSPupilSlope(i)=NaN;
        end
        % Listening
        if ~isempty(Listen)
            for j=1:size(Listen,1)
                AvgLPupilSize(i)=AvgLPupilSize(i)+mean(Diameter(Listen(j,2):Listen(j,3)));
                AvgLPupilSlope(i)=AvgLPupilSlope(i)+...
                    mean(diff(Diameter(Listen(j,2):Listen(j,3)))./diff(t_Diam(Listen(j,2):Listen(j,3)))');
            end
            AvgLPupilSize(i)=AvgLPupilSize(i)./j;
            AvgLPupilSlope(i)=AvgLPupilSlope(i)./j;
        else
            AvgLPupilSize(i)=NaN;
            AvgLPupilSlope(i)=NaN;
        end

        % Average duration of utterance per file
        AvgTimeSpe(i) = mean(Speak(:,1)); % Possible NaNs
        AvgTimeLis(i)= mean(Listen(:,1)); % Possible NaNs

        % Plots
        figure
    %     subplot(2,2,[1 2])
        plot(t_Diam,Diameter,color='black');
        hold on
        xline(TimeStart,"--",'HandleVisibility','off')
        yl=ylim();
        ylim(yl);
        % Plot rectangles (Utterance and listening time windows)
        arrayfun(@(i)rectangle('Position', [startStopU(i,1),yl(1),widthU(i),range(yl)], ...
        'EdgeColor', 'none', 'FaceColor', [0 1 0 .2]), 1:size(startStopU,1))
        arrayfun(@(i)rectangle('Position', [startStopL(i,1),yl(1),widthL(i),range(yl)], ...
        'EdgeColor', 'none', 'FaceColor', [1 0 1 .2]), 1:size(startStopL,1))
        if ~isempty(Speak)
            eline1=line(NaN,NaN,'LineWidth',2,'Color',[0 1 0 .2]);
            eline2=line(NaN,NaN,'LineWidth',2,'Color',[1 0 1 .2]);
            legend(['Baselined diameter (', eyeChosen,' Eye)'],'Utterance windows','Listening windows')
        else
            disp(['Warning: No associated Utterance/Listening data for file ', PairFiles(i).folder, '\', PairFiles(i).name, '.']);
            legend(['Baselined diameter (', eyeChosen,' Eye)'])
        end
        sgtitle(strrep(PairFiles(i).name,'_','-'))
        xticks([0:TimeStart:t_Diam(end)])
        xlabel('Time [s]')
        ylabel('Pupil baseline difference [mm]');

    %     subplot(2,2,3)
    %     % Start/End Diameter Sum (for averaging)
    %     SDWSum = zeros(1,Param.Fs*(TimeEndW+TimeStartW)); % Sum Speak-Diam-Window
    %     k=0;
    %     if ~isempty(Speak)
    %         hold on
    %         for j=1:size(Speak,1)
    %             if SWSpeakIdx(j,3)-1 <= length(BLDiam)
    %                 plot(linspace(-TimeStartW,TimeEndW,Param.Fs*(TimeEndW+TimeStartW)),...
    %                     BLDiam(SWSpeakIdx(j,1):SWSpeakIdx(j,3)-1),color=[0 1 0 .2],LineWidth=0.3)
    %                 k=k+1;
    %                 SDWSum=SDWSum+BLDiam(SWSpeakIdx(j,1):SWSpeakIdx(j,3)-1);
    %             end
    %         end
    %         plot(linspace(-TimeStartW,TimeEndW,Param.Fs*(TimeEndW+TimeStartW)),...
    %              SDWSum./k,color=[1 0 0 0.8],LineWidth=1.5)
    %         xline(0,"--")
    %         xticks([-TimeStartW:1:TimeEndW])
    %         title('Baselined utterance-evoked pupil response')
    %         xlabel('Time [s]')
    %         ylabel('Pupil baseline difference [mm]');
    %     end
    %     
    %     subplot(2,2,4)
    %     LDWSum = zeros(1,Param.Fs*(TimeEndW+TimeStartW)); % Sum Listen-Diam-Window
    %     k=0;
    %     if ~isempty(Listen)
    %         hold on
    %         for j=1:size(Listen,1)
    %             if SWListenIdx(j,3)-1 <= length(BLDiam)
    %                 plot(linspace(-TimeStartW,TimeEndW,Param.Fs*(TimeEndW+TimeStartW)),...
    %                      BLDiam(SWListenIdx(j,1):SWListenIdx(j,3)-1),color=[1 0 1 0.2],LineWidth=0.3)
    %                 k=k+1;
    %                 LDWSum=LDWSum+BLDiam(SWListenIdx(j,1):SWListenIdx(j,3)-1);
    %             end
    %         end
    %         plot(linspace(-TimeStartW,TimeEndW,Param.Fs*(TimeEndW+TimeStartW)),...
    %              LDWSum./k,color=[1 0 0 0.8],LineWidth=1.5)
    %         xline(0,"--")
    %         xticks([-TimeStartW:1:TimeEndW])
    %         title('Baselined listening-evoked pupil response')
    %         xlabel('Time [s]')
    %         ylabel('Pupil baseline difference [mm]');
    %     end

    end
end

% Global parameters from Utterance windows 
GlobalAvgLPupilSize = mean(AvgLPupilSize,'omitnan'); % Avg Listening Diameter
GlobalAvgLPupilSlope = mean(AvgLPupilSlope,'omitnan'); % Avg Listening Slope
GlobalAvgSPupilSize = mean(AvgSPupilSize,'omitnan'); % Avg Speaking Diameter
GlobalAvgSPupilSlope = mean(AvgSPupilSlope,'omitnan'); % Avg Speaking Slope
GlobalAvgTimeUtt = mean(AvgTimeSpe,'omitnan'); % Avg duration of Utt
GlobalAvgTimeLis = mean(AvgTimeLis,'omitnan'); % Avg duration of Lis
