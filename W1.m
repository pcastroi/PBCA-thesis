%% PBCA-Thesis-W1
% Pathing
% clear all; clc; close all;
BPath = strsplit(pwd,'PBCA-thesis');
addpath(genpath('data\'));
addpath([BPath{1} 'Pupil-preprocessing-tools\tools']) % For preprocessing

% Ask for pair number (pair of test subjects)
while true
    PairIn = str2double(input('Enter participant pair (from 1 to 12): ','s'));
    if fix(PairIn) == PairIn && PairIn >= 1 && PairIn <= 12
        break;
    end
    disp('Number must be an integer between 1 and 12.');
end

% Files and Utterances: different conditions
PairFiles=dir(['data\Main',sprintf('%d',PairIn),'\*.mat']);
PairUtt=load('data\utterances1110.mat');
PairUtt=PairUtt.Utterances(PairIn,:);

% Parameters
Param.Fs = 50; % Sampling frequency of pupil data
Param.RemoveBeforeAndAfter = [35 100]*1e-3; % Samples within the time range before and after NaNs will set NaNs as well.
Param.MinLengthNaNRepair = 5; % Drop values (i.e., change to NaN) before and after NaNs only for contiguous NaNs of at least __ samples. 
LPWinSize = 0.75; % [s]: Window size of hamming-window for low-pass filtering
LPWindow = hamming(round(LPWinSize*Param.Fs));
LPWindow = LPWindow/sum(LPWindow); % Hamming-window

% Run once per file
for i=1:numel(PairFiles)
    alldata(i) = load([PairFiles(i).folder, '\', PairFiles(i).name]);
    alldata_mat = cell2mat(alldata(i).data);
    
    % Preprocessing - Interpolating NaNs
    LDiamRaw = [alldata_mat.diameterLeft];
    RDiamRaw = [alldata_mat.diameterRight];
    [LDiam,LMetadata] = preprocpupil(LDiamRaw,Param);
    [RDiam,RMetadata] = preprocpupil(RDiamRaw,Param);
    
    % Low-Pass Filtering
    LDiam = conv(LDiam,LPWindow,'same');
    RDiam = conv(RDiam,LPWindow,'same');
    
    % Decide 'better' eye results
    [Min idx_decision] = min([sum(LMetadata.Isnan) sum(RMetadata.Isnan)]);
    if idx_decision == 1
        Diameter = LDiam;
        DiameterRaw = LDiamRaw';
        eyeChosen = 'Left';
    elseif idx_decision == 2
        Diameter = RDiam;
        DiameterRaw = RDiamRaw';
        eyeChosen = 'Right';
    end
    
    % Retrieve Utterances
    if contains(PairFiles(i).name,'P1')
        SpeakKey = 'utteranceCH1';
        ListenKey = 'utteranceCH2';
    elseif contains(PairFiles(i).name,'P2')
        SpeakKey = 'utteranceCH2';
        ListenKey = 'utteranceCH1';
    end

    if contains(PairFiles(i).name,'B1')
        UttB = 0;
    elseif contains(PairFiles(i).name,'B2')
        UttB = 1;
    end
    
    if contains(PairFiles(i).name,'Noise60')
        UttCond = UttB + 1;
    elseif contains(PairFiles(i).name,'Noise70')
        UttCond = UttB + 3;
    elseif contains(PairFiles(i).name,'Quiet')
        UttCond = UttB + 5;
    elseif contains(PairFiles(i).name,'SHL')
        UttCond = UttB + 7;
    end
    
    Speak = PairUtt{1,UttCond}.(SpeakKey);
    Listen = PairUtt{1,UttCond}.(ListenKey);
    binRes = PairUtt{1,UttCond}.binRes;
    
    % Downsample (rounding) Utt from 250 Hz (1/binRes) to 50 Hz
    Speak(:,2:3)=round(Speak(:,2:3)*binRes*Param.Fs);
    Listen(:,2:3)=round(Listen(:,2:3)*binRes*Param.Fs);
    
    % Time-locked indexes (based on Start or End of events)
    TimeStartW=1; % [S], time previous to Utt/Lis starts
    TimeEndW=2; % [S], time after Utt/Lis starts
    SWSpeakIdx=[Speak(:,2)-TimeStartW*Param.Fs,Speak(:,2),Speak(:,2)+TimeEndW*Param.Fs];
    SWListenIdx=[Listen(:,2)-TimeStartW*Param.Fs,Listen(:,2),Listen(:,2)+TimeEndW*Param.Fs];
    EWSpeakIdx=[Speak(:,3)-TimeStartW*Param.Fs,Speak(:,3),Speak(:,3)+TimeEndW*Param.Fs];
    EWListenIdx=[Listen(:,3)-TimeStartW*Param.Fs,Listen(:,3),Listen(:,3)+TimeEndW*Param.Fs];
    
    % Time vectors for plotting
    t_LDiam = linspace(1,length(LDiam)./Param.Fs,length(LDiam));
    t_RDiam = linspace(1,length(RDiam)./Param.Fs,length(RDiam));
    
    % Plots
    
    figure
    subplot(2,2,[1 2])
    Lplot=plot(t_LDiam,LDiam,color='blue');
    hold on
    Rplot=plot(t_RDiam,RDiam,color='red');
    yl=ylim();
    ylim(yl);
    startStopU = t_LDiam(Speak(:,2:3)); 
    widthU = startStopU(:,2)-startStopU(:,1);
    startStopL = t_LDiam(Listen(:,2:3)); 
    widthL = startStopL(:,2)-startStopL(:,1);
    % Plot rectangles (Utterance and listening time windows)
    arrayfun(@(i)rectangle('Position', [startStopU(i,1),yl(1),widthU(i),range(yl)], ...
    'EdgeColor', 'none', 'FaceColor', [0 1 0 .2]), 1:size(startStopU,1))
    arrayfun(@(i)rectangle('Position', [startStopL(i,1),yl(1),widthL(i),range(yl)], ...
    'EdgeColor', 'none', 'FaceColor', [1 0 1 .2]), 1:size(startStopL,1))
    if ~isempty(Speak)
        eline1=line(NaN,NaN,'LineWidth',2,'Color',[0 1 0 .2]);
        eline2=line(NaN,NaN,'LineWidth',2,'Color',[1 0 1 .2]);
        legend([Lplot Rplot eline1 eline2],{'Diameter Left','Diameter Right','Utterance windows','Listening windows'})
    else
        disp(['Warning: No associated Utterance/Listening data for file ', PairFiles(i).name, '.']);
        legend([Lplot Rplot],{'Diameter Left','Diameter Right'})
    end
    sgtitle(strrep(PairFiles(i).name,'_','-'))
    xticks([0:30:t_LDiam(end)])
    xlabel('Time [s]')
    ylabel('Pupil diameter [mm]')
    
    subplot(2,2,3)
    % Start/End Diameter Sum (for averaging)
    SDSum=zeros(1,Param.Fs*(TimeEndW+TimeStartW)); 
    EDSum=SDSum;
    if ~isempty(Speak)
        hold on
        for j=1:size(Speak,1)
            plot(linspace(-TimeStartW,TimeEndW,Param.Fs*(TimeEndW+TimeStartW)),...
                 Diameter(SWSpeakIdx(j,1):SWSpeakIdx(j,3)-1),color=[0 1 0 .2],LineWidth=0.3)
            SDSum=SDSum+Diameter(SWSpeakIdx(j,1):SWSpeakIdx(j,3)-1);
        end
        plot(linspace(-TimeStartW,TimeEndW,Param.Fs*(TimeEndW+TimeStartW)),...
             SDSum./size(Speak,1),color=[1 0 0 0.8],LineWidth=1)
        xline(0,"--")
        xticks([-TimeStartW:1:TimeEndW])
        title('Utterance-evoked pupil response')
        xlabel('Time [s]')
        ylabel('Pupil diameter [mm]')
    end
    
    subplot(2,2,4)
    if ~isempty(Listen)
        hold on
        for j=1:size(Listen,1)
            plot(linspace(-TimeStartW,TimeEndW,Param.Fs*(TimeEndW+TimeStartW)),...
                 Diameter(SWListenIdx(j,1):SWListenIdx(j,3)-1),color=[1 0 1 0.2],LineWidth=0.3)
             EDSum=EDSum+Diameter(SWListenIdx(j,1):SWListenIdx(j,3)-1);
        end
        plot(linspace(-TimeStartW,TimeEndW,Param.Fs*(TimeEndW+TimeStartW)),...
             EDSum./size(Listen,1),color=[1 0 0 0.8],LineWidth=1)
        xline(0,"--")
        xticks([-TimeStartW:1:TimeEndW])
        title('Listening-evoked pupil response')
        xlabel('Time [s]')
        ylabel('Pupil diameter [mm]')
    end
    
end