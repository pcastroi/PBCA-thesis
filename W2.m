%% PBCA-Thesis - Week 2 - Making sure audios and Utterances are synchronized
% Pathing
clear all; clc; close all;
BPath = strsplit(pwd,'PBCA-thesis');
addpath(genpath('audio\'));
addpath([BPath{1} 'Pupil-preprocessing-tools\tools']) % For preprocessing

% Files and Utterances: different conditions
PairIn=1; % Audios from Main1
AudFiles=dir('audio\*.wav');
PairUtt=load('data\utterances1110.mat');
PairUtt=PairUtt.Utterances(PairIn,:);

% Variables
Param.Fs = 50;

% Run once per file
for i=1:numel(AudFiles)
    [AudData,AudFs]=audioread(AudFiles(i).name);
    if contains(AudFiles(i).name,'talker1')
        SpeakKey = 'utteranceCH1';
    elseif contains(AudFiles(i).name,'talker2')
        SpeakKey = 'utteranceCH2';
    end

    if contains(AudFiles(i).name,'Rep1')
        UttB = 0;
    elseif contains(AudFiles(i).name,'Rep2')
        UttB = 1;
    end
    
    if contains(AudFiles(i).name,'SpeechNH')
        if contains(AudFiles(i).name,'Quiet')
            UttCond = UttB + 1;
        elseif contains(AudFiles(i).name,'Noise60')
            UttCond = UttB + 5;
        elseif contains(AudFiles(i).name,'Noise70')
            UttCond = UttB + 7;
        end
    elseif contains(AudFiles(i).name,'SpeechSHL')
        UttCond = UttB + 3;
    end
    
    Speak = PairUtt{1,UttCond}.(SpeakKey);
    binRes = PairUtt{1,UttCond}.binRes;
    
    % Upsample (rounding) Utt from 250 Hz (1/binRes) to 48000 Hz (audio)
    Speak(:,2:3)=round(Speak(:,2:3)*binRes*AudFs);
    
    % Time vectors and idx for plotting
    t_Aud = linspace(0,length(AudData)/AudFs,length(AudData));
    startStopU = t_Aud(Speak(:,2:3)); 
    widthU = startStopU(:,2)-startStopU(:,1);
    
    % Plots
    figure
    plot(t_Aud,AudData)
    hold on
    arrayfun(@(i)rectangle('Position', [startStopU(i,1),min(AudData),...
        widthU(i),range([min(AudData), max(AudData)])],'EdgeColor','none',...
        'FaceColor', [0 1 0 .2]), 1:size(startStopU,1))
    xticks([0:30:t_Aud(end)])
    xlabel('Time [s]')
    ylabel('Amplitude')
    title(strrep(AudFiles(i).name,'_','-'))
end