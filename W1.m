%% PBCA-Thesis-W1
% Pathing
clear all; clc; close all;
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
        DiameterRaw = LDiamRaw;
        eyeChosen = 'Left';
    elseif idx_decision == 2
        Diameter = RDiam;
        DiameterRaw = RDiamRaw;
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
    
    SpeakUtt = PairUtt{1,UttCond}.(SpeakKey);
    ListenUtt = PairUtt{1,UttCond}.(ListenKey);
    binRes = PairUtt{1,UttCond}.binRes;
    
    % Downsample (rounding) Utt from 250 Hz (1/binRes) to 50 Hz
    SpeakUtt(:,2:3)=round(SpeakUtt(:,2:3)*binRes*Param.Fs);
    ListenUtt(:,2:3)=round(ListenUtt(:,2:3)*binRes*Param.Fs);
    
    % Arrays for plotting (patch)
    SU = [sort([SpeakUtt(:,2); SpeakUtt(:,3)]), sort([SpeakUtt(:,2); SpeakUtt(:,3)])];
    LU = [sort([ListenUtt(:,2); ListenUtt(:,3)]), sort([ListenUtt(:,2); ListenUtt(:,3)])];
    
    % Time vectors for plotting
    t_LDiam = linspace(1,length(LDiam)./Param.Fs,length(LDiam));
    t_RDiam = linspace(1,length(RDiam)./Param.Fs,length(RDiam));
    
    figure
    subplot(2,2,[1 2])
    plot(t_LDiam,LDiam,color='blue');
    hold on
    plot(t_RDiam,RDiam,color='red');
    patch(SU,[-mean(LDiam)*(rem(1:size(SpeakUtt,1)*2,2)'-1),mean(LDiam)*rem(1:size(SpeakUtt,1)*2,2)'],[0 0.3 0.3 0.5])
    title(strrep(PairFiles(i).name,'_','-'))
    xticks([1:60:t_LDiam(end)])
    xlabel('Time [s]')
    ylabel('Pupil diameter [mm]')
    legend('Diameter Left','Diameter Right','Speaking Utterance', 'Listening Utterance')
    
    
end