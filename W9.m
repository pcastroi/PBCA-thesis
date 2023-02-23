%% PBCA-Thesis - Week 9 - Fixation duration? Not using diameter anymore
% Pathing
clear all; clc; close all;
BPath = strsplit(pwd,'PBCA-thesis');
addpath('tools')
addpath([BPath{1} 'Pupil-preprocessing-tools\tools']) % For preprocessing

% Files and Utterances: different conditions
TobiiFs=50;
AudFs=48000;
BLDelay=20; % s
TDelay=zeros(12,16)*NaN;
NSampDel=TDelay;
x=1;
[subDirs] = GetSubDirsFirstLevelOnly('data');
FileNames={'P1_Quiet_B1.mat','P1_Quiet_B2.mat','P1_SHL_B1.mat','P1_SHL_B2.mat','P1_Noise60_B1.mat','P1_Noise60_B2.mat','P1_Noise70_B1.mat','P1_Noise70_B2.mat','P2_Quiet_B1.mat','P2_Quiet_B2.mat','P2_SHL_B1.mat','P2_SHL_B2.mat','P2_Noise60_B1.mat','P2_Noise60_B2.mat','P2_Noise70_B1.mat','P2_Noise70_B2.mat'};

for q=1:numel(subDirs)
    PairIn = q;
    PairFiles=dir(['data\Main',sprintf('%d',PairIn),'\*.mat']);
    AudFiles=dir(['audio\Main',sprintf('%d',PairIn),'\*.mat']);
    if isempty(AudFiles)
        disp(['Warning: No associated Audio for file ', PairFiles(i).folder, '\', PairFiles(i).name, '.']);
       continue 
    end
    PairUtt=load('data\utterances1110.mat');
    PairUtt=PairUtt.Utterances(PairIn,:);
    
    for i=1:numel(FileNames)
        try
            alldata = load([PairFiles(1).folder, '\', cell2mat(FileNames(i))]);
        catch ME
            disp(['Warning: No associated Delay data for file ', PairFiles(1).folder, '\', cell2mat(FileNames(i)), '.']);
            continue
        end
        alldata_mat = cell2mat(alldata.data);
        
        [alldata_mat(:,cellfun(@(xd) ~any(isnan(xd)),{alldata_mat.gaze2d})).gaze2d] % all non-nans in gaze2d
    end
end
