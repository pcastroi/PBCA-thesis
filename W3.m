%% PBCA-Thesis - Week 3 - Delay between Tobii and Audio?
% Pathing
clear all; clc; close all;
BPath = strsplit(pwd,'PBCA-thesis');
addpath([BPath{1} 'Pupil-preprocessing-tools\tools']) % For preprocessing

% Files and Utterances: different conditions
TobiiFs=50;
BLDelay=20; % s
x=0;
[subDirs] = GetSubDirsFirstLevelOnly('data');

for q=1:numel(subDirs)
    PairIn = q;
    PairFiles=dir(['data\Main',sprintf('%d',PairIn),'\*.mat']);
    AudFiles=dir(['audio\Main',sprintf('%d',PairIn),'\*.wav']);
    if isempty(AudFiles)
       continue 
    end
    PairUtt=load('data\utterances1110.mat');
    PairUtt=PairUtt.Utterances(PairIn,:);
    
    for i=1:numel(PairFiles)
        AudKey=[];
        alldata = load([PairFiles(i).folder, '\', PairFiles(i).name]);
        alldata_mat = cell2mat(alldata.data);

        if contains(PairFiles(i).name,'P1')
            TalkerKey = 'talker1';
        elseif contains(PairFiles(i).name,'P2')
            TalkerKey = 'talker2';
        end

        if contains(PairFiles(i).name,'B1')
            RepKey='Rep1';
        elseif contains(PairFiles(i).name,'B2')
            RepKey='Rep2';
        end

        if contains(PairFiles(i).name,'Quiet')
            CondKey='SpeechNH-Quiet';
        elseif contains(PairFiles(i).name,'Noise60')
            CondKey='SpeechNH-Noise60';
        elseif contains(PairFiles(i).name,'Noise70')
            CondKey='SpeechNH-Noise70';
        elseif contains(PairFiles(i).name,'SHL')
            CondKey='SpeechSHL-Quiet';
        end

        AudKey = [CondKey,'_',RepKey,'_',TalkerKey,'.wav'];
        
        [AudData,AudFs]=audioread(['audio\Main',sprintf('%d',PairIn),'\',AudKey]);
        FolderFile = strsplit(PairFiles(i).folder,'\');
        EyeAudDelay=alldata_mat(end).timeStamp-alldata_mat(1).timeStamp-length(alldata_mat)/TobiiFs;
        if ~isempty(AudKey)
            disp(['FILE ',int2str(i),' - ',FolderFile{end},' | Tobii duration = ',sprintf('%0.2f',numel(alldata_mat)/TobiiFs),' s | Audio duration = ',sprintf('%0.2f',numel(AudData)/AudFs),' s | DIFF = ',sprintf('%0.5f',numel(alldata_mat)/TobiiFs-numel(AudData)/AudFs),' s | Delay timeStamp = ',sprintf('%0.5f',EyeAudDelay),'s | Intrinsic delay = ',sprintf('%0.5f',numel(alldata_mat)/TobiiFs-numel(AudData)/AudFs-EyeAudDelay-BLDelay),'s.']);
        else
            disp(['FILE ',int2str(i),' - ',FolderFile{end},' | No data'])
        end
    end
end