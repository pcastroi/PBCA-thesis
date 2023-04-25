%% PBCA-Thesis - Week 3 - Delay between Tobii and Audio
% Pathing
clear all; clc; close all;
BPath = strsplit(pwd,'PBCA-thesis');
addpath('tools')
addpath([BPath{1} 'Pupil-preprocessing-tools\tools']) % For preprocessing

% Files and Utterances: different conditions
TobiiFs=50;
AudFs=48000;
BLDelay=20; % s
xd=zeros(1,96);
x=1;
[subDirs] = GetSubDirsFirstLevelOnly('data');

for q=1:numel(subDirs)
    PairIn = q;
    PairFiles=dir(['data\AMEND_I\Main',sprintf('%d',PairIn),'\*.mat']);
    AudFiles=dir(['audio\AMEND_I\Main',sprintf('%d',PairIn),'\*.mat']);
    if isempty(AudFiles)
       continue 
    end
    PairUtt=load('data\AMEND_I\utterances1110.mat');
    PairUtt=PairUtt.Utterances(PairIn,:);
    
    for i=1:numel(PairFiles)
        alldata = load([PairFiles(i).folder, '\', PairFiles(i).name]);
        alldata_mat = cell2mat(alldata.data);

        if contains(PairFiles(i).name,'P2')
            DelKey = 'delayCH1';
        elseif contains(PairFiles(i).name,'P1')
            DelKey = 'delayCH2';
        end

        if contains(PairFiles(i).name,'B1')
            RepKey='Rep1';
            SpeB = 0;
        elseif contains(PairFiles(i).name,'B2')
            RepKey='Rep2';
            SpeB = 1;
        end

        if contains(PairFiles(i).name,'Quiet')
            CondKey='NH-Quiet';
            SpeCond = SpeB + 1;
        elseif contains(PairFiles(i).name,'SHL')
            CondKey='SHL-Quiet';
            SpeCond = SpeB + 3;
        elseif contains(PairFiles(i).name,'Noise60')
            CondKey='NH-Noise60';
            SpeCond = SpeB + 5;
        elseif contains(PairFiles(i).name,'Noise70')
            CondKey='NH-Noise70';
            SpeCond = SpeB + 7;
        end

        AudKey = [CondKey,'_',RepKey,'.mat'];
        AudData=load(['audio\Main',sprintf('%d',PairIn),'\',AudKey]);
        FolderFile = strsplit(PairFiles(i).folder,'\');
        
        TobAudDelay{q,SpeCond}.(DelKey)(:,1)=numel(alldata_mat)-(size(AudData.outputAll.speech(:,1),1)/AudFs+20)*TobiiFs;
        TobAudDelay{q,SpeCond}.(DelKey)(:,2)=numel(alldata_mat)/TobiiFs-size(AudData.outputAll.speech(:,1),1)/AudFs-20;
        TobAudDelay{q,SpeCond}.binRes=1/TobiiFs;
        xd(:,x)=TobAudDelay{q,SpeCond}.(DelKey)(:,2);
        x=x+1;
    end
end
save("data\delays1110.mat",'TobAudDelay');
