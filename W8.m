%% PBCA-Thesis - Week 8 - Visualize delay for each trial
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
[subDirs] = GetSubDirsFirstLevelOnly('data\AMEND_I');
FileNames={'P1_Quiet_B1.mat','P1_Quiet_B2.mat','P1_SHL_B1.mat','P1_SHL_B2.mat','P1_Noise60_B1.mat','P1_Noise60_B2.mat','P1_Noise70_B1.mat','P1_Noise70_B2.mat','P2_Quiet_B1.mat','P2_Quiet_B2.mat','P2_SHL_B1.mat','P2_SHL_B2.mat','P2_Noise60_B1.mat','P2_Noise60_B2.mat','P2_Noise70_B1.mat','P2_Noise70_B2.mat'};
LoadUtt=load('data\AMEND_I\utterances1110.mat');
LoadDelays=load('data\AMEND_I\delays1110.mat');

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

        if contains(cell2mat(FileNames(i)),'P2')
            DelKey = 'delayCH1';
        elseif contains(cell2mat(FileNames(i)),'P1')
            DelKey = 'delayCH2';
        end

        if contains(cell2mat(FileNames(i)),'B1')
            RepKey='Rep1';
            SpeB = 0;
        elseif contains(cell2mat(FileNames(i)),'B2')
            RepKey='Rep2';
            SpeB = 1;
        end

        if contains(cell2mat(FileNames(i)),'Quiet')
            CondKey='NH-Quiet';
            SpeCond = SpeB + 1;
        elseif contains(cell2mat(FileNames(i)),'SHL')
            CondKey='SHL-Quiet';
            SpeCond = SpeB + 3;
        elseif contains(cell2mat(FileNames(i)),'Noise60')
            CondKey='NH-Noise60';
            SpeCond = SpeB + 5;
        elseif contains(cell2mat(FileNames(i)),'Noise70')
            CondKey='NH-Noise70';
            SpeCond = SpeB + 7;
        end

        AudKey = [CondKey,'_',RepKey,'.mat'];
        try
            AudData=load(['audio\Main',sprintf('%d',PairIn),'\',AudKey]);
        catch ME
            disp(['Warning: No associated Audio for file ', PairFiles(1).folder, '\', cell2mat(FileNames(i)), '.']);
            continue
        end
        
        FolderFile = strsplit(PairFiles(1).folder,'\');
        
        TobAudDelay{q,SpeCond}.(DelKey)(:,1)=numel(alldata_mat)-(size(AudData.outputAll.speech(:,1),1)/AudFs+20)*TobiiFs;
        TobAudDelay{q,SpeCond}.(DelKey)(:,2)=numel(alldata_mat)/TobiiFs-size(AudData.outputAll.speech(:,1),1)/AudFs-20;
        TobAudDelay{q,SpeCond}.binRes=1/TobiiFs;
        NSampDel(q,i)=TobAudDelay{q,SpeCond}.(DelKey)(:,1);
        TDelay(q,i)=TobAudDelay{q,SpeCond}.(DelKey)(:,2);
        x=x+1;
    end
end

%% Plots
figure; hold on; ax1=gca;
h1=image(ax1,NSampDel,'CDataMapping','scaled');colormap autumn;colorbar;
set(h1, 'AlphaData', ~isnan(NSampDel))
for i = 1:numel(subDirs)
  for j = 1:numel(FileNames)
      nu = NSampDel(i,j);
      val = num2str(round(nu));
      text(j,i,val,'HorizontalAlignment','center')
  end
end
hold off;

figure; hold on; ax2=gca;
h2=image(ax2,TDelay,'CDataMapping','scaled');colormap autumn;colorbar;
set(h2, 'AlphaData', ~isnan(TDelay))
for i = 1:numel(subDirs)
  for j = 1:numel(FileNames)
      nu = TDelay(i,j);
      val = num2str(round(nu,3));
      text(j,i,val,'HorizontalAlignment','center')
  end
end
hold off;

xticks([ax1 ax2],1:1:size(NSampDel,2));xticklabels([ax1 ax2],strsplit(strrep(cell2mat(FileNames),'_','-'),'.mat'))
yticks([ax1 ax2],1:1:size(NSampDel,1));yticklabels([ax1 ax2],{'Main1','Main2','Main3','Main4','Main5','Main6','Main7','Main8','Main9','Main10','Main11','Main12'})

title(ax1,'Delay in Samples')
title(ax2,'Delay in Seconds')
