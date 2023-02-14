%% PBCA-Thesis - Week 6 - Plot histograms with lenght of gaps, pauses, gap-between and gap-within
% Pathing
clear all; clc; close all;
BPath = strsplit(pwd,'PBCA-thesis');
addpath('tools')
addpath([BPath{1} 'Pupil-preprocessing-tools\tools']) % For preprocessing

% Colors
SpeakColor=[0, 255, 0, .3*255]./255;
ListenColor=[255, 0, 0, .3*255]./255;
GapColor=[136, 247, 226, .3*255]./255;
OLWColor=[255, 190, 0, .3*255]./255;
OLBColor=[255, 238, 88, .3*255]./255;
TurnColor=[217, 217, 217, .3*255]./255;
PauseColor=[117, 117, 117, .3*255]./255;

[subDirs] = GetSubDirsFirstLevelOnly('data');
LoadDelays=load('data\delays1110.mat');
LoadUtt=load('data\utterances1110.mat');
GapDur=zeros(200,1000);
OLWDur=GapDur;
OLBDur=GapDur;
TurnDur=GapDur;
PauseDur=GapDur;
SpeakDur=GapDur;
ListenDur=GapDur;
nbins=50;
x=1;

for q=1:numel(subDirs)
    PairIn = q;
    PairFiles=dir(['data\Main',sprintf('%d',PairIn),'\*.mat']);
    PairUtt=LoadUtt.Utterances(PairIn,:);
    PairDelay=LoadDelays.TobAudDelay(PairIn,:);
    
    for i=1:numel(PairFiles)
        % Retrieve Utterances
        if contains(PairFiles(i).name,'P2')
            SpeakKey = 'utteranceCH1';
            ListenKey = 'utteranceCH2';
            SDelayKey = 'delayCH1';
            LDelayKey = 'delayCH2';
        elseif contains(PairFiles(i).name,'P1')
            SpeakKey = 'utteranceCH2';
            ListenKey = 'utteranceCH1';
            SDelayKey = 'delayCH2';
            LDelayKey = 'delayCH1';
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

        SpeakRaw = PairUtt{1,SpeCond}.(SpeakKey);
        ListenRaw = PairUtt{1,SpeCond}.(ListenKey);
        GapRaw = [PairUtt{1,SpeCond}.('gapCH1'); PairUtt{1,SpeCond}.('gapCH2')];
        OLWRaw = [PairUtt{1,SpeCond}.('overlapWCH1');PairUtt{1,SpeCond}.('overlapWCH2')];
        OLBRaw = [PairUtt{1,SpeCond}.('overlapBCH1');PairUtt{1,SpeCond}.('overlapBCH2')];
        TurnRaw = [PairUtt{1,SpeCond}.('turnCH1');PairUtt{1,SpeCond}.('turnCH2')];
        PauseRaw = [PairUtt{1,SpeCond}.('pauseCH1');PairUtt{1,SpeCond}.('pauseCH2')];
        binResUtt = PairUtt{1,SpeCond}.binRes;
%         try
%             SDelayRaw = PairDelay{1,SpeCond}.(SDelayKey);
%             LDelayRaw = PairDelay{1,SpeCond}.(LDelayKey);
%         catch ME
%             disp(['Warning: No associated Delay data for file ', PairFiles(i).folder, '\', PairFiles(i).name, '.']);
%             continue
%         end
% 
%         if or(SDelayRaw < 0,LDelayRaw < 0)
%             SDelayRaw=[0,0];
%             LDelayRaw=[0,0];
%         end
% 
        if isempty(SpeakRaw) && isempty(ListenRaw)
            disp(['Warning: No associated Utterance/Listening data for file ', PairFiles(i).folder, '\', PairFiles(i).name, '.']);
            continue
        end

%         SpeakRaw(:,2:3)=round((SpeakRaw(:,2:3)*binResUtt+BLPeriod(2)+0.36)*Param.Fs+SDelayRaw(1)/2);
%         ListenRaw(:,2:3)=round((ListenRaw(:,2:3)*binResUtt+BLPeriod(2)+0.58)*Param.Fs+LDelayRaw(1)/2);

        % Merge windows if gap <= TimeMergeGap
%         [SpeakM] = MergeWin(SpeakRaw, Param.Fs, TimeMergeGap);
%         [ListenM] = MergeWin(ListenRaw, Param.Fs, TimeMergeGap);

        GapDur(x,1:size(GapRaw,1))=GapRaw(:,1)';
        OLWDur(x,1:size(OLWRaw,1))=OLWRaw(:,1)';
        OLBDur(x,1:size(OLBRaw,1))=OLBRaw(:,1)';
        TurnDur(x,1:size(TurnRaw,1))=TurnRaw(:,1)';
        PauseDur(x,1:size(PauseRaw,1))=PauseRaw(:,1)';
        SpeakDur(x,1:size(SpeakRaw,1))=SpeakRaw(:,1)';
        ListenDur(x,1:size(ListenRaw,1))=ListenRaw(:,1)';
        x=x+1;

        %Plots
%         figure;
%         startStop = SpeakRaw(:,2:3);width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),0,width(i),1],'EdgeColor', 'none', 'FaceColor', SpeakColor), 1:size(startStop,1));
%         startStop = ListenRaw(:,2:3);width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),0,width(i),1],'EdgeColor', 'none', 'FaceColor', ListenColor), 1:size(startStop,1));
%         startStop = GapRaw(:,2:3);width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),1,width(i),1],'EdgeColor', 'none', 'FaceColor', GapColor), 1:size(startStop,1))
%         startStop = OLWRaw(:,2:3);width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),1,width(i),1],'EdgeColor', 'none', 'FaceColor', OLWColor), 1:size(startStop,1))
%         startStop = OLBRaw(:,2:3);width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),1,width(i),1],'EdgeColor', 'none', 'FaceColor', OLBColor), 1:size(startStop,1))
%         startStop = TurnRaw(:,2:3);width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),1,width(i),1],'EdgeColor', 'none', 'FaceColor', TurnColor), 1:size(startStop,1))
%         startStop = PauseRaw(:,2:3);width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),1,width(i),1],'EdgeColor', 'none', 'FaceColor', PauseColor), 1:size(startStop,1))
%         grid on
%         line(NaN,NaN,'linewidth',5,'Color',SpeakColor);
%         line(NaN,NaN,'linewidth',5,'Color',ListenColor);
%         line(NaN,NaN,'linewidth',5,'Color',GapColor);
%         line(NaN,NaN,'linewidth',5,'Color',OLWColor);
%         line(NaN,NaN,'linewidth',5,'Color',OLBColor);
%         line(NaN,NaN,'linewidth',5,'Color',TurnColor);
%         line(NaN,NaN,'linewidth',5,'Color',PauseColor);
%         legend('Speak','Listen','Gap','Overlap-Within','Overlap-Between','Turn','Pause')
        
        if contains([PairFiles(i).folder, '\', PairFiles(i).name],'\Main12\P2_SHL_B2.mat')
            GapDur(~any(GapDur,2),:)=[];
            OLWDur(~any(OLWDur,2),:)=[];
            OLBDur(~any(OLBDur,2),:)=[];
            TurnDur(~any(TurnDur,2),:)=[];
            PauseDur(~any(PauseDur,2),:)=[];
            figure;tiledlayout(2,4);
            ax1=nexttile;ax2=nexttile;ax3=nexttile;ax4=nexttile;ax5=nexttile;ax6=nexttile;ax7=nexttile;
            h1=histogram(ax1,nonzeros(SpeakDur),nbins,'FaceColor',SpeakColor(1:3));title(ax1,['Speak duration - Avg: ',sprintf('%0.5f',mean(nonzeros(SpeakDur))),' s']);
            h2=histogram(ax2,nonzeros(ListenDur),nbins,'FaceColor',ListenColor(1:3));title(ax2,['Listen duration - Avg: ',sprintf('%0.5f',mean(nonzeros(ListenDur))),' s']);
            h3=histogram(ax3,nonzeros(GapDur),nbins,'FaceColor',GapColor(1:3));title(ax3,['Gap duration - Avg: ',sprintf('%0.5f',mean(nonzeros(GapDur))),' s']);
            h4=histogram(ax4,nonzeros(PauseDur),nbins,'FaceColor',PauseColor(1:3));title(ax4,['Pause duration - Avg: ',sprintf('%0.5f',mean(nonzeros(PauseDur))),' s']);
            h5=histogram(ax5,nonzeros(OLWDur),nbins,'FaceColor',OLWColor(1:3));title(ax5,['Overlap-Within duration - Avg: ',sprintf('%0.5f',mean(nonzeros(OLWDur))),' s']);
            h6=histogram(ax6,nonzeros(OLBDur),nbins,'FaceColor',OLBColor(1:3));title(ax6,['Overlap-Between duration - Avg: ',sprintf('%0.5f',mean(nonzeros(OLBDur))),' s']);
            h7=histogram(ax7,nonzeros(TurnDur),nbins,'FaceColor',TurnColor(1:3));title(ax7,['Turn duration - Avg: ',sprintf('%0.5f',mean(nonzeros(TurnDur))),' s']);
            xlabel([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'Time [s]')
            ylabel([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'Bin Count')
            grid([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'on')
            lgd=legend([h1,h2,h3,h4,h5,h6,h7],'Speak','Listen','Gap','Pause','Overlap-Within','Overlap-Between','Turn');
            lgd.Layout.Tile = 8;
        end
    end
end
