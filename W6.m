%% PBCA-Thesis - Week 6 & Week 7 - Plot histograms with lenght of gaps, turns, pauses, overlap-between and overlap-within - Multiple Durations
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
GapDurRaw=zeros(200,1000);
GapDur=GapDurRaw;
OLWDurRaw=GapDurRaw;
OLWDur=GapDurRaw;
OLBDurRaw=GapDurRaw;
OLBDur=GapDurRaw;
TurnDurRaw=GapDurRaw;
TurnDur=GapDurRaw;
TurnDiff=GapDurRaw;
TurnDurQuiet=GapDurRaw;
TurnDurSHL=GapDurRaw;
TurnDurN60=GapDurRaw;
TurnDurN70=GapDurRaw;
PauseDurRaw=GapDurRaw;
PauseDur=GapDurRaw;
SpeakDurRaw=GapDurRaw;
ListenDurRaw=GapDurRaw;
SpeakDur=GapDurRaw;
ListenDur=GapDurRaw;

TimeMinWin = 0; % [s], Minimum time of a window
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
            GapKey = 'gapCH1';
            OLWKey = 'overlapWCH1';
            OLBKey = 'overlapBCH1';
            TurnKey = 'turnCH1';
            PauseKey = 'pauseCH1';
        elseif contains(PairFiles(i).name,'P1')
            SpeakKey = 'utteranceCH2';
            ListenKey = 'utteranceCH1';
            GapKey = 'gapCH2';
            OLWKey = 'overlapWCH2';
            OLBKey = 'overlapBCH2';
            TurnKey = 'turnCH2';
            PauseKey = 'pauseCH2';
        end

        if contains(PairFiles(i).name,'B1')
            SpeB = 0;
        elseif contains(PairFiles(i).name,'B2')
            SpeB = 1;
        end

        if contains(PairFiles(i).name,'Quiet')
            SpeCond = SpeB + 1;
            TurnDurQuiet(x,1:size(PairUtt{1,SpeCond}.(TurnKey)(PairUtt{1,SpeCond}.(TurnKey)(:,1)'>TimeMinWin),2)) = PairUtt{1,SpeCond}.(TurnKey)(PairUtt{1,SpeCond}.(TurnKey)(:,1)'>TimeMinWin);
        elseif contains(PairFiles(i).name,'SHL')
            SpeCond = SpeB + 3;
            TurnDurSHL(x,1:size(PairUtt{1,SpeCond}.(TurnKey)(PairUtt{1,SpeCond}.(TurnKey)(:,1)'>TimeMinWin),2)) = PairUtt{1,SpeCond}.(TurnKey)(PairUtt{1,SpeCond}.(TurnKey)(:,1)'>TimeMinWin);
        elseif contains(PairFiles(i).name,'Noise60')
            SpeCond = SpeB + 5;
            TurnDurN60(x,1:size(PairUtt{1,SpeCond}.(TurnKey)(PairUtt{1,SpeCond}.(TurnKey)(:,1)'>TimeMinWin),2)) = PairUtt{1,SpeCond}.(TurnKey)(PairUtt{1,SpeCond}.(TurnKey)(:,1)'>TimeMinWin);
        elseif contains(PairFiles(i).name,'Noise70')
            SpeCond = SpeB + 7;
            TurnDurN70(x,1:size(PairUtt{1,SpeCond}.(TurnKey)(PairUtt{1,SpeCond}.(TurnKey)(:,1)'>TimeMinWin),2)) = PairUtt{1,SpeCond}.(TurnKey)(PairUtt{1,SpeCond}.(TurnKey)(:,1)'>TimeMinWin);
        end

        SpeakRaw = PairUtt{1,SpeCond}.(SpeakKey);
        ListenRaw = PairUtt{1,SpeCond}.(ListenKey);
        GapRaw = PairUtt{1,SpeCond}.(GapKey);
        OLWRaw = PairUtt{1,SpeCond}.(OLWKey);
        OLBRaw = PairUtt{1,SpeCond}.(OLBKey);
        TurnRaw = PairUtt{1,SpeCond}.(TurnKey);
        PauseRaw = PairUtt{1,SpeCond}.(PauseKey);
        binResUtt = PairUtt{1,SpeCond}.binRes;

        if isempty(SpeakRaw) && isempty(ListenRaw)
            disp(['Warning: No associated Utterance/Listening data for file ', PairFiles(i).folder, '\', PairFiles(i).name, '.']);
            x=x+1;
            continue
        end

        % Discard windows if duration is < TimeMinWin
        Speak = SpeakRaw(SpeakRaw(:,1)>TimeMinWin,:);
        Listen = ListenRaw(ListenRaw(:,1)>TimeMinWin,:);
        Gap = GapRaw(GapRaw(:,1)>TimeMinWin,:);
        OLW = OLWRaw(OLWRaw(:,1)>TimeMinWin,:);
        OLB = OLBRaw(OLBRaw(:,1)>TimeMinWin,:);
        Turn = TurnRaw(TurnRaw(:,1)>TimeMinWin,:);
        Pause = PauseRaw(PauseRaw(:,1)>TimeMinWin,:);
        
%         SpeakRaw(:,2:3)=round((SpeakRaw(:,2:3)*binResUtt+BLPeriod(2)+0.36)*Param.Fs+SDelayRaw(1)/2);
%         ListenRaw(:,2:3)=round((ListenRaw(:,2:3)*binResUtt+BLPeriod(2)+0.58)*Param.Fs+LDelayRaw(1)/2);

        % Merge windows if gap <= TimeMergeGap
%         [SpeakM] = MergeWin(SpeakRaw, Param.Fs, TimeMergeGap);
%         [ListenM] = MergeWin(ListenRaw, Param.Fs, TimeMergeGap);

        GapDurRaw(x,1:size(GapRaw,1))=GapRaw(:,1)';
        OLWDurRaw(x,1:size(OLWRaw,1))=OLWRaw(:,1)';
        OLBDurRaw(x,1:size(OLBRaw,1))=OLBRaw(:,1)';
        TurnDurRaw(x,1:size(TurnRaw,1))=TurnRaw(:,1)';
        PauseDurRaw(x,1:size(PauseRaw,1))=PauseRaw(:,1)';
        SpeakDurRaw(x,1:size(SpeakRaw,1))=SpeakRaw(:,1)';
        ListenDurRaw(x,1:size(ListenRaw,1))=ListenRaw(:,1)';
        
        GapDur(x,1:size(Gap,1))=Gap(:,1)';
        OLWDur(x,1:size(OLW,1))=OLW(:,1)';
        OLBDur(x,1:size(OLB,1))=OLB(:,1)';
        TurnDur(x,1:size(Turn,1))=Turn(:,1)';
        PauseDur(x,1:size(Pause,1))=Pause(:,1)';
        SpeakDur(x,1:size(Speak,1))=Speak(:,1)';
        ListenDur(x,1:size(Listen,1))=Listen(:,1)';
        
        % Find time diff between turns (same talker)
        for j=1:size(Turn,1)-1
            TurnDiff(x,j)=(Turn(j+1,2)-Turn(j,3))*binResUtt;
        end
        
        % Median
        TurnMed=[median(nonzeros(TurnDurQuiet));median(nonzeros(TurnDurSHL));median(nonzeros(TurnDurN60));median(nonzeros(TurnDurN70))];
        
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
%         title(strrep([PairFiles(i).folder(25:end),'\',PairFiles(i).name],'_','-'))
%         legend('Speak','Listen','Gap','Overlap-Within','Overlap-Between','Turn','Pause')
        
        if contains([PairFiles(i).folder, '\', PairFiles(i).name],'\Main12\P2_SHL_B2.mat')
%             GapDur1(~any(GapDur1,2),:)=[];
%             OLWDur1(~any(OLWDur1,2),:)=[];
%             OLBDur1(~any(OLBDur1,2),:)=[];
%             TurnDur1(~any(TurnDur1,2),:)=[];
%             PauseDur1(~any(PauseDur1,2),:)=[];
            GroupRaw = struct('SpeakDur',SpeakDurRaw,'ListenDur',ListenDurRaw,'GapDur',GapDurRaw,'PauseDur',PauseDurRaw,'OLWDur',OLWDurRaw,'OLBDur',OLBDurRaw,'TurnDur',TurnDurRaw);
            Group = struct('SpeakDur',SpeakDur,'ListenDur',ListenDur,'GapDur',GapDur,'PauseDur',PauseDur,'OLWDur',OLWDur,'OLBDur',OLBDur,'TurnDur',TurnDur);
            GroupColor = struct('SpeakColor',SpeakColor,'ListenColor',ListenColor,'GapColor',GapColor,'PauseColor',PauseColor,'OLWColor',OLWColor,'OLBColor',OLBColor,'TurnColor',TurnColor);
            histplot(GroupRaw,GroupColor,nbins);
            histplot(Group,GroupColor,nbins);
            
            figure;histogram(nonzeros(TurnDiff),nbins,'FaceColor',TurnColor(1:3));title(['Time between turns - Avg: ',sprintf('%0.5f',mean(nonzeros(TurnDiff))),' s']);
            xlabel('Time [s]'),ylabel('Bin Count');grid on;
            
            figure;tiledlayout(1,4);ax1=nexttile;ax2=nexttile;ax3=nexttile;ax4=nexttile;histogram(ax1,nonzeros(TurnDurQuiet),nbins,'FaceColor',TurnColor(1:3)),histogram(ax2,nonzeros(TurnDurSHL),nbins,'FaceColor',TurnColor(1:3));
            histogram(ax3,nonzeros(TurnDurN60),nbins,'FaceColor',TurnColor(1:3));histogram(ax4,nonzeros(TurnDurN70),nbins,'FaceColor',TurnColor(1:3))
            title(ax1,['Quiet turn duration - Avg: ',sprintf('%0.5f',mean(nonzeros(TurnDurQuiet))),' s']);title(ax2,['SHL turn duration - Avg: ',sprintf('%0.5f',mean(nonzeros(TurnDurSHL))),' s']);
            title(ax3,['Noise60 turn duration - Avg: ',sprintf('%0.5f',mean(nonzeros(TurnDurN60))),' s']);title(ax4,['Noise70 turn duration - Avg: ',sprintf('%0.5f',mean(nonzeros(TurnDurN70))),' s']);
            xlabel([ax1 ax2 ax3 ax4],'Time [s]');ylabel([ax1 ax2 ax3 ax4],'Bin Count');grid([ax1 ax2 ax3 ax4],'on');
        end
    end
end


function histplot(G,C,nbins)
    figure;tiledlayout(2,4);
    ax1=nexttile;ax2=nexttile;ax3=nexttile;ax4=nexttile;ax5=nexttile;ax6=nexttile;ax7=nexttile;
    h1=histogram(ax1,nonzeros(G.SpeakDur),nbins,'FaceColor',C.SpeakColor(1:3));title(ax1,['Speak duration - Avg: ',sprintf('%0.5f',mean(nonzeros(G.SpeakDur))),' s']);
    h2=histogram(ax2,nonzeros(G.ListenDur),nbins,'FaceColor',C.ListenColor(1:3));title(ax2,['Listen duration - Avg: ',sprintf('%0.5f',mean(nonzeros(G.ListenDur))),' s']);
    h3=histogram(ax3,nonzeros(G.GapDur),nbins,'FaceColor',C.GapColor(1:3));title(ax3,['Gap duration - Avg: ',sprintf('%0.5f',mean(nonzeros(G.GapDur))),' s']);
    h4=histogram(ax4,nonzeros(G.PauseDur),nbins,'FaceColor',C.PauseColor(1:3));title(ax4,['Pause duration - Avg: ',sprintf('%0.5f',mean(nonzeros(G.PauseDur))),' s']);
    h5=histogram(ax5,nonzeros(G.OLWDur),nbins,'FaceColor',C.OLWColor(1:3));title(ax5,['Overlap-Within duration - Avg: ',sprintf('%0.5f',mean(nonzeros(G.OLWDur))),' s']);
    h6=histogram(ax6,nonzeros(G.OLBDur),nbins,'FaceColor',C.OLBColor(1:3));title(ax6,['Overlap-Between duration - Avg: ',sprintf('%0.5f',mean(nonzeros(G.OLBDur))),' s']);
    h7=histogram(ax7,nonzeros(G.TurnDur),nbins*10,'FaceColor',C.TurnColor(1:3));title(ax7,['Turn duration - Avg: ',sprintf('%0.5f',mean(nonzeros(G.TurnDur))),' s']);
    xlabel([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'Time [s]')
    ylabel([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'Bin Count')
    xlim([ax1,ax2,ax3,ax4,ax5,ax6,ax7],[-0.1 5])
    grid([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'on')
    lgd=legend([h1,h2,h3,h4,h5,h6,h7],'Speak','Listen','Gap','Pause','Overlap-Within','Overlap-Between','Turn');
    lgd.Layout.Tile = 8;
end
