clear all; clc; close all;
addpath('tools')

% Colors
SpeakColor=[0, 255, 0, .3*255]./255;
ListenColor=[255, 0, 0, .3*255]./255;

SpeakT = [0.3 1 10;0.3 18 20;0.3 27 32;35+([0.3 5 12;0.3 15 23;0.3 25 30]);0.3 68 78;12+([0.3 69 71;0.3 72 75]);0.1 91 120; 35+([0.1 92 94;0.1 96 99;0.1 100 105; 0.1 106 125])];
ListenT = [0.3 5 12;0.3 15 23;0.3 25 30;35+([0.3 1 10;0.3 18 20;0.3 27 32]);0.3 69 71;0.3 72 75;12+([0.3 68 78]);0.1 92 94;0.1 96 99;0.1 100 105; 0.1 106 115;35+([0.1 91 120])];

% SpeakT = [0.1 1 10; 0.1 20 35];
% ListenT = [0.1 2 4; 0.1 8 13];
[SpeakTO,ListenTO] = overlap_windows(SpeakT, ListenT, 50);
figure;
startStop = SpeakT(:,2:3);width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),0,width(i),1],'EdgeColor', 'none', 'FaceColor', SpeakColor), 1:size(startStop,1));
startStop = ListenT(:,2:3);width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),0,width(i),1],'EdgeColor', 'none', 'FaceColor', ListenColor), 1:size(startStop,1));
startStop = SpeakTO(:,2:3);width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),1,width(i),1],'EdgeColor', 'none', 'FaceColor', SpeakColor), 1:size(startStop,1));
startStop = ListenTO(:,2:3);width = startStop(:,2)-startStop(:,1);hold on;arrayfun(@(i)rectangle('Position', [startStop(i,1),1,width(i),1],'EdgeColor', 'none', 'FaceColor', ListenColor), 1:size(startStop,1));
grid on
yline(1,"--",'HandleVisibility','off')