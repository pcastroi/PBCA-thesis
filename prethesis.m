%% Master's thesis pre-project -- AMEND 1 data -- Pablo Castro (PBCA)
clear all; clc; close all;
BPath = strsplit(pwd,'PBCA-thesis');
addpath(genpath('data\'));
addpath([BPath{1} 'Pupil-preprocessing-tools\tools']) % For preprocessing
datadir=dir('data\*.mat');
% Parameters
Param.Fs = 250;
Param.RemoveBeforeAndAfter = [35 100]*1e-3;
Param.MinLengthNaNRepair = 5;
LPWinSize = 0.5; % [s]: Window size of hamming-window for low-pass filtering
LPWindow = hamming(round(LPWinSize*Param.Fs));
LPWindow = LPWindow/sum(LPWindow); % Hamming-window

% Initialization
LDiamRaw = zeros(16,13300);
RDiamRaw = LDiamRaw;

for i=1:length(datadir)
    alldata(i) = load(datadir(i).name);
    alldata_mat = cell2mat(alldata(i).data);
    
    % Preprocessing - Interpolating NaNs
    LDiamRaw(i,1:length([alldata_mat.diameterLeft])) = [alldata_mat.diameterLeft];
    RDiamRaw(i,1:length([alldata_mat.diameterRight])) = [alldata_mat.diameterRight];
    [LDiam(i,:),LMetadata] = preprocpupil(LDiamRaw(i,:),Param);
    [RDiam(i,:),RMetadata] = preprocpupil(RDiamRaw(i,:),Param);
    
    % Low-Pass Filtering
    LDiam(i,:) = conv(LDiam(i,:),LPWindow,'same');
    RDiam(i,:) = conv(RDiam(i,:),LPWindow,'same');
    
    % Time vectors for plotting
    t_LDiam(i,:) = linspace(1,length(LDiam(i,:)/Param.Fs),length(LDiam(i,:)));
    t_RDiam(i,:) = linspace(1,length(RDiam(i,:)/Param.Fs),length(RDiam(i,:)));
    
    figure
    plot(t_LDiam(i,:),LDiam(i,:),color='blue');
    hold on
    plot(t_RDiam(i,:),RDiam(i,:),color='red');
    title(strrep(datadir(i).name,'_','-'))
    legend('Diameter Left','Diameter Right')
end