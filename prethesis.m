%% Master's thesis pre-project -- AMEND 1 data -- Pablo Castro (PBCA)
% Main X - X pairs of participants
% <=16 files: Talkers P1 and P2, 2 repetitions B1 and B2.
clear all; clc; close all;
BPath = strsplit(pwd,'PBCA-thesis');
addpath(genpath('data\'));
addpath([BPath{1} 'Pupil-preprocessing-tools\tools']) % For preprocessing

dirlist=dir('data\');
datadir=dir('data\Main*\*.mat');

% Parameters
Param.Fs = 250;
Param.RemoveBeforeAndAfter = [35 100]*1e-3;
Param.MinLengthNaNRepair = 5;
LPWinSize = 0.75; % [s]: Window size of hamming-window for low-pass filtering
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
    [LDiam(i,:),LMetadata(i,:)] = preprocpupil(LDiamRaw(i,:),Param);
    [RDiam(i,:),RMetadata(i,:)] = preprocpupil(RDiamRaw(i,:),Param);
    
    % Low-Pass Filtering
    LDiam(i,:) = conv(LDiam(i,:),LPWindow,'same');
    RDiam(i,:) = conv(RDiam(i,:),LPWindow,'same');
    
    % Decide 'better' eye results
    [Min idx_decision] = min([sum(LMetadata(i).Isnan) sum(RMetadata(i).Isnan)]);
    if idx_decision == 1
        Diameter(i,:) = LDiam(i,:);
        DiameterRaw(i,:) = LDiamRaw(i,:);
        eyeChosen(i,:) = 'Left ';
    elseif idx_decision == 2
        Diameter(i,:) = RDiam(i,:);
        DiameterRaw(i,:) = RDiamRaw(i,:);
        eyeChosen(i,:) = 'Right';
    end
    
    % Time vectors for plotting
    t_LDiam(i,:) = linspace(1,length(LDiam(i,:))/Param.Fs,length(LDiam(i,:)));
    t_RDiam(i,:) = linspace(1,length(RDiam(i,:))./Param.Fs,length(RDiam(i,:)));
    
    figure
    plot(t_LDiam(i,:),LDiam(i,:),color='blue');
    hold on
    plot(t_RDiam(i,:),RDiam(i,:),color='red');
    title(strrep(datadir(i).name,'_','-'))
    xlabel('Time [s]')
    ylabel('Pupil diameter [mm]')
    legend('Diameter Left','Diameter Right')
end
%% List of unique filenames (dirlist(3) = Main1: constains 16 unique names)
namelist=struct2cell(dir(['data\' dirlist(3).name]));
namelist=erase(namelist(1,3:end),'.mat');

% Obtain idx and group data by filename (better eye)
for i=1:length(namelist)
    name_idx=find(contains(datadirl,namelist(i)));
    % Store (Diameter) grouped data in temp variable
    for j=1:length(name_idx)
        datagroup(j,:)=Diameter(name_idx(j),:);
    end
    % Store everything in struct
    sdata.(namelist{i})=datagroup;
end

%%
figure
plot(t_LDiam(1,:),mean(Diameter,1),color='red')
hold on
for i=1:length(datadir)
    plot(t_LDiam(i,:),Diameter(i,:),color=[0 0.3 0.3 0.4]);
end
