%% Master's thesis pre-project -- Pablo Castro (PBCA)
clear all; clc; close all;
addpath(genpath('\data'));
datadir=dir('data\*.mat');

for i=1:length(datadir)
    alldata(i) = load(datadir(i).name);
    s_ad(i,:) = [cell2struct(alldata(i).data,{'data'},1)];
    for j=1:size(s_ad(i,:),2)
        DiamLeft(j,:) = s_ad(j).data.diameterLeft;
    end
    plot(1:length(DiamLeft),DiamLeft);
end