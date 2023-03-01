%% PBCA-Thesis - Week 9 - Fixation duration? Not using diameter anymore
% Pathing
clear all; clc; close all;
BPath = strsplit(pwd,'PBCA-thesis');
addpath('tools')
addpath([BPath{1} 'Pupil-preprocessing-tools\tools']) % For preprocessing

[subDirs] = GetSubDirsFirstLevelOnly('data');
FileNames={'P1_Quiet_B1.mat','P1_Quiet_B2.mat','P1_SHL_B1.mat','P1_SHL_B2.mat','P1_Noise60_B1.mat','P1_Noise60_B2.mat','P1_Noise70_B1.mat','P1_Noise70_B2.mat','P2_Quiet_B1.mat','P2_Quiet_B2.mat','P2_SHL_B1.mat','P2_SHL_B2.mat','P2_Noise60_B1.mat','P2_Noise60_B2.mat','P2_Noise70_B1.mat','P2_Noise70_B2.mat'};

% Parameters for processing
Param.Fs = 50; % Sampling frequency of pupil data
Param.Preblink = 0.1; % [s], set to NaN, time before blink
Param.Postblink = 0.2; % [s], set to NaN, time after blink
Param.BlinkThresh = 3; % [samples], threshold of samples in between artifacts or blinks

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
        
        % In case merge_windows is transposed from origin [x;y;z] -> transpose to [x,y,z]
        if size([alldata_mat(:,cellfun(@(xd) ~any(isnan(xd)),{alldata_mat.gaze3d})).gaze3d],1) == 3
            for k=1:size(alldata_mat,2)
               alldata_mat(k).gaze3d=alldata_mat(k).gaze3d'; 
            end
        end
        
        % Replace 'NaN' to '[NaN,NaN,NaN]'
        [alldata_mat(:,cellfun(@(xd) any(isnan(xd)),{alldata_mat.gaze3d})).gaze3d] = deal([NaN,NaN,NaN]);
        
        gaze3d = vertcat(alldata_mat.gaze3d);
        
        GazeXRaw = gaze3d(:,1);
        GazeYRaw = gaze3d(:,2);
        GazeZRaw = gaze3d(:,3);
        
        % Preprocessing - 100 ms before and 200 ms after a blink -> NaN
        % NaN's appear at the same time at [X,Y,Z], we only look at X.
        XNan = find(isnan(GazeXRaw));
        
        % Pre-blink
            Pre = [XNan(1);XNan(find(diff(XNan,1) > 1) + 1)];
            for h=1:length(Pre)
                if Pre(h) >= Param.Preblink*Param.Fs % Check start delimiter
                   for k=1:Param.Preblink*Param.Fs
                       GazeXRaw(Pre(h)-k)=NaN;
                       GazeYRaw(Pre(h)-k)=NaN;
                       GazeZRaw(Pre(h)-k)=NaN;
                   end
                end
            end
        
        % Post-blink
        Post = XNan(diff(XNan,1) > Param.BlinkThresh);
        for h=1:length(Post)
            if Post(h) <= length(GazeXRaw) + Param.Postblink*Param.Fs % Check end delimiter
               for k=1:Param.Postblink*Param.Fs
                   GazeXRaw(Post(h)+k)=NaN;
                   GazeYRaw(Post(h)+k)=NaN;
                   GazeZRaw(Post(h)+k)=NaN;
               end
            end
        end
               
        GazeX = GazeXRaw;
        GazeY = GazeYRaw;
        GazeZ = GazeZRaw;
        
        % Plots
%         backimag = uint8(zeros(600,600,3));         % backimag is the background image we have; let's put black
%         [dimY, dimX, ~] = size(backimag);
%         X = round(rescale(GazeX,1,dimX));           % GazeX from [-1,1] -> convert to [1,600]
%         Y = round(rescale(GazeY,1,dimY));           % GazeY from [-1,1] -> convert to [1,600]
%         W = ones(size(GazeX));                      % we could have weights on each gaze (eg to normalise for different number of trials per subject)
%         heatmap = mat2cell(hot,256,[1 1 1]);        % we choose a heatmap colouring, eg "hot", and convert it to "cell"
%         sigma = 20;                                  % the variance parameter for the gaussian kernel
%         % Create "mask"
%         origmask = ones(dimX, dimY)*0.1;
%         for k = 1:size(GazeX,1)
%             origmask(Y(k), X(k)) = origmask(Y(k), X(k)) + W(k);
%         end
%         % Filter using a gaussian kernel
%         mask = imgaussfilt(origmask, sigma);
%         % Normalise total mass of heatmap
%         mask = rescale(mask);
%         % Colour the background image with the heatmap
%         newImage = backimag;
%         for rgbInd = 1:3
%             thisHeat = heatmap{rgbInd}( floor(mask*255) + 1 );
%             newImage(:,:,rgbInd) = (newImage(:,:,rgbInd) + uint8(thisHeat*255));
%         end
%         figure; imshow(newImage); set(gca, 'ydir', 'normal')
%         figure; imshow(origmask); set(gca, 'ydir', 'normal')
%         figure; plot(GazeX,GazeY); xlim([-1 1]); ylim([-1 1])
        
%         [alldata_mat(:,cellfun(@(xd) any(isnan(xd)),{alldata_mat.gaze3d})).gaze3d] % all non-nans in gaze3d
               
    end
end
