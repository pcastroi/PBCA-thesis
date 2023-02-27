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
Param.RemoveBeforeAndAfter = [35 100]*1e-3; % Samples within the time range before and after NaNs will set NaNs as well.
Param.MinLengthNaNRepair = 5; % Drop values (i.e., change to NaN) before and after NaNs only for contiguous NaNs of at least __ samples. 
LPWinSize = 1; % [s]: Window size of hamming-window for low-pass filtering
LPWindow = hamming(round(LPWinSize*Param.Fs));
LPWindow = LPWindow/sum(LPWindow); % Hamming-window

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
        
        % In case gaze2d is transposed from origin [x;y] -> transpose to [x,y]
        if size([alldata_mat(:,cellfun(@(xd) ~any(isnan(xd)),{alldata_mat.gaze2d})).gaze2d],1) == 2
            for k=1:size(alldata_mat,2)
               alldata_mat(k).gaze2d=alldata_mat(k).gaze2d'; 
            end
        end
        
        % Replace 'NaN' to '[NaN,NaN]'
        [alldata_mat(:,cellfun(@(xd) any(isnan(xd)),{alldata_mat.gaze2d})).gaze2d] = deal([NaN,NaN]);
        
        gaze2d = vertcat(alldata_mat.gaze2d);
        
        GazeXRaw = gaze2d(:,1);
        GazeYRaw = gaze2d(:,2);
        
        % Preprocessing - Setting outliers as NaNs (remove artifacts)
        XThreshOut = [mean(GazeXRaw,'omitnan')-2*std(GazeXRaw,'omitnan'), mean(GazeXRaw,'omitnan')+2*std(GazeXRaw,'omitnan')];
        YThreshOut = [mean(GazeYRaw,'omitnan')-2*std(GazeYRaw,'omitnan'), mean(GazeYRaw,'omitnan')+2*std(GazeYRaw,'omitnan')];
        for s=1:length(alldata_mat)
            if GazeXRaw(s,1) < XThreshOut(1) || GazeXRaw(s,1) > XThreshOut(2)
                GazeXRaw(s,1)=NaN;
            elseif GazeYRaw(s,1) < YThreshOut(1) || GazeYRaw(s,1) > YThreshOut(2)
                GazeYRaw(s,1)=NaN;
            end
        end
        
        % Processing - Interpolating NaNs
        [GazeX,XMetadata] = preprocpupil(GazeXRaw,Param);
        [GazeY,YMetadata] = preprocpupil(GazeYRaw,Param);

        % Low-Pass Filtering
        GazeXConv = conv(GazeX,LPWindow,'same'); 
        GazeYConv = conv(GazeY,LPWindow,'same');
        
        % Remove start/end artifacts (peak/dip) originated from Low-Pass Filtering
        GazeXConv(1:round(length(LPWindow)/2-1)) = mean(GazeX(1:round(length(LPWindow)/2-1)));
        GazeXConv(end-round(length(LPWindow)/2-1):end) = mean(GazeX(end-round(length(LPWindow)/2-1):end));
        GazeYConv(1:round(length(LPWindow)/2-1)) = mean(GazeY(1:round(length(LPWindow)/2-1)));
        GazeYConv(end-round(length(LPWindow)/2-1):end) = mean(GazeY(end-round(length(LPWindow)/2-1):end));
        
        GazeX = GazeXConv;
        GazeY = GazeYConv;
        
        % Limit/Clip to range [-1,1]
        GazeX(GazeX>1)=1;GazeX(GazeX<-1)=-1;
        GazeY(GazeY>1)=1;GazeY(GazeY<-1)=-1;
        
        % Plots
        backimag = uint8(zeros(600,600,3));         % backimag is the background image we have; let's put black
        [dimY, dimX, ~] = size(backimag);
        X = round(rescale(GazeX,1,dimX));           % GazeX from [-1,1] -> convert to [1,600]
        Y = round(rescale(GazeY,1,dimY));           % GazeY from [-1,1] -> convert to [1,600]
        W = ones(size(GazeX));                      % we could have weights on each gaze (eg to normalise for different number of trials per subject)
        heatmap = mat2cell(hot,256,[1 1 1]);        % we choose a heatmap colouring, eg "hot", and convert it to "cell"
        sigma = 20;                                  % the variance parameter for the gaussian kernel
        % Create "mask"
        origmask = ones(dimX, dimY)*0.1;
        for k = 1:size(GazeX,1)
            origmask(Y(k), X(k)) = origmask(Y(k), X(k)) + W(k);
        end
        % Filter using a gaussian kernel
        mask = imgaussfilt(origmask, sigma);
        % Normalise total mass of heatmap
        mask = rescale(mask);
        % Colour the background image with the heatmap
        newImage = backimag;
        for rgbInd = 1:3
            thisHeat = heatmap{rgbInd}( floor(mask*255) + 1 );
            newImage(:,:,rgbInd) = (newImage(:,:,rgbInd) + uint8(thisHeat*255));
        end
        figure; imshow(newImage); set(gca, 'ydir', 'normal')
        figure; imshow(origmask); set(gca, 'ydir', 'normal')
        figure; plot(GazeX,GazeY); xlim([-1 1]); ylim([-1 1])
        
%         [alldata_mat(:,cellfun(@(xd) any(isnan(xd)),{alldata_mat.gaze2d})).gaze2d] % all non-nans in gaze2d
        
        
%          for efi = fixationduration:nPoints
%             if sum([sum(GazeData.left_validity(1+efi-fixationduration:efi)),sum((GazeData.right_validity(1+efi-fixationduration:efi)))]) == 0 %if all rows of tracked data have Tobii validity indicator of 0
%                 if max(gazePoints(1+efi-fixationduration:efi,1)) - min(gazePoints(1+efi-fixationduration:efi,1))<= fixationhrange %if gaze x remains within spatial range
%                     if max(gazePoints(1+efi-fixationduration:efi,2)) - min(gazePoints(1+efi-fixationduration:efi,2))<= fixationvrange %if gaze y remains within spatial range
%                         fixIndicator(1+efi-fixationduration:efi) = 1;
%                     end
%                 end
%             end
%          end
        
        
    end
end
