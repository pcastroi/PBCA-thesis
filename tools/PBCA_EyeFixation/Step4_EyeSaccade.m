
%%%%%%%%%%%%%%%%%%%%%%% Gaze Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%
% In row 21(Main11,P1), we have invalid data for 3 files,
% those just have 800 samples

% In some few cases the duration was less than 4 min
% If it is less than 4 min, it is mainly finished after
% 225s
% Just in two cases it finishes around 215
% So I concider 220 sec for final analysis

clear all
close all
clc

Path= 'C:\Users\sulb\OneDrive - Demant\Documents\IFD_Project\ProjectDescription\Programs\Data\';
savePath= 'C:\Users\sulb\OneDrive - Demant\Documents\IFD_Project\ProjectDescription\Programs\MainPrograms\Version 29072022\Results\';
groupNames= {'Main1','Main2','Main3',...
    'Main4','Main5','Main6','Main7','Main8',...
    'Main9','Main10','Main11','Main12'};

fileNames= {'Quiet_B1','Quiet_B2','SHL_B1','SHL_B2',...
    'Noise60_B1','Noise60_B2','Noise70_B1','Noise70_B2'};

talkerID= {'talker1','talker2'};
fileNamesSpeech= {'NH-Quiet_Rep1','NH-Quiet_Rep2','SHL-Quiet_Rep1','SHL-Quiet_Rep2','NH-Noise60_Rep1','NH-Noise60_Rep2','NH-Noise70_Rep1','NH-Noise70_Rep2'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters
% analysis parameter
param                   = struct();
param.fs                = 100;
param.fsr               = 20;
param.debug             = 0;
param.order             = 30;
param.factor            = 2;
param.offset            = 0;
param.startT            = 2;
param.endT              = 110;
param.window            = 0.05;         % time window before and after artifact
param.threshold         = 30;
param.hp.fc             = 1;
param.hp.n              = 106;
param.lp.fc             = 9;
param.lp.n              = 40;
param.qualityThr        = 0.40;
param.limLength         = 40;
param.weakQuality       = 0;
param.betterIdx         = 0;
param.baselineDur       = 20;    %sec
param.fixThreshold      = [0.4 50];
param.angleThreshold    = [5 80];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% windowLen= round(1*param.fs);  % window length of 1 sec
% cleanPupil= zeros(length(groupNames)*length(subNames)*length(fileNames),270*param.fs);

for iGroup= 1:length(groupNames)
%     figure
    for iFile= 1:length(fileNames)
         fileDir= fullfile(Path,[groupNames{iGroup},'\Vicon\',fileNames{iFile},'_tobii.mat']);

        if isfile(fileDir)

            data= [];  data= load(fileDir).data(param.baselineDur*param.fs:end,:);

             % % Extract Gaze3D
            gazeIdxTobii1 = [18 19 20];
            gazeIdxTobii2 = [35 36 37];

            idxGaze= [gazeIdxTobii1,gazeIdxTobii2];

            %%%% Due to problem we have for Main2_part2
            if size(data,2) > 20
                gazeInfoTmp= []; gazeInfoTmp= data(:,idxGaze);
            else

                gazeInfoTmp= []; gazeInfoTmp= data(:,[idxGaze(1),idxGaze(2),idxGaze(3)]);
                gazeInfoTmp(:,[4,5,6])= 0;
            end


            % Read speech data and check the length of speech and
                % gaze data

                tmp1= []; y = []; 
                tmp1= strcat(Path,groupNames{iGroup},'\Speech\Speech',fileNamesSpeech{iFile},'_',talkerID{1},'.wav');
                [y,fs]= audioread(tmp1);

                diffTime = (length(y)/48000)- (length(gazeInfoTmp)/param.fs);
                if (diffTime < 0)
                    gazeInfo = [];
                    gazeInfo = gazeInfoTmp(abs(diffTime*param.fs):end,:);
                else
                    % Zero padding to the end as it is possible that data
                    % recording stopped in between of test
                    gazeInfo = [];
                    gazeInfo = [gazeInfoTmp;zeros(round(diffTime*param.fs),6)];
%                     remainPitch = pitchTmpCom((param.baselineDur-diffTime)*param.fs:(param.baselineDur*param.fs),:);
%                     pitch= [remainPitch;pitchTmpNew];

                end

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % compute the percentage of missing data
            for iCol= 1: size(gazeInfo,2)
                gazeInfo(find(isnan(gazeInfo(:,iCol))),iCol)= 0;

                if rem(iCol,3)==1
                    missingPer(iGroup,iFile,round(iCol/3)+1)= length(find(~gazeInfo(:,iCol)))/size(gazeInfo,1);
                end
            end

            %%%%% Check missing data

            % Session with less than 60% missing data will be kept for
            % analysis
            tmpMiss= []; tmpMiss= squeeze(missingPer(iGroup,iFile,:));

            for iSub= 1:length(tmpMiss)
                if tmpMiss(iSub) <= param.qualityThr

                %%% Remove eye artifacts
                gazeInfoNew= [];
                [gazeInfoNew,param]= eyeArtifact(gazeInfo(:,(iSub-1)*3+1:iSub*3),param);
                gazeData{(iGroup-1)*2+iSub,iFile} = gazeInfoNew; 

                % construct time vector
                param.time = linspace(0,size(gazeInfoNew,1)/param.fs,size(gazeInfoNew,1));

                [gazeInfoNew,param] = eyeResample(gazeInfoNew,param);
                [param]             = eyeAngle(gazeInfoNew,param);
                [param]             = eyeBandpass(param.pitch,param);
                [velocity,param]    = eyeVelocity(param.pitch,param); % represents the speed of changes in pitch
                [param]             = eyeFilter(velocity,param);
                [saccade]           = eyeSaccade(param.velocityFilt,param.pitch,param);
                [fixation]          = eyeFixation(saccade.vector,param);

            else  % for invalid data

                display(['Invalid Data of :',groupNames{iGroup},'_',...
                    fileNames{iFile},'_',num2str(rem(iGroup,2))])
                gazeData{(iGroup-1)*2+iSub,iFile} = [];

            end

            saccadeProperties{(iGroup-1)*2+iSub,iFile}  = saccade;
            fixationProperties{(iGroup-1)*2+iSub,iFile} = fixation;

            end

        else % for existance of files
            display(['No File for analysis:',groupNames{iGroup},'_',...
                fileNames{iFile},'_',num2str(rem(iGroup,2))])
            gazeData{(iGroup-1)*2+iSub,iFile} = [];
        end % For existance of files
    end % for ifile
end  % for length of datafile_Rep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Statistical analysis of saccade
for iSub= 1:length(groupNames)*2
    for iFile= 1:length(fileNames)
        if ~isempty(saccadeProperties{iSub,iFile})
            
            saccadeNum(iSub,iFile)         = saccadeProperties{iSub,iFile}.saccadeNum;
            saccadeStart{iSub,iFile}       = saccadeProperties{iSub,iFile}.saccadeStart;
            saccadeEnd{iSub,iFile}         = saccadeProperties{iSub,iFile}.saccadeEnd;
            saccadeDur(iSub,iFile)         = sum(saccadeProperties{iSub,iFile}.saccadeDur);
            saccadeDurMean(iSub,iFile)     = mean(saccadeProperties{iSub,iFile}.saccadeDur);
            velocityAmpMeanSac(iSub,iFile) = mean(saccadeProperties{iSub,iFile}.vA);
            velocityAmpStdSac(iSub,iFile)  = std(saccadeProperties{iSub,iFile}.vA);
            pitchAmpMeanSac(iSub,iFile)    = mean(saccadeProperties{iSub,iFile}.sA);
            velocityAmpStdSac(iSub,iFile)  = std(saccadeProperties{iSub,iFile}.sA);
        else
            saccadeNum(iSub,iFile)         = 0;
            saccadeDur(iSub,iFile)         = 0;
            saccadeDurMean(iSub,iFile)     = 0;
            velocityAmpMeanSac(iSub,iFile) = 0;
            velocityAmpStdSac(iSub,iFile)  = 0;
            pitchAmpMeanSac(iSub,iFile)    = 0;
            velocityAmpStdSac(iSub,iFile)  = 0;
        end
    end
end

%% Saccade analysis based on speaking and non-speaking time (no-overlap) 

load('C:\Users\sulb\OneDrive - Demant\Documents\IFD_Project\ProjectDescription\Programs\MainPrograms\Version 29072022\Results\utterancesSaccade.mat')

for iGroup= 2:length(groupNames)
    for iFile= 1:length(fileNames)
        tmpCH1 = []; idxCH1= []; turnCH1 = [];
        startTurnCH1 = []; endTurnCH1 = [];
        tmpCH1       = Utterances{iGroup,iFile}.turnCH1;
        % Consider turns with duration more than 500 ms
        idxCH1       = find(tmpCH1(:,1) > 1);
        turnCH1      = tmpCH1(idxCH1,:);
        startTurnCH1 = turnCH1(:,2)*param.fsr*0.004;
        endTurnCH1   = turnCH1(:,3)*param.fsr*0.004;

        tmpCH2 = []; idxCH2= []; turnCH2 = [];
        startTurnCH2 = []; endTurnCH2 = [];

        tmpCH2       = Utterances{iGroup,iFile}.turnCH2;
        % Consider turns with duration more than 500 ms
        idxCH2       = find(tmpCH2(:,1) > 1);
        turnCH2      = tmpCH2(idxCH2,:);
        startTurnCH2 = turnCH2(:,2)*param.fsr*0.004;
        endTurnCH2   = turnCH2(:,3)*param.fsr*0.004;

         if ~isempty(gazeData{(iGroup-1)*2+1,iFile})
            timeCH1 = [];  timeCH2 = [];
            timeCH1      = zeros(1,length(gazeData{(iGroup-1)*2+1,iFile}));
            
            for i=1:length(startTurnCH1); timeCH1(round(startTurnCH1(i))+1:round(endTurnCH1(i))+1) = 1; end

            % Saccade analysis
            % find when saccade start and end for each talker
            saccadeCH1 = [];  saccCH1Start = []; saccCH1End = [];
            saccadeCH1  = zeros(1,length(gazeData{(iGroup-1)*2+1,iFile}));
            saccCH1Start= saccadeStart{(iGroup-1)*2+1,iFile};
            saccCH1End  = saccadeEnd{(iGroup-1)*2+1,iFile};

            % Merge saccades if the onset difference of two saccade in sequence is less than 300 ms
            % (saccade duration is normally between 50-250 ms)
            % Then update saccade start and end matrix

            for i=1:length(saccCH1Start) 
                saccadeCH1(saccCH1Start(i,1):saccCH1End(i,1)) = 1; 
                if i<(length(saccCH1Start)-1)
                    if saccCH1Start(i+1,1)-saccCH1Start(i,1)< (0.3*param.fsr)
                       saccadeCH1(saccCH1End(i,1):saccCH1Start(i+1,1)) = 1;
                    end
                end
            end  % saccade talker 1

            % Update start and end saccade time
            saccCH1StartUpdate = []; saccCH1EndUpdate = [];
            dsaccadeCH1 = []; 
            dsaccadeCH1 = diff(saccadeCH1);
            saccCH1StartUpdate = find(dsaccadeCH1 ==1)+1;
            saccCH1EndUpdate   = find(dsaccadeCH1 ==-1);

            % Plot saccade for turn taking
%             figure('Name','Saccade talker1'),
%             plot(pitchData{(iGroup-1)*2+1,iFile});  % pitch of talker1
%             hold on
%             area(timeCH2*max(pitchData{(iGroup-1)*2+1,iFile}),'FaceColor', 'g',  'FaceAlpha', 0.25)  % plot turn talker 2
%             area(saccadeCH1*param.fsr,'FaceColor', 'r',  'FaceAlpha', 0.25)                          % plot saccade talker 1 in talker 2 turn
% 
%             figure('Name','Saccade talker2'),
%             plot(pitchData{(iGroup-1)*2+2,iFile});  % pitch of talker1
%             hold on
%             area(timeCH1*max(pitchData{(iGroup-1)*2+1,iFile}),'FaceColor', 'g',  'FaceAlpha', 0.25)  % plot turn talker 2
%             area(saccadeCH2*param.fsr,'FaceColor', 'r',  'FaceAlpha', 0.25)                          % plot saccade talker 1 in talker 2 turn

             % Find out the number of talker' saccades in speaking time of each
             % talker; means saccade of talker 1 when the same talker is
             % speaking while interlosutor speaking

             idxSacT1                               = find(timeCH1>0 & saccadeCH1>0);  %Indexes of T1 saccades when T1 is speaking
             saccadeSelfTurn((iGroup-1)*2+1,iFile)  = length(find(diff(idxSacT1)>1));
             

            % Find Onset/offset saccade delay of talker1 referring to speech
            % start of talker2
            onsetDiff = [];  offsetDiff = []; utterDurCH2 = [];
            for i=1:length(saccCH1StartUpdate)
                tmp = []; tmp = startTurnCH2-saccCH1StartUpdate(i);
                idx = []; idx = find(abs(tmp)== min (abs(tmp)));
                onsetDiff(i)  = -1*(tmp(idx(1))/param.fsr);
                utterDurCH2(i)= turnCH2(idx(1),1);
%                 tmp= [];  tmp = startTurnCH2-saccCH1End(i);
%                 idx = []; idx = find(abs(tmp)== min (abs(tmp)));
                offsetDiff(i) = -1*((startTurnCH2(idx(1))-saccCH1EndUpdate(i))/param.fsr);

            end

            distrOnset{(iGroup-1)*2+1,iFile}  = onsetDiff(onsetDiff<3 & onsetDiff>-3);
            distrOffset{(iGroup-1)*2+1,iFile} = offsetDiff(onsetDiff<3 & onsetDiff>-3);

            utterSacDur{(iGroup-1)*2+2,iFile}  = utterDurCH2;
          
            % Plot histogram of distribution
            %         [N,edges] = histcounts(distrOnset{3, 8}, 'Normalization','count');
            %         edges = edges(2:end) - (edges(2)-edges(1))/2;
            %         plot(edges, N);
            % d=pdist2(c1,c2);
            % ComparisonHistN({distrOnset{3,1},distrOnset{3,8}},edges(2:end),{[],[]},[0.05,0.05])



        else
            display(['No pitch data of :',num2str((iGroup-1)*2+1),...
                num2str(iFile)])
 end

  if ~isempty(gazeData{(iGroup-1)*2+2,iFile})
            timeCH2 = [];
            timeCH2      = zeros(1,length(gazeData{(iGroup-1)*2+2,iFile}));
            for i=1:length(startTurnCH2); timeCH2(round(startTurnCH2(i))+1:round(endTurnCH2(i))+1) = 1; end

            % Saccade analysis
            % find when saccade start and end for each talker
            saccadeCH2 = [];  saccCH2Start = []; saccCH2End = [];
            saccadeCH2  = zeros(1,length(gazeData{(iGroup-1)*2+2,iFile}));
            saccCH2Start= saccadeStart{(iGroup-1)*2+2,iFile};
            saccCH2End  = saccadeEnd{(iGroup-1)*2+2,iFile};

            % Merge saccades if the onset difference of two saccade in sequence is less than 300 ms
            % (saccade duration is normally between 50-250 ms)
            % Then update saccade start and end matrix
            for i=1:length(saccCH2Start) 
                saccadeCH2(saccCH2Start(i,1):saccCH2End(i,1)) = 1; 
                if i<(length(saccCH2Start)-1)
                    if saccCH2Start(i+1,1)-saccCH2Start(i,1)< (0.3*param.fsr)
                       saccadeCH2(saccCH2End(i,1):saccCH2Start(i+1,1)) = 1;
                    end
                end
            end  % saccade talker 2

            % Update start and end saccade time
            saccCH2StartUpdate = []; saccCH2EndUpdate = [];
            dsaccadeCH2 = []; 
            dsaccadeCH2 = diff(saccadeCH2);
            saccCH2StartUpdate = find(dsaccadeCH2 ==1)+1;
            saccCH2EndUpdate   = find(dsaccadeCH2 ==-1);
            
         
             % Find out the number of talker' saccades in speaking time of each
             % talker; means saccade of talker 1 when the same talker is
             % speaking while interlosutor speaking

             
             idxSacT2                               = find(timeCH2>0 & saccadeCH2>0);  %Indexes of T2 saccades when T2 is speaking
             saccadeSelfTurn((iGroup-1)*2+2,iFile)  = length(find(diff(idxSacT2)>1));

       

            % Find Onset/offset saccade delay of talker2 referring to speech
            % start of talker1
            onsetDiff = [];  offsetDiff = []; utterDurCH1 = [];
            for i=1:length(saccCH2StartUpdate)
                tmp = []; tmp = startTurnCH1-saccCH2StartUpdate(i);
                idx = []; idx = find(abs(tmp)== min (abs(tmp)));
                onsetDiff(i)  = -1*(tmp(idx(1))/param.fsr);
                utterDurCH1(i)= turnCH1(idx(1),1);
%                 tmp= [];  tmp = startTurnCH2-saccCH2End(i);
%                 idx = []; idx = find(abs(tmp)== min (abs(tmp)));
                offsetDiff(i) = -1*((startTurnCH1(idx(1))-saccCH2EndUpdate(i))/param.fsr);

            end


            distrOnset{(iGroup-1)*2+2,iFile}  = onsetDiff(onsetDiff<3 & onsetDiff>-3);
            distrOffset{(iGroup-1)*2+2,iFile} = offsetDiff(onsetDiff<3 & onsetDiff>-3);

            utterSacDur{(iGroup-1)*2+1,iFile}  = utterDurCH1;

        else
            display(['No pitch data of :',num2str((iGroup-1)*2+2),...
                num2str(iFile)])
  end
   if ~isempty(gazeData{(iGroup-1)*2+1,iFile}) &&  ~isempty(gazeData{(iGroup-1)*2+2,iFile})
           
           idxSacT1inT2                           = find(timeCH2>0 & saccadeCH1>0);  %Indexes of talker1 (T1) saccades when T2 is speaking
             saccadeOtherTurn((iGroup-1)*2+1,iFile) = length(find(diff(idxSacT1inT2)>1));
             idxSacT2inT1                           = find(timeCH1>0 & saccadeCH2>0);  %Indexes of T2 saccades when T1 is speaking
             saccadeOtherTurn((iGroup-1)*2+2,iFile) = length(find(diff(idxSacT2inT1)>1));
   end
    end
end

%% Plot distribution of head pitch in different noise conditions
onsetN0  = []; offsetN0  = [];
onsetN60 = []; offsetN60 = [];
onsetN70 = []; offsetN70 = [];
onsetSHL = []; offsetSHL = [];

for iGroup = 2:length(groupNames)
    for iFile= 1:length(fileNames)

        if ((iFile == 1) )     
                 onsetN0 = [onsetN0,distrOnset{(iGroup-1)*2+1,iFile},...
                     distrOnset{(iGroup-1)*2+2,iFile}];
                 offsetN0 = [offsetN0,distrOffset{(iGroup-1)*2+1,iFile},...
                     distrOffset{(iGroup-1)*2+2,iFile}];
        elseif ( (iFile == 4))
                onsetSHL = [onsetSHL,distrOnset{(iGroup-1)*2+1,iFile},...
                     distrOnset{(iGroup-1)*2+2,iFile}];
                offsetSHL = [offsetSHL,distrOffset{(iGroup-1)*2+1,iFile},...
                     distrOffset{(iGroup-1)*2+2,iFile}];
        elseif ((iFile == 6))
                onsetN60 = [onsetN60,distrOnset{(iGroup-1)*2+1,iFile},...
                     distrOnset{(iGroup-1)*2+2,iFile}];
                offsetN60 = [offsetN60,distrOffset{(iGroup-1)*2+1,iFile},...
                     distrOffset{(iGroup-1)*2+2,iFile}];
        elseif ((iFile == 7))
                onsetN70 = [onsetN70,distrOnset{(iGroup-1)*2+1,iFile},...
                     distrOnset{(iGroup-1)*2+2,iFile}];
                offsetN70 = [offsetN70,distrOffset{(iGroup-1)*2+1,iFile},...
                     distrOffset{(iGroup-1)*2+2,iFile}];
        end
    end
end

saccadeSpeech = struct('onsetN0',{onsetN0},'offsetN0',{offsetN0},...
    'onsetN60',{onsetN60},'offsetN60',{offsetN60},...
    'onsetN70',{onsetN70},'offsetN70',{offsetN70},...
    'onsetSHL',{onsetSHL},'offsetSHL',{offsetSHL});

%% Interpolation of data
for iSub= 1:length(groupNames)*2
    for iFile= 1:length(fileNames)
        if ~isempty(fixationProperties{iSub,iFile})
            
            fixationNum(iSub,iFile)  = fixationProperties{iSub,iFile}.fixationNum;
            fixationDur(iSub,iFile)  = sum(fixationProperties{iSub,iFile}.fixationDur);
            fixationDurMean(iSub,iFile)  = mean(fixationProperties{iSub,iFile}.fixationDur);
            
        else
            fixationNumFix(iSub,iFile)  = 0;
            fixationDurFix(iSub,iFile)  = 0;
            fixationDurMean(iSub,iFile) = 0;
            
        end
    end
end

%%%% Interpolation for missing data

saccadeNumInter = zeros(size(saccadeNum));
for col = 1:size(saccadeNumInter,2)
    b= []; b = saccadeNum(:,col)~=0;
    Y= []; Y = cumsum(b-diff([1;b])/2);
    saccadeNumInter(:,col) = interp1(1:nnz(b),saccadeNum(b,col),Y,'linear');
    if (saccadeNum(end,col) == 0)
        saccadeNumInter(end,col)= (saccadeNumInter(end-1,col)+ saccadeNumInter(end-2,col))/2;
    end
end

saccadeDurInter = zeros(size(saccadeDur));
for col = 1:size(saccadeDurInter,2)
    b= []; b = saccadeDur(:,col)~=0;
    Y= []; Y = cumsum(b-diff([1;b])/2);
    saccadeDurInter(:,col) = interp1(1:nnz(b),saccadeDur(b,col),Y,'linear');
    if (saccadeDur(end,col) == 0)
        saccadeDurInter(end,col)= (saccadeDurInter(end-1,col)+ saccadeDurInter(end-2,col))/2;
    end
end

saccadeDurMeanInter = zeros(size(saccadeDurMean));
for col = 1:size(saccadeDurMeanInter,2)
    b= []; b = saccadeDurMean(:,col)~=0;
    Y= []; Y = cumsum(b-diff([1;b])/2);
    saccadeDurMeanInter(:,col) = interp1(1:nnz(b),saccadeDurMean(b,col),Y,'linear');
    if (saccadeDurMean(end,col) == 0)
        saccadeDurMeanInter(end,col)= (saccadeDurMeanInter(end-1,col)+ saccadeDurMeanInter(end-2,col))/2;
    end
end

fixationNumInter = zeros(size(fixationNum));
for col = 1:size(fixationNumInter,2)
    b= []; b = fixationNum(:,col)~=0;
    Y= []; Y = cumsum(b-diff([1;b])/2);
    fixationNumInter(:,col) = interp1(1:nnz(b),fixationNum(b,col),Y,'linear');
    if (fixationNum(end,col) == 0)
        fixationNumInter(end,col)= (fixationNumInter(end-1,col)+ fixationNumInter(end-2,col))/2;
    end
end

fixationDurInter = zeros(size(fixationDur));
for col = 1:size(fixationDurInter,2)
    b= []; b = fixationDur(:,col)~=0;
    Y= []; Y = cumsum(b-diff([1;b])/2);
    fixationDurInter(:,col) = interp1(1:nnz(b),fixationDur(b,col),Y,'linear');
    if (fixationDur(end,col) == 0)
        fixationDurInter(end,col)= (fixationDurInter(end-1,col)+ fixationDurInter(end-2,col))/2;
    end
end

fixationDurMeanInter = zeros(size(fixationDurMean));
for col = 1:size(fixationDurMeanInter,2)
    b= []; b = fixationDurMean(:,col)~=0;
    Y= []; Y = cumsum(b-diff([1;b])/2);
    fixationDurMeanInter(:,col) = interp1(1:nnz(b),fixationDurMean(b,col),Y,'linear');
    if (fixationDurMean(end,col) == 0)
        fixationDurMeanInter(end,col)= (fixationDurMeanInter(end-1,col)+ fixationDurMeanInter(end-2,col))/2;
    end
end


% save(fullfile(savePath,['fixationNum.mat']),'fixationNumInter')
% save(fullfile(savePath,['fixationDur.mat']),'fixationDurInter')
% save(fullfile(savePath,['fixationDurMean.mat']),'fixationDurMeanInter')
save(fullfile(savePath,['EyesaccadeNum.mat']),'saccadeNumInter')
save(fullfile(savePath,['EyesaccadeDur.mat']),'saccadeDurInter')
save(fullfile(savePath,['EyesaccadeDurMean.mat']),'saccadeDurMeanInter')

% save(fullfile(savePath,['saccadeSpeech_eyeTest.mat']),'saccadeSpeech')
