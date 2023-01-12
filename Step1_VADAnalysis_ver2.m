
%%%%%% MAIN PROGRAM FOR SPEECH ANALYSIS %%%%%%%%
%%%% The difference between these two versions is recent version contain
%%%% correct formant analysis (only include voiced parts)
%%%% Main 1: issues in file 3 and 4 (Please refer to DataOverview.xlsx)
%%%% Main 10: Issue in file 3 and 4 (Please refer to DataOverview.xlsx)
%%%% This program compute overlap duration, dB SPL and formant analysis
%%%% Written by SULB
%%%% Last modified: 22.08.2022
%%% To find ou the effect of background noise, we need to compute splSilence 
%%% in both noise conditions of 60 and 70 dBA
%%% for noise 60, the SPL increased with 4 dB from the quiet conditions and for noise 70 it was around 8 dB 


clear all
close all
clc

Path=pwd;
% savePath= 'C:\Users\sulb\OneDrive - Demant\Documents\IFD_Project\ProjectDescription\Programs\MainPrograms\Version 29072022\Results\';
groupName= {'Main1','Main2','Main3','Main4','Main5','Main6',...
    'Main7','Main8','Main9','Main10','Main11','Main12'};
fileNames= {'NH-Quiet_Rep1','NH-Quiet_Rep2','SHL-Quiet_Rep1','SHL-Quiet_Rep2','NH-Noise60_Rep1','NH-Noise60_Rep2','NH-Noise70_Rep1','NH-Noise70_Rep2'};
talkerID= {'talker1','talker2'};

% Parameters
fs= 48000;
param                   = struct();
param.utteranceCH1      = [];
param.utteranceCH2      = [];
param.binRes            = 0.0040;



for iGroup= 2
    for iFile= 1:length(fileNames)
        
        display([groupName{iGroup},'_',fileNames{iFile}])

        param= struct('utteranceCH1',[],'utteranceCH2',[],...
            'binRes',0.0040);
        tmp1= []; tmp1= strcat(Path,'\audio\',groupName{iGroup},'\Speech',fileNames{iFile},'_',talkerID{1},'.wav');
        [y1,fs]= audioread(tmp1);

        % %%% Added by SULB %%%%%
        % [y1,fs]= audioread('C:\Users\sulb\OneDrive - Demant\Documents\IFD_Project\ProjectDescription\Programs\Pilot_Tst\convexchange_classification\convexchange_classification\talker1_noise.wav');
        % % data= audioread('C:\Users\sulb\OneDrive - Demant\Documents\IFD_Project\ProjectDescription\Programs\Pilot_Tst\Recordings_2channel_22050Hz_16bit\pair1\pair1_L1-noise_rep1.wav');
        %%%%%%%%%%%%%%%%%%%%%%%
        [actArr1, t1,buffSig1] = voiceActivityDetection(y1, fs);
%         thereshold(iGroup,iFile,1)= THF2;
        % plot(d1)
        % hold on
        % plot(t.*fs, actArr1)

        % %% d2 from the second talker
        % [y2,fs]= audioread('C:\Users\sulb\OneDrive - Demant\Documents\IFD_Project\ProjectDescription\Programs\Pilot_Tst\convexchange_classification\convexchange_classification\talker2_noise.wav');
        % % data= audioread('C:\Users\sulb\OneDrive - Demant\Documents\IFD_Project\ProjectDescription\Programs\Pilot_Tst\Recordings_2channel_22050Hz_16bit\pair1\pair1_L1-noise_rep1.wav');
        %%%%%%%%%%%%%%%%%%%%%%%
        tmp2= []; tmp2= strcat(Path,'\audio\',groupName{iGroup},'\Speech',fileNames{iFile},'_',talkerID{2},'.wav');
        [y2,fs]= audioread(tmp2);

        fs= 48000;
        [actArr2, t2,buffSig2] = voiceActivityDetection(y2, fs);
%         thereshold(iGroup,iFile,2)= THF2;
        % plot(d1)
        % hold on
        % plot(t.*fs, actArr2)

        %%%%%%%%%%%%%%%%%%%% Analysis of overlap, gap, pause and utterances %%%%%%%%%
        actArr= [actArr1 ; actArr2];
        [overlapWCH1, overlapWCH2, overlapBCH1, overlapBCH2, gapCH1, ...
            gapCH2, pauseCH1, pauseCH2, utteranceCH1, utteranceCH2, turnCH1, turnCH2] =communicativeStateClassification(actArr, param.binRes);
        overLapLength((iGroup-1)*2+1,iFile)= sum(overlapWCH1(:,1))+ sum(overlapBCH1(:,1));
        overLapLength((iGroup-1)*2+2,iFile)= sum(overlapWCH2(:,1))+ sum(overlapBCH2(:,1));

        %%% Seperate overlap between and within

        overLapWithinLength((iGroup-1)*2+1,iFile) = sum(overlapWCH1(:,1));
        overLapBetweenLength((iGroup-1)*2+1,iFile)= sum(overlapBCH1(:,1));

        overLapWithinLength((iGroup-1)*2+2,iFile) = sum(overlapWCH2(:,1));
        overLapBetweenLength((iGroup-1)*2+2,iFile)= sum(overlapWCH2(:,1));

        param.utteranceCH1= utteranceCH1;
        param.utteranceCH2= utteranceCH2;
        Utterances{iGroup,iFile}= param;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%% Start VAD with play
        play_sound = 0;
        gain_sound = 10;

        i = 1;
        idx_left = 1;
        idx_right = 2;

        %% VAD analysis
        sampleBin = param.binRes/(1/fs);
        len_trial = round(size(y1,1)/sampleBin);
        rest_trial = size(y1,1)-len_trial*sampleBin;
        t_vect_aud = 1/fs:1/fs:(len_trial*sampleBin)/fs;

        if rest_trial<0
            t_vect_aud_final= []; t_vect_aud_final= t_vect_aud(1:length(t_vect_aud)+rest_trial);
        else
            t_vect_aud_final= []; t_vect_aud_final= t_vect_aud;
        end

        t_vect_act = param.binRes:param.binRes:len_trial*param.binRes;
        vect_first = y1(1:size(t_vect_aud_final,2));

        %%%%%%%%%%%%%%%% Plot the results of VAD
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                h_seperate= figure('Position', [10 10 2000 550]);
        
                subplot(2,1,1)
        
                plot(t_vect_aud_final,vect_first, 'b')
                hold on
                area(t_vect_act(1:length(actArr1)), max(vect_first)*double(actArr1),'FaceColor', 'b',  'FaceAlpha', 0.25)
                ylabel('Amplitude','FontSize',12,'FontWeight','bold')
                xlabel('Time [s]','FontSize',12,'FontWeight','bold')
                title('talker1')
                set(gca,'YLim', [-1.25*max(vect_first), 1.25*max(vect_first)]);
                set(gca,'FontSize',14,'FontWeight','bold')
                box off
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %         h_both= figure;
        %         hold on
        %         area(t_vect_act(1:length(actArr1)), max(vect_first)*double(actArr1),'FaceColor', 'b',  'FaceAlpha', 0.25)
        %         ylabel('Amplitude')
        %         xlabel('Time [s]')
        %         set(gca,'YLim', [-1.25*max(vect_first), 1.25*max(vect_first)]);

        %%%%%%% Check if the signal cleaned properly by thresholding on
        %%%%%%% power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        idx= []; idx= find(actArr1);
        startPoint= []; startPoint= (idx-1)*4+1;  % 4 is the real duration of window
        endPoint= []; endPoint= startPoint+5;

        yfirst= []; yFirst= zeros(length(y1),1);
        for i= 1:length(idx)
            yFirst(round((startPoint(i)/1000)*fs):round((endPoint(i)/1000)*fs),1)= 1;
        end

        tmp= []; tmp= yFirst(1:length(y1),1).*y1;
        yyClean= [];  yyClean= y1(find(tmp~=0));
        yySilence= []; yySilence= y1(find(tmp==0));


        %%%%%%%%%%%%% Formant Frequency %%%%%%%%%%%%%

        [freqMean(iGroup,iFile,1),...
            freqStd(iGroup,iFile,1),...
            freqIQ(iGroup,iFile,1)]= formantAnalysis(yyClean,fs);
        [freqMeanSilence(iGroup,iFile,1),...
            freqStdSilence(iGroup,iFile,1),...
            freqIQSilence(iGroup,iFile,1)]= formantAnalysis(yySilence,fs);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%% Sound Pressure Level (dB SPL) %%%%%%%%%%

            pow=[]; pow= sqrt(mean(yyClean.^2));
            spl=[]; spl= 20.*log10(pow/(20e-6));
            powSilence=[]; powSilence= sqrt(mean(yySilence.^2));
            splSilence=[]; splSilence= 20.*log10(powSilence/(20e-6));
            splAll(iGroup,iFile,1)= spl;
            splAllSilence(iGroup,iFile,1)= splSilence;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %         for i= 1:length(idx)
            %             yyClean= [yyClean; y1(round((startPoint(i)/1000))*fs:round((endPoint(i)/1000))*fs)];
            %         end
            %         display('It is done')
            %
            yy= []; yy = gain_sound*yyClean;
            player = audioplayer(yy, fs);
%                     play(player)
            %         clear player
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%
            %calculate timeline

            %%%%%%%%%%%%%%% The second talker
            len_trial= []; len_trial = round(size(y2,1)/sampleBin);
            rest_trial= []; rest_trial = size(y2,1)-len_trial*sampleBin;
            t_vect_aud= []; t_vect_aud = 1/fs:1/fs:(len_trial*sampleBin)/fs;
            if rest_trial<0
            t_vect_aud_final= []; t_vect_aud_final= t_vect_aud(1:length(t_vect_aud)+rest_trial);
            else
                t_vect_aud_final= []; t_vect_aud_final= t_vect_aud;
            end
            t_vect_act= []; t_vect_act = param.binRes:param.binRes:len_trial*param.binRes;
            vect_first= []; vect_first = y2(1:size(t_vect_aud_final,2));

            %%%%%%%%%%%%%%% Plot the results of VAD %%%%%%%%%%%%%%%%%%%%%%

%                     figure(h_seperate)
%                     subplot(2,1,2)
%                     plot(t_vect_aud_final,vect_first, 'r')
%                     hold on
%                     area(t_vect_act(1:length(actArr2)), max(vect_first)*double(actArr2),'FaceColor', 'r',  'FaceAlpha', 0.25)
%                     ylabel('Amplitude','FontSize',12,'FontWeight','bold')
%                     xlabel('Time [s]','FontSize',12,'FontWeight','bold')
%                     title('talker2')
%                     set(gca,'YLim', [-1.25*max(vect_first), 1.25*max(vect_first)]);
%                     box off
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %         figure(h_both)
            %         area(t_vect_act(1:length(actArr2)), max(vect_first)*double(actArr2),'FaceColor', 'r',  'FaceAlpha', 0.25)
            %         ylabel('Amplitude')
            %         xlabel('Time [s]')
            %         set(gca,'YLim', [-1.25*max(vect_first), 1.25*max(vect_first)]);

            %%%%%%% Check if the signal cleaned properly by thresholding on
            %%%%%%% power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            idx= []; idx= find(actArr2);
            startPoint= []; startPoint= (idx-1)*4+1;
            endPoint= []; endPoint= startPoint+5;

            yfirst= []; yFirst= zeros(length(y2),1);
            for i= 1:length(idx)
                yFirst(round((startPoint(i)/1000)*fs):round((endPoint(i)/1000)*fs),1)= 1;
            end

            tmp= []; tmp= yFirst(1:length(y2),1).*y2;
            yyClean= [];  yyClean= y2(find(tmp~=0));
            yySilence= []; yySilence= y2(find(tmp==0));

            %%%%%%%%%%%%% Formant Frequency %%%%%%%%%%%%%

            [freqMean(iGroup,iFile,2),...
                freqStd(iGroup,iFile,2),...
                freqIQ(iGroup,iFile,2)]= formantAnalysis(yyClean,fs);
             [freqMeanSilence(iGroup,iFile,2),...
                freqStdSilence(iGroup,iFile,2),...
                freqIQSilence(iGroup,iFile,2)]= formantAnalysis(yySilence,fs);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%% Sound Pressure Level (dB SPL) %%%%%%%%%%

                pow=[]; pow= sqrt(mean(yyClean.^2));
                spl=[]; spl= 20.*log10(pow/(20e-6));
                powSilence=[]; powSilence= sqrt(mean(yySilence.^2));
                splSilence=[]; splSilence= 20.*log10(powSilence/(20e-6));
                splAll(iGroup,iFile,2)= spl;
                splAllSilence(iGroup,iFile,2)= splSilence;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                %         for i= 1:length(idx)
                %             yyClean= [yyClean; y1(round((startPoint(i)/1000))*fs:round((endPoint(i)/1000))*fs)];
                %         end
                %         display('It is done')
                %
                yy= []; yy = gain_sound*yyClean;
                        player = audioplayer(yy, fs);
%                         play(player)
                %         clear player
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




                %%
                % %%%%%%%%%%% Only overlap duration (both within and between) %%%%%%%%%
                %
                overlapCnt= 0; startPoint= []; endPoint= [];
                overlapidx= find(and(actArr1,actArr2));
                startPointOverlap= (overlapidx-1)*4+1;
                overLapLengthTest(iGroup,iFile)= (length(startPointOverlap)*5)/1000;  % s
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                end
end
%%
% Bar plot based on noise level

noiseCond(1,1)=mean(overLapLength(4,1)); % Repetition 1 of Quiet
noiseCond(2,1)=mean(overLapLength(4,2));
noiseCond(1,2)=mean(overLapLength(4,3));
noiseCond(2,2)=mean(overLapLength(4,4));
noiseCond(1,3)=mean(overLapLength(4,5));
noiseCond(2,3)=mean(overLapLength(4,6));
noiseCond(1,4)=mean(overLapLength(4,7));
noiseCond(2,4)=mean(overLapLength(4,8));


X = categorical({'Quiet','SHL','Noise60','Noise70'});
X = reordercats(X,{'Quiet','SHL','Noise60','Noise70'});
figure,
b= bar(X,noiseCond);
legend('Rep1','Rep2')
%%
% Seperate subjects in overlapLength matrix for LMM statistics
overLapLengthStat(1:12,:)  = overLapLength(1:2:end,:);
overLapLengthStat(13:24,:) = overLapLength(2:2:end,:);

splStat (1:12,:)  = splAll(:,:,1);
splStat (13:24,:) = splAll(:,:,2);

freqMeanStat(1:12,:)  = freqMean(:,:,1);
freqMeanStat(13:24,:) = freqMean(:,:,2);

freqStdStat(1:12,:)  = freqStd(:,:,1);
freqStdStat(13:24,:) = freqStd(:,:,2);

freqIQStat(1:12,:)  = freqIQ(:,:,1);
freqIQStat(13:24,:) = freqIQ(:,:,2);

% save(fullfile(savePath,['Overlap1110.mat']),'overLapLengthStat')
% save(fullfile(savePath,['SPL1110.mat']),'splStat')
% save(fullfile(savePath,['formantMean1110.mat']),'freqMeanStat')
% save(fullfile(savePath,['formantStd1110.mat']),'freqStdStat')
% save(fullfile(savePath,['formantIQ1110.mat']),'freqIQStat')
% save(fullfile(savePath,['utterances1110.mat']),'Utterances')



