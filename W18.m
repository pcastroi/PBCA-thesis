%% PBCA-Thesis - Week 18,20,21 - AMEND II data - Pupillometry
% Pathing
clear all; clc; close all;
BPath = strsplit(pwd,'PBCA-thesis');
addpath('tools')
addpath([BPath{1} 'Pupil-preprocessing-tools\tools']) % For preprocessing
addpath([BPath{1} 'PUPILS-preprocessing-pipeline']) % For preprocessing

% Colors
SColor = [53, 155, 67]./255;
LColor = [204, 36, 0]./255;
NHColor = [204, 62, 0]./255;
HIColor = [164, 0, 204]./255;
UNColor = [1, 35, 33]./255;
AAColor = [4, 105, 100]./255;
ABColor = [7, 192, 182]./255;
QuietColor = [204, 152, 0]./255;
SHLColor = [123, 31, 162]./255;
N60Color = [0, 196, 215]./255;
N70Color = [2, 36, 223]./255;

% Important variables
[subDirs_I] = GetSubDirsFirstLevelOnly('data\AMEND_I');
[subDirs_II] = GetSubDirsFirstLevelOnly('data\AMEND_II');
subDirs_II(contains(subDirs_II,{'Pilot 1','Pair01','Pair14'}))=[]; % Only using data from Pair01 to Pair13, others removed
FileNames_I={'P1_N0_B1.mat','P1_N0_B2.mat','P1_N60_B1.mat','P1_N60_B2.mat','P1_Noise60_B1.mat','P1_Noise60_B2.mat','P1_Noise70_B1.mat','P1_Noise70_B2.mat','P2_N0_B1.mat','P2_N0_B2.mat','P2_N60_B1.mat','P2_N60_B2.mat','P2_Noise60_B1.mat','P2_Noise60_B2.mat','P2_Noise70_B1.mat','P2_Noise70_B2.mat'};
FileNames_II={'UNHI_N0.mat','UNHI_N60.mat','UNHI_N70.mat','AAHI_N0.mat','AAHI_N60.mat','AAHI_N70.mat','ABHI_N0.mat','ABHI_N60.mat','ABHI_N70.mat'};
LoadUtt_I=load('data\AMEND_I\utterances1110.mat');
LoadUtt_II=load('data\AMEND_II\Utterances0805.mat');
LoadDelays_I=load('data\AMEND_I\delays1110.mat');
% LoadDelays_II=load('data\AMEND_II\delays1110.mat'); % Should be no delay
Param.Fs = 50; % Sampling frequency of pupil data
Param.RemoveBeforeAndAfter = [0 0]*1e-3; % Samples within the time range before and after NaNs will set NaNs as well.
Param.MinLengthNaNRepair = 0; % Drop values (i.e., change to NaN) before and after NaNs only for contiguous NaNs of at least __ samples.
Param.Preblink = 0.1; % [s], set to NaN, time before blink
Param.Postblink = 0.2; % [s], set to NaN, time after blink
Param.BlinkThresh = 3; % [samples], threshold of samples in between artifacts or blinks
MaxVisualAngle = 1.5; % [degrees], critical visual angle for fixation definition
LPWinSize = 0.5; % [s]: Window size of hamming-window for low-pass filtering
FilterWidth = round((LPWinSize*Param.Fs)/2); % [samples]: Width of hamming filter used for fixation duration
LPWindow = hamming(round(LPWinSize*Param.Fs));
LPWindow = LPWindow/sum(LPWindow); % Hamming-window
AudFs = 48000; % [Hz], sampling frequency of the audio files
AdapBL = 0.3; % [s], Duration of baseline prior to event
TimeMinWin = 0.5; % [s], Minimum time of a window
TimeInitialMerge = 0.3; % [s], Time threshold for merging windows initially
TimeMerge = 2; % [s], Time threshold for merging windows after rejecting small windows
VisRatio = 0.8; % Percentage of data to visualize in plots
RejectRatio = 0.4; % Rejection threshold based on the ratio of NaNs in data
RejectDelay = 0.5; % [s], Rejection threshold based on delay between timestamps and n-samples
TimeStartW = 0.5; % [s], time before Utt/Lis starts
TimeEndW = 0; % [s], time after Utt/Lis starts
x = 1; % idx to store global values
NCond_II = numel(FileNames_II); % N of conditions - AMEND II
NTPs_II = 2*numel(subDirs_II); % Total N of TPs - AMEND II
TPsOrder_II = zeros(NTPs_II,NCond_II); % Vector that will contain the indexes of trials/file for each TP - AMEND II

NCols=60*Param.Fs; % Duration (samples) of each window
NRows=50; % Number of windows per trial
NLayers=numel(FileNames_II)*numel(subDirs_II); % Number of trials

DGSW = zeros(NLayers,NRows,NCols); DGLW = DGSW; % Global Diameter Speaking/Listening Windows
DSW_NH = DGSW; DLW_NH = DGSW; % NH Diameter Speaking/Listening Windows
DSW_HI = DGSW; DLW_HI = DGSW; % HI Diameter Speaking/Listening Windows
DSW_N0 = DGSW; DLW_N0 = DGSW; % N0 Diameter Speaking/Listening Windows
DSW_N60 = DGSW; DLW_N60 = DGSW; % N60 Diameter Speaking/Listening Windows
DSW_N70 = DGSW; DLW_N70 = DGSW; % N70 Diameter Speaking/Listening Windows
DSW_NH_N0 = DGSW; DLW_NH_N0 = DGSW; % NH N0 Diameter Speaking/Listening Windows
DSW_NH_N60 = DGSW; DLW_NH_N60 = DGSW; % NH N60 Diameter Speaking/Listening Windows
DSW_NH_N70 = DGSW; DLW_NH_N70 = DGSW; % NH N70 Diameter Speaking/Listening Windows
DSW_HI_N0 = DGSW; DLW_HI_N0 = DGSW; % HI N0 Diameter Speaking/Listening Windows
DSW_HI_N60 = DGSW; DLW_HI_N60 = DGSW; % HI N60 Diameter Speaking/Listening Windows
DSW_HI_N70 = DGSW; DLW_HI_N70 = DGSW; % HI N70 Diameter Speaking/Listening Windows
DSW_NH_UN = DGSW; DLW_NH_UN = DGSW; % NH UN Diameter Speaking/Listening Windows
DSW_NH_AA = DGSW; DLW_NH_AA = DGSW; % NH AA Diameter Speaking/Listening Windows
DSW_NH_AB = DGSW; DLW_NH_AB = DGSW; % NH AB Diameter Speaking/Listening Windows
DSW_HI_UN = DGSW; DLW_HI_UN = DGSW; % HI UN Diameter Speaking/Listening Windows
DSW_HI_AA = DGSW; DLW_HI_AA = DGSW; % HI AA Diameter Speaking/Listening Windows
DSW_HI_AB = DGSW; DLW_HI_AB = DGSW; % HI AB Diameter Speaking/Listening Windows
DSW_NH_UN_N0 = DGSW; DLW_NH_UN_N0 = DGSW; % NH UN N0 Diameter Speaking/Listening Windows
DSW_NH_UN_N60 = DGSW; DLW_NH_UN_N60 = DGSW; % NH UN N60 Diameter Speaking/Listening Windows
DSW_NH_UN_N70 = DGSW; DLW_NH_UN_N70 = DGSW; % NH UN N70 Diameter Speaking/Listening Windows
DSW_NH_AA_N0 = DGSW; DLW_NH_AA_N0 = DGSW; % NH AA N0 Diameter Speaking/Listening Windows
DSW_NH_AA_N60 = DGSW; DLW_NH_AA_N60 = DGSW; % NH AA N60 Diameter Speaking/Listening Windows
DSW_NH_AA_N70 = DGSW; DLW_NH_AA_N70 = DGSW; % NH AA N70 Diameter Speaking/Listening Windows
DSW_NH_AB_N0 = DGSW; DLW_NH_AB_N0 = DGSW; % NH AB N0 Diameter Speaking/Listening Windows
DSW_NH_AB_N60 = DGSW; DLW_NH_AB_N60 = DGSW; % NH AB N60 Diameter Speaking/Listening Windows
DSW_NH_AB_N70 = DGSW; DLW_NH_AB_N70 = DGSW; % NH AB N70 Diameter Speaking/Listening Windows
DSW_HI_UN_N0 = DGSW; DLW_HI_UN_N0 = DGSW; % HI UN N0 Diameter Speaking/Listening Windows
DSW_HI_UN_N60 = DGSW; DLW_HI_UN_N60 = DGSW; % HI UN N60 Diameter Speaking/Listening Windows
DSW_HI_UN_N70 = DGSW; DLW_HI_UN_N70 = DGSW; % HI UN N70 Diameter Speaking/Listening Windows
DSW_HI_AA_N0 = DGSW; DLW_HI_AA_N0 = DGSW; % HI AA N0 Diameter Speaking/Listening Windows
DSW_HI_AA_N60 = DGSW; DLW_HI_AA_N60 = DGSW; % HI AA N60 Diameter Speaking/Listening Windows
DSW_HI_AA_N70 = DGSW; DLW_HI_AA_N70 = DGSW; % HI AA N70 Diameter Speaking/Listening Windows
DSW_HI_AB_N0 = DGSW; DLW_HI_AB_N0 = DGSW; % HI AB N0 Diameter Speaking/Listening Windows
DSW_HI_AB_N60 = DGSW; DLW_HI_AB_N60 = DGSW; % HI AB N60 Diameter Speaking/Listening Windows
DSW_HI_AB_N70 = DGSW; DLW_HI_AB_N70 = DGSW; % HI AB N70 Diameter Speaking/Listening Windows

DGSW_B = DGSW; DGLW_B = DGSW; % Global Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_NH_B = DGSW; DLW_NH_B = DGSW; % NH Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_HI_B = DGSW; DLW_HI_B = DGSW; % HI Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_N0_B = DGSW; DLW_N0_B = DGSW; % N0 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_N60_B = DGSW; DLW_N60_B = DGSW; % N60 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_N70_B = DGSW; DLW_N70_B = DGSW; % N70 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_NH_N0_B = DGSW; DLW_NH_N0_B = DGSW; % NH N0 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_NH_N60_B = DGSW; DLW_NH_N60_B = DGSW; % NH N60 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_NH_N70_B = DGSW; DLW_NH_N70_B = DGSW; % NH N70 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_HI_N0_B = DGSW; DLW_HI_N0_B = DGSW; % HI N0 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_HI_N60_B = DGSW; DLW_HI_N60_B = DGSW; % HI N60 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_HI_N70_B = DGSW; DLW_HI_N70_B = DGSW; % HI N70 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_NH_UN_B = DGSW; DLW_NH_UN_B = DGSW; % NH UN Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_NH_AA_B = DGSW; DLW_NH_AA_B = DGSW; % NH AA Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_NH_AB_B = DGSW; DLW_NH_AB_B = DGSW; % NH AB Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_HI_UN_B = DGSW; DLW_HI_UN_B = DGSW; % HI UN Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_HI_AA_B = DGSW; DLW_HI_AA_B = DGSW; % HI AA Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_HI_AB_B = DGSW; DLW_HI_AB_B = DGSW; % HI AB Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_NH_UN_N0_B = DGSW; DLW_NH_UN_N0_B = DGSW; % NH UN N0 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_NH_UN_N60_B = DGSW; DLW_NH_UN_N60_B = DGSW; % NH UN N60 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_NH_UN_N70_B = DGSW; DLW_NH_UN_N70_B = DGSW; % NH UN N70 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_NH_AA_N0_B = DGSW; DLW_NH_AA_N0_B = DGSW; % NH AA N0 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_NH_AA_N60_B = DGSW; DLW_NH_AA_N60_B = DGSW; % NH AA N60 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_NH_AA_N70_B = DGSW; DLW_NH_AA_N70_B = DGSW; % NH AA N70 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_NH_AB_N0_B = DGSW; DLW_NH_AB_N0_B = DGSW; % NH AB N0 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_NH_AB_N60_B = DGSW; DLW_NH_AB_N60_B = DGSW; % NH AB N60 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_NH_AB_N70_B = DGSW; DLW_NH_AB_N70_B = DGSW; % NH AB N70 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_HI_UN_N0_B = DGSW; DLW_HI_UN_N0_B = DGSW; % HI UN N0 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_HI_UN_N60_B = DGSW; DLW_HI_UN_N60_B = DGSW; % HI UN N60 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_HI_UN_N70_B = DGSW; DLW_HI_UN_N70_B = DGSW; % HI UN N70 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_HI_AA_N0_B = DGSW; DLW_HI_AA_N0_B = DGSW; % HI AA N0 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_HI_AA_N60_B = DGSW; DLW_HI_AA_N60_B = DGSW; % HI AA N60 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_HI_AA_N70_B = DGSW; DLW_HI_AA_N70_B = DGSW; % HI AA N70 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_HI_AB_N0_B = DGSW; DLW_HI_AB_N0_B = DGSW; % HI AB N0 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_HI_AB_N60_B = DGSW; DLW_HI_AB_N60_B = DGSW; % HI AB N60 Diameter Speaking/Listening Windows Adaptive-Baseline corrected
DSW_HI_AB_N70_B = DGSW; DLW_HI_AB_N70_B = DGSW; % HI AB N70 Diameter Speaking/Listening Windows Adaptive-Baseline corrected

FGSW = DGSW; FGLW = DGSW; % Global Fixation Speaking/Listening Windows
FSW_NH = DGSW; FLW_NH = DGSW; % NH Fixation Speaking/Listening Windows
FSW_HI = DGSW; FLW_HI = DGSW; % HI Fixation Speaking/Listening Windows
FSW_N0 = DGSW; FLW_N0 = DGSW; % N0 Fixation Speaking/Listening Windows
FSW_N60 = DGSW; FLW_N60 = DGSW; % N60 Fixation Speaking/Listening Windows
FSW_N70 = DGSW; FLW_N70 = DGSW; % N70 Fixation Speaking/Listening Windows
FSW_NH_N0 = DGSW; FLW_NH_N0 = DGSW; % NH N0 Fixation Speaking/Listening Windows
FSW_NH_N60 = DGSW; FLW_NH_N60 = DGSW; % NH N60 Fixation Speaking/Listening Windows
FSW_NH_N70 = DGSW; FLW_NH_N70 = DGSW; % NH N70 Fixation Speaking/Listening Windows
FSW_HI_N0 = DGSW; FLW_HI_N0 = DGSW; % HI N0 Fixation Speaking/Listening Windows
FSW_HI_N60 = DGSW; FLW_HI_N60 = DGSW; % HI N60 Fixation Speaking/Listening Windows
FSW_HI_N70 = DGSW; FLW_HI_N70 = DGSW; % HI N70 Fixation Speaking/Listening Windows
FSW_NH_UN = DGSW; FLW_NH_UN = DGSW; % NH UN Fixation Speaking/Listening Windows
FSW_NH_AA = DGSW; FLW_NH_AA = DGSW; % NH AA Fixation Speaking/Listening Windows
FSW_NH_AB = DGSW; FLW_NH_AB = DGSW; % NH AB Fixation Speaking/Listening Windows
FSW_HI_UN = DGSW; FLW_HI_UN = DGSW; % HI UN Fixation Speaking/Listening Windows
FSW_HI_AA = DGSW; FLW_HI_AA = DGSW; % HI AA Fixation Speaking/Listening Windows
FSW_HI_AB = DGSW; FLW_HI_AB = DGSW; % HI AB Fixation Speaking/Listening Windows
FSW_NH_UN_N0 = DGSW; FLW_NH_UN_N0 = DGSW; % NH UN N0 Fixation Speaking/Listening Windows
FSW_NH_UN_N60 = DGSW; FLW_NH_UN_N60 = DGSW; % NH UN N60 Fixation Speaking/Listening Windows
FSW_NH_UN_N70 = DGSW; FLW_NH_UN_N70 = DGSW; % NH UN N70 Fixation Speaking/Listening Windows
FSW_NH_AA_N0 = DGSW; FLW_NH_AA_N0 = DGSW; % NH AA N0 Fixation Speaking/Listening Windows
FSW_NH_AA_N60 = DGSW; FLW_NH_AA_N60 = DGSW; % NH AA N60 Fixation Speaking/Listening Windows
FSW_NH_AA_N70 = DGSW; FLW_NH_AA_N70 = DGSW; % NH AA N70 Fixation Speaking/Listening Windows
FSW_NH_AB_N0 = DGSW; FLW_NH_AB_N0 = DGSW; % NH AB N0 Fixation Speaking/Listening Windows
FSW_NH_AB_N60 = DGSW; FLW_NH_AB_N60 = DGSW; % NH AB N60 Fixation Speaking/Listening Windows
FSW_NH_AB_N70 = DGSW; FLW_NH_AB_N70 = DGSW; % NH AB N70 Fixation Speaking/Listening Windows
FSW_HI_UN_N0 = DGSW; FLW_HI_UN_N0 = DGSW; % HI UN N0 Fixation Speaking/Listening Windows
FSW_HI_UN_N60 = DGSW; FLW_HI_UN_N60 = DGSW; % HI UN N60 Fixation Speaking/Listening Windows
FSW_HI_UN_N70 = DGSW; FLW_HI_UN_N70 = DGSW; % HI UN N70 Fixation Speaking/Listening Windows
FSW_HI_AA_N0 = DGSW; FLW_HI_AA_N0 = DGSW; % HI AA N0 Fixation Speaking/Listening Windows
FSW_HI_AA_N60 = DGSW; FLW_HI_AA_N60 = DGSW; % HI AA N60 Fixation Speaking/Listening Windows
FSW_HI_AA_N70 = DGSW; FLW_HI_AA_N70 = DGSW; % HI AA N70 Fixation Speaking/Listening Windows
FSW_HI_AB_N0 = DGSW; FLW_HI_AB_N0 = DGSW; % HI AB N0 Fixation Speaking/Listening Windows
FSW_HI_AB_N60 = DGSW; FLW_HI_AB_N60 = DGSW; % HI AB N60 Fixation Speaking/Listening Windows
FSW_HI_AB_N70 = DGSW; FLW_HI_AB_N70 = DGSW; % HI AB N70 Fixation Speaking/Listening Windows

GSDur = zeros(NLayers,NRows); % Global Speaking Windows Duration
GLDur = zeros(NLayers,NRows); % Global Listening Windows Duration
SDur_N0 = zeros(NLayers,NRows); % N0 Speaking Windows Duration
LDur_N0 = zeros(NLayers,NRows); % N0 Listening Windows Duration
SDur_N60 = zeros(NLayers,NRows); % N60 Speaking Windows Duration
LDur_N60 = zeros(NLayers,NRows); % N60 Listening Windows Duration
SDur_N70 = zeros(NLayers,NRows); % N70 Speaking Windows Duration
LDur_N70 = zeros(NLayers,NRows); % N70 Listening Windows Duration

% AMEND II Loop
for q=1:numel(subDirs_II)
    PairIn_II = q;
    for ChosenFolder = {'\NH\','\HI\'}
        PairFolder_II=[pwd,'\data\AMEND_II\',cell2mat(subDirs_II(q)),cell2mat(ChosenFolder)]; % Folder naming changed
        PairFiles_II=dir(PairFolder_II); % Folder naming changed
        PairUtt_II=LoadUtt_II.Utterances(PairIn_II,:);
        
        if isempty(PairUtt_II{1})
            disp(['Warning: No Utterance found for folder "',cell2mat(ChosenFolder),cell2mat(subDirs_II(q)),'"'])
            continue
        end
%         PairDelay_II=LoadDelays_II.TobAudDelay(PairIn_II,:);
    
        for i=1:numel(FileNames_II)
            try
                alldata = load([PairFiles_II(1).folder, '\', cell2mat(FileNames_II(i))]);
            catch ME
                disp(['Warning: File ', PairFiles_II(1).folder, '\', cell2mat(FileNames_II(i)), ' not found (no Gaze data).']);
                continue
            end
            alldata_mat = cell2mat(alldata.data);

            % Extract Timestamps [s]
            timestamps = vertcat(alldata_mat.timeStamp);
            
            % DIAMETER
            % Replace blanks '[]' for 'NaN' in fields diameterLeft and diameterRight
            [alldata_mat(cellfun(@isempty,{alldata_mat.diameterLeft})).diameterLeft] = deal(NaN);
            [alldata_mat(cellfun(@isempty,{alldata_mat.diameterRight})).diameterRight] = deal(NaN);

            LDiamRaw = [alldata_mat.diameterLeft];
            RDiamRaw = [alldata_mat.diameterRight];

            % Preprocessing - Setting outliers as NaNs (remove artifacts)
            % New artifact-removal method
            standardRawSettings = rawDataFilter();
            [LvalOut,LspeedFiltData,LdevFiltData] = rawDataFilter(linspace(0,length(LDiamRaw)./Param.Fs,length(LDiamRaw))',LDiamRaw',standardRawSettings);
            [RvalOut,RspeedFiltData,RdevFiltData] = rawDataFilter(linspace(0,length(RDiamRaw)./Param.Fs,length(RDiamRaw))',RDiamRaw',standardRawSettings);
            
            LDiamRaw(~LvalOut)=NaN;
            RDiamRaw(~RvalOut)=NaN;
            
            % Processing - Interpolating NaNs
            [LDiam,LMetadata] = preprocpupil(LDiamRaw,Param);
            [RDiam,RMetadata] = preprocpupil(RDiamRaw,Param);

            % Low-Pass Filtering
            LDiamConv = conv(LDiam,LPWindow,'same'); 
            RDiamConv = conv(RDiam,LPWindow,'same');

            % Remove start/end artifacts (peak/dip) originated from Low-Pass Filtering
            LDiamConv(1:round(length(LPWindow)/2-1)) = mean(LDiam(1:round(length(LPWindow)/2-1)));
            LDiamConv(end-round(length(LPWindow)/2-1):end) = mean(LDiam(end-round(length(LPWindow)/2-1):end));
            RDiamConv(1:round(length(LPWindow)/2-1)) = mean(RDiam(1:round(length(LPWindow)/2-1)));
            RDiamConv(end-round(length(LPWindow)/2-1):end) = mean(RDiam(end-round(length(LPWindow)/2-1):end));

            % Comment below if No LP filtering wanted
            LDiam = LDiamConv;
            RDiam = RDiamConv;

            % Decide 'better' eye results
            [Min, idx_decision] = min([sum(LMetadata.Isnan) sum(RMetadata.Isnan)]);
            if idx_decision == 1
                Diameter = LDiam;
                DiameterRaw = LDiamRaw';
                eyeChosen = 'Left';
                DiamNaN = sum(LMetadata.Isnan);
            elseif idx_decision == 2
                Diameter = RDiam;
                DiameterRaw = RDiamRaw';
                eyeChosen = 'Right';
                DiamNaN = sum(RMetadata.Isnan);
            end

            if DiamNaN/length(Diameter) >= RejectRatio
                disp(['Warning: File ', PairFiles_II(1).folder, '\', cell2mat(FileNames_II(i)), ' was rejected because pupil size contains too many NaNs (',sprintf('%0.2f',100*DiamNaN/length(Diameter)),'%).'])
                continue
            end
            
            
            % GAZE 3D
            % In case gaze3d is transposed from origin [x;y;z] -> transpose to [x,y,z]
            if size([alldata_mat(:,cellfun(@(xd) ~any(isnan(xd)),{alldata_mat.gaze3d})).gaze3d],1) == 3
                for k=1:size(alldata_mat,2)
                   alldata_mat(k).gaze3d=alldata_mat(k).gaze3d'; 
                end
            end

            % Replace blanks '[]' for 'NaN' in gaze3d
            [alldata_mat(cellfun(@isempty,{alldata_mat.gaze3d})).gaze3d] = deal(NaN);

            % Replace 'NaN' to '[NaN,NaN,NaN]' in gaze3d
            [alldata_mat(:,cellfun(@(xd) any(isnan(xd)),{alldata_mat.gaze3d})).gaze3d] = deal([NaN,NaN,NaN]);
            
            gaze3d = vertcat(alldata_mat.gaze3d);
        
            GazeX = gaze3d(:,1);
            GazeY = gaze3d(:,2);
            GazeZ = gaze3d(:,3);

            % NaN's appear at the same time at [X,Y,Z], we only look at X.
            XNan = find(isnan(GazeX));

            % Reject a file if gaze3d has too many NaNs
            if length(XNan)/length(gaze3d) >= RejectRatio
                disp(['Warning: File ', PairFiles_II(1).folder, '\', cell2mat(FileNames_II(i)), ' was rejected because gaze3d contains too many NaNs (',sprintf('%0.2f',100*length(XNan)/length(gaze3d)),'%).'])
                continue
            end
            
            % Preprocessing - PUPILS pipeline
            % 1) Calculate angular velocity
            [az, el] = calcAngular(GazeX,GazeY,GazeZ);
            gazeAz_velocity = [0;diff(az)/(1/Param.Fs)];
            gazeEl_velocity = [0;diff(el)/(1/Param.Fs)];
            gaze_vel=sqrt(gazeAz_velocity.^2.*cosd(el).^2+gazeEl_velocity.^2);

            % 2) Select preprocessing options
            options = struct;
            options.fs = Param.Fs;            % sampling frequency (Hz)
            options.blink_rule = 'std';       % Rule for blink detection 'std' / 'vel'
            options.pre_blink_t   = 100;      % region to interpolate before blink (ms)
            options.post_blink_t  = 200;      % region to interpolate after blink (ms)
            options.xy_units = 'mm';          % xy coordinate units 'px' / 'mm' / 'cm'
            options.vel_threshold =  30;      % velocity threshold for saccade detection
            options.min_sacc_duration = 10;   % minimum saccade duration (ms)
            options.interpolate_saccades = 0; % Specify whether saccadic distortions should be interpolated 1-yes 0-noB
            options.pre_sacc_t   = 50;        % Region to interpolate before saccade (ms)
            options.post_sacc_t  = 100;       % Region to interpolate after saccade (ms)
            options.low_pass_fc   = 10;       % Low-pass filter cut-off frequency (Hz)

            % 3) Using PUPILS toolbox 
            %   -IN [samples]: [timestamps, X, Y, Z, Diameter, gaze_vel]
            %   -OUT[samples]: [timestamps, X, Y, Z, Diameter, gaze_vel, blinks, saccades, fixations, interpolated, denoised]
            data = [timestamps*Param.Fs GazeX GazeY DiameterRaw gaze_vel];
            [proc_data, proc_info] = processPupilData(data, options);

            % 4) Fixation duration - from fixation samples (from PUPILS)
            GazeX_fix = GazeX; GazeY_fix = GazeY; GazeZ_fix = GazeZ;
            GazeX_fix(proc_data(:,8)==0)=NaN; GazeY_fix(proc_data(:,8)==0)=NaN; GazeZ_fix(proc_data(:,8)==0)=NaN;
            Fixation = fix_duration([GazeX_fix,GazeY_fix,GazeZ_fix],MaxVisualAngle,Param.Fs);        

            % 5) Filtering - Fixation duration with hamming window (F = 3)
            Fixation = ndnanfilter(Fixation,'hamming',FilterWidth);
            
            % SPEAKING/LISTENING WINDOWS
            % W18 - Add nan padding (for valid baseline)
            NaNPadS = 2*TimeStartW*Param.Fs;
            Diameter = [nan*ones(NaNPadS,1);Diameter];
            Fixation = [nan*ones(NaNPadS,1);Fixation];
            
            % Retrieve Utterances
            if contains(ChosenFolder,'NH')
                SpeakKey = 'utteranceCH1';
                ListenKey = 'utteranceCH2';
            elseif contains(ChosenFolder,'HI')
                SpeakKey = 'utteranceCH2';
                ListenKey = 'utteranceCH1';
            end

            if contains(cell2mat(FileNames_II(i)),'UN')
                SpeAid = 0;
            elseif contains(cell2mat(FileNames_II(i)),'AA')
                SpeAid = 3;
            elseif contains(cell2mat(FileNames_II(i)),'AB')
                SpeAid = 6;
            end

            if contains(cell2mat(FileNames_II(i)),'N0')
                SpeCond = SpeAid + 1;
            elseif contains(cell2mat(FileNames_II(i)),'N60')
                SpeCond = SpeAid + 2;
            elseif contains(cell2mat(FileNames_II(i)),'N70')
                SpeCond = SpeAid + 3;
            end
            
            SpeakRaw = PairUtt_II{1,SpeCond}.(SpeakKey);
            ListenRaw = PairUtt_II{1,SpeCond}.(ListenKey);
            binResUtt = PairUtt_II{1,SpeCond}.binRes;

            if isempty(SpeakRaw) && isempty(ListenRaw)
                disp(['Warning: File ',PairFiles_II(1).folder, '\', cell2mat(FileNames_II(i)),' was rejected for not having associated Utterance windows.']);
                continue
            end
            
            % DIFFERENT PROCESSING FOR AMEND II (NanPadS)
            % Downsample (rounding) Utt from 250 Hz (1/binRes) to 50 Hz
            SpeakRaw(:,2:3)=round((SpeakRaw(:,2:3)*binResUtt)*Param.Fs)+NaNPadS;
            ListenRaw(:,2:3)=round((ListenRaw(:,2:3)*binResUtt)*Param.Fs)+NaNPadS;

            % Merge windows if duration between windows <= TimeInitialMerge (300 ms)
            SpeakMI = merge_windows(SpeakRaw, Param.Fs, TimeInitialMerge);
            ListenMI = merge_windows(ListenRaw, Param.Fs, TimeInitialMerge);

            % Discard windows if duration is < TimeMinWin (500 ms)
            SpeakDI = SpeakMI(SpeakMI(:,1)>TimeMinWin,:);
            ListenDI = ListenMI(ListenMI(:,1)>TimeMinWin,:);

            % Merge again if duration between windows <= TimeMerge (2 s)
            SpeakM = merge_windows(SpeakDI, Param.Fs, TimeMerge);
            ListenM = merge_windows(ListenDI, Param.Fs, TimeMerge);

            % Discard windows if duration is < 2*TimeMinWin (1 s)
            SpeakD = SpeakM(SpeakM(:,1)>2*TimeMinWin,:);
            ListenD = ListenM(ListenM(:,1)>2*TimeMinWin,:);

            % Added from W6.m -> Cut/Split overlaps
            [Speak,Listen] = overlap_windows(SpeakD,ListenD,Param.Fs);

            % Time-locked indexes (based on Start or End of events)
            WSpeakIdx=[Speak(:,2)-TimeStartW*Param.Fs,Speak(:,2),Speak(:,3)+TimeEndW*Param.Fs];
            WListenIdx=[Listen(:,2)-TimeStartW*Param.Fs,Listen(:,2),Listen(:,3)+TimeEndW*Param.Fs];
            
        %         EWSpeakIdx=[Speak(:,3)-TimeStartW*Param.Fs,Speak(:,3),Speak(:,3)+TimeEndW*Param.Fs];
        %         EWListenIdx=[Listen(:,3)-TimeStartW*Param.Fs,Listen(:,3),Listen(:,3)+TimeEndW*Param.Fs];

        %         t_Diam = linspace(0,length(Diameter)./Param.Fs,length(Diameter));

            % FIXATION+DIAMETER GROUPS
            % Storing Speaking/Listening by conditions
            if contains(cell2mat(FileNames_II(i)),'N0')
                for j=1:size(Speak,1)
                    DSW_N0(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    DSW_N0_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_N0_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    FSW_N0(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    SDur_N0(i,j) = Speak(j,1);
                end
                for j=1:size(Listen,1)
                    DLW_N0(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                    DLW_N0_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_N0_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    FLW_N0(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                    LDur_N0(i,j) = Listen(j,1);
                end
            elseif contains(cell2mat(FileNames_II(i)),'N60')
                for j=1:size(Speak,1)
                    DSW_N60(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    DSW_N60_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_N60_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    FSW_N60(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    SDur_N60(i,j) = Speak(j,1);
                end
                for j=1:size(Listen,1)
                    DLW_N60(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                    DLW_N60_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_N60_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    FLW_N60(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                    LDur_N60(i,j) = Listen(j,1);
                end
            elseif contains(cell2mat(FileNames_II(i)),'N70')
                for j=1:size(Speak,1)
                    DSW_N70(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    DSW_N70_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_N70_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    FSW_N70(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    SDur_N70(i,j) = Speak(j,1);
                end
                for j=1:size(Listen,1)
                    DLW_N70(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                    DLW_N70_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_N70_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    FLW_N70(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                    LDur_N70(i,j) = Listen(j,1);
                end
            end
            % Either NH or HI
            if contains(ChosenFolder,'NH')
                for j=1:size(Speak,1)
                    DSW_NH(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_NH)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    DSW_NH_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_NH_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                	FSW_NH(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_NH)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                end
                for j=1:size(Listen,1)
                    DLW_NH(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_NH)-length(WListenIdx(j,1):Listen(j,3)))'];
                    DLW_NH_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_NH_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                	FLW_NH(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_NH)-length(WListenIdx(j,1):Listen(j,3)))'];
                end
                if contains(cell2mat(FileNames_II(i)),'N0')
                    for j=1:size(Speak,1)
                        DSW_NH_N0(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_NH_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        DSW_NH_N0_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_NH_N0_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    	FSW_NH_N0(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_NH_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    end
                    for j=1:size(Listen,1)
                        DLW_NH_N0(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_NH_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                        DLW_NH_N0_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_NH_N0_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    	FLW_NH_N0(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_NH_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                    end
                elseif contains(cell2mat(FileNames_II(i)),'N60')
                    for j=1:size(Speak,1)
                        DSW_NH_N60(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_NH_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        DSW_NH_N60_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_NH_N60_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    	FSW_NH_N60(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_NH_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    end
                    for j=1:size(Listen,1)
                        DLW_NH_N60(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_NH_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                        DLW_NH_N60_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_NH_N60_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    	FLW_NH_N60(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_NH_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                    end
                elseif contains(cell2mat(FileNames_II(i)),'N70')
                    for j=1:size(Speak,1)
                        DSW_NH_N70(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_NH_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        DSW_NH_N70_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_NH_N70_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    	FSW_NH_N70(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_NH_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    end
                    for j=1:size(Listen,1)
                        DLW_NH_N70(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_NH_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                        DLW_NH_N70_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_NH_N70_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    	FLW_NH_N70(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_NH_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                    end
                end
                if contains(cell2mat(FileNames_II(i)),'UN')
                    for j=1:size(Speak,1)
                        DSW_NH_UN(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_NH_UN)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        DSW_NH_UN_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_NH_UN_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    	FSW_NH_UN(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_NH_UN)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    end
                    for j=1:size(Listen,1)
                        DLW_NH_UN(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_NH_UN)-length(WListenIdx(j,1):Listen(j,3)))'];
                        DLW_NH_UN_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_NH_UN_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    	FLW_NH_UN(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_NH_UN)-length(WListenIdx(j,1):Listen(j,3)))'];
                    end
                    if contains(cell2mat(FileNames_II(i)),'N0')
                        for j=1:size(Speak,1)
                            DSW_NH_UN_N0(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_NH_UN_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            DSW_NH_UN_N0_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_NH_UN_N0_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        	FSW_NH_UN_N0(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_NH_UN_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            DLW_NH_UN_N0(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_NH_UN_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                            DLW_NH_UN_N0_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_NH_UN_N0_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                            FLW_NH_UN_N0(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_NH_UN_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    elseif contains(cell2mat(FileNames_II(i)),'N60')
                        for j=1:size(Speak,1)
                            DSW_NH_UN_N60(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_NH_UN_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            DSW_NH_UN_N60_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_NH_UN_N60_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        	FSW_NH_UN_N60(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_NH_UN_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            DLW_NH_UN_N60(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_NH_UN_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                            DLW_NH_UN_N60_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_NH_UN_N60_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        	FLW_NH_UN_N60(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_NH_UN_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    elseif contains(cell2mat(FileNames_II(i)),'N70')
                        for j=1:size(Speak,1)
                            DSW_NH_UN_N70(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_NH_UN_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            DSW_NH_UN_N70_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_NH_UN_N70_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        	FSW_NH_UN_N70(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_NH_UN_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            DLW_NH_UN_N70(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_NH_UN_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                            DLW_NH_UN_N70_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_NH_UN_N70_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        	FLW_NH_UN_N70(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_NH_UN_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    end
                elseif contains(cell2mat(FileNames_II(i)),'AA')
                    for j=1:size(Speak,1)
                        DSW_NH_AA(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_NH_AA)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        DSW_NH_AA_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_NH_AA_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    	FSW_NH_AA(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_NH_AA)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    end
                    for j=1:size(Listen,1)
                        DLW_NH_AA(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_NH_AA)-length(WListenIdx(j,1):Listen(j,3)))'];
                        DLW_NH_AA_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_NH_AA_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    	FLW_NH_AA(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_NH_AA)-length(WListenIdx(j,1):Listen(j,3)))'];
                    end
                    if contains(cell2mat(FileNames_II(i)),'N0')
                        for j=1:size(Speak,1)
                            DSW_NH_AA_N0(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_NH_AA_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            DSW_NH_AA_N0_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_NH_AA_N0_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        	FSW_NH_AA_N0(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_NH_AA_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            DLW_NH_AA_N0(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_NH_AA_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                            DLW_NH_AA_N0_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_NH_AA_N0_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        	FLW_NH_AA_N0(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_NH_AA_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    elseif contains(cell2mat(FileNames_II(i)),'N60')
                        for j=1:size(Speak,1)
                            DSW_NH_AA_N60(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_NH_AA_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            DSW_NH_AA_N60_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_NH_AA_N60_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        	FSW_NH_AA_N60(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_NH_AA_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            DLW_NH_AA_N60(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_NH_AA_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                            DLW_NH_AA_N60_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_NH_AA_N60_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        	FLW_NH_AA_N60(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_NH_AA_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    elseif contains(cell2mat(FileNames_II(i)),'N70')
                        for j=1:size(Speak,1)
                            DSW_NH_AA_N70(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_NH_AA_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            DSW_NH_AA_N70_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_NH_AA_N70_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        	FSW_NH_AA_N70(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_NH_AA_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            DLW_NH_AA_N70(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_NH_AA_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                            DLW_NH_AA_N70_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_NH_AA_N70_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        	FLW_NH_AA_N70(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_NH_AA_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    end
                elseif contains(cell2mat(FileNames_II(i)),'AB')
                    for j=1:size(Speak,1)
                        DSW_NH_AB(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_NH_AB)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        DSW_NH_AB_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_NH_AB_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    	FSW_NH_AB(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_NH_AB)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    end
                    for j=1:size(Listen,1)
                        DLW_NH_AB(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_NH_AB)-length(WListenIdx(j,1):Listen(j,3)))'];
                        DLW_NH_AB_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_NH_AB_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    	FLW_NH_AB(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_NH_AB)-length(WListenIdx(j,1):Listen(j,3)))'];
                    end
                    if contains(cell2mat(FileNames_II(i)),'N0')
                        for j=1:size(Speak,1)
                            DSW_NH_AB_N0(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_NH_AB_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            DSW_NH_AB_N0_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_NH_AB_N0_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        	FSW_NH_AB_N0(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_NH_AB_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            DLW_NH_AB_N0(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_NH_AB_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                            DLW_NH_AB_N0_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_NH_AB_N0_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        	FLW_NH_AB_N0(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_NH_AB_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    elseif contains(cell2mat(FileNames_II(i)),'N60')
                        for j=1:size(Speak,1)
                            DSW_NH_AB_N60(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_NH_AB_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            DSW_NH_AB_N60_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_NH_AB_N60_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        	FSW_NH_AB_N60(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_NH_AB_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            DLW_NH_AB_N60(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_NH_AB_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                            DLW_NH_AB_N60_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_NH_AB_N60_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        	FLW_NH_AB_N60(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_NH_AB_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    elseif contains(cell2mat(FileNames_II(i)),'N70')
                        for j=1:size(Speak,1)
                            DSW_NH_AB_N70(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_NH_AB_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            DSW_NH_AB_N70_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_NH_AB_N70_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        	FSW_NH_AB_N70(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_NH_AB_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            DLW_NH_AB_N70(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_NH_AB_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                            DLW_NH_AB_N70_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_NH_AB_N70_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        	FLW_NH_AB_N70(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_NH_AB_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    end
                end
            elseif contains(ChosenFolder,'HI')
                for j=1:size(Speak,1)
                    DSW_HI(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_HI)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    DSW_HI_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_HI_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                	FSW_HI(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_HI)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                end
                for j=1:size(Listen,1)
                    DLW_HI(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_HI)-length(WListenIdx(j,1):Listen(j,3)))'];
                    DLW_HI_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_HI_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                	FLW_HI(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_HI)-length(WListenIdx(j,1):Listen(j,3)))'];
                end
                if contains(cell2mat(FileNames_II(i)),'N0')
                    for j=1:size(Speak,1)
                        DSW_HI_N0(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_HI_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        DSW_HI_N0_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_HI_N0_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    	FSW_HI_N0(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_HI_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    end
                    for j=1:size(Listen,1)
                        DLW_HI_N0(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_HI_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                        DLW_HI_N0_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_HI_N0_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    	FLW_HI_N0(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_HI_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                    end
                elseif contains(cell2mat(FileNames_II(i)),'N60')
                    for j=1:size(Speak,1)
                        DSW_HI_N60(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_HI_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        DSW_HI_N60_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_HI_N60_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    	FSW_HI_N60(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_HI_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    end
                    for j=1:size(Listen,1)
                        DLW_HI_N60(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_HI_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                        DLW_HI_N60_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_HI_N60_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    	FLW_HI_N60(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_HI_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                    end
                elseif contains(cell2mat(FileNames_II(i)),'N70')
                    for j=1:size(Speak,1)
                        DSW_HI_N70(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_HI_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        DSW_HI_N70_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_HI_N70_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    	FSW_HI_N70(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_HI_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    end
                    for j=1:size(Listen,1)
                        DLW_HI_N70(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_HI_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                        DLW_HI_N70_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_HI_N70_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    	FLW_HI_N70(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_HI_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                    end
                end
                if contains(cell2mat(FileNames_II(i)),'UN')
                    for j=1:size(Speak,1)
                        DSW_HI_UN(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_HI_UN)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        DSW_HI_UN_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_HI_UN_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    	FSW_HI_UN(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_HI_UN)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    end
                    for j=1:size(Listen,1)
                        DLW_HI_UN(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_HI_UN)-length(WListenIdx(j,1):Listen(j,3)))'];
                        DLW_HI_UN_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_HI_UN_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    	FLW_HI_UN(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_HI_UN)-length(WListenIdx(j,1):Listen(j,3)))'];
                    end
                    if contains(cell2mat(FileNames_II(i)),'N0')
                        for j=1:size(Speak,1)
                            DSW_HI_UN_N0(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_HI_UN_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            DSW_HI_UN_N0_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_HI_UN_N0_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        	FSW_HI_UN_N0(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_HI_UN_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            DLW_HI_UN_N0(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_HI_UN_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                            DLW_HI_UN_N0_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_HI_UN_N0_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                            FLW_HI_UN_N0(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_HI_UN_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    elseif contains(cell2mat(FileNames_II(i)),'N60')
                        for j=1:size(Speak,1)
                            DSW_HI_UN_N60(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_HI_UN_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            DSW_HI_UN_N60_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_HI_UN_N60_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        	FSW_HI_UN_N60(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_HI_UN_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            DLW_HI_UN_N60(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_HI_UN_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                            DLW_HI_UN_N60_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_HI_UN_N60_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        	FLW_HI_UN_N60(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_HI_UN_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    elseif contains(cell2mat(FileNames_II(i)),'N70')
                        for j=1:size(Speak,1)
                            DSW_HI_UN_N70(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_HI_UN_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            DSW_HI_UN_N70_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_HI_UN_N70_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        	FSW_HI_UN_N70(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_HI_UN_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            DLW_HI_UN_N70(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_HI_UN_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                            DLW_HI_UN_N70_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_HI_UN_N70_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        	FLW_HI_UN_N70(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_HI_UN_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    end
                elseif contains(cell2mat(FileNames_II(i)),'AA')
                    for j=1:size(Speak,1)
                        DSW_HI_AA(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_HI_AA)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        DSW_HI_AA_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_HI_AA_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    	FSW_HI_AA(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_HI_AA)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    end
                    for j=1:size(Listen,1)
                        DLW_HI_AA(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_HI_AA)-length(WListenIdx(j,1):Listen(j,3)))'];
                        DLW_HI_AA_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_HI_AA_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    	FLW_HI_AA(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_HI_AA)-length(WListenIdx(j,1):Listen(j,3)))'];
                    end
                    if contains(cell2mat(FileNames_II(i)),'N0')
                        for j=1:size(Speak,1)
                            DSW_HI_AA_N0(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_HI_AA_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            DSW_HI_AA_N0_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_HI_AA_N0_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        	FSW_HI_AA_N0(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_HI_AA_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            DLW_HI_AA_N0(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_HI_AA_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                            DLW_HI_AA_N0_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_HI_AA_N0_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        	FLW_HI_AA_N0(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_HI_AA_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    elseif contains(cell2mat(FileNames_II(i)),'N60')
                        for j=1:size(Speak,1)
                            DSW_HI_AA_N60(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_HI_AA_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            DSW_HI_AA_N60_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_HI_AA_N60_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        	FSW_HI_AA_N60(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_HI_AA_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            DLW_HI_AA_N60(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_HI_AA_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                            DLW_HI_AA_N60_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_HI_AA_N60_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        	FLW_HI_AA_N60(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_HI_AA_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    elseif contains(cell2mat(FileNames_II(i)),'N70')
                        for j=1:size(Speak,1)
                            DSW_HI_AA_N70(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_HI_AA_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            DSW_HI_AA_N70_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_HI_AA_N70_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        	FSW_HI_AA_N70(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_HI_AA_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            DLW_HI_AA_N70(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_HI_AA_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                            DLW_HI_AA_N70_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_HI_AA_N70_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        	FLW_HI_AA_N70(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_HI_AA_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    end
                elseif contains(cell2mat(FileNames_II(i)),'AB')
                    for j=1:size(Speak,1)
                        DSW_HI_AB(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_HI_AB)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        DSW_HI_AB_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_HI_AB_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    	FSW_HI_AB(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_HI_AB)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                    end
                    for j=1:size(Listen,1)
                        DLW_HI_AB(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_HI_AB)-length(WListenIdx(j,1):Listen(j,3)))'];
                        DLW_HI_AB_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_HI_AB_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                    	FLW_HI_AB(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_HI_AB)-length(WListenIdx(j,1):Listen(j,3)))'];
                    end
                    if contains(cell2mat(FileNames_II(i)),'N0')
                        for j=1:size(Speak,1)
                            DSW_HI_AB_N0(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_HI_AB_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            DSW_HI_AB_N0_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_HI_AB_N0_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        	FSW_HI_AB_N0(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_HI_AB_N0)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            DLW_HI_AB_N0(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_HI_AB_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                            DLW_HI_AB_N0_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_HI_AB_N0_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        	FLW_HI_AB_N0(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_HI_AB_N0)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    elseif contains(cell2mat(FileNames_II(i)),'N60')
                        for j=1:size(Speak,1)
                            DSW_HI_AB_N60(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_HI_AB_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            DSW_HI_AB_N60_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_HI_AB_N60_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        	FSW_HI_AB_N60(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_HI_AB_N60)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            DLW_HI_AB_N60(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_HI_AB_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                            DLW_HI_AB_N60_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_HI_AB_N60_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        	FLW_HI_AB_N60(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_HI_AB_N60)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    elseif contains(cell2mat(FileNames_II(i)),'N70')
                        for j=1:size(Speak,1)
                            DSW_HI_AB_N70(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DSW_HI_AB_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                            DSW_HI_AB_N70_B(x,j,:)=[Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DSW_HI_AB_N70_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        	FSW_HI_AB_N70(x,j,:)=[Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FSW_HI_AB_N70)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                        end
                        for j=1:size(Listen,1)
                            DLW_HI_AB_N70(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DLW_HI_AB_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                            DLW_HI_AB_N70_B(x,j,:)=[Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DLW_HI_AB_N70_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                        	FLW_HI_AB_N70(x,j,:)=[Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FLW_HI_AB_N70)-length(WListenIdx(j,1):Listen(j,3)))'];
                        end
                    end
                end
            end

            % Storing Global Speaking/Listening windows
            for j=1:size(Speak,1)
                % Add nan-padding when necessary
                DGSW(x,j,:) = [Diameter(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(DGSW)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                DGSW_B(x,j,:) = [Diameter(WSpeakIdx(j,1):Speak(j,3))-mean(Diameter(Speak(j,2)-AdapBL*Param.Fs:Speak(j,2)));NaN*ones(1,length(DGSW_B)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                FGSW(x,j,:) = [Fixation(WSpeakIdx(j,1):Speak(j,3));NaN*ones(1,length(FGSW)-length(WSpeakIdx(j,1):Speak(j,3)))'];
                GSDur(x,j) = Speak(j,1);
            end

            for j=1:size(Listen,1)
                % Add nan-padding when necessary
                DGLW(x,j,:) = [Diameter(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(DGLW)-length(WListenIdx(j,1):Listen(j,3)))'];
                DGLW_B(x,j,:) = [Diameter(WListenIdx(j,1):Listen(j,3))-mean(Diameter(Listen(j,2)-AdapBL*Param.Fs:Listen(j,2)));NaN*ones(1,length(DGLW_B)-length(WListenIdx(j,1):Listen(j,3)))'];
                FGLW(x,j,:) = [Fixation(WListenIdx(j,1):Listen(j,3));NaN*ones(1,length(FGLW)-length(WListenIdx(j,1):Listen(j,3)))'];
                GLDur(x,j) = Listen(j,1);
            end
            
            % Store idx of non-rejected files associated with TPs
            if contains(ChosenFolder,'HI')
                TPsOrder_II(2*q,i) = x;
            elseif contains(ChosenFolder,'NH')
                TPsOrder_II(2*q-1,i) = x;
            end

            % Increase index of num of files used
            x=x+1;
            
        end
    end
end
%% Out-of-files-loop calculations
% Clean empty rows and layers, set 0's to NaN
DGSW(~any(DGSW,[2 3]),:,:)=[];DGSW(:,~any(DGSW,[1 3]),:)=[];DGSW(DGSW==0)=NaN;
DGSW_B(~any(DGSW_B,[2 3]),:,:)=[];DGSW_B(:,~any(DGSW_B,[1 3]),:)=[];DGSW_B(DGSW_B==0)=NaN;
DGLW(~any(DGLW,[2 3]),:,:)=[];DGLW(:,~any(DGLW,[1 3]),:)=[];DGLW(DGLW==0)=NaN;
DGLW_B(~any(DGLW_B,[2 3]),:,:)=[];DGLW_B(:,~any(DGLW_B,[1 3]),:)=[];DGLW_B(DGLW_B==0)=NaN;

DSW_NH(~any(DSW_NH,[2 3]),:,:)=[];DSW_NH(:,~any(DSW_NH,[1 3]),:)=[];DSW_NH(DSW_NH==0)=NaN;
DSW_NH_B(~any(DSW_NH_B,[2 3]),:,:)=[];DSW_NH_B(:,~any(DSW_NH_B,[1 3]),:)=[];DSW_NH_B(DSW_NH_B==0)=NaN;
DLW_NH(~any(DLW_NH,[2 3]),:,:)=[];DLW_NH(:,~any(DLW_NH,[1 3]),:)=[];DLW_NH(DLW_NH==0)=NaN;
DLW_NH_B(~any(DLW_NH_B,[2 3]),:,:)=[];DLW_NH_B(:,~any(DLW_NH_B,[1 3]),:)=[];DLW_NH_B(DLW_NH_B==0)=NaN;
DSW_HI(~any(DSW_HI,[2 3]),:,:)=[];DSW_HI(:,~any(DSW_HI,[1 3]),:)=[];DSW_HI(DSW_HI==0)=NaN;
DSW_HI_B(~any(DSW_HI_B,[2 3]),:,:)=[];DSW_HI_B(:,~any(DSW_HI_B,[1 3]),:)=[];DSW_HI_B(DSW_HI_B==0)=NaN;
DLW_HI(~any(DLW_HI,[2 3]),:,:)=[];DLW_HI(:,~any(DLW_HI,[1 3]),:)=[];DLW_HI(DLW_HI==0)=NaN;
DLW_HI_B(~any(DLW_HI_B,[2 3]),:,:)=[];DLW_HI_B(:,~any(DLW_HI_B,[1 3]),:)=[];DLW_HI_B(DLW_HI_B==0)=NaN;

DSW_N0(~any(DSW_N0,[2 3]),:,:)=[];DSW_N0(:,~any(DSW_N0,[1 3]),:)=[];DSW_N0(DSW_N0==0)=NaN;
DSW_N0_B(~any(DSW_N0_B,[2 3]),:,:)=[];DSW_N0_B(:,~any(DSW_N0_B,[1 3]),:)=[];DSW_N0_B(DSW_N0_B==0)=NaN;
DLW_N0(~any(DLW_N0,[2 3]),:,:)=[];DLW_N0(:,~any(DLW_N0,[1 3]),:)=[];DLW_N0(DLW_N0==0)=NaN;
DLW_N0_B(~any(DLW_N0_B,[2 3]),:,:)=[];DLW_N0_B(:,~any(DLW_N0_B,[1 3]),:)=[];DLW_N0_B(DLW_N0_B==0)=NaN;
DSW_N60(~any(DSW_N60,[2 3]),:,:)=[];DSW_N60(:,~any(DSW_N60,[1 3]),:)=[];DSW_N60(DSW_N60==0)=NaN;
DSW_N60_B(~any(DSW_N60_B,[2 3]),:,:)=[];DSW_N60_B(:,~any(DSW_N60_B,[1 3]),:)=[];DSW_N60_B(DSW_N60_B==0)=NaN;
DLW_N60(~any(DLW_N60,[2 3]),:,:)=[];DLW_N60(:,~any(DLW_N60,[1 3]),:)=[];DLW_N60(DLW_N60==0)=NaN;
DLW_N60_B(~any(DLW_N60_B,[2 3]),:,:)=[];DLW_N60_B(:,~any(DLW_N60_B,[1 3]),:)=[];DLW_N60_B(DLW_N60_B==0)=NaN;
DSW_N70(~any(DSW_N70,[2 3]),:,:)=[];DSW_N70(:,~any(DSW_N70,[1 3]),:)=[];DSW_N70(DSW_N70==0)=NaN;
DSW_N70_B(~any(DSW_N70_B,[2 3]),:,:)=[];DSW_N70_B(:,~any(DSW_N70_B,[1 3]),:)=[];DSW_N70_B(DSW_N70_B==0)=NaN;
DLW_N70(~any(DLW_N70,[2 3]),:,:)=[];DLW_N70(:,~any(DLW_N70,[1 3]),:)=[];DLW_N70(DLW_N70==0)=NaN;
DLW_N70_B(~any(DLW_N70_B,[2 3]),:,:)=[];DLW_N70_B(:,~any(DLW_N70_B,[1 3]),:)=[];DLW_N70_B(DLW_N70_B==0)=NaN;

DSW_NH_N0(~any(DSW_NH_N0,[2 3]),:,:)=[];DSW_NH_N0(:,~any(DSW_NH_N0,[1 3]),:)=[];DSW_NH_N0(DSW_NH_N0==0)=NaN;
DSW_NH_N0_B(~any(DSW_NH_N0_B,[2 3]),:,:)=[];DSW_NH_N0_B(:,~any(DSW_NH_N0_B,[1 3]),:)=[];DSW_NH_N0_B(DSW_NH_N0_B==0)=NaN;
DLW_NH_N0(~any(DLW_NH_N0,[2 3]),:,:)=[];DLW_NH_N0(:,~any(DLW_NH_N0,[1 3]),:)=[];DLW_NH_N0(DLW_NH_N0==0)=NaN;
DLW_NH_N0_B(~any(DLW_NH_N0_B,[2 3]),:,:)=[];DLW_NH_N0_B(:,~any(DLW_NH_N0_B,[1 3]),:)=[];DLW_NH_N0_B(DLW_NH_N0_B==0)=NaN;
DSW_NH_N60(~any(DSW_NH_N60,[2 3]),:,:)=[];DSW_NH_N60(:,~any(DSW_NH_N60,[1 3]),:)=[];DSW_NH_N60(DSW_NH_N60==0)=NaN;
DSW_NH_N60_B(~any(DSW_NH_N60_B,[2 3]),:,:)=[];DSW_NH_N60_B(:,~any(DSW_NH_N60_B,[1 3]),:)=[];DSW_NH_N60_B(DSW_NH_N60_B==0)=NaN;
DLW_NH_N60(~any(DLW_NH_N60,[2 3]),:,:)=[];DLW_NH_N60(:,~any(DLW_NH_N60,[1 3]),:)=[];DLW_NH_N60(DLW_NH_N60==0)=NaN;
DLW_NH_N60_B(~any(DLW_NH_N60_B,[2 3]),:,:)=[];DLW_NH_N60_B(:,~any(DLW_NH_N60_B,[1 3]),:)=[];DLW_NH_N60_B(DLW_NH_N60_B==0)=NaN;
DSW_NH_N70(~any(DSW_NH_N70,[2 3]),:,:)=[];DSW_NH_N70(:,~any(DSW_NH_N70,[1 3]),:)=[];DSW_NH_N70(DSW_NH_N70==0)=NaN;
DSW_NH_N70_B(~any(DSW_NH_N70_B,[2 3]),:,:)=[];DSW_NH_N70_B(:,~any(DSW_NH_N70_B,[1 3]),:)=[];DSW_NH_N70_B(DSW_NH_N70_B==0)=NaN;
DLW_NH_N70(~any(DLW_NH_N70,[2 3]),:,:)=[];DLW_NH_N70(:,~any(DLW_NH_N70,[1 3]),:)=[];DLW_NH_N70(DLW_NH_N70==0)=NaN;
DLW_NH_N70_B(~any(DLW_NH_N70_B,[2 3]),:,:)=[];DLW_NH_N70_B(:,~any(DLW_NH_N70_B,[1 3]),:)=[];DLW_NH_N70_B(DLW_NH_N70_B==0)=NaN;

DSW_HI_N0(~any(DSW_HI_N0,[2 3]),:,:)=[];DSW_HI_N0(:,~any(DSW_HI_N0,[1 3]),:)=[];DSW_HI_N0(DSW_HI_N0==0)=NaN;
DSW_HI_N0_B(~any(DSW_HI_N0_B,[2 3]),:,:)=[];DSW_HI_N0_B(:,~any(DSW_HI_N0_B,[1 3]),:)=[];DSW_HI_N0_B(DSW_HI_N0_B==0)=NaN;
DLW_HI_N0(~any(DLW_HI_N0,[2 3]),:,:)=[];DLW_HI_N0(:,~any(DLW_HI_N0,[1 3]),:)=[];DLW_HI_N0(DLW_HI_N0==0)=NaN;
DLW_HI_N0_B(~any(DLW_HI_N0_B,[2 3]),:,:)=[];DLW_HI_N0_B(:,~any(DLW_HI_N0_B,[1 3]),:)=[];DLW_HI_N0_B(DLW_HI_N0_B==0)=NaN;
DSW_HI_N60(~any(DSW_HI_N60,[2 3]),:,:)=[];DSW_HI_N60(:,~any(DSW_HI_N60,[1 3]),:)=[];DSW_HI_N60(DSW_HI_N60==0)=NaN;
DSW_HI_N60_B(~any(DSW_HI_N60_B,[2 3]),:,:)=[];DSW_HI_N60_B(:,~any(DSW_HI_N60_B,[1 3]),:)=[];DSW_HI_N60_B(DSW_HI_N60_B==0)=NaN;
DLW_HI_N60(~any(DLW_HI_N60,[2 3]),:,:)=[];DLW_HI_N60(:,~any(DLW_HI_N60,[1 3]),:)=[];DLW_HI_N60(DLW_HI_N60==0)=NaN;
DLW_HI_N60_B(~any(DLW_HI_N60_B,[2 3]),:,:)=[];DLW_HI_N60_B(:,~any(DLW_HI_N60_B,[1 3]),:)=[];DLW_HI_N60_B(DLW_HI_N60_B==0)=NaN;
DSW_HI_N70(~any(DSW_HI_N70,[2 3]),:,:)=[];DSW_HI_N70(:,~any(DSW_HI_N70,[1 3]),:)=[];DSW_HI_N70(DSW_HI_N70==0)=NaN;
DSW_HI_N70_B(~any(DSW_HI_N70_B,[2 3]),:,:)=[];DSW_HI_N70_B(:,~any(DSW_HI_N70_B,[1 3]),:)=[];DSW_HI_N70_B(DSW_HI_N70_B==0)=NaN;
DLW_HI_N70(~any(DLW_HI_N70,[2 3]),:,:)=[];DLW_HI_N70(:,~any(DLW_HI_N70,[1 3]),:)=[];DLW_HI_N70(DLW_HI_N70==0)=NaN;
DLW_HI_N70_B(~any(DLW_HI_N70_B,[2 3]),:,:)=[];DLW_HI_N70_B(:,~any(DLW_HI_N70_B,[1 3]),:)=[];DLW_HI_N70_B(DLW_HI_N70_B==0)=NaN;

DSW_NH_UN(~any(DSW_NH_UN,[2 3]),:,:)=[];DSW_NH_UN(:,~any(DSW_NH_UN,[1 3]),:)=[];DSW_NH_UN(DSW_NH_UN==0)=NaN;
DSW_NH_UN_B(~any(DSW_NH_UN_B,[2 3]),:,:)=[];DSW_NH_UN_B(:,~any(DSW_NH_UN_B,[1 3]),:)=[];DSW_NH_UN_B(DSW_NH_UN_B==0)=NaN;
DLW_NH_UN(~any(DLW_NH_UN,[2 3]),:,:)=[];DLW_NH_UN(:,~any(DLW_NH_UN,[1 3]),:)=[];DLW_NH_UN(DLW_NH_UN==0)=NaN;
DLW_NH_UN_B(~any(DLW_NH_UN_B,[2 3]),:,:)=[];DLW_NH_UN_B(:,~any(DLW_NH_UN_B,[1 3]),:)=[];DLW_NH_UN_B(DLW_NH_UN_B==0)=NaN;
DSW_NH_AA(~any(DSW_NH_AA,[2 3]),:,:)=[];DSW_NH_AA(:,~any(DSW_NH_AA,[1 3]),:)=[];DSW_NH_AA(DSW_NH_AA==0)=NaN;
DSW_NH_AA_B(~any(DSW_NH_AA_B,[2 3]),:,:)=[];DSW_NH_AA_B(:,~any(DSW_NH_AA_B,[1 3]),:)=[];DSW_NH_AA_B(DSW_NH_AA_B==0)=NaN;
DLW_NH_AA(~any(DLW_NH_AA,[2 3]),:,:)=[];DLW_NH_AA(:,~any(DLW_NH_AA,[1 3]),:)=[];DLW_NH_AA(DLW_NH_AA==0)=NaN;
DLW_NH_AA_B(~any(DLW_NH_AA_B,[2 3]),:,:)=[];DLW_NH_AA_B(:,~any(DLW_NH_AA_B,[1 3]),:)=[];DLW_NH_AA_B(DLW_NH_AA_B==0)=NaN;
DSW_NH_AB(~any(DSW_NH_AB,[2 3]),:,:)=[];DSW_NH_AB(:,~any(DSW_NH_AB,[1 3]),:)=[];DSW_NH_AB(DSW_NH_AB==0)=NaN;
DSW_NH_AB_B(~any(DSW_NH_AB_B,[2 3]),:,:)=[];DSW_NH_AB_B(:,~any(DSW_NH_AB_B,[1 3]),:)=[];DSW_NH_AB_B(DSW_NH_AB_B==0)=NaN;
DLW_NH_AB(~any(DLW_NH_AB,[2 3]),:,:)=[];DLW_NH_AB(:,~any(DLW_NH_AB,[1 3]),:)=[];DLW_NH_AB(DLW_NH_AB==0)=NaN;
DLW_NH_AB_B(~any(DLW_NH_AB_B,[2 3]),:,:)=[];DLW_NH_AB_B(:,~any(DLW_NH_AB_B,[1 3]),:)=[];DLW_NH_AB_B(DLW_NH_AB_B==0)=NaN;

DSW_HI_UN(~any(DSW_HI_UN,[2 3]),:,:)=[];DSW_HI_UN(:,~any(DSW_HI_UN,[1 3]),:)=[];DSW_HI_UN(DSW_HI_UN==0)=NaN;
DSW_HI_UN_B(~any(DSW_HI_UN_B,[2 3]),:,:)=[];DSW_HI_UN_B(:,~any(DSW_HI_UN_B,[1 3]),:)=[];DSW_HI_UN_B(DSW_HI_UN_B==0)=NaN;
DLW_HI_UN(~any(DLW_HI_UN,[2 3]),:,:)=[];DLW_HI_UN(:,~any(DLW_HI_UN,[1 3]),:)=[];DLW_HI_UN(DLW_HI_UN==0)=NaN;
DLW_HI_UN_B(~any(DLW_HI_UN_B,[2 3]),:,:)=[];DLW_HI_UN_B(:,~any(DLW_HI_UN_B,[1 3]),:)=[];DLW_HI_UN_B(DLW_HI_UN_B==0)=NaN;
DSW_HI_AA(~any(DSW_HI_AA,[2 3]),:,:)=[];DSW_HI_AA(:,~any(DSW_HI_AA,[1 3]),:)=[];DSW_HI_AA(DSW_HI_AA==0)=NaN;
DSW_HI_AA_B(~any(DSW_HI_AA_B,[2 3]),:,:)=[];DSW_HI_AA_B(:,~any(DSW_HI_AA_B,[1 3]),:)=[];DSW_HI_AA_B(DSW_HI_AA_B==0)=NaN;
DLW_HI_AA(~any(DLW_HI_AA,[2 3]),:,:)=[];DLW_HI_AA(:,~any(DLW_HI_AA,[1 3]),:)=[];DLW_HI_AA(DLW_HI_AA==0)=NaN;
DLW_HI_AA_B(~any(DLW_HI_AA_B,[2 3]),:,:)=[];DLW_HI_AA_B(:,~any(DLW_HI_AA_B,[1 3]),:)=[];DLW_HI_AA_B(DLW_HI_AA_B==0)=NaN;
DSW_HI_AB(~any(DSW_HI_AB,[2 3]),:,:)=[];DSW_HI_AB(:,~any(DSW_HI_AB,[1 3]),:)=[];DSW_HI_AB(DSW_HI_AB==0)=NaN;
DSW_HI_AB_B(~any(DSW_HI_AB_B,[2 3]),:,:)=[];DSW_HI_AB_B(:,~any(DSW_HI_AB_B,[1 3]),:)=[];DSW_HI_AB_B(DSW_HI_AB_B==0)=NaN;
DLW_HI_AB(~any(DLW_HI_AB,[2 3]),:,:)=[];DLW_HI_AB(:,~any(DLW_HI_AB,[1 3]),:)=[];DLW_HI_AB(DLW_HI_AB==0)=NaN;
DLW_HI_AB_B(~any(DLW_HI_AB_B,[2 3]),:,:)=[];DLW_HI_AB_B(:,~any(DLW_HI_AB_B,[1 3]),:)=[];DLW_HI_AB_B(DLW_HI_AB_B==0)=NaN;

DSW_NH_UN_N0(~any(DSW_NH_UN_N0,[2 3]),:,:)=[];DSW_NH_UN_N0(:,~any(DSW_NH_UN_N0,[1 3]),:)=[];DSW_NH_UN_N0(DSW_NH_UN_N0==0)=NaN;
DSW_NH_UN_N0_B(~any(DSW_NH_UN_N0_B,[2 3]),:,:)=[];DSW_NH_UN_N0_B(:,~any(DSW_NH_UN_N0_B,[1 3]),:)=[];DSW_NH_UN_N0_B(DSW_NH_UN_N0_B==0)=NaN;
DLW_NH_UN_N0(~any(DLW_NH_UN_N0,[2 3]),:,:)=[];DLW_NH_UN_N0(:,~any(DLW_NH_UN_N0,[1 3]),:)=[];DLW_NH_UN_N0(DLW_NH_UN_N0==0)=NaN;
DLW_NH_UN_N0_B(~any(DLW_NH_UN_N0_B,[2 3]),:,:)=[];DLW_NH_UN_N0_B(:,~any(DLW_NH_UN_N0_B,[1 3]),:)=[];DLW_NH_UN_N0_B(DLW_NH_UN_N0_B==0)=NaN;
DSW_NH_UN_N60(~any(DSW_NH_UN_N60,[2 3]),:,:)=[];DSW_NH_UN_N60(:,~any(DSW_NH_UN_N60,[1 3]),:)=[];DSW_NH_UN_N60(DSW_NH_UN_N60==0)=NaN;
DSW_NH_UN_N60_B(~any(DSW_NH_UN_N60_B,[2 3]),:,:)=[];DSW_NH_UN_N60_B(:,~any(DSW_NH_UN_N60_B,[1 3]),:)=[];DSW_NH_UN_N60_B(DSW_NH_UN_N60_B==0)=NaN;
DLW_NH_UN_N60(~any(DLW_NH_UN_N60,[2 3]),:,:)=[];DLW_NH_UN_N60(:,~any(DLW_NH_UN_N60,[1 3]),:)=[];DLW_NH_UN_N60(DLW_NH_UN_N60==0)=NaN;
DLW_NH_UN_N60_B(~any(DLW_NH_UN_N60_B,[2 3]),:,:)=[];DLW_NH_UN_N60_B(:,~any(DLW_NH_UN_N60_B,[1 3]),:)=[];DLW_NH_UN_N60_B(DLW_NH_UN_N60_B==0)=NaN;
DSW_NH_UN_N70(~any(DSW_NH_UN_N70,[2 3]),:,:)=[];DSW_NH_UN_N70(:,~any(DSW_NH_UN_N70,[1 3]),:)=[];DSW_NH_UN_N70(DSW_NH_UN_N70==0)=NaN;
DSW_NH_UN_N70_B(~any(DSW_NH_UN_N70_B,[2 3]),:,:)=[];DSW_NH_UN_N70_B(:,~any(DSW_NH_UN_N70_B,[1 3]),:)=[];DSW_NH_UN_N70_B(DSW_NH_UN_N70_B==0)=NaN;
DLW_NH_UN_N70(~any(DLW_NH_UN_N70,[2 3]),:,:)=[];DLW_NH_UN_N70(:,~any(DLW_NH_UN_N70,[1 3]),:)=[];DLW_NH_UN_N70(DLW_NH_UN_N70==0)=NaN;
DLW_NH_UN_N70_B(~any(DLW_NH_UN_N70_B,[2 3]),:,:)=[];DLW_NH_UN_N70_B(:,~any(DLW_NH_UN_N70_B,[1 3]),:)=[];DLW_NH_UN_N70_B(DLW_NH_UN_N70_B==0)=NaN;
DSW_NH_AA_N0(~any(DSW_NH_AA_N0,[2 3]),:,:)=[];DSW_NH_AA_N0(:,~any(DSW_NH_AA_N0,[1 3]),:)=[];DSW_NH_AA_N0(DSW_NH_AA_N0==0)=NaN;
DSW_NH_AA_N0_B(~any(DSW_NH_AA_N0_B,[2 3]),:,:)=[];DSW_NH_AA_N0_B(:,~any(DSW_NH_AA_N0_B,[1 3]),:)=[];DSW_NH_AA_N0_B(DSW_NH_AA_N0_B==0)=NaN;
DLW_NH_AA_N0(~any(DLW_NH_AA_N0,[2 3]),:,:)=[];DLW_NH_AA_N0(:,~any(DLW_NH_AA_N0,[1 3]),:)=[];DLW_NH_AA_N0(DLW_NH_AA_N0==0)=NaN;
DLW_NH_AA_N0_B(~any(DLW_NH_AA_N0_B,[2 3]),:,:)=[];DLW_NH_AA_N0_B(:,~any(DLW_NH_AA_N0_B,[1 3]),:)=[];DLW_NH_AA_N0_B(DLW_NH_AA_N0_B==0)=NaN;
DSW_NH_AA_N60(~any(DSW_NH_AA_N60,[2 3]),:,:)=[];DSW_NH_AA_N60(:,~any(DSW_NH_AA_N60,[1 3]),:)=[];DSW_NH_AA_N60(DSW_NH_AA_N60==0)=NaN;
DSW_NH_AA_N60_B(~any(DSW_NH_AA_N60_B,[2 3]),:,:)=[];DSW_NH_AA_N60_B(:,~any(DSW_NH_AA_N60_B,[1 3]),:)=[];DSW_NH_AA_N60_B(DSW_NH_AA_N60_B==0)=NaN;
DLW_NH_AA_N60(~any(DLW_NH_AA_N60,[2 3]),:,:)=[];DLW_NH_AA_N60(:,~any(DLW_NH_AA_N60,[1 3]),:)=[];DLW_NH_AA_N60(DLW_NH_AA_N60==0)=NaN;
DLW_NH_AA_N60_B(~any(DLW_NH_AA_N60_B,[2 3]),:,:)=[];DLW_NH_AA_N60_B(:,~any(DLW_NH_AA_N60_B,[1 3]),:)=[];DLW_NH_AA_N60_B(DLW_NH_AA_N60_B==0)=NaN;
DSW_NH_AA_N70(~any(DSW_NH_AA_N70,[2 3]),:,:)=[];DSW_NH_AA_N70(:,~any(DSW_NH_AA_N70,[1 3]),:)=[];DSW_NH_AA_N70(DSW_NH_AA_N70==0)=NaN;
DSW_NH_AA_N70_B(~any(DSW_NH_AA_N70_B,[2 3]),:,:)=[];DSW_NH_AA_N70_B(:,~any(DSW_NH_AA_N70_B,[1 3]),:)=[];DSW_NH_AA_N70_B(DSW_NH_AA_N70_B==0)=NaN;
DLW_NH_AA_N70(~any(DLW_NH_AA_N70,[2 3]),:,:)=[];DLW_NH_AA_N70(:,~any(DLW_NH_AA_N70,[1 3]),:)=[];DLW_NH_AA_N70(DLW_NH_AA_N70==0)=NaN;
DLW_NH_AA_N70_B(~any(DLW_NH_AA_N70_B,[2 3]),:,:)=[];DLW_NH_AA_N70_B(:,~any(DLW_NH_AA_N70_B,[1 3]),:)=[];DLW_NH_AA_N70_B(DLW_NH_AA_N70_B==0)=NaN;
DSW_NH_AB_N0(~any(DSW_NH_AB_N0,[2 3]),:,:)=[];DSW_NH_AB_N0(:,~any(DSW_NH_AB_N0,[1 3]),:)=[];DSW_NH_AB_N0(DSW_NH_AB_N0==0)=NaN;
DSW_NH_AB_N0_B(~any(DSW_NH_AB_N0_B,[2 3]),:,:)=[];DSW_NH_AB_N0_B(:,~any(DSW_NH_AB_N0_B,[1 3]),:)=[];DSW_NH_AB_N0_B(DSW_NH_AB_N0_B==0)=NaN;
DLW_NH_AB_N0(~any(DLW_NH_AB_N0,[2 3]),:,:)=[];DLW_NH_AB_N0(:,~any(DLW_NH_AB_N0,[1 3]),:)=[];DLW_NH_AB_N0(DLW_NH_AB_N0==0)=NaN;
DLW_NH_AB_N0_B(~any(DLW_NH_AB_N0_B,[2 3]),:,:)=[];DLW_NH_AB_N0_B(:,~any(DLW_NH_AB_N0_B,[1 3]),:)=[];DLW_NH_AB_N0_B(DLW_NH_AB_N0_B==0)=NaN;
DSW_NH_AB_N60(~any(DSW_NH_AB_N60,[2 3]),:,:)=[];DSW_NH_AB_N60(:,~any(DSW_NH_AB_N60,[1 3]),:)=[];DSW_NH_AB_N60(DSW_NH_AB_N60==0)=NaN;
DSW_NH_AB_N60_B(~any(DSW_NH_AB_N60_B,[2 3]),:,:)=[];DSW_NH_AB_N60_B(:,~any(DSW_NH_AB_N60_B,[1 3]),:)=[];DSW_NH_AB_N60_B(DSW_NH_AB_N60_B==0)=NaN;
DLW_NH_AB_N60(~any(DLW_NH_AB_N60,[2 3]),:,:)=[];DLW_NH_AB_N60(:,~any(DLW_NH_AB_N60,[1 3]),:)=[];DLW_NH_AB_N60(DLW_NH_AB_N60==0)=NaN;
DLW_NH_AB_N60_B(~any(DLW_NH_AB_N60_B,[2 3]),:,:)=[];DLW_NH_AB_N60_B(:,~any(DLW_NH_AB_N60_B,[1 3]),:)=[];DLW_NH_AB_N60_B(DLW_NH_AB_N60_B==0)=NaN;
DSW_NH_AB_N70(~any(DSW_NH_AB_N70,[2 3]),:,:)=[];DSW_NH_AB_N70(:,~any(DSW_NH_AB_N70,[1 3]),:)=[];DSW_NH_AB_N70(DSW_NH_AB_N70==0)=NaN;
DSW_NH_AB_N70_B(~any(DSW_NH_AB_N70_B,[2 3]),:,:)=[];DSW_NH_AB_N70_B(:,~any(DSW_NH_AB_N70_B,[1 3]),:)=[];DSW_NH_AB_N70_B(DSW_NH_AB_N70_B==0)=NaN;
DLW_NH_AB_N70(~any(DLW_NH_AB_N70,[2 3]),:,:)=[];DLW_NH_AB_N70(:,~any(DLW_NH_AB_N70,[1 3]),:)=[];DLW_NH_AB_N70(DLW_NH_AB_N70==0)=NaN;
DLW_NH_AB_N70_B(~any(DLW_NH_AB_N70_B,[2 3]),:,:)=[];DLW_NH_AB_N70_B(:,~any(DLW_NH_AB_N70_B,[1 3]),:)=[];DLW_NH_AB_N70_B(DLW_NH_AB_N70_B==0)=NaN;

DSW_HI_UN_N0(~any(DSW_HI_UN_N0,[2 3]),:,:)=[];DSW_HI_UN_N0(:,~any(DSW_HI_UN_N0,[1 3]),:)=[];DSW_HI_UN_N0(DSW_HI_UN_N0==0)=NaN;
DSW_HI_UN_N0_B(~any(DSW_HI_UN_N0_B,[2 3]),:,:)=[];DSW_HI_UN_N0_B(:,~any(DSW_HI_UN_N0_B,[1 3]),:)=[];DSW_HI_UN_N0_B(DSW_HI_UN_N0_B==0)=NaN;
DLW_HI_UN_N0(~any(DLW_HI_UN_N0,[2 3]),:,:)=[];DLW_HI_UN_N0(:,~any(DLW_HI_UN_N0,[1 3]),:)=[];DLW_HI_UN_N0(DLW_HI_UN_N0==0)=NaN;
DLW_HI_UN_N0_B(~any(DLW_HI_UN_N0_B,[2 3]),:,:)=[];DLW_HI_UN_N0_B(:,~any(DLW_HI_UN_N0_B,[1 3]),:)=[];DLW_HI_UN_N0_B(DLW_HI_UN_N0_B==0)=NaN;
DSW_HI_UN_N60(~any(DSW_HI_UN_N60,[2 3]),:,:)=[];DSW_HI_UN_N60(:,~any(DSW_HI_UN_N60,[1 3]),:)=[];DSW_HI_UN_N60(DSW_HI_UN_N60==0)=NaN;
DSW_HI_UN_N60_B(~any(DSW_HI_UN_N60_B,[2 3]),:,:)=[];DSW_HI_UN_N60_B(:,~any(DSW_HI_UN_N60_B,[1 3]),:)=[];DSW_HI_UN_N60_B(DSW_HI_UN_N60_B==0)=NaN;
DLW_HI_UN_N60(~any(DLW_HI_UN_N60,[2 3]),:,:)=[];DLW_HI_UN_N60(:,~any(DLW_HI_UN_N60,[1 3]),:)=[];DLW_HI_UN_N60(DLW_HI_UN_N60==0)=NaN;
DLW_HI_UN_N60_B(~any(DLW_HI_UN_N60_B,[2 3]),:,:)=[];DLW_HI_UN_N60_B(:,~any(DLW_HI_UN_N60_B,[1 3]),:)=[];DLW_HI_UN_N60_B(DLW_HI_UN_N60_B==0)=NaN;
DSW_HI_UN_N70(~any(DSW_HI_UN_N70,[2 3]),:,:)=[];DSW_HI_UN_N70(:,~any(DSW_HI_UN_N70,[1 3]),:)=[];DSW_HI_UN_N70(DSW_HI_UN_N70==0)=NaN;
DSW_HI_UN_N70_B(~any(DSW_HI_UN_N70_B,[2 3]),:,:)=[];DSW_HI_UN_N70_B(:,~any(DSW_HI_UN_N70_B,[1 3]),:)=[];DSW_HI_UN_N70_B(DSW_HI_UN_N70_B==0)=NaN;
DLW_HI_UN_N70(~any(DLW_HI_UN_N70,[2 3]),:,:)=[];DLW_HI_UN_N70(:,~any(DLW_HI_UN_N70,[1 3]),:)=[];DLW_HI_UN_N70(DLW_HI_UN_N70==0)=NaN;
DLW_HI_UN_N70_B(~any(DLW_HI_UN_N70_B,[2 3]),:,:)=[];DLW_HI_UN_N70_B(:,~any(DLW_HI_UN_N70_B,[1 3]),:)=[];DLW_HI_UN_N70_B(DLW_HI_UN_N70_B==0)=NaN;
DSW_HI_AA_N0(~any(DSW_HI_AA_N0,[2 3]),:,:)=[];DSW_HI_AA_N0(:,~any(DSW_HI_AA_N0,[1 3]),:)=[];DSW_HI_AA_N0(DSW_HI_AA_N0==0)=NaN;
DSW_HI_AA_N0_B(~any(DSW_HI_AA_N0_B,[2 3]),:,:)=[];DSW_HI_AA_N0_B(:,~any(DSW_HI_AA_N0_B,[1 3]),:)=[];DSW_HI_AA_N0_B(DSW_HI_AA_N0_B==0)=NaN;
DLW_HI_AA_N0(~any(DLW_HI_AA_N0,[2 3]),:,:)=[];DLW_HI_AA_N0(:,~any(DLW_HI_AA_N0,[1 3]),:)=[];DLW_HI_AA_N0(DLW_HI_AA_N0==0)=NaN;
DLW_HI_AA_N0_B(~any(DLW_HI_AA_N0_B,[2 3]),:,:)=[];DLW_HI_AA_N0_B(:,~any(DLW_HI_AA_N0_B,[1 3]),:)=[];DLW_HI_AA_N0_B(DLW_HI_AA_N0_B==0)=NaN;
DSW_HI_AA_N60(~any(DSW_HI_AA_N60,[2 3]),:,:)=[];DSW_HI_AA_N60(:,~any(DSW_HI_AA_N60,[1 3]),:)=[];DSW_HI_AA_N60(DSW_HI_AA_N60==0)=NaN;
DSW_HI_AA_N60_B(~any(DSW_HI_AA_N60_B,[2 3]),:,:)=[];DSW_HI_AA_N60_B(:,~any(DSW_HI_AA_N60_B,[1 3]),:)=[];DSW_HI_AA_N60_B(DSW_HI_AA_N60_B==0)=NaN;
DLW_HI_AA_N60(~any(DLW_HI_AA_N60,[2 3]),:,:)=[];DLW_HI_AA_N60(:,~any(DLW_HI_AA_N60,[1 3]),:)=[];DLW_HI_AA_N60(DLW_HI_AA_N60==0)=NaN;
DLW_HI_AA_N60_B(~any(DLW_HI_AA_N60_B,[2 3]),:,:)=[];DLW_HI_AA_N60_B(:,~any(DLW_HI_AA_N60_B,[1 3]),:)=[];DLW_HI_AA_N60_B(DLW_HI_AA_N60_B==0)=NaN;
DSW_HI_AA_N70(~any(DSW_HI_AA_N70,[2 3]),:,:)=[];DSW_HI_AA_N70(:,~any(DSW_HI_AA_N70,[1 3]),:)=[];DSW_HI_AA_N70(DSW_HI_AA_N70==0)=NaN;
DSW_HI_AA_N70_B(~any(DSW_HI_AA_N70_B,[2 3]),:,:)=[];DSW_HI_AA_N70_B(:,~any(DSW_HI_AA_N70_B,[1 3]),:)=[];DSW_HI_AA_N70_B(DSW_HI_AA_N70_B==0)=NaN;
DLW_HI_AA_N70(~any(DLW_HI_AA_N70,[2 3]),:,:)=[];DLW_HI_AA_N70(:,~any(DLW_HI_AA_N70,[1 3]),:)=[];DLW_HI_AA_N70(DLW_HI_AA_N70==0)=NaN;
DLW_HI_AA_N70_B(~any(DLW_HI_AA_N70_B,[2 3]),:,:)=[];DLW_HI_AA_N70_B(:,~any(DLW_HI_AA_N70_B,[1 3]),:)=[];DLW_HI_AA_N70_B(DLW_HI_AA_N70_B==0)=NaN;
DSW_HI_AB_N0(~any(DSW_HI_AB_N0,[2 3]),:,:)=[];DSW_HI_AB_N0(:,~any(DSW_HI_AB_N0,[1 3]),:)=[];DSW_HI_AB_N0(DSW_HI_AB_N0==0)=NaN;
DSW_HI_AB_N0_B(~any(DSW_HI_AB_N0_B,[2 3]),:,:)=[];DSW_HI_AB_N0_B(:,~any(DSW_HI_AB_N0_B,[1 3]),:)=[];DSW_HI_AB_N0_B(DSW_HI_AB_N0_B==0)=NaN;
DLW_HI_AB_N0(~any(DLW_HI_AB_N0,[2 3]),:,:)=[];DLW_HI_AB_N0(:,~any(DLW_HI_AB_N0,[1 3]),:)=[];DLW_HI_AB_N0(DLW_HI_AB_N0==0)=NaN;
DLW_HI_AB_N0_B(~any(DLW_HI_AB_N0_B,[2 3]),:,:)=[];DLW_HI_AB_N0_B(:,~any(DLW_HI_AB_N0_B,[1 3]),:)=[];DLW_HI_AB_N0_B(DLW_HI_AB_N0_B==0)=NaN;
DSW_HI_AB_N60(~any(DSW_HI_AB_N60,[2 3]),:,:)=[];DSW_HI_AB_N60(:,~any(DSW_HI_AB_N60,[1 3]),:)=[];DSW_HI_AB_N60(DSW_HI_AB_N60==0)=NaN;
DSW_HI_AB_N60_B(~any(DSW_HI_AB_N60_B,[2 3]),:,:)=[];DSW_HI_AB_N60_B(:,~any(DSW_HI_AB_N60_B,[1 3]),:)=[];DSW_HI_AB_N60_B(DSW_HI_AB_N60_B==0)=NaN;
DLW_HI_AB_N60(~any(DLW_HI_AB_N60,[2 3]),:,:)=[];DLW_HI_AB_N60(:,~any(DLW_HI_AB_N60,[1 3]),:)=[];DLW_HI_AB_N60(DLW_HI_AB_N60==0)=NaN;
DLW_HI_AB_N60_B(~any(DLW_HI_AB_N60_B,[2 3]),:,:)=[];DLW_HI_AB_N60_B(:,~any(DLW_HI_AB_N60_B,[1 3]),:)=[];DLW_HI_AB_N60_B(DLW_HI_AB_N60_B==0)=NaN;
DSW_HI_AB_N70(~any(DSW_HI_AB_N70,[2 3]),:,:)=[];DSW_HI_AB_N70(:,~any(DSW_HI_AB_N70,[1 3]),:)=[];DSW_HI_AB_N70(DSW_HI_AB_N70==0)=NaN;
DSW_HI_AB_N70_B(~any(DSW_HI_AB_N70_B,[2 3]),:,:)=[];DSW_HI_AB_N70_B(:,~any(DSW_HI_AB_N70_B,[1 3]),:)=[];DSW_HI_AB_N70_B(DSW_HI_AB_N70_B==0)=NaN;
DLW_HI_AB_N70(~any(DLW_HI_AB_N70,[2 3]),:,:)=[];DLW_HI_AB_N70(:,~any(DLW_HI_AB_N70,[1 3]),:)=[];DLW_HI_AB_N70(DLW_HI_AB_N70==0)=NaN;
DLW_HI_AB_N70_B(~any(DLW_HI_AB_N70_B,[2 3]),:,:)=[];DLW_HI_AB_N70_B(:,~any(DLW_HI_AB_N70_B,[1 3]),:)=[];DLW_HI_AB_N70_B(DLW_HI_AB_N70_B==0)=NaN;

FGSW(~any(FGSW,[2 3]),:,:)=[];FGSW(:,~any(FGSW,[1 3]),:)=[];FGSW(FGSW==0)=NaN;
FGLW(~any(FGLW,[2 3]),:,:)=[];FGLW(:,~any(FGLW,[1 3]),:)=[];FGLW(FGLW==0)=NaN;

FSW_NH(~any(FSW_NH,[2 3]),:,:)=[];FSW_NH(:,~any(FSW_NH,[1 3]),:)=[];FSW_NH(FSW_NH==0)=NaN;
FLW_NH(~any(FLW_NH,[2 3]),:,:)=[];FLW_NH(:,~any(FLW_NH,[1 3]),:)=[];FLW_NH(FLW_NH==0)=NaN;
FSW_HI(~any(FSW_HI,[2 3]),:,:)=[];FSW_HI(:,~any(FSW_HI,[1 3]),:)=[];FSW_HI(FSW_HI==0)=NaN;
FLW_HI(~any(FLW_HI,[2 3]),:,:)=[];FLW_HI(:,~any(FLW_HI,[1 3]),:)=[];FLW_HI(FLW_HI==0)=NaN;

FSW_N0(~any(FSW_N0,[2 3]),:,:)=[];FSW_N0(:,~any(FSW_N0,[1 3]),:)=[];FSW_N0(FSW_N0==0)=NaN;
FLW_N0(~any(FLW_N0,[2 3]),:,:)=[];FLW_N0(:,~any(FLW_N0,[1 3]),:)=[];FLW_N0(FLW_N0==0)=NaN;
FSW_N60(~any(FSW_N60,[2 3]),:,:)=[];FSW_N60(:,~any(FSW_N60,[1 3]),:)=[];FSW_N60(FSW_N60==0)=NaN;
FLW_N60(~any(FLW_N60,[2 3]),:,:)=[];FLW_N60(:,~any(FLW_N60,[1 3]),:)=[];FLW_N60(FLW_N60==0)=NaN;
FSW_N70(~any(FSW_N70,[2 3]),:,:)=[];FSW_N70(:,~any(FSW_N70,[1 3]),:)=[];FSW_N70(FSW_N70==0)=NaN;
FLW_N70(~any(FLW_N70,[2 3]),:,:)=[];FLW_N70(:,~any(FLW_N70,[1 3]),:)=[];FLW_N70(FLW_N70==0)=NaN;

FSW_NH_N0(~any(FSW_NH_N0,[2 3]),:,:)=[];FSW_NH_N0(:,~any(FSW_NH_N0,[1 3]),:)=[];FSW_NH_N0(FSW_NH_N0==0)=NaN;
FLW_NH_N0(~any(FLW_NH_N0,[2 3]),:,:)=[];FLW_NH_N0(:,~any(FLW_NH_N0,[1 3]),:)=[];FLW_NH_N0(FLW_NH_N0==0)=NaN;
FSW_NH_N60(~any(FSW_NH_N60,[2 3]),:,:)=[];FSW_NH_N60(:,~any(FSW_NH_N60,[1 3]),:)=[];FSW_NH_N60(FSW_NH_N60==0)=NaN;
FLW_NH_N60(~any(FLW_NH_N60,[2 3]),:,:)=[];FLW_NH_N60(:,~any(FLW_NH_N60,[1 3]),:)=[];FLW_NH_N60(FLW_NH_N60==0)=NaN;
FSW_NH_N70(~any(FSW_NH_N70,[2 3]),:,:)=[];FSW_NH_N70(:,~any(FSW_NH_N70,[1 3]),:)=[];FSW_NH_N70(FSW_NH_N70==0)=NaN;
FLW_NH_N70(~any(FLW_NH_N70,[2 3]),:,:)=[];FLW_NH_N70(:,~any(FLW_NH_N70,[1 3]),:)=[];FLW_NH_N70(FLW_NH_N70==0)=NaN;

FSW_HI_N0(~any(FSW_HI_N0,[2 3]),:,:)=[];FSW_HI_N0(:,~any(FSW_HI_N0,[1 3]),:)=[];FSW_HI_N0(FSW_HI_N0==0)=NaN;
FLW_HI_N0(~any(FLW_HI_N0,[2 3]),:,:)=[];FLW_HI_N0(:,~any(FLW_HI_N0,[1 3]),:)=[];FLW_HI_N0(FLW_HI_N0==0)=NaN;
FSW_HI_N60(~any(FSW_HI_N60,[2 3]),:,:)=[];FSW_HI_N60(:,~any(FSW_HI_N60,[1 3]),:)=[];FSW_HI_N60(FSW_HI_N60==0)=NaN;
FLW_HI_N60(~any(FLW_HI_N60,[2 3]),:,:)=[];FLW_HI_N60(:,~any(FLW_HI_N60,[1 3]),:)=[];FLW_HI_N60(FLW_HI_N60==0)=NaN;
FSW_HI_N70(~any(FSW_HI_N70,[2 3]),:,:)=[];FSW_HI_N70(:,~any(FSW_HI_N70,[1 3]),:)=[];FSW_HI_N70(FSW_HI_N70==0)=NaN;
FLW_HI_N70(~any(FLW_HI_N70,[2 3]),:,:)=[];FLW_HI_N70(:,~any(FLW_HI_N70,[1 3]),:)=[];FLW_HI_N70(FLW_HI_N70==0)=NaN;

FSW_NH_UN(~any(FSW_NH_UN,[2 3]),:,:)=[];FSW_NH_UN(:,~any(FSW_NH_UN,[1 3]),:)=[];FSW_NH_UN(FSW_NH_UN==0)=NaN;
FLW_NH_UN(~any(FLW_NH_UN,[2 3]),:,:)=[];FLW_NH_UN(:,~any(FLW_NH_UN,[1 3]),:)=[];FLW_NH_UN(FLW_NH_UN==0)=NaN;
FSW_NH_AA(~any(FSW_NH_AA,[2 3]),:,:)=[];FSW_NH_AA(:,~any(FSW_NH_AA,[1 3]),:)=[];FSW_NH_AA(FSW_NH_AA==0)=NaN;
FLW_NH_AA(~any(FLW_NH_AA,[2 3]),:,:)=[];FLW_NH_AA(:,~any(FLW_NH_AA,[1 3]),:)=[];FLW_NH_AA(FLW_NH_AA==0)=NaN;
FSW_NH_AB(~any(FSW_NH_AB,[2 3]),:,:)=[];FSW_NH_AB(:,~any(FSW_NH_AB,[1 3]),:)=[];FSW_NH_AB(FSW_NH_AB==0)=NaN;
FLW_NH_AB(~any(FLW_NH_AB,[2 3]),:,:)=[];FLW_NH_AB(:,~any(FLW_NH_AB,[1 3]),:)=[];FLW_NH_AB(FLW_NH_AB==0)=NaN;

FSW_HI_UN(~any(FSW_HI_UN,[2 3]),:,:)=[];FSW_HI_UN(:,~any(FSW_HI_UN,[1 3]),:)=[];FSW_HI_UN(FSW_HI_UN==0)=NaN;
FLW_HI_UN(~any(FLW_HI_UN,[2 3]),:,:)=[];FLW_HI_UN(:,~any(FLW_HI_UN,[1 3]),:)=[];FLW_HI_UN(FLW_HI_UN==0)=NaN;
FSW_HI_AA(~any(FSW_HI_AA,[2 3]),:,:)=[];FSW_HI_AA(:,~any(FSW_HI_AA,[1 3]),:)=[];FSW_HI_AA(FSW_HI_AA==0)=NaN;
FLW_HI_AA(~any(FLW_HI_AA,[2 3]),:,:)=[];FLW_HI_AA(:,~any(FLW_HI_AA,[1 3]),:)=[];FLW_HI_AA(FLW_HI_AA==0)=NaN;
FSW_HI_AB(~any(FSW_HI_AB,[2 3]),:,:)=[];FSW_HI_AB(:,~any(FSW_HI_AB,[1 3]),:)=[];FSW_HI_AB(FSW_HI_AB==0)=NaN;
FLW_HI_AB(~any(FLW_HI_AB,[2 3]),:,:)=[];FLW_HI_AB(:,~any(FLW_HI_AB,[1 3]),:)=[];FLW_HI_AB(FLW_HI_AB==0)=NaN;

FSW_NH_UN_N0(~any(FSW_NH_UN_N0,[2 3]),:,:)=[];FSW_NH_UN_N0(:,~any(FSW_NH_UN_N0,[1 3]),:)=[];FSW_NH_UN_N0(FSW_NH_UN_N0==0)=NaN;
FLW_NH_UN_N0(~any(FLW_NH_UN_N0,[2 3]),:,:)=[];FLW_NH_UN_N0(:,~any(FLW_NH_UN_N0,[1 3]),:)=[];FLW_NH_UN_N0(FLW_NH_UN_N0==0)=NaN;
FSW_NH_UN_N60(~any(FSW_NH_UN_N60,[2 3]),:,:)=[];FSW_NH_UN_N60(:,~any(FSW_NH_UN_N60,[1 3]),:)=[];FSW_NH_UN_N60(FSW_NH_UN_N60==0)=NaN;
FLW_NH_UN_N60(~any(FLW_NH_UN_N60,[2 3]),:,:)=[];FLW_NH_UN_N60(:,~any(FLW_NH_UN_N60,[1 3]),:)=[];FLW_NH_UN_N60(FLW_NH_UN_N60==0)=NaN;
FSW_NH_UN_N70(~any(FSW_NH_UN_N70,[2 3]),:,:)=[];FSW_NH_UN_N70(:,~any(FSW_NH_UN_N70,[1 3]),:)=[];FSW_NH_UN_N70(FSW_NH_UN_N70==0)=NaN;
FLW_NH_UN_N70(~any(FLW_NH_UN_N70,[2 3]),:,:)=[];FLW_NH_UN_N70(:,~any(FLW_NH_UN_N70,[1 3]),:)=[];FLW_NH_UN_N70(FLW_NH_UN_N70==0)=NaN;
FSW_NH_AA_N0(~any(FSW_NH_AA_N0,[2 3]),:,:)=[];FSW_NH_AA_N0(:,~any(FSW_NH_AA_N0,[1 3]),:)=[];FSW_NH_AA_N0(FSW_NH_AA_N0==0)=NaN;
FLW_NH_AA_N0(~any(FLW_NH_AA_N0,[2 3]),:,:)=[];FLW_NH_AA_N0(:,~any(FLW_NH_AA_N0,[1 3]),:)=[];FLW_NH_AA_N0(FLW_NH_AA_N0==0)=NaN;
FSW_NH_AA_N60(~any(FSW_NH_AA_N60,[2 3]),:,:)=[];FSW_NH_AA_N60(:,~any(FSW_NH_AA_N60,[1 3]),:)=[];FSW_NH_AA_N60(FSW_NH_AA_N60==0)=NaN;
FLW_NH_AA_N60(~any(FLW_NH_AA_N60,[2 3]),:,:)=[];FLW_NH_AA_N60(:,~any(FLW_NH_AA_N60,[1 3]),:)=[];FLW_NH_AA_N60(FLW_NH_AA_N60==0)=NaN;
FSW_NH_AA_N70(~any(FSW_NH_AA_N70,[2 3]),:,:)=[];FSW_NH_AA_N70(:,~any(FSW_NH_AA_N70,[1 3]),:)=[];FSW_NH_AA_N70(FSW_NH_AA_N70==0)=NaN;
FLW_NH_AA_N70(~any(FLW_NH_AA_N70,[2 3]),:,:)=[];FLW_NH_AA_N70(:,~any(FLW_NH_AA_N70,[1 3]),:)=[];FLW_NH_AA_N70(FLW_NH_AA_N70==0)=NaN;
FSW_NH_AB_N0(~any(FSW_NH_AB_N0,[2 3]),:,:)=[];FSW_NH_AB_N0(:,~any(FSW_NH_AB_N0,[1 3]),:)=[];FSW_NH_AB_N0(FSW_NH_AB_N0==0)=NaN;
FLW_NH_AB_N0(~any(FLW_NH_AB_N0,[2 3]),:,:)=[];FLW_NH_AB_N0(:,~any(FLW_NH_AB_N0,[1 3]),:)=[];FLW_NH_AB_N0(FLW_NH_AB_N0==0)=NaN;
FSW_NH_AB_N60(~any(FSW_NH_AB_N60,[2 3]),:,:)=[];FSW_NH_AB_N60(:,~any(FSW_NH_AB_N60,[1 3]),:)=[];FSW_NH_AB_N60(FSW_NH_AB_N60==0)=NaN;
FLW_NH_AB_N60(~any(FLW_NH_AB_N60,[2 3]),:,:)=[];FLW_NH_AB_N60(:,~any(FLW_NH_AB_N60,[1 3]),:)=[];FLW_NH_AB_N60(FLW_NH_AB_N60==0)=NaN;
FSW_NH_AB_N70(~any(FSW_NH_AB_N70,[2 3]),:,:)=[];FSW_NH_AB_N70(:,~any(FSW_NH_AB_N70,[1 3]),:)=[];FSW_NH_AB_N70(FSW_NH_AB_N70==0)=NaN;
FLW_NH_AB_N70(~any(FLW_NH_AB_N70,[2 3]),:,:)=[];FLW_NH_AB_N70(:,~any(FLW_NH_AB_N70,[1 3]),:)=[];FLW_NH_AB_N70(FLW_NH_AB_N70==0)=NaN;

FSW_HI_UN_N0(~any(FSW_HI_UN_N0,[2 3]),:,:)=[];FSW_HI_UN_N0(:,~any(FSW_HI_UN_N0,[1 3]),:)=[];FSW_HI_UN_N0(FSW_HI_UN_N0==0)=NaN;
FLW_HI_UN_N0(~any(FLW_HI_UN_N0,[2 3]),:,:)=[];FLW_HI_UN_N0(:,~any(FLW_HI_UN_N0,[1 3]),:)=[];FLW_HI_UN_N0(FLW_HI_UN_N0==0)=NaN;
FSW_HI_UN_N60(~any(FSW_HI_UN_N60,[2 3]),:,:)=[];FSW_HI_UN_N60(:,~any(FSW_HI_UN_N60,[1 3]),:)=[];FSW_HI_UN_N60(FSW_HI_UN_N60==0)=NaN;
FLW_HI_UN_N60(~any(FLW_HI_UN_N60,[2 3]),:,:)=[];FLW_HI_UN_N60(:,~any(FLW_HI_UN_N60,[1 3]),:)=[];FLW_HI_UN_N60(FLW_HI_UN_N60==0)=NaN;
FSW_HI_UN_N70(~any(FSW_HI_UN_N70,[2 3]),:,:)=[];FSW_HI_UN_N70(:,~any(FSW_HI_UN_N70,[1 3]),:)=[];FSW_HI_UN_N70(FSW_HI_UN_N70==0)=NaN;
FLW_HI_UN_N70(~any(FLW_HI_UN_N70,[2 3]),:,:)=[];FLW_HI_UN_N70(:,~any(FLW_HI_UN_N70,[1 3]),:)=[];FLW_HI_UN_N70(FLW_HI_UN_N70==0)=NaN;
FSW_HI_AA_N0(~any(FSW_HI_AA_N0,[2 3]),:,:)=[];FSW_HI_AA_N0(:,~any(FSW_HI_AA_N0,[1 3]),:)=[];FSW_HI_AA_N0(FSW_HI_AA_N0==0)=NaN;
FLW_HI_AA_N0(~any(FLW_HI_AA_N0,[2 3]),:,:)=[];FLW_HI_AA_N0(:,~any(FLW_HI_AA_N0,[1 3]),:)=[];FLW_HI_AA_N0(FLW_HI_AA_N0==0)=NaN;
FSW_HI_AA_N60(~any(FSW_HI_AA_N60,[2 3]),:,:)=[];FSW_HI_AA_N60(:,~any(FSW_HI_AA_N60,[1 3]),:)=[];FSW_HI_AA_N60(FSW_HI_AA_N60==0)=NaN;
FLW_HI_AA_N60(~any(FLW_HI_AA_N60,[2 3]),:,:)=[];FLW_HI_AA_N60(:,~any(FLW_HI_AA_N60,[1 3]),:)=[];FLW_HI_AA_N60(FLW_HI_AA_N60==0)=NaN;
FSW_HI_AA_N70(~any(FSW_HI_AA_N70,[2 3]),:,:)=[];FSW_HI_AA_N70(:,~any(FSW_HI_AA_N70,[1 3]),:)=[];FSW_HI_AA_N70(FSW_HI_AA_N70==0)=NaN;
FLW_HI_AA_N70(~any(FLW_HI_AA_N70,[2 3]),:,:)=[];FLW_HI_AA_N70(:,~any(FLW_HI_AA_N70,[1 3]),:)=[];FLW_HI_AA_N70(FLW_HI_AA_N70==0)=NaN;
FSW_HI_AB_N0(~any(FSW_HI_AB_N0,[2 3]),:,:)=[];FSW_HI_AB_N0(:,~any(FSW_HI_AB_N0,[1 3]),:)=[];FSW_HI_AB_N0(FSW_HI_AB_N0==0)=NaN;
FLW_HI_AB_N0(~any(FLW_HI_AB_N0,[2 3]),:,:)=[];FLW_HI_AB_N0(:,~any(FLW_HI_AB_N0,[1 3]),:)=[];FLW_HI_AB_N0(FLW_HI_AB_N0==0)=NaN;
FSW_HI_AB_N60(~any(FSW_HI_AB_N60,[2 3]),:,:)=[];FSW_HI_AB_N60(:,~any(FSW_HI_AB_N60,[1 3]),:)=[];FSW_HI_AB_N60(FSW_HI_AB_N60==0)=NaN;
FLW_HI_AB_N60(~any(FLW_HI_AB_N60,[2 3]),:,:)=[];FLW_HI_AB_N60(:,~any(FLW_HI_AB_N60,[1 3]),:)=[];FLW_HI_AB_N60(FLW_HI_AB_N60==0)=NaN;
FSW_HI_AB_N70(~any(FSW_HI_AB_N70,[2 3]),:,:)=[];FSW_HI_AB_N70(:,~any(FSW_HI_AB_N70,[1 3]),:)=[];FSW_HI_AB_N70(FSW_HI_AB_N70==0)=NaN;
FLW_HI_AB_N70(~any(FLW_HI_AB_N70,[2 3]),:,:)=[];FLW_HI_AB_N70(:,~any(FLW_HI_AB_N70,[1 3]),:)=[];FLW_HI_AB_N70(FLW_HI_AB_N70==0)=NaN;

% Calculate global means omitting NaNs, LP filtering with mean-padding at 
% start/end of each group of Speaking/Listening windows
DGSW_Mean = ndnanfilter(reshape(mean(DGSW,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DGSW_B_Mean = ndnanfilter(reshape(mean(DGSW_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DGLW_Mean = ndnanfilter(reshape(mean(DGLW,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DGLW_B_Mean = ndnanfilter(reshape(mean(DGLW_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

DSW_NH_Mean = ndnanfilter(reshape(mean(DSW_NH,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_B_Mean = ndnanfilter(reshape(mean(DSW_NH_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_Mean = ndnanfilter(reshape(mean(DLW_NH,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_B_Mean = ndnanfilter(reshape(mean(DLW_NH_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_Mean = ndnanfilter(reshape(mean(DSW_HI,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_B_Mean = ndnanfilter(reshape(mean(DSW_HI_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_Mean = ndnanfilter(reshape(mean(DLW_HI,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_B_Mean = ndnanfilter(reshape(mean(DLW_HI_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

DSW_N0_Mean = ndnanfilter(reshape(mean(DSW_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_N0_B_Mean = ndnanfilter(reshape(mean(DSW_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_N0_Mean = ndnanfilter(reshape(mean(DLW_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_N0_B_Mean = ndnanfilter(reshape(mean(DLW_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_N60_Mean = ndnanfilter(reshape(mean(DSW_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_N60_B_Mean = ndnanfilter(reshape(mean(DSW_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_N60_Mean = ndnanfilter(reshape(mean(DLW_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_N60_B_Mean = ndnanfilter(reshape(mean(DLW_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_N70_Mean = ndnanfilter(reshape(mean(DSW_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_N70_B_Mean = ndnanfilter(reshape(mean(DSW_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_N70_Mean = ndnanfilter(reshape(mean(DLW_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_N70_B_Mean = ndnanfilter(reshape(mean(DLW_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

DSW_NH_N0_Mean = ndnanfilter(reshape(mean(DSW_NH_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_N0_B_Mean = ndnanfilter(reshape(mean(DSW_NH_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_N0_Mean = ndnanfilter(reshape(mean(DLW_NH_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_N0_B_Mean = ndnanfilter(reshape(mean(DLW_NH_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_N60_Mean = ndnanfilter(reshape(mean(DSW_NH_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_N60_B_Mean = ndnanfilter(reshape(mean(DSW_NH_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_N60_Mean = ndnanfilter(reshape(mean(DLW_NH_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_N60_B_Mean = ndnanfilter(reshape(mean(DLW_NH_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_N70_Mean = ndnanfilter(reshape(mean(DSW_NH_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_N70_B_Mean = ndnanfilter(reshape(mean(DSW_NH_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_N70_Mean = ndnanfilter(reshape(mean(DLW_NH_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_N70_B_Mean = ndnanfilter(reshape(mean(DLW_NH_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

DSW_HI_N0_Mean = ndnanfilter(reshape(mean(DSW_HI_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_N0_B_Mean = ndnanfilter(reshape(mean(DSW_HI_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_N0_Mean = ndnanfilter(reshape(mean(DLW_HI_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_N0_B_Mean = ndnanfilter(reshape(mean(DLW_HI_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_N60_Mean = ndnanfilter(reshape(mean(DSW_HI_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_N60_B_Mean = ndnanfilter(reshape(mean(DSW_HI_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_N60_Mean = ndnanfilter(reshape(mean(DLW_HI_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_N60_B_Mean = ndnanfilter(reshape(mean(DLW_HI_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_N70_Mean = ndnanfilter(reshape(mean(DSW_HI_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_N70_B_Mean = ndnanfilter(reshape(mean(DSW_HI_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_N70_Mean = ndnanfilter(reshape(mean(DLW_HI_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_N70_B_Mean = ndnanfilter(reshape(mean(DLW_HI_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

DSW_HI_UN_Mean = ndnanfilter(reshape(mean(DSW_HI_UN,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_UN_B_Mean = ndnanfilter(reshape(mean(DSW_HI_UN_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_UN_Mean = ndnanfilter(reshape(mean(DLW_HI_UN,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_UN_B_Mean = ndnanfilter(reshape(mean(DLW_HI_UN_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_AA_Mean = ndnanfilter(reshape(mean(DSW_HI_AA,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_AA_B_Mean = ndnanfilter(reshape(mean(DSW_HI_AA_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_AA_Mean = ndnanfilter(reshape(mean(DLW_HI_AA,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_AA_B_Mean = ndnanfilter(reshape(mean(DLW_HI_AA_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_AB_Mean = ndnanfilter(reshape(mean(DSW_HI_AB,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_AB_B_Mean = ndnanfilter(reshape(mean(DSW_HI_AB_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_AB_Mean = ndnanfilter(reshape(mean(DLW_HI_AB,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_AB_B_Mean = ndnanfilter(reshape(mean(DLW_HI_AB_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

DSW_HI_UN_Mean = ndnanfilter(reshape(mean(DSW_HI_UN,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_UN_B_Mean = ndnanfilter(reshape(mean(DSW_HI_UN_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_UN_Mean = ndnanfilter(reshape(mean(DLW_HI_UN,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_UN_B_Mean = ndnanfilter(reshape(mean(DLW_HI_UN_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_AA_Mean = ndnanfilter(reshape(mean(DSW_HI_AA,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_AA_B_Mean = ndnanfilter(reshape(mean(DSW_HI_AA_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_AA_Mean = ndnanfilter(reshape(mean(DLW_HI_AA,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_AA_B_Mean = ndnanfilter(reshape(mean(DLW_HI_AA_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_AB_Mean = ndnanfilter(reshape(mean(DSW_HI_AB,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_AB_B_Mean = ndnanfilter(reshape(mean(DSW_HI_AB_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_AB_Mean = ndnanfilter(reshape(mean(DLW_HI_AB,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_AB_B_Mean = ndnanfilter(reshape(mean(DLW_HI_AB_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

DSW_NH_UN_N0_Mean = ndnanfilter(reshape(mean(DSW_NH_UN_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_UN_N0_B_Mean = ndnanfilter(reshape(mean(DSW_NH_UN_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_UN_N0_Mean = ndnanfilter(reshape(mean(DLW_NH_UN_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_UN_N0_B_Mean = ndnanfilter(reshape(mean(DLW_NH_UN_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_UN_N60_Mean = ndnanfilter(reshape(mean(DSW_NH_UN_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_UN_N60_B_Mean = ndnanfilter(reshape(mean(DSW_NH_UN_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_UN_N60_Mean = ndnanfilter(reshape(mean(DLW_NH_UN_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_UN_N60_B_Mean = ndnanfilter(reshape(mean(DLW_NH_UN_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_UN_N70_Mean = ndnanfilter(reshape(mean(DSW_NH_UN_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_UN_N70_B_Mean = ndnanfilter(reshape(mean(DSW_NH_UN_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_UN_N70_Mean = ndnanfilter(reshape(mean(DLW_NH_UN_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_UN_N70_B_Mean = ndnanfilter(reshape(mean(DLW_NH_UN_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_AA_N0_Mean = ndnanfilter(reshape(mean(DSW_NH_AA_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_AA_N0_B_Mean = ndnanfilter(reshape(mean(DSW_NH_AA_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_AA_N0_Mean = ndnanfilter(reshape(mean(DLW_NH_AA_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_AA_N0_B_Mean = ndnanfilter(reshape(mean(DLW_NH_AA_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_AA_N60_Mean = ndnanfilter(reshape(mean(DSW_NH_AA_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_AA_N60_B_Mean = ndnanfilter(reshape(mean(DSW_NH_AA_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_AA_N60_Mean = ndnanfilter(reshape(mean(DLW_NH_AA_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_AA_N60_B_Mean = ndnanfilter(reshape(mean(DLW_NH_AA_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_AA_N70_Mean = ndnanfilter(reshape(mean(DSW_NH_AA_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_AA_N70_B_Mean = ndnanfilter(reshape(mean(DSW_NH_AA_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_AA_N70_Mean = ndnanfilter(reshape(mean(DLW_NH_AA_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_AA_N70_B_Mean = ndnanfilter(reshape(mean(DLW_NH_AA_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_AB_N0_Mean = ndnanfilter(reshape(mean(DSW_NH_AB_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_AB_N0_B_Mean = ndnanfilter(reshape(mean(DSW_NH_AB_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_AB_N0_Mean = ndnanfilter(reshape(mean(DLW_NH_AB_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_AB_N0_B_Mean = ndnanfilter(reshape(mean(DLW_NH_AB_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_AB_N60_Mean = ndnanfilter(reshape(mean(DSW_NH_AB_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_AB_N60_B_Mean = ndnanfilter(reshape(mean(DSW_NH_AB_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_AB_N60_Mean = ndnanfilter(reshape(mean(DLW_NH_AB_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_AB_N60_B_Mean = ndnanfilter(reshape(mean(DLW_NH_AB_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_AB_N70_Mean = ndnanfilter(reshape(mean(DSW_NH_AB_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_NH_AB_N70_B_Mean = ndnanfilter(reshape(mean(DSW_NH_AB_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_AB_N70_Mean = ndnanfilter(reshape(mean(DLW_NH_AB_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_NH_AB_N70_B_Mean = ndnanfilter(reshape(mean(DLW_NH_AB_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

DSW_HI_UN_N0_Mean = ndnanfilter(reshape(mean(DSW_HI_UN_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_UN_N0_B_Mean = ndnanfilter(reshape(mean(DSW_HI_UN_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_UN_N0_Mean = ndnanfilter(reshape(mean(DLW_HI_UN_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_UN_N0_B_Mean = ndnanfilter(reshape(mean(DLW_HI_UN_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_UN_N60_Mean = ndnanfilter(reshape(mean(DSW_HI_UN_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_UN_N60_B_Mean = ndnanfilter(reshape(mean(DSW_HI_UN_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_UN_N60_Mean = ndnanfilter(reshape(mean(DLW_HI_UN_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_UN_N60_B_Mean = ndnanfilter(reshape(mean(DLW_HI_UN_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_UN_N70_Mean = ndnanfilter(reshape(mean(DSW_HI_UN_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_UN_N70_B_Mean = ndnanfilter(reshape(mean(DSW_HI_UN_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_UN_N70_Mean = ndnanfilter(reshape(mean(DLW_HI_UN_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_UN_N70_B_Mean = ndnanfilter(reshape(mean(DLW_HI_UN_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_AA_N0_Mean = ndnanfilter(reshape(mean(DSW_HI_AA_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_AA_N0_B_Mean = ndnanfilter(reshape(mean(DSW_HI_AA_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_AA_N0_Mean = ndnanfilter(reshape(mean(DLW_HI_AA_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_AA_N0_B_Mean = ndnanfilter(reshape(mean(DLW_HI_AA_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_AA_N60_Mean = ndnanfilter(reshape(mean(DSW_HI_AA_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_AA_N60_B_Mean = ndnanfilter(reshape(mean(DSW_HI_AA_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_AA_N60_Mean = ndnanfilter(reshape(mean(DLW_HI_AA_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_AA_N60_B_Mean = ndnanfilter(reshape(mean(DLW_HI_AA_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_AA_N70_Mean = ndnanfilter(reshape(mean(DSW_HI_AA_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_AA_N70_B_Mean = ndnanfilter(reshape(mean(DSW_HI_AA_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_AA_N70_Mean = ndnanfilter(reshape(mean(DLW_HI_AA_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_AA_N70_B_Mean = ndnanfilter(reshape(mean(DLW_HI_AA_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_AB_N0_Mean = ndnanfilter(reshape(mean(DSW_HI_AB_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_AB_N0_B_Mean = ndnanfilter(reshape(mean(DSW_HI_AB_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_AB_N0_Mean = ndnanfilter(reshape(mean(DLW_HI_AB_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_AB_N0_B_Mean = ndnanfilter(reshape(mean(DLW_HI_AB_N0_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_AB_N60_Mean = ndnanfilter(reshape(mean(DSW_HI_AB_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_AB_N60_B_Mean = ndnanfilter(reshape(mean(DSW_HI_AB_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_AB_N60_Mean = ndnanfilter(reshape(mean(DLW_HI_AB_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_AB_N60_B_Mean = ndnanfilter(reshape(mean(DLW_HI_AB_N60_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_AB_N70_Mean = ndnanfilter(reshape(mean(DSW_HI_AB_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DSW_HI_AB_N70_B_Mean = ndnanfilter(reshape(mean(DSW_HI_AB_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_AB_N70_Mean = ndnanfilter(reshape(mean(DLW_HI_AB_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
DLW_HI_AB_N70_B_Mean = ndnanfilter(reshape(mean(DLW_HI_AB_N70_B,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

FGSW_Mean = ndnanfilter(reshape(mean(FGSW,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FGLW_Mean = ndnanfilter(reshape(mean(FGLW,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

FSW_NH_Mean = ndnanfilter(reshape(mean(FSW_NH,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_NH_Mean = ndnanfilter(reshape(mean(FLW_NH,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_HI_Mean = ndnanfilter(reshape(mean(FSW_HI,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_HI_Mean = ndnanfilter(reshape(mean(FLW_HI,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

FSW_N0_Mean = ndnanfilter(reshape(mean(FSW_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_N0_Mean = ndnanfilter(reshape(mean(FLW_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_N60_Mean = ndnanfilter(reshape(mean(FSW_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_N60_Mean = ndnanfilter(reshape(mean(FLW_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_N70_Mean = ndnanfilter(reshape(mean(FSW_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_N70_Mean = ndnanfilter(reshape(mean(FLW_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

FSW_NH_N0_Mean = ndnanfilter(reshape(mean(FSW_NH_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_NH_N0_Mean = ndnanfilter(reshape(mean(FLW_NH_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_NH_N60_Mean = ndnanfilter(reshape(mean(FSW_NH_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_NH_N60_Mean = ndnanfilter(reshape(mean(FLW_NH_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_NH_N70_Mean = ndnanfilter(reshape(mean(FSW_NH_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_NH_N70_Mean = ndnanfilter(reshape(mean(FLW_NH_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

FSW_HI_N0_Mean = ndnanfilter(reshape(mean(FSW_HI_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_HI_N0_Mean = ndnanfilter(reshape(mean(FLW_HI_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_HI_N60_Mean = ndnanfilter(reshape(mean(FSW_HI_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_HI_N60_Mean = ndnanfilter(reshape(mean(FLW_HI_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_HI_N70_Mean = ndnanfilter(reshape(mean(FSW_HI_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_HI_N70_Mean = ndnanfilter(reshape(mean(FLW_HI_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

FSW_NH_UN_Mean = ndnanfilter(reshape(mean(FSW_NH_UN,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_NH_UN_Mean = ndnanfilter(reshape(mean(FLW_NH_UN,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_NH_AA_Mean = ndnanfilter(reshape(mean(FSW_NH_AA,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_NH_AA_Mean = ndnanfilter(reshape(mean(FLW_NH_AA,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_NH_AB_Mean = ndnanfilter(reshape(mean(FSW_NH_AB,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_NH_AB_Mean = ndnanfilter(reshape(mean(FLW_NH_AB,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

FSW_HI_UN_Mean = ndnanfilter(reshape(mean(FSW_HI_UN,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_HI_UN_Mean = ndnanfilter(reshape(mean(FLW_HI_UN,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_HI_AA_Mean = ndnanfilter(reshape(mean(FSW_HI_AA,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_HI_AA_Mean = ndnanfilter(reshape(mean(FLW_HI_AA,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_HI_AB_Mean = ndnanfilter(reshape(mean(FSW_HI_AB,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_HI_AB_Mean = ndnanfilter(reshape(mean(FLW_HI_AB,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

FSW_NH_UN_N0_Mean = ndnanfilter(reshape(mean(FSW_NH_UN_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_NH_UN_N0_Mean = ndnanfilter(reshape(mean(FLW_NH_UN_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_NH_UN_N60_Mean = ndnanfilter(reshape(mean(FSW_NH_UN_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_NH_UN_N60_Mean = ndnanfilter(reshape(mean(FLW_NH_UN_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_NH_UN_N70_Mean = ndnanfilter(reshape(mean(FSW_NH_UN_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_NH_UN_N70_Mean = ndnanfilter(reshape(mean(FLW_NH_UN_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_NH_AA_N0_Mean = ndnanfilter(reshape(mean(FSW_NH_AA_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_NH_AA_N0_Mean = ndnanfilter(reshape(mean(FLW_NH_AA_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_NH_AA_N60_Mean = ndnanfilter(reshape(mean(FSW_NH_AA_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_NH_AA_N60_Mean = ndnanfilter(reshape(mean(FLW_NH_AA_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_NH_AA_N70_Mean = ndnanfilter(reshape(mean(FSW_NH_AA_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_NH_AA_N70_Mean = ndnanfilter(reshape(mean(FLW_NH_AA_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_NH_AB_N0_Mean = ndnanfilter(reshape(mean(FSW_NH_AB_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_NH_AB_N0_Mean = ndnanfilter(reshape(mean(FLW_NH_AB_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_NH_AB_N60_Mean = ndnanfilter(reshape(mean(FSW_NH_AB_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_NH_AB_N60_Mean = ndnanfilter(reshape(mean(FLW_NH_AB_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_NH_AB_N70_Mean = ndnanfilter(reshape(mean(FSW_NH_AB_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_NH_AB_N70_Mean = ndnanfilter(reshape(mean(FLW_NH_AB_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

FSW_HI_UN_N0_Mean = ndnanfilter(reshape(mean(FSW_HI_UN_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_HI_UN_N0_Mean = ndnanfilter(reshape(mean(FLW_HI_UN_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_HI_UN_N60_Mean = ndnanfilter(reshape(mean(FSW_HI_UN_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_HI_UN_N60_Mean = ndnanfilter(reshape(mean(FLW_HI_UN_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_HI_UN_N70_Mean = ndnanfilter(reshape(mean(FSW_HI_UN_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_HI_UN_N70_Mean = ndnanfilter(reshape(mean(FLW_HI_UN_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_HI_AA_N0_Mean = ndnanfilter(reshape(mean(FSW_HI_AA_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_HI_AA_N0_Mean = ndnanfilter(reshape(mean(FLW_HI_AA_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_HI_AA_N60_Mean = ndnanfilter(reshape(mean(FSW_HI_AA_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_HI_AA_N60_Mean = ndnanfilter(reshape(mean(FLW_HI_AA_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_HI_AA_N70_Mean = ndnanfilter(reshape(mean(FSW_HI_AA_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_HI_AA_N70_Mean = ndnanfilter(reshape(mean(FLW_HI_AA_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_HI_AB_N0_Mean = ndnanfilter(reshape(mean(FSW_HI_AB_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_HI_AB_N0_Mean = ndnanfilter(reshape(mean(FLW_HI_AB_N0,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_HI_AB_N60_Mean = ndnanfilter(reshape(mean(FSW_HI_AB_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_HI_AB_N60_Mean = ndnanfilter(reshape(mean(FLW_HI_AB_N60,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FSW_HI_AB_N70_Mean = ndnanfilter(reshape(mean(FSW_HI_AB_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);
FLW_HI_AB_N70_Mean = ndnanfilter(reshape(mean(FLW_HI_AB_N70,[1 2],'omitnan'),[],1)','hamming',FilterWidth);

% Calculate SEM as: std(X)/sqrt(squeeze(sum(~isnan(X),[1 2]))))
DGSW_SEM = (squeeze(std(DGSW,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DGSW),[1 2]))))';
DGSW_B_SEM = (squeeze(std(DGSW_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DGSW_B),[1 2]))))';
DGLW_SEM = (squeeze(std(DGLW,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DGLW),[1 2]))))';
DGLW_B_SEM = (squeeze(std(DGLW_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DGLW_B),[1 2]))))';

DSW_NH_SEM = (squeeze(std(DSW_NH,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH),[1 2]))))';
DSW_NH_B_SEM = (squeeze(std(DSW_NH_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_B),[1 2]))))';
DLW_NH_SEM = (squeeze(std(DLW_NH,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH),[1 2]))))';
DLW_NH_B_SEM = (squeeze(std(DLW_NH_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_B),[1 2]))))';
DSW_HI_SEM = (squeeze(std(DSW_HI,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI),[1 2]))))';
DSW_HI_B_SEM = (squeeze(std(DSW_HI_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_B),[1 2]))))';
DLW_HI_SEM = (squeeze(std(DLW_HI,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI),[1 2]))))';
DLW_HI_B_SEM = (squeeze(std(DLW_HI_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_B),[1 2]))))';

DSW_N0_SEM = (squeeze(std(DSW_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_N0),[1 2]))))';
DSW_N0_B_SEM = (squeeze(std(DSW_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_N0_B),[1 2]))))';
DLW_N0_SEM = (squeeze(std(DLW_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_N0),[1 2]))))';
DLW_N0_B_SEM = (squeeze(std(DLW_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_N0_B),[1 2]))))';
DSW_N60_SEM = (squeeze(std(DSW_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_N60),[1 2]))))';
DSW_N60_B_SEM = (squeeze(std(DSW_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_N60_B),[1 2]))))';
DLW_N60_SEM = (squeeze(std(DLW_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_N60),[1 2]))))';
DLW_N60_B_SEM = (squeeze(std(DLW_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_N60_B),[1 2]))))';
DSW_N70_SEM = (squeeze(std(DSW_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_N70),[1 2]))))';
DSW_N70_B_SEM = (squeeze(std(DSW_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_N70_B),[1 2]))))';
DLW_N70_SEM = (squeeze(std(DLW_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_N70),[1 2]))))';
DLW_N70_B_SEM = (squeeze(std(DLW_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_N70_B),[1 2]))))';

DSW_NH_N0_SEM = (squeeze(std(DSW_NH_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_N0),[1 2]))))';
DSW_NH_N0_B_SEM = (squeeze(std(DSW_NH_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_N0_B),[1 2]))))';
DLW_NH_N0_SEM = (squeeze(std(DLW_NH_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_N0),[1 2]))))';
DLW_NH_N0_B_SEM = (squeeze(std(DLW_NH_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_N0_B),[1 2]))))';
DSW_NH_N60_SEM = (squeeze(std(DSW_NH_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_N60),[1 2]))))';
DSW_NH_N60_B_SEM = (squeeze(std(DSW_NH_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_N60_B),[1 2]))))';
DLW_NH_N60_SEM = (squeeze(std(DLW_NH_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_N60),[1 2]))))';
DLW_NH_N60_B_SEM = (squeeze(std(DLW_NH_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_N60_B),[1 2]))))';
DSW_NH_N70_SEM = (squeeze(std(DSW_NH_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_N70),[1 2]))))';
DSW_NH_N70_B_SEM = (squeeze(std(DSW_NH_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_N70_B),[1 2]))))';
DLW_NH_N70_SEM = (squeeze(std(DLW_NH_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_N70),[1 2]))))';
DLW_NH_N70_B_SEM = (squeeze(std(DLW_NH_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_N70_B),[1 2]))))';

DSW_HI_N0_SEM = (squeeze(std(DSW_HI_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_N0),[1 2]))))';
DSW_HI_N0_B_SEM = (squeeze(std(DSW_HI_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_N0_B),[1 2]))))';
DLW_HI_N0_SEM = (squeeze(std(DLW_HI_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_N0),[1 2]))))';
DLW_HI_N0_B_SEM = (squeeze(std(DLW_HI_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_N0_B),[1 2]))))';
DSW_HI_N60_SEM = (squeeze(std(DSW_HI_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_N60),[1 2]))))';
DSW_HI_N60_B_SEM = (squeeze(std(DSW_HI_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_N60_B),[1 2]))))';
DLW_HI_N60_SEM = (squeeze(std(DLW_HI_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_N60),[1 2]))))';
DLW_HI_N60_B_SEM = (squeeze(std(DLW_HI_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_N60_B),[1 2]))))';
DSW_HI_N70_SEM = (squeeze(std(DSW_HI_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_N70),[1 2]))))';
DSW_HI_N70_B_SEM = (squeeze(std(DSW_HI_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_N70_B),[1 2]))))';
DLW_HI_N70_SEM = (squeeze(std(DLW_HI_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_N70),[1 2]))))';
DLW_HI_N70_B_SEM = (squeeze(std(DLW_HI_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_N70_B),[1 2]))))';

DSW_NH_UN_SEM = (squeeze(std(DSW_NH_UN,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_UN),[1 2]))))';
DSW_NH_UN_B_SEM = (squeeze(std(DSW_NH_UN_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_UN_B),[1 2]))))';
DLW_NH_UN_SEM = (squeeze(std(DLW_NH_UN,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_UN),[1 2]))))';
DLW_NH_UN_B_SEM = (squeeze(std(DLW_NH_UN_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_UN_B),[1 2]))))';
DSW_NH_AA_SEM = (squeeze(std(DSW_NH_AA,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_AA),[1 2]))))';
DSW_NH_AA_B_SEM = (squeeze(std(DSW_NH_AA_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_AA_B),[1 2]))))';
DLW_NH_AA_SEM = (squeeze(std(DLW_NH_AA,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_AA),[1 2]))))';
DLW_NH_AA_B_SEM = (squeeze(std(DLW_NH_AA_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_AA_B),[1 2]))))';
DSW_NH_AB_SEM = (squeeze(std(DSW_NH_AB,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_AB),[1 2]))))';
DSW_NH_AB_B_SEM = (squeeze(std(DSW_NH_AB_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_AB_B),[1 2]))))';
DLW_NH_AB_SEM = (squeeze(std(DLW_NH_AB,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_AB),[1 2]))))';
DLW_NH_AB_B_SEM = (squeeze(std(DLW_NH_AB_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_AB_B),[1 2]))))';

DSW_HI_UN_SEM = (squeeze(std(DSW_HI_UN,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_UN),[1 2]))))';
DSW_HI_UN_B_SEM = (squeeze(std(DSW_HI_UN_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_UN_B),[1 2]))))';
DLW_HI_UN_SEM = (squeeze(std(DLW_HI_UN,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_UN),[1 2]))))';
DLW_HI_UN_B_SEM = (squeeze(std(DLW_HI_UN_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_UN_B),[1 2]))))';
DSW_HI_AA_SEM = (squeeze(std(DSW_HI_AA,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_AA),[1 2]))))';
DSW_HI_AA_B_SEM = (squeeze(std(DSW_HI_AA_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_AA_B),[1 2]))))';
DLW_HI_AA_SEM = (squeeze(std(DLW_HI_AA,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_AA),[1 2]))))';
DLW_HI_AA_B_SEM = (squeeze(std(DLW_HI_AA_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_AA_B),[1 2]))))';
DSW_HI_AB_SEM = (squeeze(std(DSW_HI_AB,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_AB),[1 2]))))';
DSW_HI_AB_B_SEM = (squeeze(std(DSW_HI_AB_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_AB_B),[1 2]))))';
DLW_HI_AB_SEM = (squeeze(std(DLW_HI_AB,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_AB),[1 2]))))';
DLW_HI_AB_B_SEM = (squeeze(std(DLW_HI_AB_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_AB_B),[1 2]))))';

DSW_NH_UN_N0_SEM = (squeeze(std(DSW_NH_UN_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_UN_N0),[1 2]))))';
DSW_NH_UN_N0_B_SEM = (squeeze(std(DSW_NH_UN_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_UN_N0_B),[1 2]))))';
DLW_NH_UN_N0_SEM = (squeeze(std(DLW_NH_UN_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_UN_N0),[1 2]))))';
DLW_NH_UN_N0_B_SEM = (squeeze(std(DLW_NH_UN_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_UN_N0_B),[1 2]))))';
DSW_NH_UN_N60_SEM = (squeeze(std(DSW_NH_UN_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_UN_N60),[1 2]))))';
DSW_NH_UN_N60_B_SEM = (squeeze(std(DSW_NH_UN_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_UN_N60_B),[1 2]))))';
DLW_NH_UN_N60_SEM = (squeeze(std(DLW_NH_UN_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_UN_N60),[1 2]))))';
DLW_NH_UN_N60_B_SEM = (squeeze(std(DLW_NH_UN_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_UN_N60_B),[1 2]))))';
DSW_NH_UN_N70_SEM = (squeeze(std(DSW_NH_UN_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_UN_N70),[1 2]))))';
DSW_NH_UN_N70_B_SEM = (squeeze(std(DSW_NH_UN_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_UN_N70_B),[1 2]))))';
DLW_NH_UN_N70_SEM = (squeeze(std(DLW_NH_UN_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_UN_N70),[1 2]))))';
DLW_NH_UN_N70_B_SEM = (squeeze(std(DLW_NH_UN_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_UN_N70_B),[1 2]))))';
DSW_NH_AA_N0_SEM = (squeeze(std(DSW_NH_AA_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_AA_N0),[1 2]))))';
DSW_NH_AA_N0_B_SEM = (squeeze(std(DSW_NH_AA_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_AA_N0_B),[1 2]))))';
DLW_NH_AA_N0_SEM = (squeeze(std(DLW_NH_AA_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_AA_N0),[1 2]))))';
DLW_NH_AA_N0_B_SEM = (squeeze(std(DLW_NH_AA_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_AA_N0_B),[1 2]))))';
DSW_NH_AA_N60_SEM = (squeeze(std(DSW_NH_AA_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_AA_N60),[1 2]))))';
DSW_NH_AA_N60_B_SEM = (squeeze(std(DSW_NH_AA_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_AA_N60_B),[1 2]))))';
DLW_NH_AA_N60_SEM = (squeeze(std(DLW_NH_AA_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_AA_N60),[1 2]))))';
DLW_NH_AA_N60_B_SEM = (squeeze(std(DLW_NH_AA_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_AA_N60_B),[1 2]))))';
DSW_NH_AA_N70_SEM = (squeeze(std(DSW_NH_AA_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_AA_N70),[1 2]))))';
DSW_NH_AA_N70_B_SEM = (squeeze(std(DSW_NH_AA_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_AA_N70_B),[1 2]))))';
DLW_NH_AA_N70_SEM = (squeeze(std(DLW_NH_AA_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_AA_N70),[1 2]))))';
DLW_NH_AA_N70_B_SEM = (squeeze(std(DLW_NH_AA_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_AA_N70_B),[1 2]))))';
DSW_NH_AB_N0_SEM = (squeeze(std(DSW_NH_AB_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_AB_N0),[1 2]))))';
DSW_NH_AB_N0_B_SEM = (squeeze(std(DSW_NH_AB_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_AB_N0_B),[1 2]))))';
DLW_NH_AB_N0_SEM = (squeeze(std(DLW_NH_AB_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_AB_N0),[1 2]))))';
DLW_NH_AB_N0_B_SEM = (squeeze(std(DLW_NH_AB_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_AB_N0_B),[1 2]))))';
DSW_NH_AB_N60_SEM = (squeeze(std(DSW_NH_AB_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_AB_N60),[1 2]))))';
DSW_NH_AB_N60_B_SEM = (squeeze(std(DSW_NH_AB_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_AB_N60_B),[1 2]))))';
DLW_NH_AB_N60_SEM = (squeeze(std(DLW_NH_AB_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_AB_N60),[1 2]))))';
DLW_NH_AB_N60_B_SEM = (squeeze(std(DLW_NH_AB_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_AB_N60_B),[1 2]))))';
DSW_NH_AB_N70_SEM = (squeeze(std(DSW_NH_AB_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_AB_N70),[1 2]))))';
DSW_NH_AB_N70_B_SEM = (squeeze(std(DSW_NH_AB_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_NH_AB_N70_B),[1 2]))))';
DLW_NH_AB_N70_SEM = (squeeze(std(DLW_NH_AB_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_AB_N70),[1 2]))))';
DLW_NH_AB_N70_B_SEM = (squeeze(std(DLW_NH_AB_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_NH_AB_N70_B),[1 2]))))';

DSW_HI_UN_N0_SEM = (squeeze(std(DSW_HI_UN_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_UN_N0),[1 2]))))';
DSW_HI_UN_N0_B_SEM = (squeeze(std(DSW_HI_UN_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_UN_N0_B),[1 2]))))';
DLW_HI_UN_N0_SEM = (squeeze(std(DLW_HI_UN_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_UN_N0),[1 2]))))';
DLW_HI_UN_N0_B_SEM = (squeeze(std(DLW_HI_UN_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_UN_N0_B),[1 2]))))';
DSW_HI_UN_N60_SEM = (squeeze(std(DSW_HI_UN_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_UN_N60),[1 2]))))';
DSW_HI_UN_N60_B_SEM = (squeeze(std(DSW_HI_UN_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_UN_N60_B),[1 2]))))';
DLW_HI_UN_N60_SEM = (squeeze(std(DLW_HI_UN_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_UN_N60),[1 2]))))';
DLW_HI_UN_N60_B_SEM = (squeeze(std(DLW_HI_UN_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_UN_N60_B),[1 2]))))';
DSW_HI_UN_N70_SEM = (squeeze(std(DSW_HI_UN_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_UN_N70),[1 2]))))';
DSW_HI_UN_N70_B_SEM = (squeeze(std(DSW_HI_UN_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_UN_N70_B),[1 2]))))';
DLW_HI_UN_N70_SEM = (squeeze(std(DLW_HI_UN_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_UN_N70),[1 2]))))';
DLW_HI_UN_N70_B_SEM = (squeeze(std(DLW_HI_UN_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_UN_N70_B),[1 2]))))';
DSW_HI_AA_N0_SEM = (squeeze(std(DSW_HI_AA_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_AA_N0),[1 2]))))';
DSW_HI_AA_N0_B_SEM = (squeeze(std(DSW_HI_AA_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_AA_N0_B),[1 2]))))';
DLW_HI_AA_N0_SEM = (squeeze(std(DLW_HI_AA_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_AA_N0),[1 2]))))';
DLW_HI_AA_N0_B_SEM = (squeeze(std(DLW_HI_AA_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_AA_N0_B),[1 2]))))';
DSW_HI_AA_N60_SEM = (squeeze(std(DSW_HI_AA_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_AA_N60),[1 2]))))';
DSW_HI_AA_N60_B_SEM = (squeeze(std(DSW_HI_AA_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_AA_N60_B),[1 2]))))';
DLW_HI_AA_N60_SEM = (squeeze(std(DLW_HI_AA_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_AA_N60),[1 2]))))';
DLW_HI_AA_N60_B_SEM = (squeeze(std(DLW_HI_AA_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_AA_N60_B),[1 2]))))';
DSW_HI_AA_N70_SEM = (squeeze(std(DSW_HI_AA_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_AA_N70),[1 2]))))';
DSW_HI_AA_N70_B_SEM = (squeeze(std(DSW_HI_AA_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_AA_N70_B),[1 2]))))';
DLW_HI_AA_N70_SEM = (squeeze(std(DLW_HI_AA_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_AA_N70),[1 2]))))';
DLW_HI_AA_N70_B_SEM = (squeeze(std(DLW_HI_AA_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_AA_N70_B),[1 2]))))';
DSW_HI_AB_N0_SEM = (squeeze(std(DSW_HI_AB_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_AB_N0),[1 2]))))';
DSW_HI_AB_N0_B_SEM = (squeeze(std(DSW_HI_AB_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_AB_N0_B),[1 2]))))';
DLW_HI_AB_N0_SEM = (squeeze(std(DLW_HI_AB_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_AB_N0),[1 2]))))';
DLW_HI_AB_N0_B_SEM = (squeeze(std(DLW_HI_AB_N0_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_AB_N0_B),[1 2]))))';
DSW_HI_AB_N60_SEM = (squeeze(std(DSW_HI_AB_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_AB_N60),[1 2]))))';
DSW_HI_AB_N60_B_SEM = (squeeze(std(DSW_HI_AB_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_AB_N60_B),[1 2]))))';
DLW_HI_AB_N60_SEM = (squeeze(std(DLW_HI_AB_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_AB_N60),[1 2]))))';
DLW_HI_AB_N60_B_SEM = (squeeze(std(DLW_HI_AB_N60_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_AB_N60_B),[1 2]))))';
DSW_HI_AB_N70_SEM = (squeeze(std(DSW_HI_AB_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_AB_N70),[1 2]))))';
DSW_HI_AB_N70_B_SEM = (squeeze(std(DSW_HI_AB_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DSW_HI_AB_N70_B),[1 2]))))';
DLW_HI_AB_N70_SEM = (squeeze(std(DLW_HI_AB_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_AB_N70),[1 2]))))';
DLW_HI_AB_N70_B_SEM = (squeeze(std(DLW_HI_AB_N70_B,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(DLW_HI_AB_N70_B),[1 2]))))';

FGSW_SEM = (squeeze(std(FGSW,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FGSW),[1 2]))))';
FGLW_SEM = (squeeze(std(FGLW,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FGLW),[1 2]))))';

FSW_NH_SEM = (squeeze(std(FSW_NH,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_NH),[1 2]))))';
FLW_NH_SEM = (squeeze(std(FLW_NH,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_NH),[1 2]))))';
FSW_HI_SEM = (squeeze(std(FSW_HI,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_HI),[1 2]))))';
FLW_HI_SEM = (squeeze(std(FLW_HI,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_HI),[1 2]))))';

FSW_N0_SEM = (squeeze(std(FSW_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_N0),[1 2]))))';
FLW_N0_SEM = (squeeze(std(FLW_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_N0),[1 2]))))';
FSW_N60_SEM = (squeeze(std(FSW_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_N60),[1 2]))))';
FLW_N60_SEM = (squeeze(std(FLW_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_N60),[1 2]))))';
FSW_N70_SEM = (squeeze(std(FSW_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_N70),[1 2]))))';
FLW_N70_SEM = (squeeze(std(FLW_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_N70),[1 2]))))';

FSW_NH_N0_SEM = (squeeze(std(FSW_NH_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_NH_N0),[1 2]))))';
FLW_NH_N0_SEM = (squeeze(std(FLW_NH_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_NH_N0),[1 2]))))';
FSW_NH_N60_SEM = (squeeze(std(FSW_NH_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_NH_N60),[1 2]))))';
FLW_NH_N60_SEM = (squeeze(std(FLW_NH_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_NH_N60),[1 2]))))';
FSW_NH_N70_SEM = (squeeze(std(FSW_NH_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_NH_N70),[1 2]))))';
FLW_NH_N70_SEM = (squeeze(std(FLW_NH_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_NH_N70),[1 2]))))';

FSW_HI_N0_SEM = (squeeze(std(FSW_HI_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_HI_N0),[1 2]))))';
FLW_HI_N0_SEM = (squeeze(std(FLW_HI_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_HI_N0),[1 2]))))';
FSW_HI_N60_SEM = (squeeze(std(FSW_HI_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_HI_N60),[1 2]))))';
FLW_HI_N60_SEM = (squeeze(std(FLW_HI_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_HI_N60),[1 2]))))';
FSW_HI_N70_SEM = (squeeze(std(FSW_HI_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_HI_N70),[1 2]))))';
FLW_HI_N70_SEM = (squeeze(std(FLW_HI_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_HI_N70),[1 2]))))';

FSW_NH_UN_SEM = (squeeze(std(FSW_NH_UN,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_NH_UN),[1 2]))))';
FLW_NH_UN_SEM = (squeeze(std(FLW_NH_UN,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_NH_UN),[1 2]))))';
FSW_NH_AA_SEM = (squeeze(std(FSW_NH_AA,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_NH_AA),[1 2]))))';
FLW_NH_AA_SEM = (squeeze(std(FLW_NH_AA,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_NH_AA),[1 2]))))';
FSW_NH_AB_SEM = (squeeze(std(FSW_NH_AB,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_NH_AB),[1 2]))))';
FLW_NH_AB_SEM = (squeeze(std(FLW_NH_AB,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_NH_AB),[1 2]))))';

FSW_HI_UN_SEM = (squeeze(std(FSW_HI_UN,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_HI_UN),[1 2]))))';
FLW_HI_UN_SEM = (squeeze(std(FLW_HI_UN,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_HI_UN),[1 2]))))';
FSW_HI_AA_SEM = (squeeze(std(FSW_HI_AA,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_HI_AA),[1 2]))))';
FLW_HI_AA_SEM = (squeeze(std(FLW_HI_AA,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_HI_AA),[1 2]))))';
FSW_HI_AB_SEM = (squeeze(std(FSW_HI_AB,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_HI_AB),[1 2]))))';
FLW_HI_AB_SEM = (squeeze(std(FLW_HI_AB,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_HI_AB),[1 2]))))';

FSW_NH_UN_N0_SEM = (squeeze(std(FSW_NH_UN_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_NH_UN_N0),[1 2]))))';
FLW_NH_UN_N0_SEM = (squeeze(std(FLW_NH_UN_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_NH_UN_N0),[1 2]))))';
FSW_NH_UN_N60_SEM = (squeeze(std(FSW_NH_UN_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_NH_UN_N60),[1 2]))))';
FLW_NH_UN_N60_SEM = (squeeze(std(FLW_NH_UN_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_NH_UN_N60),[1 2]))))';
FSW_NH_UN_N70_SEM = (squeeze(std(FSW_NH_UN_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_NH_UN_N70),[1 2]))))';
FLW_NH_UN_N70_SEM = (squeeze(std(FLW_NH_UN_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_NH_UN_N70),[1 2]))))';
FSW_NH_AA_N0_SEM = (squeeze(std(FSW_NH_AA_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_NH_AA_N0),[1 2]))))';
FLW_NH_AA_N0_SEM = (squeeze(std(FLW_NH_AA_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_NH_AA_N0),[1 2]))))';
FSW_NH_AA_N60_SEM = (squeeze(std(FSW_NH_AA_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_NH_AA_N60),[1 2]))))';
FLW_NH_AA_N60_SEM = (squeeze(std(FLW_NH_AA_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_NH_AA_N60),[1 2]))))';
FSW_NH_AA_N70_SEM = (squeeze(std(FSW_NH_AA_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_NH_AA_N70),[1 2]))))';
FLW_NH_AA_N70_SEM = (squeeze(std(FLW_NH_AA_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_NH_AA_N70),[1 2]))))';
FSW_NH_AB_N0_SEM = (squeeze(std(FSW_NH_AB_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_NH_AB_N0),[1 2]))))';
FLW_NH_AB_N0_SEM = (squeeze(std(FLW_NH_AB_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_NH_AB_N0),[1 2]))))';
FSW_NH_AB_N60_SEM = (squeeze(std(FSW_NH_AB_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_NH_AB_N60),[1 2]))))';
FLW_NH_AB_N60_SEM = (squeeze(std(FLW_NH_AB_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_NH_AB_N60),[1 2]))))';
FSW_NH_AB_N70_SEM = (squeeze(std(FSW_NH_AB_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_NH_AB_N70),[1 2]))))';
FLW_NH_AB_N70_SEM = (squeeze(std(FLW_NH_AB_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_NH_AB_N70),[1 2]))))';

FSW_HI_UN_N0_SEM = (squeeze(std(FSW_HI_UN_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_HI_UN_N0),[1 2]))))';
FLW_HI_UN_N0_SEM = (squeeze(std(FLW_HI_UN_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_HI_UN_N0),[1 2]))))';
FSW_HI_UN_N60_SEM = (squeeze(std(FSW_HI_UN_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_HI_UN_N60),[1 2]))))';
FLW_HI_UN_N60_SEM = (squeeze(std(FLW_HI_UN_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_HI_UN_N60),[1 2]))))';
FSW_HI_UN_N70_SEM = (squeeze(std(FSW_HI_UN_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_HI_UN_N70),[1 2]))))';
FLW_HI_UN_N70_SEM = (squeeze(std(FLW_HI_UN_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_HI_UN_N70),[1 2]))))';
FSW_HI_AA_N0_SEM = (squeeze(std(FSW_HI_AA_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_HI_AA_N0),[1 2]))))';
FLW_HI_AA_N0_SEM = (squeeze(std(FLW_HI_AA_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_HI_AA_N0),[1 2]))))';
FSW_HI_AA_N60_SEM = (squeeze(std(FSW_HI_AA_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_HI_AA_N60),[1 2]))))';
FLW_HI_AA_N60_SEM = (squeeze(std(FLW_HI_AA_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_HI_AA_N60),[1 2]))))';
FSW_HI_AA_N70_SEM = (squeeze(std(FSW_HI_AA_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_HI_AA_N70),[1 2]))))';
FLW_HI_AA_N70_SEM = (squeeze(std(FLW_HI_AA_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_HI_AA_N70),[1 2]))))';
FSW_HI_AB_N0_SEM = (squeeze(std(FSW_HI_AB_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_HI_AB_N0),[1 2]))))';
FLW_HI_AB_N0_SEM = (squeeze(std(FLW_HI_AB_N0,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_HI_AB_N0),[1 2]))))';
FSW_HI_AB_N60_SEM = (squeeze(std(FSW_HI_AB_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_HI_AB_N60),[1 2]))))';
FLW_HI_AB_N60_SEM = (squeeze(std(FLW_HI_AB_N60,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_HI_AB_N60),[1 2]))))';
FSW_HI_AB_N70_SEM = (squeeze(std(FSW_HI_AB_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FSW_HI_AB_N70),[1 2]))))';
FLW_HI_AB_N70_SEM = (squeeze(std(FLW_HI_AB_N70,0,[1 2],'omitnan'))./sqrt(squeeze(sum(~isnan(FLW_HI_AB_N70),[1 2]))))';
%% Calculate Grand Average over the mean of the participants
TPsG = nonzeros(TPsOrder_II);
TPsNH = nonzeros(TPsOrder_II(1:2:NTPs_II,:));
TPsHI = nonzeros(TPsOrder_II(2:2:NTPs_II,:));
TPsUN = nonzeros(TPsOrder_II(:,1:3));
TPsAA = nonzeros(TPsOrder_II(:,4:6));
TPsAB = nonzeros(TPsOrder_II(:,7:9));
TPsN0 = nonzeros(TPsOrder_II(:,1:3:NCond_II));
TPsN60 = nonzeros(TPsOrder_II(:,2:3:NCond_II));
TPsN70 = nonzeros(TPsOrder_II(:,3:3:NCond_II));

TP_DGSW = zeros(NTPs_II,NCols); TP_DGLW = TP_DGSW; % Global TP_Diameter Speaking/Listening Windows
TP_DSW_NH = TP_DGSW; TP_DLW_NH = TP_DGSW; % NH TP_Diameter Speaking/Listening Windows
TP_DSW_HI = TP_DGSW; TP_DLW_HI = TP_DGSW; % HI TP_Diameter Speaking/Listening Windows
TP_DSW_N0 = TP_DGSW; TP_DLW_N0 = TP_DGSW; % N0 TP_Diameter Speaking/Listening Windows
TP_DSW_N60 = TP_DGSW; TP_DLW_N60 = TP_DGSW; % N60 TP_Diameter Speaking/Listening Windows
TP_DSW_N70 = TP_DGSW; TP_DLW_N70 = TP_DGSW; % N70 TP_Diameter Speaking/Listening Windows
TP_DSW_NH_N0 = TP_DGSW; TP_DLW_NH_N0 = TP_DGSW; % NH N0 TP_Diameter Speaking/Listening Windows
TP_DSW_NH_N60 = TP_DGSW; TP_DLW_NH_N60 = TP_DGSW; % NH N60 TP_Diameter Speaking/Listening Windows
TP_DSW_NH_N70 = TP_DGSW; TP_DLW_NH_N70 = TP_DGSW; % NH N70 TP_Diameter Speaking/Listening Windows
TP_DSW_HI_N0 = TP_DGSW; TP_DLW_HI_N0 = TP_DGSW; % HI N0 TP_Diameter Speaking/Listening Windows
TP_DSW_HI_N60 = TP_DGSW; TP_DLW_HI_N60 = TP_DGSW; % HI N60 TP_Diameter Speaking/Listening Windows
TP_DSW_HI_N70 = TP_DGSW; TP_DLW_HI_N70 = TP_DGSW; % HI N70 TP_Diameter Speaking/Listening Windows
TP_DSW_NH_UN = TP_DGSW; TP_DLW_NH_UN = TP_DGSW; % NH UN TP_Diameter Speaking/Listening Windows
TP_DSW_NH_AA = TP_DGSW; TP_DLW_NH_AA = TP_DGSW; % NH AA TP_Diameter Speaking/Listening Windows
TP_DSW_NH_AB = TP_DGSW; TP_DLW_NH_AB = TP_DGSW; % NH AB TP_Diameter Speaking/Listening Windows
TP_DSW_NH_UN_N0 = TP_DGSW; TP_DLW_NH_UN_N0 = TP_DGSW; % NH UN N0 TP_Diameter Speaking/Listening Windows
TP_DSW_NH_UN_N60 = TP_DGSW; TP_DLW_NH_UN_N60 = TP_DGSW; % NH UN N60 TP_Diameter Speaking/Listening Windows
TP_DSW_NH_UN_N70 = TP_DGSW; TP_DLW_NH_UN_N70 = TP_DGSW; % NH UN N70 TP_Diameter Speaking/Listening Windows
TP_DSW_NH_AA_N0 = TP_DGSW; TP_DLW_NH_AA_N0 = TP_DGSW; % NH AA N0 TP_Diameter Speaking/Listening Windows
TP_DSW_NH_AA_N60 = TP_DGSW; TP_DLW_NH_AA_N60 = TP_DGSW; % NH AA N60 TP_Diameter Speaking/Listening Windows
TP_DSW_NH_AA_N70 = TP_DGSW; TP_DLW_NH_AA_N70 = TP_DGSW; % NH AA N70 TP_Diameter Speaking/Listening Windows
TP_DSW_NH_AB_N0 = TP_DGSW; TP_DLW_NH_AB_N0 = TP_DGSW; % NH AB N0 TP_Diameter Speaking/Listening Windows
TP_DSW_NH_AB_N60 = TP_DGSW; TP_DLW_NH_AB_N60 = TP_DGSW; % NH AB N60 TP_Diameter Speaking/Listening Windows
TP_DSW_NH_AB_N70 = TP_DGSW; TP_DLW_NH_AB_N70 = TP_DGSW; % NH AB N70 TP_Diameter Speaking/Listening Windows
TP_DSW_HI_UN = TP_DGSW; TP_DLW_HI_UN = TP_DGSW; % HI UN TP_Diameter Speaking/Listening Windows
TP_DSW_HI_AA = TP_DGSW; TP_DLW_HI_AA = TP_DGSW; % HI AA TP_Diameter Speaking/Listening Windows
TP_DSW_HI_AB = TP_DGSW; TP_DLW_HI_AB = TP_DGSW; % HI AB TP_Diameter Speaking/Listening Windows
TP_DSW_HI_UN_N0 = TP_DGSW; TP_DLW_HI_UN_N0 = TP_DGSW; % HI UN N0 TP_Diameter Speaking/Listening Windows
TP_DSW_HI_UN_N60 = TP_DGSW; TP_DLW_HI_UN_N60 = TP_DGSW; % HI UN N60 TP_Diameter Speaking/Listening Windows
TP_DSW_HI_UN_N70 = TP_DGSW; TP_DLW_HI_UN_N70 = TP_DGSW; % HI UN N70 TP_Diameter Speaking/Listening Windows
TP_DSW_HI_AA_N0 = TP_DGSW; TP_DLW_HI_AA_N0 = TP_DGSW; % HI AA N0 TP_Diameter Speaking/Listening Windows
TP_DSW_HI_AA_N60 = TP_DGSW; TP_DLW_HI_AA_N60 = TP_DGSW; % HI AA N60 TP_Diameter Speaking/Listening Windows
TP_DSW_HI_AA_N70 = TP_DGSW; TP_DLW_HI_AA_N70 = TP_DGSW; % HI AA N70 TP_Diameter Speaking/Listening Windows
TP_DSW_HI_AB_N0 = TP_DGSW; TP_DLW_HI_AB_N0 = TP_DGSW; % HI AB N0 TP_Diameter Speaking/Listening Windows
TP_DSW_HI_AB_N60 = TP_DGSW; TP_DLW_HI_AB_N60 = TP_DGSW; % HI AB N60 TP_Diameter Speaking/Listening Windows
TP_DSW_HI_AB_N70 = TP_DGSW; TP_DLW_HI_AB_N70 = TP_DGSW; % HI AB N70 TP_Diameter Speaking/Listening Windows

TP_DGSW_B = TP_DGSW; TP_DGLW_B = TP_DGSW; % Global TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_NH_B = TP_DGSW; TP_DLW_NH_B = TP_DGSW; % NH TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_HI_B = TP_DGSW; TP_DLW_HI_B = TP_DGSW; % HI TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_N0_B = TP_DGSW; TP_DLW_N0_B = TP_DGSW; % N0 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_N60_B = TP_DGSW; TP_DLW_N60_B = TP_DGSW; % N60 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_N70_B = TP_DGSW; TP_DLW_N70_B = TP_DGSW; % N70 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_NH_N0_B = TP_DGSW; TP_DLW_NH_N0_B = TP_DGSW; % NH N0 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_NH_N60_B = TP_DGSW; TP_DLW_NH_N60_B = TP_DGSW; % NH N60 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_NH_N70_B = TP_DGSW; TP_DLW_NH_N70_B = TP_DGSW; % NH N70 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_HI_N0_B = TP_DGSW; TP_DLW_HI_N0_B = TP_DGSW; % HI N0 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_HI_N60_B = TP_DGSW; TP_DLW_HI_N60_B = TP_DGSW; % HI N60 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_HI_N70_B = TP_DGSW; TP_DLW_HI_N70_B = TP_DGSW; % HI N70 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_NH_UN_B = TP_DGSW; TP_DLW_NH_UN_B = TP_DGSW; % NH UN TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_NH_AA_B = TP_DGSW; TP_DLW_NH_AA_B = TP_DGSW; % NH AA TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_NH_AB_B = TP_DGSW; TP_DLW_NH_AB_B = TP_DGSW; % NH AB TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_NH_UN_N0_B = TP_DGSW; TP_DLW_NH_UN_N0_B = TP_DGSW; % NH UN N0 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_NH_UN_N60_B = TP_DGSW; TP_DLW_NH_UN_N60_B = TP_DGSW; % NH UN N60 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_NH_UN_N70_B = TP_DGSW; TP_DLW_NH_UN_N70_B = TP_DGSW; % NH UN N70 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_NH_AA_N0_B = TP_DGSW; TP_DLW_NH_AA_N0_B = TP_DGSW; % NH AA N0 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_NH_AA_N60_B = TP_DGSW; TP_DLW_NH_AA_N60_B = TP_DGSW; % NH AA N60 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_NH_AA_N70_B = TP_DGSW; TP_DLW_NH_AA_N70_B = TP_DGSW; % NH AA N70 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_NH_AB_N0_B = TP_DGSW; TP_DLW_NH_AB_N0_B = TP_DGSW; % NH AB N0 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_NH_AB_N60_B = TP_DGSW; TP_DLW_NH_AB_N60_B = TP_DGSW; % NH AB N60 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_NH_AB_N70_B = TP_DGSW; TP_DLW_NH_AB_N70_B = TP_DGSW; % NH AB N70 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_HI_UN_B = TP_DGSW; TP_DLW_HI_UN_B = TP_DGSW; % HI UN TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_HI_AA_B = TP_DGSW; TP_DLW_HI_AA_B = TP_DGSW; % HI AA TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_HI_AB_B = TP_DGSW; TP_DLW_HI_AB_B = TP_DGSW; % HI AB TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_HI_UN_N0_B = TP_DGSW; TP_DLW_HI_UN_N0_B = TP_DGSW; % HI UN N0 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_HI_UN_N60_B = TP_DGSW; TP_DLW_HI_UN_N60_B = TP_DGSW; % HI UN N60 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_HI_UN_N70_B = TP_DGSW; TP_DLW_HI_UN_N70_B = TP_DGSW; % HI UN N70 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_HI_AA_N0_B = TP_DGSW; TP_DLW_HI_AA_N0_B = TP_DGSW; % HI AA N0 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_HI_AA_N60_B = TP_DGSW; TP_DLW_HI_AA_N60_B = TP_DGSW; % HI AA N60 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_HI_AA_N70_B = TP_DGSW; TP_DLW_HI_AA_N70_B = TP_DGSW; % HI AA N70 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_HI_AB_N0_B = TP_DGSW; TP_DLW_HI_AB_N0_B = TP_DGSW; % HI AB N0 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_HI_AB_N60_B = TP_DGSW; TP_DLW_HI_AB_N60_B = TP_DGSW; % HI AB N60 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected
TP_DSW_HI_AB_N70_B = TP_DGSW; TP_DLW_HI_AB_N70_B = TP_DGSW; % HI AB N70 TP_Diameter Speaking/Listening Windows Adaptive_Baseline corrected

TP_FGSW = TP_DGSW; TP_FGLW = TP_DGSW; % Global TP_Fixation Speaking/Listening Windows
TP_FSW_NH = TP_DGSW; TP_FLW_NH = TP_DGSW; % NH TP_Fixation Speaking/Listening Windows
TP_FSW_HI = TP_DGSW; TP_FLW_HI = TP_DGSW; % HI TP_Fixation Speaking/Listening Windows
TP_FSW_N0 = TP_DGSW; TP_FLW_N0 = TP_DGSW; % N0 TP_Fixation Speaking/Listening Windows
TP_FSW_N60 = TP_DGSW; TP_FLW_N60 = TP_DGSW; % N60 TP_Fixation Speaking/Listening Windows
TP_FSW_N70 = TP_DGSW; TP_FLW_N70 = TP_DGSW; % N70 TP_Fixation Speaking/Listening Windows
TP_FSW_NH_N0 = TP_DGSW; TP_FLW_NH_N0 = TP_DGSW; % NH N0 TP_Fixation Speaking/Listening Windows
TP_FSW_NH_N60 = TP_DGSW; TP_FLW_NH_N60 = TP_DGSW; % NH N60 TP_Fixation Speaking/Listening Windows
TP_FSW_NH_N70 = TP_DGSW; TP_FLW_NH_N70 = TP_DGSW; % NH N70 TP_Fixation Speaking/Listening Windows
TP_FSW_HI_N0 = TP_DGSW; TP_FLW_HI_N0 = TP_DGSW; % HI N0 TP_Fixation Speaking/Listening Windows
TP_FSW_HI_N60 = TP_DGSW; TP_FLW_HI_N60 = TP_DGSW; % HI N60 TP_Fixation Speaking/Listening Windows
TP_FSW_HI_N70 = TP_DGSW; TP_FLW_HI_N70 = TP_DGSW; % HI N70 TP_Fixation Speaking/Listening Windows
TP_FSW_NH_UN = TP_DGSW; TP_FLW_NH_UN = TP_DGSW; % NH UN TP_Fixation Speaking/Listening Windows
TP_FSW_NH_AA = TP_DGSW; TP_FLW_NH_AA = TP_DGSW; % NH AA TP_Fixation Speaking/Listening Windows
TP_FSW_NH_AB = TP_DGSW; TP_FLW_NH_AB = TP_DGSW; % NH AB TP_Fixation Speaking/Listening Windows
TP_FSW_NH_UN_N0 = TP_DGSW; TP_FLW_NH_UN_N0 = TP_DGSW; % NH UN N0 TP_Fixation Speaking/Listening Windows
TP_FSW_NH_UN_N60 = TP_DGSW; TP_FLW_NH_UN_N60 = TP_DGSW; % NH UN N60 TP_Fixation Speaking/Listening Windows
TP_FSW_NH_UN_N70 = TP_DGSW; TP_FLW_NH_UN_N70 = TP_DGSW; % NH UN N70 TP_Fixation Speaking/Listening Windows
TP_FSW_NH_AA_N0 = TP_DGSW; TP_FLW_NH_AA_N0 = TP_DGSW; % NH AA N0 TP_Fixation Speaking/Listening Windows
TP_FSW_NH_AA_N60 = TP_DGSW; TP_FLW_NH_AA_N60 = TP_DGSW; % NH AA N60 TP_Fixation Speaking/Listening Windows
TP_FSW_NH_AA_N70 = TP_DGSW; TP_FLW_NH_AA_N70 = TP_DGSW; % NH AA N70 TP_Fixation Speaking/Listening Windows
TP_FSW_NH_AB_N0 = TP_DGSW; TP_FLW_NH_AB_N0 = TP_DGSW; % NH AB N0 TP_Fixation Speaking/Listening Windows
TP_FSW_NH_AB_N60 = TP_DGSW; TP_FLW_NH_AB_N60 = TP_DGSW; % NH AB N60 TP_Fixation Speaking/Listening Windows
TP_FSW_NH_AB_N70 = TP_DGSW; TP_FLW_NH_AB_N70 = TP_DGSW; % NH AB N70 TP_Fixation Speaking/Listening Windows
TP_FSW_HI_UN = TP_DGSW; TP_FLW_HI_UN = TP_DGSW; % HI UN TP_Fixation Speaking/Listening Windows
TP_FSW_HI_AA = TP_DGSW; TP_FLW_HI_AA = TP_DGSW; % HI AA TP_Fixation Speaking/Listening Windows
TP_FSW_HI_AB = TP_DGSW; TP_FLW_HI_AB = TP_DGSW; % HI AB TP_Fixation Speaking/Listening Windows
TP_FSW_HI_UN_N0 = TP_DGSW; TP_FLW_HI_UN_N0 = TP_DGSW; % HI UN N0 TP_Fixation Speaking/Listening Windows
TP_FSW_HI_UN_N60 = TP_DGSW; TP_FLW_HI_UN_N60 = TP_DGSW; % HI UN N60 TP_Fixation Speaking/Listening Windows
TP_FSW_HI_UN_N70 = TP_DGSW; TP_FLW_HI_UN_N70 = TP_DGSW; % HI UN N70 TP_Fixation Speaking/Listening Windows
TP_FSW_HI_AA_N0 = TP_DGSW; TP_FLW_HI_AA_N0 = TP_DGSW; % HI AA N0 TP_Fixation Speaking/Listening Windows
TP_FSW_HI_AA_N60 = TP_DGSW; TP_FLW_HI_AA_N60 = TP_DGSW; % HI AA N60 TP_Fixation Speaking/Listening Windows
TP_FSW_HI_AA_N70 = TP_DGSW; TP_FLW_HI_AA_N70 = TP_DGSW; % HI AA N70 TP_Fixation Speaking/Listening Windows
TP_FSW_HI_AB_N0 = TP_DGSW; TP_FLW_HI_AB_N0 = TP_DGSW; % HI AB N0 TP_Fixation Speaking/Listening Windows
TP_FSW_HI_AB_N60 = TP_DGSW; TP_FLW_HI_AB_N60 = TP_DGSW; % HI AB N60 TP_Fixation Speaking/Listening Windows
TP_FSW_HI_AB_N70 = TP_DGSW; TP_FLW_HI_AB_N70 = TP_DGSW; % HI AB N70 TP_Fixation Speaking/Listening Windows

for i=1:NTPs_II
    TPsF = nonzeros(TPsOrder_II(i,:));
    TP_DGSW(i,:) = (reshape(mean(DGSW(intersect(TPsF,TPsG),:,:),[1 2],'omitnan'),[],1))';
    TP_DGSW_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,TPsG),:,:),[1 2],'omitnan'),[],1))';
    TP_DGLW(i,:) = (reshape(mean(DGLW(intersect(TPsF,TPsG),:,:),[1 2],'omitnan'),[],1))';
    TP_DGLW_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,TPsG),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH(i,:) = (reshape(mean(DGSW(intersect(TPsF,TPsNH),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,TPsNH),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH(i,:) = (reshape(mean(DGLW(intersect(TPsF,TPsNH),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,TPsNH),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI(i,:) = (reshape(mean(DGSW(intersect(TPsF,TPsHI),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,TPsHI),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI(i,:) = (reshape(mean(DGLW(intersect(TPsF,TPsHI),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,TPsHI),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_N0(i,:) = (reshape(mean(DGSW(intersect(TPsF,TPsN0),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_N0_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,TPsN0),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_N0(i,:) = (reshape(mean(DGLW(intersect(TPsF,TPsN0),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_N0_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,TPsN0),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_N60(i,:) = (reshape(mean(DGSW(intersect(TPsF,TPsN60),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_N60_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,TPsN60),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_N60(i,:) = (reshape(mean(DGLW(intersect(TPsF,TPsN60),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_N60_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,TPsN60),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_N70(i,:) = (reshape(mean(DGSW(intersect(TPsF,TPsN70),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_N70_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,TPsN70),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_N70(i,:) = (reshape(mean(DGLW(intersect(TPsF,TPsN70),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_N70_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,TPsN70),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_N0(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(TPsNH,TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_N0_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(TPsNH,TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_N0(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(TPsNH,TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_N0_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(TPsNH,TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_N60(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(TPsNH,TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_N60_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(TPsNH,TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_N60(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(TPsNH,TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_N60_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(TPsNH,TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_N70(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(TPsNH,TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_N70_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(TPsNH,TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_N70(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(TPsNH,TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_N70_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(TPsNH,TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_N0(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(TPsHI,TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_N0_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(TPsHI,TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_N0(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(TPsHI,TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_N0_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(TPsHI,TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_N60(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(TPsHI,TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_N60_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(TPsHI,TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_N60(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(TPsHI,TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_N60_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(TPsHI,TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_N70(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(TPsHI,TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_N70_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(TPsHI,TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_N70(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(TPsHI,TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_N70_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(TPsHI,TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_UN(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(TPsNH,TPsUN)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_UN_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(TPsNH,TPsUN)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_UN(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(TPsNH,TPsUN)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_UN_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(TPsNH,TPsUN)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_AA(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(TPsNH,TPsAA)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_AA_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(TPsNH,TPsAA)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_AA(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(TPsNH,TPsAA)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_AA_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(TPsNH,TPsAA)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_AB(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(TPsNH,TPsAB)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_AB_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(TPsNH,TPsAB)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_AB(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(TPsNH,TPsAB)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_AB_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(TPsNH,TPsAB)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_UN_N0(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(intersect(TPsNH,TPsUN),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_UN_N0_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(intersect(TPsNH,TPsUN),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_UN_N0(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(intersect(TPsNH,TPsUN),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_UN_N0_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(intersect(TPsNH,TPsUN),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_UN_N60(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(intersect(TPsNH,TPsUN),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_UN_N60_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(intersect(TPsNH,TPsUN),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_UN_N60(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(intersect(TPsNH,TPsUN),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_UN_N60_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(intersect(TPsNH,TPsUN),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_UN_N70(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(intersect(TPsNH,TPsUN),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_UN_N70_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(intersect(TPsNH,TPsUN),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_UN_N70(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(intersect(TPsNH,TPsUN),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_UN_N70_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(intersect(TPsNH,TPsUN),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_AA_N0(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(intersect(TPsNH,TPsAA),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_AA_N0_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(intersect(TPsNH,TPsAA),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_AA_N0(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(intersect(TPsNH,TPsAA),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_AA_N0_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(intersect(TPsNH,TPsAA),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_AA_N60(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(intersect(TPsNH,TPsAA),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_AA_N60_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(intersect(TPsNH,TPsAA),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_AA_N60(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(intersect(TPsNH,TPsAA),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_AA_N60_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(intersect(TPsNH,TPsAA),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_AA_N70(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(intersect(TPsNH,TPsAA),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_AA_N70_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(intersect(TPsNH,TPsAA),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_AA_N70(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(intersect(TPsNH,TPsAA),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_AA_N70_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(intersect(TPsNH,TPsAA),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_AB_N0(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(intersect(TPsNH,TPsAB),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_AB_N0_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(intersect(TPsNH,TPsAB),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_AB_N0(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(intersect(TPsNH,TPsAB),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_AB_N0_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(intersect(TPsNH,TPsAB),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_AB_N60(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(intersect(TPsNH,TPsAB),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_AB_N60_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(intersect(TPsNH,TPsAB),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_AB_N60(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(intersect(TPsNH,TPsAB),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_AB_N60_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(intersect(TPsNH,TPsAB),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_AB_N70(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(intersect(TPsNH,TPsAB),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_NH_AB_N70_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(intersect(TPsNH,TPsAB),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_AB_N70(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(intersect(TPsNH,TPsAB),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_NH_AB_N70_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(intersect(TPsNH,TPsAB),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_UN(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(TPsHI,TPsUN)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_UN_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(TPsHI,TPsUN)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_UN(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(TPsHI,TPsUN)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_UN_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(TPsHI,TPsUN)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_AA(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(TPsHI,TPsAA)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_AA_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(TPsHI,TPsAA)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_AA(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(TPsHI,TPsAA)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_AA_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(TPsHI,TPsAA)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_AB(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(TPsHI,TPsAB)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_AB_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(TPsHI,TPsAB)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_AB(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(TPsHI,TPsAB)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_AB_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(TPsHI,TPsAB)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_UN_N0(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(intersect(TPsHI,TPsUN),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_UN_N0_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(intersect(TPsHI,TPsUN),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_UN_N0(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(intersect(TPsHI,TPsUN),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_UN_N0_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(intersect(TPsHI,TPsUN),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_UN_N60(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(intersect(TPsHI,TPsUN),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_UN_N60_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(intersect(TPsHI,TPsUN),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_UN_N60(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(intersect(TPsHI,TPsUN),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_UN_N60_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(intersect(TPsHI,TPsUN),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_UN_N70(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(intersect(TPsHI,TPsUN),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_UN_N70_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(intersect(TPsHI,TPsUN),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_UN_N70(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(intersect(TPsHI,TPsUN),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_UN_N70_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(intersect(TPsHI,TPsUN),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_AA_N0(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(intersect(TPsHI,TPsAA),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_AA_N0_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(intersect(TPsHI,TPsAA),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_AA_N0(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(intersect(TPsHI,TPsAA),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_AA_N0_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(intersect(TPsHI,TPsAA),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_AA_N60(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(intersect(TPsHI,TPsAA),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_AA_N60_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(intersect(TPsHI,TPsAA),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_AA_N60(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(intersect(TPsHI,TPsAA),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_AA_N60_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(intersect(TPsHI,TPsAA),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_AA_N70(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(intersect(TPsHI,TPsAA),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_AA_N70_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(intersect(TPsHI,TPsAA),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_AA_N70(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(intersect(TPsHI,TPsAA),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_AA_N70_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(intersect(TPsHI,TPsAA),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_AB_N0(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(intersect(TPsHI,TPsAB),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_AB_N0_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(intersect(TPsHI,TPsAB),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_AB_N0(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(intersect(TPsHI,TPsAB),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_AB_N0_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(intersect(TPsHI,TPsAB),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_AB_N60(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(intersect(TPsHI,TPsAB),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_AB_N60_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(intersect(TPsHI,TPsAB),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_AB_N60(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(intersect(TPsHI,TPsAB),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_AB_N60_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(intersect(TPsHI,TPsAB),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_AB_N70(i,:) = (reshape(mean(DGSW(intersect(TPsF,intersect(intersect(TPsHI,TPsAB),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DSW_HI_AB_N70_B(i,:) = (reshape(mean(DGSW_B(intersect(TPsF,intersect(intersect(TPsHI,TPsAB),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_AB_N70(i,:) = (reshape(mean(DGLW(intersect(TPsF,intersect(intersect(TPsHI,TPsAB),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_DLW_HI_AB_N70_B(i,:) = (reshape(mean(DGLW_B(intersect(TPsF,intersect(intersect(TPsHI,TPsAB),TPsN70)),:,:),[1 2],'omitnan'),[],1))';

    TP_FGSW(i,:) = (reshape(mean(FGSW(intersect(TPsF,TPsG),:,:),[1 2],'omitnan'),[],1))';
    TP_FGLW(i,:) = (reshape(mean(FGLW(intersect(TPsF,TPsG),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_NH(i,:) = (reshape(mean(FGSW(intersect(TPsF,TPsNH),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_NH(i,:) = (reshape(mean(FGLW(intersect(TPsF,TPsNH),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_HI(i,:) = (reshape(mean(FGSW(intersect(TPsF,TPsHI),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_HI(i,:) = (reshape(mean(FGLW(intersect(TPsF,TPsHI),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_N0(i,:) = (reshape(mean(FGSW(intersect(TPsF,TPsN0),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_N0(i,:) = (reshape(mean(FGLW(intersect(TPsF,TPsN0),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_N60(i,:) = (reshape(mean(FGSW(intersect(TPsF,TPsN60),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_N60(i,:) = (reshape(mean(FGLW(intersect(TPsF,TPsN60),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_N70(i,:) = (reshape(mean(FGSW(intersect(TPsF,TPsN70),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_N70(i,:) = (reshape(mean(FGLW(intersect(TPsF,TPsN70),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_NH_N0(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(TPsNH,TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_NH_N0(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(TPsNH,TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_NH_N60(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(TPsNH,TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_NH_N60(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(TPsNH,TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_NH_N70(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(TPsNH,TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_NH_N70(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(TPsNH,TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_HI_N0(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(TPsHI,TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_HI_N0(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(TPsHI,TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_HI_N60(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(TPsHI,TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_HI_N60(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(TPsHI,TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_HI_N70(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(TPsHI,TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_HI_N70(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(TPsHI,TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_NH_UN(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(TPsNH,TPsUN)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_NH_UN(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(TPsNH,TPsUN)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_NH_AA(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(TPsNH,TPsAA)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_NH_AA(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(TPsNH,TPsAA)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_NH_AB(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(TPsNH,TPsAB)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_NH_AB(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(TPsNH,TPsAB)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_NH_UN_N0(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(intersect(TPsNH,TPsUN),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_NH_UN_N0(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(intersect(TPsNH,TPsUN),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_NH_UN_N60(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(intersect(TPsNH,TPsUN),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_NH_UN_N60(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(intersect(TPsNH,TPsUN),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_NH_UN_N70(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(intersect(TPsNH,TPsUN),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_NH_UN_N70(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(intersect(TPsNH,TPsUN),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_NH_AA_N0(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(intersect(TPsNH,TPsAA),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_NH_AA_N0(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(intersect(TPsNH,TPsAA),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_NH_AA_N60(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(intersect(TPsNH,TPsAA),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_NH_AA_N60(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(intersect(TPsNH,TPsAA),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_NH_AA_N70(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(intersect(TPsNH,TPsAA),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_NH_AA_N70(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(intersect(TPsNH,TPsAA),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_NH_AB_N0(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(intersect(TPsNH,TPsAB),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_NH_AB_N0(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(intersect(TPsNH,TPsAB),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_NH_AB_N60(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(intersect(TPsNH,TPsAB),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_NH_AB_N60(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(intersect(TPsNH,TPsAB),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_NH_AB_N70(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(intersect(TPsNH,TPsAB),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_NH_AB_N70(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(intersect(TPsNH,TPsAB),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_HI_UN(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(TPsHI,TPsUN)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_HI_UN(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(TPsHI,TPsUN)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_HI_AA(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(TPsHI,TPsAA)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_HI_AA(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(TPsHI,TPsAA)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_HI_AB(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(TPsHI,TPsAB)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_HI_AB(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(TPsHI,TPsAB)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_HI_UN_N0(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(intersect(TPsHI,TPsUN),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_HI_UN_N0(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(intersect(TPsHI,TPsUN),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_HI_UN_N60(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(intersect(TPsHI,TPsUN),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_HI_UN_N60(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(intersect(TPsHI,TPsUN),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_HI_UN_N70(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(intersect(TPsHI,TPsUN),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_HI_UN_N70(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(intersect(TPsHI,TPsUN),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_HI_AA_N0(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(intersect(TPsHI,TPsAA),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_HI_AA_N0(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(intersect(TPsHI,TPsAA),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_HI_AA_N60(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(intersect(TPsHI,TPsAA),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_HI_AA_N60(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(intersect(TPsHI,TPsAA),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_HI_AA_N70(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(intersect(TPsHI,TPsAA),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_HI_AA_N70(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(intersect(TPsHI,TPsAA),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_HI_AB_N0(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(intersect(TPsHI,TPsAB),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_HI_AB_N0(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(intersect(TPsHI,TPsAB),TPsN0)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_HI_AB_N60(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(intersect(TPsHI,TPsAB),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_HI_AB_N60(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(intersect(TPsHI,TPsAB),TPsN60)),:,:),[1 2],'omitnan'),[],1))';
    TP_FSW_HI_AB_N70(i,:) = (reshape(mean(FGSW(intersect(TPsF,intersect(intersect(TPsHI,TPsAB),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
    TP_FLW_HI_AB_N70(i,:) = (reshape(mean(FGLW(intersect(TPsF,intersect(intersect(TPsHI,TPsAB),TPsN70)),:,:),[1 2],'omitnan'),[],1))';
end

% Change 0 to NaN
TP_DGSW(TP_DGSW==0)=NaN;
TP_DGSW_B(TP_DGSW_B==0)=NaN;
TP_DGLW(TP_DGLW==0)=NaN;
TP_DGLW_B(TP_DGLW_B==0)=NaN;

TP_DSW_NH(TP_DSW_NH==0)=NaN;
TP_DSW_NH_B(TP_DSW_NH_B==0)=NaN;
TP_DLW_NH(TP_DLW_NH==0)=NaN;
TP_DLW_NH_B(TP_DLW_NH_B==0)=NaN;
TP_DSW_HI(TP_DSW_HI==0)=NaN;
TP_DSW_HI_B(TP_DSW_HI_B==0)=NaN;
TP_DLW_HI(TP_DLW_HI==0)=NaN;
TP_DLW_HI_B(TP_DLW_HI_B==0)=NaN;

TP_DSW_N0(TP_DSW_N0==0)=NaN;
TP_DSW_N0_B(TP_DSW_N0_B==0)=NaN;
TP_DLW_N0(TP_DLW_N0==0)=NaN;
TP_DLW_N0_B(TP_DLW_N0_B==0)=NaN;
TP_DSW_N60(TP_DSW_N60==0)=NaN;
TP_DSW_N60_B(TP_DSW_N60_B==0)=NaN;
TP_DLW_N60(TP_DLW_N60==0)=NaN;
TP_DLW_N60_B(TP_DLW_N60_B==0)=NaN;
TP_DSW_N70(TP_DSW_N70==0)=NaN;
TP_DSW_N70_B(TP_DSW_N70_B==0)=NaN;
TP_DLW_N70(TP_DLW_N70==0)=NaN;
TP_DLW_N70_B(TP_DLW_N70_B==0)=NaN;

TP_DSW_NH_N0(TP_DSW_NH_N0==0)=NaN;
TP_DSW_NH_N0_B(TP_DSW_NH_N0_B==0)=NaN;
TP_DLW_NH_N0(TP_DLW_NH_N0==0)=NaN;
TP_DLW_NH_N0_B(TP_DLW_NH_N0_B==0)=NaN;
TP_DSW_NH_N60(TP_DSW_NH_N60==0)=NaN;
TP_DSW_NH_N60_B(TP_DSW_NH_N60_B==0)=NaN;
TP_DLW_NH_N60(TP_DLW_NH_N60==0)=NaN;
TP_DLW_NH_N60_B(TP_DLW_NH_N60_B==0)=NaN;
TP_DSW_NH_N70(TP_DSW_NH_N70==0)=NaN;
TP_DSW_NH_N70_B(TP_DSW_NH_N70_B==0)=NaN;
TP_DLW_NH_N70(TP_DLW_NH_N70==0)=NaN;
TP_DLW_NH_N70_B(TP_DLW_NH_N70_B==0)=NaN;

TP_DSW_HI_N0(TP_DSW_HI_N0==0)=NaN;
TP_DSW_HI_N0_B(TP_DSW_HI_N0_B==0)=NaN;
TP_DLW_HI_N0(TP_DLW_HI_N0==0)=NaN;
TP_DLW_HI_N0_B(TP_DLW_HI_N0_B==0)=NaN;
TP_DSW_HI_N60(TP_DSW_HI_N60==0)=NaN;
TP_DSW_HI_N60_B(TP_DSW_HI_N60_B==0)=NaN;
TP_DLW_HI_N60(TP_DLW_HI_N60==0)=NaN;
TP_DLW_HI_N60_B(TP_DLW_HI_N60_B==0)=NaN;
TP_DSW_HI_N70(TP_DSW_HI_N70==0)=NaN;
TP_DSW_HI_N70_B(TP_DSW_HI_N70_B==0)=NaN;
TP_DLW_HI_N70(TP_DLW_HI_N70==0)=NaN;
TP_DLW_HI_N70_B(TP_DLW_HI_N70_B==0)=NaN;

TP_DSW_NH_UN(TP_DSW_NH_UN==0)=NaN;
TP_DSW_NH_UN_B(TP_DSW_NH_UN_B==0)=NaN;
TP_DLW_NH_UN(TP_DLW_NH_UN==0)=NaN;
TP_DLW_NH_UN_B(TP_DLW_NH_UN_B==0)=NaN;
TP_DSW_NH_AA(TP_DSW_NH_AA==0)=NaN;
TP_DSW_NH_AA_B(TP_DSW_NH_AA_B==0)=NaN;
TP_DLW_NH_AA(TP_DLW_NH_AA==0)=NaN;
TP_DLW_NH_AA_B(TP_DLW_NH_AA_B==0)=NaN;
TP_DSW_NH_AB(TP_DSW_NH_AB==0)=NaN;
TP_DSW_NH_AB_B(TP_DSW_NH_AB_B==0)=NaN;
TP_DLW_NH_AB(TP_DLW_NH_AB==0)=NaN;
TP_DLW_NH_AB_B(TP_DLW_NH_AB_B==0)=NaN;

TP_DSW_HI_UN(TP_DSW_HI_UN==0)=NaN;
TP_DSW_HI_UN_B(TP_DSW_HI_UN_B==0)=NaN;
TP_DLW_HI_UN(TP_DLW_HI_UN==0)=NaN;
TP_DLW_HI_UN_B(TP_DLW_HI_UN_B==0)=NaN;
TP_DSW_HI_AA(TP_DSW_HI_AA==0)=NaN;
TP_DSW_HI_AA_B(TP_DSW_HI_AA_B==0)=NaN;
TP_DLW_HI_AA(TP_DLW_HI_AA==0)=NaN;
TP_DLW_HI_AA_B(TP_DLW_HI_AA_B==0)=NaN;
TP_DSW_HI_AB(TP_DSW_HI_AB==0)=NaN;
TP_DSW_HI_AB_B(TP_DSW_HI_AB_B==0)=NaN;
TP_DLW_HI_AB(TP_DLW_HI_AB==0)=NaN;
TP_DLW_HI_AB_B(TP_DLW_HI_AB_B==0)=NaN;

TP_DSW_NH_UN_N0(TP_DSW_NH_UN_N0==0)=NaN;
TP_DSW_NH_UN_N0_B(TP_DSW_NH_UN_N0_B==0)=NaN;
TP_DLW_NH_UN_N0(TP_DLW_NH_UN_N0==0)=NaN;
TP_DLW_NH_UN_N0_B(TP_DLW_NH_UN_N0_B==0)=NaN;
TP_DSW_NH_UN_N60(TP_DSW_NH_UN_N60==0)=NaN;
TP_DSW_NH_UN_N60_B(TP_DSW_NH_UN_N60_B==0)=NaN;
TP_DLW_NH_UN_N60(TP_DLW_NH_UN_N60==0)=NaN;
TP_DLW_NH_UN_N60_B(TP_DLW_NH_UN_N60_B==0)=NaN;
TP_DSW_NH_UN_N70(TP_DSW_NH_UN_N70==0)=NaN;
TP_DSW_NH_UN_N70_B(TP_DSW_NH_UN_N70_B==0)=NaN;
TP_DLW_NH_UN_N70(TP_DLW_NH_UN_N70==0)=NaN;
TP_DLW_NH_UN_N70_B(TP_DLW_NH_UN_N70_B==0)=NaN;
TP_DSW_NH_AA_N0(TP_DSW_NH_AA_N0==0)=NaN;
TP_DSW_NH_AA_N0_B(TP_DSW_NH_AA_N0_B==0)=NaN;
TP_DLW_NH_AA_N0(TP_DLW_NH_AA_N0==0)=NaN;
TP_DLW_NH_AA_N0_B(TP_DLW_NH_AA_N0_B==0)=NaN;
TP_DSW_NH_AA_N60(TP_DSW_NH_AA_N60==0)=NaN;
TP_DSW_NH_AA_N60_B(TP_DSW_NH_AA_N60_B==0)=NaN;
TP_DLW_NH_AA_N60(TP_DLW_NH_AA_N60==0)=NaN;
TP_DLW_NH_AA_N60_B(TP_DLW_NH_AA_N60_B==0)=NaN;
TP_DSW_NH_AA_N70(TP_DSW_NH_AA_N70==0)=NaN;
TP_DSW_NH_AA_N70_B(TP_DSW_NH_AA_N70_B==0)=NaN;
TP_DLW_NH_AA_N70(TP_DLW_NH_AA_N70==0)=NaN;
TP_DLW_NH_AA_N70_B(TP_DLW_NH_AA_N70_B==0)=NaN;
TP_DSW_NH_AB_N0(TP_DSW_NH_AB_N0==0)=NaN;
TP_DSW_NH_AB_N0_B(TP_DSW_NH_AB_N0_B==0)=NaN;
TP_DLW_NH_AB_N0(TP_DLW_NH_AB_N0==0)=NaN;
TP_DLW_NH_AB_N0_B(TP_DLW_NH_AB_N0_B==0)=NaN;
TP_DSW_NH_AB_N60(TP_DSW_NH_AB_N60==0)=NaN;
TP_DSW_NH_AB_N60_B(TP_DSW_NH_AB_N60_B==0)=NaN;
TP_DLW_NH_AB_N60(TP_DLW_NH_AB_N60==0)=NaN;
TP_DLW_NH_AB_N60_B(TP_DLW_NH_AB_N60_B==0)=NaN;
TP_DSW_NH_AB_N70(TP_DSW_NH_AB_N70==0)=NaN;
TP_DSW_NH_AB_N70_B(TP_DSW_NH_AB_N70_B==0)=NaN;
TP_DLW_NH_AB_N70(TP_DLW_NH_AB_N70==0)=NaN;
TP_DLW_NH_AB_N70_B(TP_DLW_NH_AB_N70_B==0)=NaN;

TP_DSW_HI_UN_N0(TP_DSW_HI_UN_N0==0)=NaN;
TP_DSW_HI_UN_N0_B(TP_DSW_HI_UN_N0_B==0)=NaN;
TP_DLW_HI_UN_N0(TP_DLW_HI_UN_N0==0)=NaN;
TP_DLW_HI_UN_N0_B(TP_DLW_HI_UN_N0_B==0)=NaN;
TP_DSW_HI_UN_N60(TP_DSW_HI_UN_N60==0)=NaN;
TP_DSW_HI_UN_N60_B(TP_DSW_HI_UN_N60_B==0)=NaN;
TP_DLW_HI_UN_N60(TP_DLW_HI_UN_N60==0)=NaN;
TP_DLW_HI_UN_N60_B(TP_DLW_HI_UN_N60_B==0)=NaN;
TP_DSW_HI_UN_N70(TP_DSW_HI_UN_N70==0)=NaN;
TP_DSW_HI_UN_N70_B(TP_DSW_HI_UN_N70_B==0)=NaN;
TP_DLW_HI_UN_N70(TP_DLW_HI_UN_N70==0)=NaN;
TP_DLW_HI_UN_N70_B(TP_DLW_HI_UN_N70_B==0)=NaN;
TP_DSW_HI_AA_N0(TP_DSW_HI_AA_N0==0)=NaN;
TP_DSW_HI_AA_N0_B(TP_DSW_HI_AA_N0_B==0)=NaN;
TP_DLW_HI_AA_N0(TP_DLW_HI_AA_N0==0)=NaN;
TP_DLW_HI_AA_N0_B(TP_DLW_HI_AA_N0_B==0)=NaN;
TP_DSW_HI_AA_N60(TP_DSW_HI_AA_N60==0)=NaN;
TP_DSW_HI_AA_N60_B(TP_DSW_HI_AA_N60_B==0)=NaN;
TP_DLW_HI_AA_N60(TP_DLW_HI_AA_N60==0)=NaN;
TP_DLW_HI_AA_N60_B(TP_DLW_HI_AA_N60_B==0)=NaN;
TP_DSW_HI_AA_N70(TP_DSW_HI_AA_N70==0)=NaN;
TP_DSW_HI_AA_N70_B(TP_DSW_HI_AA_N70_B==0)=NaN;
TP_DLW_HI_AA_N70(TP_DLW_HI_AA_N70==0)=NaN;
TP_DLW_HI_AA_N70_B(TP_DLW_HI_AA_N70_B==0)=NaN;
TP_DSW_HI_AB_N0(TP_DSW_HI_AB_N0==0)=NaN;
TP_DSW_HI_AB_N0_B(TP_DSW_HI_AB_N0_B==0)=NaN;
TP_DLW_HI_AB_N0(TP_DLW_HI_AB_N0==0)=NaN;
TP_DLW_HI_AB_N0_B(TP_DLW_HI_AB_N0_B==0)=NaN;
TP_DSW_HI_AB_N60(TP_DSW_HI_AB_N60==0)=NaN;
TP_DSW_HI_AB_N60_B(TP_DSW_HI_AB_N60_B==0)=NaN;
TP_DLW_HI_AB_N60(TP_DLW_HI_AB_N60==0)=NaN;
TP_DLW_HI_AB_N60_B(TP_DLW_HI_AB_N60_B==0)=NaN;
TP_DSW_HI_AB_N70(TP_DSW_HI_AB_N70==0)=NaN;
TP_DSW_HI_AB_N70_B(TP_DSW_HI_AB_N70_B==0)=NaN;
TP_DLW_HI_AB_N70(TP_DLW_HI_AB_N70==0)=NaN;
TP_DLW_HI_AB_N70_B(TP_DLW_HI_AB_N70_B==0)=NaN;

TP_FGSW(TP_FGSW==0)=NaN;
TP_FGLW(TP_FGLW==0)=NaN;

TP_FSW_NH(TP_FSW_NH==0)=NaN;
TP_FLW_NH(TP_FLW_NH==0)=NaN;
TP_FSW_HI(TP_FSW_HI==0)=NaN;
TP_FLW_HI(TP_FLW_HI==0)=NaN;

TP_FSW_N0(TP_FSW_N0==0)=NaN;
TP_FLW_N0(TP_FLW_N0==0)=NaN;
TP_FSW_N60(TP_FSW_N60==0)=NaN;
TP_FLW_N60(TP_FLW_N60==0)=NaN;
TP_FSW_N70(TP_FSW_N70==0)=NaN;
TP_FLW_N70(TP_FLW_N70==0)=NaN;

TP_FSW_NH_N0(TP_FSW_NH_N0==0)=NaN;
TP_FLW_NH_N0(TP_FLW_NH_N0==0)=NaN;
TP_FSW_NH_N60(TP_FSW_NH_N60==0)=NaN;
TP_FLW_NH_N60(TP_FLW_NH_N60==0)=NaN;
TP_FSW_NH_N70(TP_FSW_NH_N70==0)=NaN;
TP_FLW_NH_N70(TP_FLW_NH_N70==0)=NaN;

TP_FSW_HI_N0(TP_FSW_HI_N0==0)=NaN;
TP_FLW_HI_N0(TP_FLW_HI_N0==0)=NaN;
TP_FSW_HI_N60(TP_FSW_HI_N60==0)=NaN;
TP_FLW_HI_N60(TP_FLW_HI_N60==0)=NaN;
TP_FSW_HI_N70(TP_FSW_HI_N70==0)=NaN;
TP_FLW_HI_N70(TP_FLW_HI_N70==0)=NaN;

TP_FSW_NH_UN(TP_FSW_NH_UN==0)=NaN;
TP_FLW_NH_UN(TP_FLW_NH_UN==0)=NaN;
TP_FSW_NH_AA(TP_FSW_NH_AA==0)=NaN;
TP_FLW_NH_AA(TP_FLW_NH_AA==0)=NaN;
TP_FSW_NH_AB(TP_FSW_NH_AB==0)=NaN;
TP_FLW_NH_AB(TP_FLW_NH_AB==0)=NaN;

TP_FSW_HI_UN(TP_FSW_HI_UN==0)=NaN;
TP_FLW_HI_UN(TP_FLW_HI_UN==0)=NaN;
TP_FSW_HI_AA(TP_FSW_HI_AA==0)=NaN;
TP_FLW_HI_AA(TP_FLW_HI_AA==0)=NaN;
TP_FSW_HI_AB(TP_FSW_HI_AB==0)=NaN;
TP_FLW_HI_AB(TP_FLW_HI_AB==0)=NaN;

TP_FSW_NH_UN_N0(TP_FSW_NH_UN_N0==0)=NaN;
TP_FLW_NH_UN_N0(TP_FLW_NH_UN_N0==0)=NaN;
TP_FSW_NH_UN_N60(TP_FSW_NH_UN_N60==0)=NaN;
TP_FLW_NH_UN_N60(TP_FLW_NH_UN_N60==0)=NaN;
TP_FSW_NH_UN_N70(TP_FSW_NH_UN_N70==0)=NaN;
TP_FLW_NH_UN_N70(TP_FLW_NH_UN_N70==0)=NaN;
TP_FSW_NH_AA_N0(TP_FSW_NH_AA_N0==0)=NaN;
TP_FLW_NH_AA_N0(TP_FLW_NH_AA_N0==0)=NaN;
TP_FSW_NH_AA_N60(TP_FSW_NH_AA_N60==0)=NaN;
TP_FLW_NH_AA_N60(TP_FLW_NH_AA_N60==0)=NaN;
TP_FSW_NH_AA_N70(TP_FSW_NH_AA_N70==0)=NaN;
TP_FLW_NH_AA_N70(TP_FLW_NH_AA_N70==0)=NaN;
TP_FSW_NH_AB_N0(TP_FSW_NH_AB_N0==0)=NaN;
TP_FLW_NH_AB_N0(TP_FLW_NH_AB_N0==0)=NaN;
TP_FSW_NH_AB_N60(TP_FSW_NH_AB_N60==0)=NaN;
TP_FLW_NH_AB_N60(TP_FLW_NH_AB_N60==0)=NaN;
TP_FSW_NH_AB_N70(TP_FSW_NH_AB_N70==0)=NaN;
TP_FLW_NH_AB_N70(TP_FLW_NH_AB_N70==0)=NaN;

TP_FSW_HI_UN_N0(TP_FSW_HI_UN_N0==0)=NaN;
TP_FLW_HI_UN_N0(TP_FLW_HI_UN_N0==0)=NaN;
TP_FSW_HI_UN_N60(TP_FSW_HI_UN_N60==0)=NaN;
TP_FLW_HI_UN_N60(TP_FLW_HI_UN_N60==0)=NaN;
TP_FSW_HI_UN_N70(TP_FSW_HI_UN_N70==0)=NaN;
TP_FLW_HI_UN_N70(TP_FLW_HI_UN_N70==0)=NaN;
TP_FSW_HI_AA_N0(TP_FSW_HI_AA_N0==0)=NaN;
TP_FLW_HI_AA_N0(TP_FLW_HI_AA_N0==0)=NaN;
TP_FSW_HI_AA_N60(TP_FSW_HI_AA_N60==0)=NaN;
TP_FLW_HI_AA_N60(TP_FLW_HI_AA_N60==0)=NaN;
TP_FSW_HI_AA_N70(TP_FSW_HI_AA_N70==0)=NaN;
TP_FLW_HI_AA_N70(TP_FLW_HI_AA_N70==0)=NaN;
TP_FSW_HI_AB_N0(TP_FSW_HI_AB_N0==0)=NaN;
TP_FLW_HI_AB_N0(TP_FLW_HI_AB_N0==0)=NaN;
TP_FSW_HI_AB_N60(TP_FSW_HI_AB_N60==0)=NaN;
TP_FLW_HI_AB_N60(TP_FLW_HI_AB_N60==0)=NaN;
TP_FSW_HI_AB_N70(TP_FSW_HI_AB_N70==0)=NaN;
TP_FLW_HI_AB_N70(TP_FLW_HI_AB_N70==0)=NaN;

% LP Filter + Mean
TP_DGSW_Mean = ndnanfilter(mean(TP_DGSW,'omitnan'),'hamming',FilterWidth);
TP_DGSW_B_Mean = ndnanfilter(mean(TP_DGSW_B,'omitnan'),'hamming',FilterWidth);
TP_DGLW_Mean = ndnanfilter(mean(TP_DGLW,'omitnan'),'hamming',FilterWidth);
TP_DGLW_B_Mean = ndnanfilter(mean(TP_DGLW_B,'omitnan'),'hamming',FilterWidth);

TP_DSW_NH_Mean = ndnanfilter(mean(TP_DSW_NH,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_B_Mean = ndnanfilter(mean(TP_DSW_NH_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_Mean = ndnanfilter(mean(TP_DLW_NH,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_B_Mean = ndnanfilter(mean(TP_DLW_NH_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_Mean = ndnanfilter(mean(TP_DSW_HI,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_B_Mean = ndnanfilter(mean(TP_DSW_HI_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_Mean = ndnanfilter(mean(TP_DLW_HI,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_B_Mean = ndnanfilter(mean(TP_DLW_HI_B,'omitnan'),'hamming',FilterWidth);

TP_DSW_N0_Mean = ndnanfilter(mean(TP_DSW_N0,'omitnan'),'hamming',FilterWidth);
TP_DSW_N0_B_Mean = ndnanfilter(mean(TP_DSW_N0_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_N0_Mean = ndnanfilter(mean(TP_DLW_N0,'omitnan'),'hamming',FilterWidth);
TP_DLW_N0_B_Mean = ndnanfilter(mean(TP_DLW_N0_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_N60_Mean = ndnanfilter(mean(TP_DSW_N60,'omitnan'),'hamming',FilterWidth);
TP_DSW_N60_B_Mean = ndnanfilter(mean(TP_DSW_N60_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_N60_Mean = ndnanfilter(mean(TP_DLW_N60,'omitnan'),'hamming',FilterWidth);
TP_DLW_N60_B_Mean = ndnanfilter(mean(TP_DLW_N60_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_N70_Mean = ndnanfilter(mean(TP_DSW_N70,'omitnan'),'hamming',FilterWidth);
TP_DSW_N70_B_Mean = ndnanfilter(mean(TP_DSW_N70_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_N70_Mean = ndnanfilter(mean(TP_DLW_N70,'omitnan'),'hamming',FilterWidth);
TP_DLW_N70_B_Mean = ndnanfilter(mean(TP_DLW_N70_B,'omitnan'),'hamming',FilterWidth);

TP_DSW_NH_N0_Mean = ndnanfilter(mean(TP_DSW_NH_N0,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_N0_B_Mean = ndnanfilter(mean(TP_DSW_NH_N0_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_N0_Mean = ndnanfilter(mean(TP_DLW_NH_N0,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_N0_B_Mean = ndnanfilter(mean(TP_DLW_NH_N0_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_N60_Mean = ndnanfilter(mean(TP_DSW_NH_N60,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_N60_B_Mean = ndnanfilter(mean(TP_DSW_NH_N60_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_N60_Mean = ndnanfilter(mean(TP_DLW_NH_N60,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_N60_B_Mean = ndnanfilter(mean(TP_DLW_NH_N60_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_N70_Mean = ndnanfilter(mean(TP_DSW_NH_N70,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_N70_B_Mean = ndnanfilter(mean(TP_DSW_NH_N70_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_N70_Mean = ndnanfilter(mean(TP_DLW_NH_N70,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_N70_B_Mean = ndnanfilter(mean(TP_DLW_NH_N70_B,'omitnan'),'hamming',FilterWidth);

TP_DSW_HI_N0_Mean = ndnanfilter(mean(TP_DSW_HI_N0,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_N0_B_Mean = ndnanfilter(mean(TP_DSW_HI_N0_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_N0_Mean = ndnanfilter(mean(TP_DLW_HI_N0,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_N0_B_Mean = ndnanfilter(mean(TP_DLW_HI_N0_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_N60_Mean = ndnanfilter(mean(TP_DSW_HI_N60,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_N60_B_Mean = ndnanfilter(mean(TP_DSW_HI_N60_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_N60_Mean = ndnanfilter(mean(TP_DLW_HI_N60,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_N60_B_Mean = ndnanfilter(mean(TP_DLW_HI_N60_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_N70_Mean = ndnanfilter(mean(TP_DSW_HI_N70,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_N70_B_Mean = ndnanfilter(mean(TP_DSW_HI_N70_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_N70_Mean = ndnanfilter(mean(TP_DLW_HI_N70,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_N70_B_Mean = ndnanfilter(mean(TP_DLW_HI_N70_B,'omitnan'),'hamming',FilterWidth);

TP_DSW_NH_UN_Mean = ndnanfilter(mean(TP_DSW_NH_UN,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_UN_B_Mean = ndnanfilter(mean(TP_DSW_NH_UN_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_UN_Mean = ndnanfilter(mean(TP_DLW_NH_UN,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_UN_B_Mean = ndnanfilter(mean(TP_DLW_NH_UN_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_AA_Mean = ndnanfilter(mean(TP_DSW_NH_AA,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_AA_B_Mean = ndnanfilter(mean(TP_DSW_NH_AA_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_AA_Mean = ndnanfilter(mean(TP_DLW_NH_AA,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_AA_B_Mean = ndnanfilter(mean(TP_DLW_NH_AA_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_AB_Mean = ndnanfilter(mean(TP_DSW_NH_AB,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_AB_B_Mean = ndnanfilter(mean(TP_DSW_NH_AB_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_AB_Mean = ndnanfilter(mean(TP_DLW_NH_AB,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_AB_B_Mean = ndnanfilter(mean(TP_DLW_NH_AB_B,'omitnan'),'hamming',FilterWidth);

TP_DSW_HI_UN_Mean = ndnanfilter(mean(TP_DSW_HI_UN,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_UN_B_Mean = ndnanfilter(mean(TP_DSW_HI_UN_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_UN_Mean = ndnanfilter(mean(TP_DLW_HI_UN,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_UN_B_Mean = ndnanfilter(mean(TP_DLW_HI_UN_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_AA_Mean = ndnanfilter(mean(TP_DSW_HI_AA,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_AA_B_Mean = ndnanfilter(mean(TP_DSW_HI_AA_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_AA_Mean = ndnanfilter(mean(TP_DLW_HI_AA,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_AA_B_Mean = ndnanfilter(mean(TP_DLW_HI_AA_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_AB_Mean = ndnanfilter(mean(TP_DSW_HI_AB,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_AB_B_Mean = ndnanfilter(mean(TP_DSW_HI_AB_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_AB_Mean = ndnanfilter(mean(TP_DLW_HI_AB,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_AB_B_Mean = ndnanfilter(mean(TP_DLW_HI_AB_B,'omitnan'),'hamming',FilterWidth);

TP_DSW_NH_UN_N0_Mean = ndnanfilter(mean(TP_DSW_NH_UN_N0,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_UN_N0_B_Mean = ndnanfilter(mean(TP_DSW_NH_UN_N0_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_UN_N0_Mean = ndnanfilter(mean(TP_DLW_NH_UN_N0,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_UN_N0_B_Mean = ndnanfilter(mean(TP_DLW_NH_UN_N0_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_UN_N60_Mean = ndnanfilter(mean(TP_DSW_NH_UN_N60,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_UN_N60_B_Mean = ndnanfilter(mean(TP_DSW_NH_UN_N60_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_UN_N60_Mean = ndnanfilter(mean(TP_DLW_NH_UN_N60,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_UN_N60_B_Mean = ndnanfilter(mean(TP_DLW_NH_UN_N60_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_UN_N70_Mean = ndnanfilter(mean(TP_DSW_NH_UN_N70,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_UN_N70_B_Mean = ndnanfilter(mean(TP_DSW_NH_UN_N70_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_UN_N70_Mean = ndnanfilter(mean(TP_DLW_NH_UN_N70,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_UN_N70_B_Mean = ndnanfilter(mean(TP_DLW_NH_UN_N70_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_AA_N0_Mean = ndnanfilter(mean(TP_DSW_NH_AA_N0,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_AA_N0_B_Mean = ndnanfilter(mean(TP_DSW_NH_AA_N0_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_AA_N0_Mean = ndnanfilter(mean(TP_DLW_NH_AA_N0,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_AA_N0_B_Mean = ndnanfilter(mean(TP_DLW_NH_AA_N0_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_AA_N60_Mean = ndnanfilter(mean(TP_DSW_NH_AA_N60,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_AA_N60_B_Mean = ndnanfilter(mean(TP_DSW_NH_AA_N60_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_AA_N60_Mean = ndnanfilter(mean(TP_DLW_NH_AA_N60,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_AA_N60_B_Mean = ndnanfilter(mean(TP_DLW_NH_AA_N60_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_AA_N70_Mean = ndnanfilter(mean(TP_DSW_NH_AA_N70,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_AA_N70_B_Mean = ndnanfilter(mean(TP_DSW_NH_AA_N70_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_AA_N70_Mean = ndnanfilter(mean(TP_DLW_NH_AA_N70,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_AA_N70_B_Mean = ndnanfilter(mean(TP_DLW_NH_AA_N70_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_AB_N0_Mean = ndnanfilter(mean(TP_DSW_NH_AB_N0,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_AB_N0_B_Mean = ndnanfilter(mean(TP_DSW_NH_AB_N0_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_AB_N0_Mean = ndnanfilter(mean(TP_DLW_NH_AB_N0,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_AB_N0_B_Mean = ndnanfilter(mean(TP_DLW_NH_AB_N0_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_AB_N60_Mean = ndnanfilter(mean(TP_DSW_NH_AB_N60,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_AB_N60_B_Mean = ndnanfilter(mean(TP_DSW_NH_AB_N60_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_AB_N60_Mean = ndnanfilter(mean(TP_DLW_NH_AB_N60,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_AB_N60_B_Mean = ndnanfilter(mean(TP_DLW_NH_AB_N60_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_AB_N70_Mean = ndnanfilter(mean(TP_DSW_NH_AB_N70,'omitnan'),'hamming',FilterWidth);
TP_DSW_NH_AB_N70_B_Mean = ndnanfilter(mean(TP_DSW_NH_AB_N70_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_AB_N70_Mean = ndnanfilter(mean(TP_DLW_NH_AB_N70,'omitnan'),'hamming',FilterWidth);
TP_DLW_NH_AB_N70_B_Mean = ndnanfilter(mean(TP_DLW_NH_AB_N70_B,'omitnan'),'hamming',FilterWidth);

TP_DSW_HI_UN_N0_Mean = ndnanfilter(mean(TP_DSW_HI_UN_N0,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_UN_N0_B_Mean = ndnanfilter(mean(TP_DSW_HI_UN_N0_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_UN_N0_Mean = ndnanfilter(mean(TP_DLW_HI_UN_N0,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_UN_N0_B_Mean = ndnanfilter(mean(TP_DLW_HI_UN_N0_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_UN_N60_Mean = ndnanfilter(mean(TP_DSW_HI_UN_N60,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_UN_N60_B_Mean = ndnanfilter(mean(TP_DSW_HI_UN_N60_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_UN_N60_Mean = ndnanfilter(mean(TP_DLW_HI_UN_N60,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_UN_N60_B_Mean = ndnanfilter(mean(TP_DLW_HI_UN_N60_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_UN_N70_Mean = ndnanfilter(mean(TP_DSW_HI_UN_N70,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_UN_N70_B_Mean = ndnanfilter(mean(TP_DSW_HI_UN_N70_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_UN_N70_Mean = ndnanfilter(mean(TP_DLW_HI_UN_N70,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_UN_N70_B_Mean = ndnanfilter(mean(TP_DLW_HI_UN_N70_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_AA_N0_Mean = ndnanfilter(mean(TP_DSW_HI_AA_N0,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_AA_N0_B_Mean = ndnanfilter(mean(TP_DSW_HI_AA_N0_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_AA_N0_Mean = ndnanfilter(mean(TP_DLW_HI_AA_N0,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_AA_N0_B_Mean = ndnanfilter(mean(TP_DLW_HI_AA_N0_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_AA_N60_Mean = ndnanfilter(mean(TP_DSW_HI_AA_N60,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_AA_N60_B_Mean = ndnanfilter(mean(TP_DSW_HI_AA_N60_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_AA_N60_Mean = ndnanfilter(mean(TP_DLW_HI_AA_N60,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_AA_N60_B_Mean = ndnanfilter(mean(TP_DLW_HI_AA_N60_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_AA_N70_Mean = ndnanfilter(mean(TP_DSW_HI_AA_N70,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_AA_N70_B_Mean = ndnanfilter(mean(TP_DSW_HI_AA_N70_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_AA_N70_Mean = ndnanfilter(mean(TP_DLW_HI_AA_N70,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_AA_N70_B_Mean = ndnanfilter(mean(TP_DLW_HI_AA_N70_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_AB_N0_Mean = ndnanfilter(mean(TP_DSW_HI_AB_N0,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_AB_N0_B_Mean = ndnanfilter(mean(TP_DSW_HI_AB_N0_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_AB_N0_Mean = ndnanfilter(mean(TP_DLW_HI_AB_N0,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_AB_N0_B_Mean = ndnanfilter(mean(TP_DLW_HI_AB_N0_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_AB_N60_Mean = ndnanfilter(mean(TP_DSW_HI_AB_N60,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_AB_N60_B_Mean = ndnanfilter(mean(TP_DSW_HI_AB_N60_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_AB_N60_Mean = ndnanfilter(mean(TP_DLW_HI_AB_N60,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_AB_N60_B_Mean = ndnanfilter(mean(TP_DLW_HI_AB_N60_B,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_AB_N70_Mean = ndnanfilter(mean(TP_DSW_HI_AB_N70,'omitnan'),'hamming',FilterWidth);
TP_DSW_HI_AB_N70_B_Mean = ndnanfilter(mean(TP_DSW_HI_AB_N70_B,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_AB_N70_Mean = ndnanfilter(mean(TP_DLW_HI_AB_N70,'omitnan'),'hamming',FilterWidth);
TP_DLW_HI_AB_N70_B_Mean = ndnanfilter(mean(TP_DLW_HI_AB_N70_B,'omitnan'),'hamming',FilterWidth);

TP_FGSW_Mean = ndnanfilter(mean(TP_FGSW,'omitnan'),'hamming',FilterWidth);
TP_FGLW_Mean = ndnanfilter(mean(TP_FGLW,'omitnan'),'hamming',FilterWidth);

TP_FSW_NH_Mean = ndnanfilter(mean(TP_FSW_NH,'omitnan'),'hamming',FilterWidth);
TP_FLW_NH_Mean = ndnanfilter(mean(TP_FLW_NH,'omitnan'),'hamming',FilterWidth);
TP_FSW_HI_Mean = ndnanfilter(mean(TP_FSW_HI,'omitnan'),'hamming',FilterWidth);
TP_FLW_HI_Mean = ndnanfilter(mean(TP_FLW_HI,'omitnan'),'hamming',FilterWidth);

TP_FSW_N0_Mean = ndnanfilter(mean(TP_FSW_N0,'omitnan'),'hamming',FilterWidth);
TP_FLW_N0_Mean = ndnanfilter(mean(TP_FLW_N0,'omitnan'),'hamming',FilterWidth);
TP_FSW_N60_Mean = ndnanfilter(mean(TP_FSW_N60,'omitnan'),'hamming',FilterWidth);
TP_FLW_N60_Mean = ndnanfilter(mean(TP_FLW_N60,'omitnan'),'hamming',FilterWidth);
TP_FSW_N70_Mean = ndnanfilter(mean(TP_FSW_N70,'omitnan'),'hamming',FilterWidth);
TP_FLW_N70_Mean = ndnanfilter(mean(TP_FLW_N70,'omitnan'),'hamming',FilterWidth);

TP_FSW_NH_N0_Mean = ndnanfilter(mean(TP_FSW_NH_N0,'omitnan'),'hamming',FilterWidth);
TP_FLW_NH_N0_Mean = ndnanfilter(mean(TP_FLW_NH_N0,'omitnan'),'hamming',FilterWidth);
TP_FSW_NH_N60_Mean = ndnanfilter(mean(TP_FSW_NH_N60,'omitnan'),'hamming',FilterWidth);
TP_FLW_NH_N60_Mean = ndnanfilter(mean(TP_FLW_NH_N60,'omitnan'),'hamming',FilterWidth);
TP_FSW_NH_N70_Mean = ndnanfilter(mean(TP_FSW_NH_N70,'omitnan'),'hamming',FilterWidth);
TP_FLW_NH_N70_Mean = ndnanfilter(mean(TP_FLW_NH_N70,'omitnan'),'hamming',FilterWidth);

TP_FSW_HI_N0_Mean = ndnanfilter(mean(TP_FSW_HI_N0,'omitnan'),'hamming',FilterWidth);
TP_FLW_HI_N0_Mean = ndnanfilter(mean(TP_FLW_HI_N0,'omitnan'),'hamming',FilterWidth);
TP_FSW_HI_N60_Mean = ndnanfilter(mean(TP_FSW_HI_N60,'omitnan'),'hamming',FilterWidth);
TP_FLW_HI_N60_Mean = ndnanfilter(mean(TP_FLW_HI_N60,'omitnan'),'hamming',FilterWidth);
TP_FSW_HI_N70_Mean = ndnanfilter(mean(TP_FSW_HI_N70,'omitnan'),'hamming',FilterWidth);
TP_FLW_HI_N70_Mean = ndnanfilter(mean(TP_FLW_HI_N70,'omitnan'),'hamming',FilterWidth);

TP_FSW_NH_UN_Mean = ndnanfilter(mean(TP_FSW_NH_UN,'omitnan'),'hamming',FilterWidth);
TP_FLW_NH_UN_Mean = ndnanfilter(mean(TP_FLW_NH_UN,'omitnan'),'hamming',FilterWidth);
TP_FSW_NH_AA_Mean = ndnanfilter(mean(TP_FSW_NH_AA,'omitnan'),'hamming',FilterWidth);
TP_FLW_NH_AA_Mean = ndnanfilter(mean(TP_FLW_NH_AA,'omitnan'),'hamming',FilterWidth);
TP_FSW_NH_AB_Mean = ndnanfilter(mean(TP_FSW_NH_AB,'omitnan'),'hamming',FilterWidth);
TP_FLW_NH_AB_Mean = ndnanfilter(mean(TP_FLW_NH_AB,'omitnan'),'hamming',FilterWidth);

TP_FSW_HI_UN_Mean = ndnanfilter(mean(TP_FSW_HI_UN,'omitnan'),'hamming',FilterWidth);
TP_FLW_HI_UN_Mean = ndnanfilter(mean(TP_FLW_HI_UN,'omitnan'),'hamming',FilterWidth);
TP_FSW_HI_AA_Mean = ndnanfilter(mean(TP_FSW_HI_AA,'omitnan'),'hamming',FilterWidth);
TP_FLW_HI_AA_Mean = ndnanfilter(mean(TP_FLW_HI_AA,'omitnan'),'hamming',FilterWidth);
TP_FSW_HI_AB_Mean = ndnanfilter(mean(TP_FSW_HI_AB,'omitnan'),'hamming',FilterWidth);
TP_FLW_HI_AB_Mean = ndnanfilter(mean(TP_FLW_HI_AB,'omitnan'),'hamming',FilterWidth);

TP_FSW_NH_UN_N0_Mean = ndnanfilter(mean(TP_FSW_NH_UN_N0,'omitnan'),'hamming',FilterWidth);
TP_FLW_NH_UN_N0_Mean = ndnanfilter(mean(TP_FLW_NH_UN_N0,'omitnan'),'hamming',FilterWidth);
TP_FSW_NH_UN_N60_Mean = ndnanfilter(mean(TP_FSW_NH_UN_N60,'omitnan'),'hamming',FilterWidth);
TP_FLW_NH_UN_N60_Mean = ndnanfilter(mean(TP_FLW_NH_UN_N60,'omitnan'),'hamming',FilterWidth);
TP_FSW_NH_UN_N70_Mean = ndnanfilter(mean(TP_FSW_NH_UN_N70,'omitnan'),'hamming',FilterWidth);
TP_FLW_NH_UN_N70_Mean = ndnanfilter(mean(TP_FLW_NH_UN_N70,'omitnan'),'hamming',FilterWidth);
TP_FSW_NH_AA_N0_Mean = ndnanfilter(mean(TP_FSW_NH_AA_N0,'omitnan'),'hamming',FilterWidth);
TP_FLW_NH_AA_N0_Mean = ndnanfilter(mean(TP_FLW_NH_AA_N0,'omitnan'),'hamming',FilterWidth);
TP_FSW_NH_AA_N60_Mean = ndnanfilter(mean(TP_FSW_NH_AA_N60,'omitnan'),'hamming',FilterWidth);
TP_FLW_NH_AA_N60_Mean = ndnanfilter(mean(TP_FLW_NH_AA_N60,'omitnan'),'hamming',FilterWidth);
TP_FSW_NH_AA_N70_Mean = ndnanfilter(mean(TP_FSW_NH_AA_N70,'omitnan'),'hamming',FilterWidth);
TP_FLW_NH_AA_N70_Mean = ndnanfilter(mean(TP_FLW_NH_AA_N70,'omitnan'),'hamming',FilterWidth);
TP_FSW_NH_AB_N0_Mean = ndnanfilter(mean(TP_FSW_NH_AB_N0,'omitnan'),'hamming',FilterWidth);
TP_FLW_NH_AB_N0_Mean = ndnanfilter(mean(TP_FLW_NH_AB_N0,'omitnan'),'hamming',FilterWidth);
TP_FSW_NH_AB_N60_Mean = ndnanfilter(mean(TP_FSW_NH_AB_N60,'omitnan'),'hamming',FilterWidth);
TP_FLW_NH_AB_N60_Mean = ndnanfilter(mean(TP_FLW_NH_AB_N60,'omitnan'),'hamming',FilterWidth);
TP_FSW_NH_AB_N70_Mean = ndnanfilter(mean(TP_FSW_NH_AB_N70,'omitnan'),'hamming',FilterWidth);
TP_FLW_NH_AB_N70_Mean = ndnanfilter(mean(TP_FLW_NH_AB_N70,'omitnan'),'hamming',FilterWidth);

TP_FSW_HI_UN_N0_Mean = ndnanfilter(mean(TP_FSW_HI_UN_N0,'omitnan'),'hamming',FilterWidth);
TP_FLW_HI_UN_N0_Mean = ndnanfilter(mean(TP_FLW_HI_UN_N0,'omitnan'),'hamming',FilterWidth);
TP_FSW_HI_UN_N60_Mean = ndnanfilter(mean(TP_FSW_HI_UN_N60,'omitnan'),'hamming',FilterWidth);
TP_FLW_HI_UN_N60_Mean = ndnanfilter(mean(TP_FLW_HI_UN_N60,'omitnan'),'hamming',FilterWidth);
TP_FSW_HI_UN_N70_Mean = ndnanfilter(mean(TP_FSW_HI_UN_N70,'omitnan'),'hamming',FilterWidth);
TP_FLW_HI_UN_N70_Mean = ndnanfilter(mean(TP_FLW_HI_UN_N70,'omitnan'),'hamming',FilterWidth);
TP_FSW_HI_AA_N0_Mean = ndnanfilter(mean(TP_FSW_HI_AA_N0,'omitnan'),'hamming',FilterWidth);
TP_FLW_HI_AA_N0_Mean = ndnanfilter(mean(TP_FLW_HI_AA_N0,'omitnan'),'hamming',FilterWidth);
TP_FSW_HI_AA_N60_Mean = ndnanfilter(mean(TP_FSW_HI_AA_N60,'omitnan'),'hamming',FilterWidth);
TP_FLW_HI_AA_N60_Mean = ndnanfilter(mean(TP_FLW_HI_AA_N60,'omitnan'),'hamming',FilterWidth);
TP_FSW_HI_AA_N70_Mean = ndnanfilter(mean(TP_FSW_HI_AA_N70,'omitnan'),'hamming',FilterWidth);
TP_FLW_HI_AA_N70_Mean = ndnanfilter(mean(TP_FLW_HI_AA_N70,'omitnan'),'hamming',FilterWidth);
TP_FSW_HI_AB_N0_Mean = ndnanfilter(mean(TP_FSW_HI_AB_N0,'omitnan'),'hamming',FilterWidth);
TP_FLW_HI_AB_N0_Mean = ndnanfilter(mean(TP_FLW_HI_AB_N0,'omitnan'),'hamming',FilterWidth);
TP_FSW_HI_AB_N60_Mean = ndnanfilter(mean(TP_FSW_HI_AB_N60,'omitnan'),'hamming',FilterWidth);
TP_FLW_HI_AB_N60_Mean = ndnanfilter(mean(TP_FLW_HI_AB_N60,'omitnan'),'hamming',FilterWidth);
TP_FSW_HI_AB_N70_Mean = ndnanfilter(mean(TP_FSW_HI_AB_N70,'omitnan'),'hamming',FilterWidth);
TP_FLW_HI_AB_N70_Mean = ndnanfilter(mean(TP_FLW_HI_AB_N70,'omitnan'),'hamming',FilterWidth);

% Calculate SEM as: std(X)/sqrt(squeeze(sum(~isnan(X),[1 2]))))
TP_DGSW_SEM = (std(TP_DGSW,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DGSW),[1 2])))');
TP_DGSW_B_SEM = (std(TP_DGSW_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DGSW_B),[1 2])))');
TP_DGLW_SEM = (std(TP_DGLW,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DGLW),[1 2])))');
TP_DGLW_B_SEM = (std(TP_DGLW_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DGLW_B),[1 2])))');

TP_DSW_NH_SEM = (std(TP_DSW_NH,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH),[1 2])))');
TP_DSW_NH_B_SEM = (std(TP_DSW_NH_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_B),[1 2])))');
TP_DLW_NH_SEM = (std(TP_DLW_NH,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH),[1 2])))');
TP_DLW_NH_B_SEM = (std(TP_DLW_NH_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_B),[1 2])))');
TP_DSW_HI_SEM = (std(TP_DSW_HI,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI),[1 2])))');
TP_DSW_HI_B_SEM = (std(TP_DSW_HI_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_B),[1 2])))');
TP_DLW_HI_SEM = (std(TP_DLW_HI,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI),[1 2])))');
TP_DLW_HI_B_SEM = (std(TP_DLW_HI_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_B),[1 2])))');

TP_DSW_N0_SEM = (std(TP_DSW_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_N0),[1 2])))');
TP_DSW_N0_B_SEM = (std(TP_DSW_N0_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_N0_B),[1 2])))');
TP_DLW_N0_SEM = (std(TP_DLW_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_N0),[1 2])))');
TP_DLW_N0_B_SEM = (std(TP_DLW_N0_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_N0_B),[1 2])))');
TP_DSW_N60_SEM = (std(TP_DSW_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_N60),[1 2])))');
TP_DSW_N60_B_SEM = (std(TP_DSW_N60_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_N60_B),[1 2])))');
TP_DLW_N60_SEM = (std(TP_DLW_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_N60),[1 2])))');
TP_DLW_N60_B_SEM = (std(TP_DLW_N60_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_N60_B),[1 2])))');
TP_DSW_N70_SEM = (std(TP_DSW_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_N70),[1 2])))');
TP_DSW_N70_B_SEM = (std(TP_DSW_N70_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_N70_B),[1 2])))');
TP_DLW_N70_SEM = (std(TP_DLW_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_N70),[1 2])))');
TP_DLW_N70_B_SEM = (std(TP_DLW_N70_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_N70_B),[1 2])))');

TP_DSW_NH_N0_SEM = (std(TP_DSW_NH_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_N0),[1 2])))');
TP_DSW_NH_N0_B_SEM = (std(TP_DSW_NH_N0_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_N0_B),[1 2])))');
TP_DLW_NH_N0_SEM = (std(TP_DLW_NH_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_N0),[1 2])))');
TP_DLW_NH_N0_B_SEM = (std(TP_DLW_NH_N0_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_N0_B),[1 2])))');
TP_DSW_NH_N60_SEM = (std(TP_DSW_NH_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_N60),[1 2])))');
TP_DSW_NH_N60_B_SEM = (std(TP_DSW_NH_N60_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_N60_B),[1 2])))');
TP_DLW_NH_N60_SEM = (std(TP_DLW_NH_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_N60),[1 2])))');
TP_DLW_NH_N60_B_SEM = (std(TP_DLW_NH_N60_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_N60_B),[1 2])))');
TP_DSW_NH_N70_SEM = (std(TP_DSW_NH_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_N70),[1 2])))');
TP_DSW_NH_N70_B_SEM = (std(TP_DSW_NH_N70_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_N70_B),[1 2])))');
TP_DLW_NH_N70_SEM = (std(TP_DLW_NH_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_N70),[1 2])))');
TP_DLW_NH_N70_B_SEM = (std(TP_DLW_NH_N70_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_N70_B),[1 2])))');

TP_DSW_HI_N0_SEM = (std(TP_DSW_HI_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_N0),[1 2])))');
TP_DSW_HI_N0_B_SEM = (std(TP_DSW_HI_N0_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_N0_B),[1 2])))');
TP_DLW_HI_N0_SEM = (std(TP_DLW_HI_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_N0),[1 2])))');
TP_DLW_HI_N0_B_SEM = (std(TP_DLW_HI_N0_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_N0_B),[1 2])))');
TP_DSW_HI_N60_SEM = (std(TP_DSW_HI_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_N60),[1 2])))');
TP_DSW_HI_N60_B_SEM = (std(TP_DSW_HI_N60_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_N60_B),[1 2])))');
TP_DLW_HI_N60_SEM = (std(TP_DLW_HI_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_N60),[1 2])))');
TP_DLW_HI_N60_B_SEM = (std(TP_DLW_HI_N60_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_N60_B),[1 2])))');
TP_DSW_HI_N70_SEM = (std(TP_DSW_HI_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_N70),[1 2])))');
TP_DSW_HI_N70_B_SEM = (std(TP_DSW_HI_N70_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_N70_B),[1 2])))');
TP_DLW_HI_N70_SEM = (std(TP_DLW_HI_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_N70),[1 2])))');
TP_DLW_HI_N70_B_SEM = (std(TP_DLW_HI_N70_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_N70_B),[1 2])))');

TP_DSW_NH_UN_SEM = (std(TP_DSW_NH_UN,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_UN),[1 2])))');
TP_DSW_NH_UN_B_SEM = (std(TP_DSW_NH_UN_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_UN_B),[1 2])))');
TP_DLW_NH_UN_SEM = (std(TP_DLW_NH_UN,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_UN),[1 2])))');
TP_DLW_NH_UN_B_SEM = (std(TP_DLW_NH_UN_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_UN_B),[1 2])))');
TP_DSW_NH_AA_SEM = (std(TP_DSW_NH_AA,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_AA),[1 2])))');
TP_DSW_NH_AA_B_SEM = (std(TP_DSW_NH_AA_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_AA_B),[1 2])))');
TP_DLW_NH_AA_SEM = (std(TP_DLW_NH_AA,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_AA),[1 2])))');
TP_DLW_NH_AA_B_SEM = (std(TP_DLW_NH_AA_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_AA_B),[1 2])))');
TP_DSW_NH_AB_SEM = (std(TP_DSW_NH_AB,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_AB),[1 2])))');
TP_DSW_NH_AB_B_SEM = (std(TP_DSW_NH_AB_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_AB_B),[1 2])))');
TP_DLW_NH_AB_SEM = (std(TP_DLW_NH_AB,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_AB),[1 2])))');
TP_DLW_NH_AB_B_SEM = (std(TP_DLW_NH_AB_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_AB_B),[1 2])))');

TP_DSW_HI_UN_SEM = (std(TP_DSW_HI_UN,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_UN),[1 2])))');
TP_DSW_HI_UN_B_SEM = (std(TP_DSW_HI_UN_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_UN_B),[1 2])))');
TP_DLW_HI_UN_SEM = (std(TP_DLW_HI_UN,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_UN),[1 2])))');
TP_DLW_HI_UN_B_SEM = (std(TP_DLW_HI_UN_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_UN_B),[1 2])))');
TP_DSW_HI_AA_SEM = (std(TP_DSW_HI_AA,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_AA),[1 2])))');
TP_DSW_HI_AA_B_SEM = (std(TP_DSW_HI_AA_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_AA_B),[1 2])))');
TP_DLW_HI_AA_SEM = (std(TP_DLW_HI_AA,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_AA),[1 2])))');
TP_DLW_HI_AA_B_SEM = (std(TP_DLW_HI_AA_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_AA_B),[1 2])))');
TP_DSW_HI_AB_SEM = (std(TP_DSW_HI_AB,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_AB),[1 2])))');
TP_DSW_HI_AB_B_SEM = (std(TP_DSW_HI_AB_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_AB_B),[1 2])))');
TP_DLW_HI_AB_SEM = (std(TP_DLW_HI_AB,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_AB),[1 2])))');
TP_DLW_HI_AB_B_SEM = (std(TP_DLW_HI_AB_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_AB_B),[1 2])))');

TP_DSW_NH_UN_N0_SEM = (std(TP_DSW_NH_UN_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_UN_N0),[1 2])))');
TP_DSW_NH_UN_N0_B_SEM = (std(TP_DSW_NH_UN_N0_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_UN_N0_B),[1 2])))');
TP_DLW_NH_UN_N0_SEM = (std(TP_DLW_NH_UN_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_UN_N0),[1 2])))');
TP_DLW_NH_UN_N0_B_SEM = (std(TP_DLW_NH_UN_N0_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_UN_N0_B),[1 2])))');
TP_DSW_NH_UN_N60_SEM = (std(TP_DSW_NH_UN_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_UN_N60),[1 2])))');
TP_DSW_NH_UN_N60_B_SEM = (std(TP_DSW_NH_UN_N60_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_UN_N60_B),[1 2])))');
TP_DLW_NH_UN_N60_SEM = (std(TP_DLW_NH_UN_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_UN_N60),[1 2])))');
TP_DLW_NH_UN_N60_B_SEM = (std(TP_DLW_NH_UN_N60_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_UN_N60_B),[1 2])))');
TP_DSW_NH_UN_N70_SEM = (std(TP_DSW_NH_UN_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_UN_N70),[1 2])))');
TP_DSW_NH_UN_N70_B_SEM = (std(TP_DSW_NH_UN_N70_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_UN_N70_B),[1 2])))');
TP_DLW_NH_UN_N70_SEM = (std(TP_DLW_NH_UN_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_UN_N70),[1 2])))');
TP_DLW_NH_UN_N70_B_SEM = (std(TP_DLW_NH_UN_N70_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_UN_N70_B),[1 2])))');
TP_DSW_NH_AA_N0_SEM = (std(TP_DSW_NH_AA_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_AA_N0),[1 2])))');
TP_DSW_NH_AA_N0_B_SEM = (std(TP_DSW_NH_AA_N0_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_AA_N0_B),[1 2])))');
TP_DLW_NH_AA_N0_SEM = (std(TP_DLW_NH_AA_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_AA_N0),[1 2])))');
TP_DLW_NH_AA_N0_B_SEM = (std(TP_DLW_NH_AA_N0_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_AA_N0_B),[1 2])))');
TP_DSW_NH_AA_N60_SEM = (std(TP_DSW_NH_AA_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_AA_N60),[1 2])))');
TP_DSW_NH_AA_N60_B_SEM = (std(TP_DSW_NH_AA_N60_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_AA_N60_B),[1 2])))');
TP_DLW_NH_AA_N60_SEM = (std(TP_DLW_NH_AA_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_AA_N60),[1 2])))');
TP_DLW_NH_AA_N60_B_SEM = (std(TP_DLW_NH_AA_N60_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_AA_N60_B),[1 2])))');
TP_DSW_NH_AA_N70_SEM = (std(TP_DSW_NH_AA_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_AA_N70),[1 2])))');
TP_DSW_NH_AA_N70_B_SEM = (std(TP_DSW_NH_AA_N70_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_AA_N70_B),[1 2])))');
TP_DLW_NH_AA_N70_SEM = (std(TP_DLW_NH_AA_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_AA_N70),[1 2])))');
TP_DLW_NH_AA_N70_B_SEM = (std(TP_DLW_NH_AA_N70_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_AA_N70_B),[1 2])))');
TP_DSW_NH_AB_N0_SEM = (std(TP_DSW_NH_AB_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_AB_N0),[1 2])))');
TP_DSW_NH_AB_N0_B_SEM = (std(TP_DSW_NH_AB_N0_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_AB_N0_B),[1 2])))');
TP_DLW_NH_AB_N0_SEM = (std(TP_DLW_NH_AB_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_AB_N0),[1 2])))');
TP_DLW_NH_AB_N0_B_SEM = (std(TP_DLW_NH_AB_N0_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_AB_N0_B),[1 2])))');
TP_DSW_NH_AB_N60_SEM = (std(TP_DSW_NH_AB_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_AB_N60),[1 2])))');
TP_DSW_NH_AB_N60_B_SEM = (std(TP_DSW_NH_AB_N60_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_AB_N60_B),[1 2])))');
TP_DLW_NH_AB_N60_SEM = (std(TP_DLW_NH_AB_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_AB_N60),[1 2])))');
TP_DLW_NH_AB_N60_B_SEM = (std(TP_DLW_NH_AB_N60_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_AB_N60_B),[1 2])))');
TP_DSW_NH_AB_N70_SEM = (std(TP_DSW_NH_AB_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_AB_N70),[1 2])))');
TP_DSW_NH_AB_N70_B_SEM = (std(TP_DSW_NH_AB_N70_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_NH_AB_N70_B),[1 2])))');
TP_DLW_NH_AB_N70_SEM = (std(TP_DLW_NH_AB_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_AB_N70),[1 2])))');
TP_DLW_NH_AB_N70_B_SEM = (std(TP_DLW_NH_AB_N70_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_NH_AB_N70_B),[1 2])))');

TP_DSW_HI_UN_N0_SEM = (std(TP_DSW_HI_UN_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_UN_N0),[1 2])))');
TP_DSW_HI_UN_N0_B_SEM = (std(TP_DSW_HI_UN_N0_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_UN_N0_B),[1 2])))');
TP_DLW_HI_UN_N0_SEM = (std(TP_DLW_HI_UN_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_UN_N0),[1 2])))');
TP_DLW_HI_UN_N0_B_SEM = (std(TP_DLW_HI_UN_N0_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_UN_N0_B),[1 2])))');
TP_DSW_HI_UN_N60_SEM = (std(TP_DSW_HI_UN_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_UN_N60),[1 2])))');
TP_DSW_HI_UN_N60_B_SEM = (std(TP_DSW_HI_UN_N60_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_UN_N60_B),[1 2])))');
TP_DLW_HI_UN_N60_SEM = (std(TP_DLW_HI_UN_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_UN_N60),[1 2])))');
TP_DLW_HI_UN_N60_B_SEM = (std(TP_DLW_HI_UN_N60_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_UN_N60_B),[1 2])))');
TP_DSW_HI_UN_N70_SEM = (std(TP_DSW_HI_UN_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_UN_N70),[1 2])))');
TP_DSW_HI_UN_N70_B_SEM = (std(TP_DSW_HI_UN_N70_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_UN_N70_B),[1 2])))');
TP_DLW_HI_UN_N70_SEM = (std(TP_DLW_HI_UN_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_UN_N70),[1 2])))');
TP_DLW_HI_UN_N70_B_SEM = (std(TP_DLW_HI_UN_N70_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_UN_N70_B),[1 2])))');
TP_DSW_HI_AA_N0_SEM = (std(TP_DSW_HI_AA_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_AA_N0),[1 2])))');
TP_DSW_HI_AA_N0_B_SEM = (std(TP_DSW_HI_AA_N0_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_AA_N0_B),[1 2])))');
TP_DLW_HI_AA_N0_SEM = (std(TP_DLW_HI_AA_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_AA_N0),[1 2])))');
TP_DLW_HI_AA_N0_B_SEM = (std(TP_DLW_HI_AA_N0_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_AA_N0_B),[1 2])))');
TP_DSW_HI_AA_N60_SEM = (std(TP_DSW_HI_AA_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_AA_N60),[1 2])))');
TP_DSW_HI_AA_N60_B_SEM = (std(TP_DSW_HI_AA_N60_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_AA_N60_B),[1 2])))');
TP_DLW_HI_AA_N60_SEM = (std(TP_DLW_HI_AA_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_AA_N60),[1 2])))');
TP_DLW_HI_AA_N60_B_SEM = (std(TP_DLW_HI_AA_N60_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_AA_N60_B),[1 2])))');
TP_DSW_HI_AA_N70_SEM = (std(TP_DSW_HI_AA_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_AA_N70),[1 2])))');
TP_DSW_HI_AA_N70_B_SEM = (std(TP_DSW_HI_AA_N70_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_AA_N70_B),[1 2])))');
TP_DLW_HI_AA_N70_SEM = (std(TP_DLW_HI_AA_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_AA_N70),[1 2])))');
TP_DLW_HI_AA_N70_B_SEM = (std(TP_DLW_HI_AA_N70_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_AA_N70_B),[1 2])))');
TP_DSW_HI_AB_N0_SEM = (std(TP_DSW_HI_AB_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_AB_N0),[1 2])))');
TP_DSW_HI_AB_N0_B_SEM = (std(TP_DSW_HI_AB_N0_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_AB_N0_B),[1 2])))');
TP_DLW_HI_AB_N0_SEM = (std(TP_DLW_HI_AB_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_AB_N0),[1 2])))');
TP_DLW_HI_AB_N0_B_SEM = (std(TP_DLW_HI_AB_N0_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_AB_N0_B),[1 2])))');
TP_DSW_HI_AB_N60_SEM = (std(TP_DSW_HI_AB_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_AB_N60),[1 2])))');
TP_DSW_HI_AB_N60_B_SEM = (std(TP_DSW_HI_AB_N60_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_AB_N60_B),[1 2])))');
TP_DLW_HI_AB_N60_SEM = (std(TP_DLW_HI_AB_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_AB_N60),[1 2])))');
TP_DLW_HI_AB_N60_B_SEM = (std(TP_DLW_HI_AB_N60_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_AB_N60_B),[1 2])))');
TP_DSW_HI_AB_N70_SEM = (std(TP_DSW_HI_AB_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_AB_N70),[1 2])))');
TP_DSW_HI_AB_N70_B_SEM = (std(TP_DSW_HI_AB_N70_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DSW_HI_AB_N70_B),[1 2])))');
TP_DLW_HI_AB_N70_SEM = (std(TP_DLW_HI_AB_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_AB_N70),[1 2])))');
TP_DLW_HI_AB_N70_B_SEM = (std(TP_DLW_HI_AB_N70_B,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(DLW_HI_AB_N70_B),[1 2])))');

TP_FGSW_SEM = (std(TP_FGSW,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FGSW),[1 2])))');
TP_FGLW_SEM = (std(TP_FGLW,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FGLW),[1 2])))');

TP_FSW_NH_SEM = (std(TP_FSW_NH,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_NH),[1 2])))');
TP_FLW_NH_SEM = (std(TP_FLW_NH,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_NH),[1 2])))');
TP_FSW_HI_SEM = (std(TP_FSW_HI,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_HI),[1 2])))');
TP_FLW_HI_SEM = (std(TP_FLW_HI,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_HI),[1 2])))');

TP_FSW_N0_SEM = (std(TP_FSW_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_N0),[1 2])))');
TP_FLW_N0_SEM = (std(TP_FLW_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_N0),[1 2])))');
TP_FSW_N60_SEM = (std(TP_FSW_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_N60),[1 2])))');
TP_FLW_N60_SEM = (std(TP_FLW_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_N60),[1 2])))');
TP_FSW_N70_SEM = (std(TP_FSW_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_N70),[1 2])))');
TP_FLW_N70_SEM = (std(TP_FLW_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_N70),[1 2])))');

TP_FSW_NH_N0_SEM = (std(TP_FSW_NH_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_NH_N0),[1 2])))');
TP_FLW_NH_N0_SEM = (std(TP_FLW_NH_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_NH_N0),[1 2])))');
TP_FSW_NH_N60_SEM = (std(TP_FSW_NH_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_NH_N60),[1 2])))');
TP_FLW_NH_N60_SEM = (std(TP_FLW_NH_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_NH_N60),[1 2])))');
TP_FSW_NH_N70_SEM = (std(TP_FSW_NH_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_NH_N70),[1 2])))');
TP_FLW_NH_N70_SEM = (std(TP_FLW_NH_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_NH_N70),[1 2])))');

TP_FSW_HI_N0_SEM = (std(TP_FSW_HI_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_HI_N0),[1 2])))');
TP_FLW_HI_N0_SEM = (std(TP_FLW_HI_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_HI_N0),[1 2])))');
TP_FSW_HI_N60_SEM = (std(TP_FSW_HI_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_HI_N60),[1 2])))');
TP_FLW_HI_N60_SEM = (std(TP_FLW_HI_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_HI_N60),[1 2])))');
TP_FSW_HI_N70_SEM = (std(TP_FSW_HI_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_HI_N70),[1 2])))');
TP_FLW_HI_N70_SEM = (std(TP_FLW_HI_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_HI_N70),[1 2])))');

TP_FSW_NH_UN_SEM = (std(TP_FSW_NH_UN,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_NH_UN),[1 2])))');
TP_FLW_NH_UN_SEM = (std(TP_FLW_NH_UN,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_NH_UN),[1 2])))');
TP_FSW_NH_AA_SEM = (std(TP_FSW_NH_AA,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_NH_AA),[1 2])))');
TP_FLW_NH_AA_SEM = (std(TP_FLW_NH_AA,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_NH_AA),[1 2])))');
TP_FSW_NH_AB_SEM = (std(TP_FSW_NH_AB,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_NH_AB),[1 2])))');
TP_FLW_NH_AB_SEM = (std(TP_FLW_NH_AB,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_NH_AB),[1 2])))');

TP_FSW_HI_UN_SEM = (std(TP_FSW_HI_UN,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_HI_UN),[1 2])))');
TP_FLW_HI_UN_SEM = (std(TP_FLW_HI_UN,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_HI_UN),[1 2])))');
TP_FSW_HI_AA_SEM = (std(TP_FSW_HI_AA,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_HI_AA),[1 2])))');
TP_FLW_HI_AA_SEM = (std(TP_FLW_HI_AA,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_HI_AA),[1 2])))');
TP_FSW_HI_AB_SEM = (std(TP_FSW_HI_AB,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_HI_AB),[1 2])))');
TP_FLW_HI_AB_SEM = (std(TP_FLW_HI_AB,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_HI_AB),[1 2])))');

TP_FSW_NH_UN_N0_SEM = (std(TP_FSW_NH_UN_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_NH_UN_N0),[1 2])))');
TP_FLW_NH_UN_N0_SEM = (std(TP_FLW_NH_UN_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_NH_UN_N0),[1 2])))');
TP_FSW_NH_UN_N60_SEM = (std(TP_FSW_NH_UN_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_NH_UN_N60),[1 2])))');
TP_FLW_NH_UN_N60_SEM = (std(TP_FLW_NH_UN_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_NH_UN_N60),[1 2])))');
TP_FSW_NH_UN_N70_SEM = (std(TP_FSW_NH_UN_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_NH_UN_N70),[1 2])))');
TP_FLW_NH_UN_N70_SEM = (std(TP_FLW_NH_UN_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_NH_UN_N70),[1 2])))');
TP_FSW_NH_AA_N0_SEM = (std(TP_FSW_NH_AA_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_NH_AA_N0),[1 2])))');
TP_FLW_NH_AA_N0_SEM = (std(TP_FLW_NH_AA_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_NH_AA_N0),[1 2])))');
TP_FSW_NH_AA_N60_SEM = (std(TP_FSW_NH_AA_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_NH_AA_N60),[1 2])))');
TP_FLW_NH_AA_N60_SEM = (std(TP_FLW_NH_AA_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_NH_AA_N60),[1 2])))');
TP_FSW_NH_AA_N70_SEM = (std(TP_FSW_NH_AA_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_NH_AA_N70),[1 2])))');
TP_FLW_NH_AA_N70_SEM = (std(TP_FLW_NH_AA_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_NH_AA_N70),[1 2])))');
TP_FSW_NH_AB_N0_SEM = (std(TP_FSW_NH_AB_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_NH_AB_N0),[1 2])))');
TP_FLW_NH_AB_N0_SEM = (std(TP_FLW_NH_AB_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_NH_AB_N0),[1 2])))');
TP_FSW_NH_AB_N60_SEM = (std(TP_FSW_NH_AB_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_NH_AB_N60),[1 2])))');
TP_FLW_NH_AB_N60_SEM = (std(TP_FLW_NH_AB_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_NH_AB_N60),[1 2])))');
TP_FSW_NH_AB_N70_SEM = (std(TP_FSW_NH_AB_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_NH_AB_N70),[1 2])))');
TP_FLW_NH_AB_N70_SEM = (std(TP_FLW_NH_AB_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_NH_AB_N70),[1 2])))');

TP_FSW_HI_UN_N0_SEM = (std(TP_FSW_HI_UN_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_HI_UN_N0),[1 2])))');
TP_FLW_HI_UN_N0_SEM = (std(TP_FLW_HI_UN_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_HI_UN_N0),[1 2])))');
TP_FSW_HI_UN_N60_SEM = (std(TP_FSW_HI_UN_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_HI_UN_N60),[1 2])))');
TP_FLW_HI_UN_N60_SEM = (std(TP_FLW_HI_UN_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_HI_UN_N60),[1 2])))');
TP_FSW_HI_UN_N70_SEM = (std(TP_FSW_HI_UN_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_HI_UN_N70),[1 2])))');
TP_FLW_HI_UN_N70_SEM = (std(TP_FLW_HI_UN_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_HI_UN_N70),[1 2])))');
TP_FSW_HI_AA_N0_SEM = (std(TP_FSW_HI_AA_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_HI_AA_N0),[1 2])))');
TP_FLW_HI_AA_N0_SEM = (std(TP_FLW_HI_AA_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_HI_AA_N0),[1 2])))');
TP_FSW_HI_AA_N60_SEM = (std(TP_FSW_HI_AA_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_HI_AA_N60),[1 2])))');
TP_FLW_HI_AA_N60_SEM = (std(TP_FLW_HI_AA_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_HI_AA_N60),[1 2])))');
TP_FSW_HI_AA_N70_SEM = (std(TP_FSW_HI_AA_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_HI_AA_N70),[1 2])))');
TP_FLW_HI_AA_N70_SEM = (std(TP_FLW_HI_AA_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_HI_AA_N70),[1 2])))');
TP_FSW_HI_AB_N0_SEM = (std(TP_FSW_HI_AB_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_HI_AB_N0),[1 2])))');
TP_FLW_HI_AB_N0_SEM = (std(TP_FLW_HI_AB_N0,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_HI_AB_N0),[1 2])))');
TP_FSW_HI_AB_N60_SEM = (std(TP_FSW_HI_AB_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_HI_AB_N60),[1 2])))');
TP_FLW_HI_AB_N60_SEM = (std(TP_FLW_HI_AB_N60,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_HI_AB_N60),[1 2])))');
TP_FSW_HI_AB_N70_SEM = (std(TP_FSW_HI_AB_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FSW_HI_AB_N70),[1 2])))');
TP_FLW_HI_AB_N70_SEM = (std(TP_FLW_HI_AB_N70,0,1,'omitnan')./sqrt(squeeze(sum(~isnan(FLW_HI_AB_N70),[1 2])))');
%% Global Plots
f1 = figure;t1=tiledlayout(1,2);ax1 = nexttile;ax2 = nexttile;
f2 = figure;t2=tiledlayout(1,2);ax3 = nexttile;ax4 = nexttile;
f3 = figure;t3=tiledlayout(1,2);ax5 = nexttile;ax6 = nexttile;
f4 = figure;t4=tiledlayout(1,2);ax7 = nexttile;ax8 = nexttile;
f5 = figure;t5=tiledlayout(1,2);ax9 = nexttile;ax10 = nexttile;
f6 = figure;t6=tiledlayout(1,2);ax11 = nexttile;ax12 = nexttile;
f7 = figure;t7=tiledlayout(1,2);ax13 = nexttile;ax14 = nexttile;
f8 = figure;t8=tiledlayout(1,2);ax15 = nexttile;ax16 = nexttile;
f9 = figure;t9=tiledlayout(1,2);ax17 = nexttile;ax18 = nexttile;
f10 = figure;t10=tiledlayout(1,2);ax19 = nexttile;ax20 = nexttile;
f11 = figure;t11=tiledlayout(1,2);ax21 = nexttile;ax22 = nexttile;
f12 = figure;t12=tiledlayout(1,2);ax23 = nexttile;ax24 = nexttile;
f13 = figure;t13=tiledlayout(1,2);ax25 = nexttile;ax26 = nexttile;
f14 = figure;t14=tiledlayout(1,2);ax27 = nexttile;ax28 = nexttile;
f15 = figure;t15=tiledlayout(1,2);ax29 = nexttile;ax30 = nexttile;
f16 = figure;t16=tiledlayout(1,2);ax31 = nexttile;ax32 = nexttile;
f17 = figure;t17=tiledlayout(1,2);ax33 = nexttile;ax34 = nexttile;
f18 = figure;t18=tiledlayout(1,2);ax35 = nexttile;ax36 = nexttile;

hold([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10 ax11 ax12 ax13 ax14 ax15 ax16 ax17 ax18 ax19 ax20 ax21 ax22 ax23 ax24 ax25 ax26 ax27 ax28 ax29 ax30 ax31 ax32 ax33 ax34 ax35 ax36],'on')
grid([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10 ax11 ax12 ax13 ax14 ax15 ax16 ax17 ax18 ax19 ax20 ax21 ax22 ax23 ax24 ax25 ax26 ax27 ax28 ax29 ax30 ax31 ax32 ax33 ax34 ax35 ax36],'on')

% Plot event onset and baseline markers
xline(ax1,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax2,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax3,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax4,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax5,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax6,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax7,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax8,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax9,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax10,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax11,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax12,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax13,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax14,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax15,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax16,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax17,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax18,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax19,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax20,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax21,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax22,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax23,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax24,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax25,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax26,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax27,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax28,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax29,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax30,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax31,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax32,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax33,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax34,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax35,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')
xline(ax36,0,'--','Event onset','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','handlevisibility','off')

xline(ax3,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
xline(ax4,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
xline(ax7,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
xline(ax8,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
xline(ax11,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
xline(ax12,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
xline(ax15,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
xline(ax16,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
xline(ax19,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
xline(ax20,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
xline(ax23,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')
xline(ax24,-AdapBL,'--','Baseline','LabelVerticalAlignment','top','LabelOrientation','horizontal','handlevisibility','off')

% plot(linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3)),GSW_Mean,color=SColor,linewidth=2);hold on;
% GSW_Mean(isnan(GSW_Mean))=0;GSW_SEM(isnan(GSW_SEM))=0;
% fill([linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3)), flipud(linspace(-TimeStartW,size(GSW,3)/Param.Fs,size(GSW,3))')'],[(GSW_Mean+GSW_SEM), flipud((GSW_Mean-GSW_SEM)')'],SColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')

% Define time vector
t = linspace(-TimeStartW,size(DGSW,3)/Param.Fs,size(DGSW,3));

% Diam
plot(ax1,t,TP_DSW_NH_UN_N0_Mean,color=QuietColor)
plot(ax1,t,TP_DSW_NH_UN_N60_Mean,color=N60Color)
plot(ax1,t,TP_DSW_NH_UN_N70_Mean,color=N70Color)
plot(ax2,t,TP_DLW_NH_UN_N0_Mean,color=QuietColor)
plot(ax2,t,TP_DLW_NH_UN_N60_Mean,color=N60Color)
plot(ax2,t,TP_DLW_NH_UN_N70_Mean,color=N70Color)
plot(ax3,t,TP_DSW_NH_UN_N0_B_Mean,color=QuietColor)
plot(ax3,t,TP_DSW_NH_UN_N60_B_Mean,color=N60Color)
plot(ax3,t,TP_DSW_NH_UN_N70_B_Mean,color=N70Color)
plot(ax4,t,TP_DLW_NH_UN_N0_B_Mean,color=QuietColor)
plot(ax4,t,TP_DLW_NH_UN_N60_B_Mean,color=N60Color)
plot(ax4,t,TP_DLW_NH_UN_N70_B_Mean,color=N70Color)
plot(ax5,t,TP_DSW_HI_UN_N0_Mean,color=QuietColor)
plot(ax5,t,TP_DSW_HI_UN_N60_Mean,color=N60Color)
plot(ax5,t,TP_DSW_HI_UN_N70_Mean,color=N70Color)
plot(ax6,t,TP_DLW_HI_UN_N0_Mean,color=QuietColor)
plot(ax6,t,TP_DLW_HI_UN_N60_Mean,color=N60Color)
plot(ax6,t,TP_DLW_HI_UN_N70_Mean,color=N70Color)
plot(ax7,t,TP_DSW_HI_UN_N0_B_Mean,color=QuietColor)
plot(ax7,t,TP_DSW_HI_UN_N60_B_Mean,color=N60Color)
plot(ax7,t,TP_DSW_HI_UN_N70_B_Mean,color=N70Color)
plot(ax8,t,TP_DLW_HI_UN_N0_B_Mean,color=QuietColor)
plot(ax8,t,TP_DLW_HI_UN_N60_B_Mean,color=N60Color)
plot(ax8,t,TP_DLW_HI_UN_N70_B_Mean,color=N70Color)

plot(ax9,t,TP_DSW_NH_AA_N0_Mean,color=QuietColor)
plot(ax9,t,TP_DSW_NH_AA_N60_Mean,color=N60Color)
plot(ax9,t,TP_DSW_NH_AA_N70_Mean,color=N70Color)
plot(ax10,t,TP_DLW_NH_AA_N0_Mean,color=QuietColor)
plot(ax10,t,TP_DLW_NH_AA_N60_Mean,color=N60Color)
plot(ax10,t,TP_DLW_NH_AA_N70_Mean,color=N70Color)
plot(ax11,t,TP_DSW_NH_AA_N0_B_Mean,color=QuietColor)
plot(ax11,t,TP_DSW_NH_AA_N60_B_Mean,color=N60Color)
plot(ax11,t,TP_DSW_NH_AA_N70_B_Mean,color=N70Color)
plot(ax12,t,TP_DLW_NH_AA_N0_B_Mean,color=QuietColor)
plot(ax12,t,TP_DLW_NH_AA_N60_B_Mean,color=N60Color)
plot(ax12,t,TP_DLW_NH_AA_N70_B_Mean,color=N70Color)
plot(ax13,t,TP_DSW_HI_AA_N0_Mean,color=QuietColor)
plot(ax13,t,TP_DSW_HI_AA_N60_Mean,color=N60Color)
plot(ax13,t,TP_DSW_HI_AA_N70_Mean,color=N70Color)
plot(ax14,t,TP_DLW_HI_AA_N0_Mean,color=QuietColor)
plot(ax14,t,TP_DLW_HI_AA_N60_Mean,color=N60Color)
plot(ax14,t,TP_DLW_HI_AA_N70_Mean,color=N70Color)
plot(ax15,t,TP_DSW_HI_AA_N0_B_Mean,color=QuietColor)
plot(ax15,t,TP_DSW_HI_AA_N60_B_Mean,color=N60Color)
plot(ax15,t,TP_DSW_HI_AA_N70_B_Mean,color=N70Color)
plot(ax16,t,TP_DLW_HI_AA_N0_B_Mean,color=QuietColor)
plot(ax16,t,TP_DLW_HI_AA_N60_B_Mean,color=N60Color)
plot(ax16,t,TP_DLW_HI_AA_N70_B_Mean,color=N70Color)

plot(ax17,t,TP_DSW_NH_AB_N0_Mean,color=QuietColor)
plot(ax17,t,TP_DSW_NH_AB_N60_Mean,color=N60Color)
plot(ax17,t,TP_DSW_NH_AB_N70_Mean,color=N70Color)
plot(ax18,t,TP_DLW_NH_AB_N0_Mean,color=QuietColor)
plot(ax18,t,TP_DLW_NH_AB_N60_Mean,color=N60Color)
plot(ax18,t,TP_DLW_NH_AB_N70_Mean,color=N70Color)
plot(ax19,t,TP_DSW_NH_AB_N0_B_Mean,color=QuietColor)
plot(ax19,t,TP_DSW_NH_AB_N60_B_Mean,color=N60Color)
plot(ax19,t,TP_DSW_NH_AB_N70_B_Mean,color=N70Color)
plot(ax20,t,TP_DLW_NH_AB_N0_B_Mean,color=QuietColor)
plot(ax20,t,TP_DLW_NH_AB_N60_B_Mean,color=N60Color)
plot(ax20,t,TP_DLW_NH_AB_N70_B_Mean,color=N70Color)
plot(ax21,t,TP_DSW_HI_AB_N0_Mean,color=QuietColor)
plot(ax21,t,TP_DSW_HI_AB_N60_Mean,color=N60Color)
plot(ax21,t,TP_DSW_HI_AB_N70_Mean,color=N70Color)
plot(ax22,t,TP_DLW_HI_AB_N0_Mean,color=QuietColor)
plot(ax22,t,TP_DLW_HI_AB_N60_Mean,color=N60Color)
plot(ax22,t,TP_DLW_HI_AB_N70_Mean,color=N70Color)
plot(ax23,t,TP_DSW_HI_AB_N0_B_Mean,color=QuietColor)
plot(ax23,t,TP_DSW_HI_AB_N60_B_Mean,color=N60Color)
plot(ax23,t,TP_DSW_HI_AB_N70_B_Mean,color=N70Color)
plot(ax24,t,TP_DLW_HI_AB_N0_B_Mean,color=QuietColor)
plot(ax24,t,TP_DLW_HI_AB_N60_B_Mean,color=N60Color)
plot(ax24,t,TP_DLW_HI_AB_N70_B_Mean,color=N70Color)

% Fix
plot(ax25,t,TP_FSW_NH_UN_N0_Mean,color=QuietColor)
plot(ax25,t,TP_FSW_NH_UN_N60_Mean,color=N60Color)
plot(ax25,t,TP_FSW_NH_UN_N70_Mean,color=N70Color)
plot(ax26,t,TP_FLW_NH_UN_N0_Mean,color=QuietColor)
plot(ax26,t,TP_FLW_NH_UN_N60_Mean,color=N60Color)
plot(ax26,t,TP_FLW_NH_UN_N70_Mean,color=N70Color)
plot(ax27,t,TP_FSW_HI_UN_N0_Mean,color=QuietColor)
plot(ax27,t,TP_FSW_HI_UN_N60_Mean,color=N60Color)
plot(ax27,t,TP_FSW_HI_UN_N70_Mean,color=N70Color)
plot(ax28,t,TP_FLW_HI_UN_N0_Mean,color=QuietColor)
plot(ax28,t,TP_FLW_HI_UN_N60_Mean,color=N60Color)
plot(ax28,t,TP_FLW_HI_UN_N70_Mean,color=N70Color)

plot(ax29,t,TP_FSW_NH_AA_N0_Mean,color=QuietColor)
plot(ax29,t,TP_FSW_NH_AA_N60_Mean,color=N60Color)
plot(ax29,t,TP_FSW_NH_AA_N70_Mean,color=N70Color)
plot(ax30,t,TP_FLW_NH_AA_N0_Mean,color=QuietColor)
plot(ax30,t,TP_FLW_NH_AA_N60_Mean,color=N60Color)
plot(ax30,t,TP_FLW_NH_AA_N70_Mean,color=N70Color)
plot(ax31,t,TP_FSW_HI_AA_N0_Mean,color=QuietColor)
plot(ax31,t,TP_FSW_HI_AA_N60_Mean,color=N60Color)
plot(ax31,t,TP_FSW_HI_AA_N70_Mean,color=N70Color)
plot(ax32,t,TP_FLW_HI_AA_N0_Mean,color=QuietColor)
plot(ax32,t,TP_FLW_HI_AA_N60_Mean,color=N60Color)
plot(ax32,t,TP_FLW_HI_AA_N70_Mean,color=N70Color)

plot(ax33,t,TP_FSW_NH_AB_N0_Mean,color=QuietColor)
plot(ax33,t,TP_FSW_NH_AB_N60_Mean,color=N60Color)
plot(ax33,t,TP_FSW_NH_AB_N70_Mean,color=N70Color)
plot(ax34,t,TP_FLW_NH_AB_N0_Mean,color=QuietColor)
plot(ax34,t,TP_FLW_NH_AB_N60_Mean,color=N60Color)
plot(ax34,t,TP_FLW_NH_AB_N70_Mean,color=N70Color)
plot(ax35,t,TP_FSW_HI_AB_N0_Mean,color=QuietColor)
plot(ax35,t,TP_FSW_HI_AB_N60_Mean,color=N60Color)
plot(ax35,t,TP_FSW_HI_AB_N70_Mean,color=N70Color)
plot(ax36,t,TP_FLW_HI_AB_N0_Mean,color=QuietColor)
plot(ax36,t,TP_FLW_HI_AB_N60_Mean,color=N60Color)
plot(ax36,t,TP_FLW_HI_AB_N70_Mean,color=N70Color)

% Mean and SEM: Set NaNs to 0 in order to use "fill" correctly
TP_DGSW_Mean(isnan(TP_DGSW_Mean))=0;TP_DGSW_SEM(isnan(TP_DGSW_SEM))=0;
TP_DGLW_Mean(isnan(TP_DGLW_Mean))=0;TP_DGLW_SEM(isnan(TP_DGLW_SEM))=0;
TP_DSW_NH_Mean(isnan(TP_DSW_NH_Mean))=0;TP_DSW_NH_SEM(isnan(TP_DSW_NH_SEM))=0;
TP_DLW_NH_Mean(isnan(TP_DLW_NH_Mean))=0;TP_DLW_NH_SEM(isnan(TP_DLW_NH_SEM))=0;
TP_DSW_N0_Mean(isnan(TP_DSW_N0_Mean))=0;TP_DSW_N0_SEM(isnan(TP_DSW_N0_SEM))=0;
TP_DSW_N60_Mean(isnan(TP_DSW_N60_Mean))=0;TP_DSW_N60_SEM(isnan(TP_DSW_N60_SEM))=0;
TP_DSW_N70_Mean(isnan(TP_DSW_N70_Mean))=0;TP_DSW_N70_SEM(isnan(TP_DSW_N70_SEM))=0;
TP_DLW_N0_Mean(isnan(TP_DLW_N0_Mean))=0;TP_DLW_N0_SEM(isnan(TP_DLW_N0_SEM))=0;
TP_DLW_N60_Mean(isnan(TP_DLW_N60_Mean))=0;TP_DLW_N60_SEM(isnan(TP_DLW_N60_SEM))=0;
TP_DLW_N70_Mean(isnan(TP_DLW_N70_Mean))=0;TP_DLW_N70_SEM(isnan(TP_DLW_N70_SEM))=0;
TP_DSW_NH_N0_Mean(isnan(TP_DSW_NH_N0_Mean))=0;TP_DSW_NH_N0_SEM(isnan(TP_DSW_NH_N0_SEM))=0;
TP_DSW_NH_N60_Mean(isnan(TP_DSW_NH_N60_Mean))=0;TP_DSW_NH_N60_SEM(isnan(TP_DSW_NH_N60_SEM))=0;
TP_DSW_NH_N70_Mean(isnan(TP_DSW_NH_N70_Mean))=0;TP_DSW_NH_N70_SEM(isnan(TP_DSW_NH_N70_SEM))=0;
TP_DLW_NH_N0_Mean(isnan(TP_DLW_NH_N0_Mean))=0;TP_DLW_NH_N0_SEM(isnan(TP_DLW_NH_N0_SEM))=0;
TP_DLW_NH_N60_Mean(isnan(TP_DLW_NH_N60_Mean))=0;TP_DLW_NH_N60_SEM(isnan(TP_DLW_NH_N60_SEM))=0;
TP_DLW_NH_N70_Mean(isnan(TP_DLW_NH_N70_Mean))=0;TP_DLW_NH_N70_SEM(isnan(TP_DLW_NH_N70_SEM))=0;
TP_DSW_HI_Mean(isnan(TP_DSW_HI_Mean))=0;TP_DSW_HI_SEM(isnan(TP_DSW_HI_SEM))=0;
TP_DLW_HI_Mean(isnan(TP_DLW_HI_Mean))=0;TP_DLW_HI_SEM(isnan(TP_DLW_HI_SEM))=0;
TP_DSW_HI_N0_Mean(isnan(TP_DSW_HI_N0_Mean))=0;TP_DSW_HI_N0_SEM(isnan(TP_DSW_HI_N0_SEM))=0;
TP_DSW_HI_N60_Mean(isnan(TP_DSW_HI_N60_Mean))=0;TP_DSW_HI_N60_SEM(isnan(TP_DSW_HI_N60_SEM))=0;
TP_DSW_HI_N70_Mean(isnan(TP_DSW_HI_N70_Mean))=0;TP_DSW_HI_N70_SEM(isnan(TP_DSW_HI_N70_SEM))=0;
TP_DLW_HI_N0_Mean(isnan(TP_DLW_HI_N0_Mean))=0;TP_DLW_HI_N0_SEM(isnan(TP_DLW_HI_N0_SEM))=0;
TP_DLW_HI_N60_Mean(isnan(TP_DLW_HI_N60_Mean))=0;TP_DLW_HI_N60_SEM(isnan(TP_DLW_HI_N60_SEM))=0;
TP_DLW_HI_N70_Mean(isnan(TP_DLW_HI_N70_Mean))=0;TP_DLW_HI_N70_SEM(isnan(TP_DLW_HI_N70_SEM))=0;
TP_DSW_NH_UN_Mean(isnan(TP_DSW_NH_UN_Mean))=0;TP_DSW_NH_UN_SEM(isnan(TP_DSW_NH_UN_SEM))=0;
TP_DLW_NH_UN_Mean(isnan(TP_DLW_NH_UN_Mean))=0;TP_DLW_NH_UN_SEM(isnan(TP_DLW_NH_UN_SEM))=0;
TP_DSW_NH_AA_Mean(isnan(TP_DSW_NH_AA_Mean))=0;TP_DSW_NH_AA_SEM(isnan(TP_DSW_NH_AA_SEM))=0;
TP_DLW_NH_AA_Mean(isnan(TP_DLW_NH_AA_Mean))=0;TP_DLW_NH_AA_SEM(isnan(TP_DLW_NH_AA_SEM))=0;
TP_DSW_NH_AB_Mean(isnan(TP_DSW_NH_AB_Mean))=0;TP_DSW_NH_AB_SEM(isnan(TP_DSW_NH_AB_SEM))=0;
TP_DLW_NH_AB_Mean(isnan(TP_DLW_NH_AB_Mean))=0;TP_DLW_NH_AB_SEM(isnan(TP_DLW_NH_AB_SEM))=0;
TP_DSW_HI_UN_Mean(isnan(TP_DSW_HI_UN_Mean))=0;TP_DSW_HI_UN_SEM(isnan(TP_DSW_HI_UN_SEM))=0;
TP_DLW_HI_UN_Mean(isnan(TP_DLW_HI_UN_Mean))=0;TP_DLW_HI_UN_SEM(isnan(TP_DLW_HI_UN_SEM))=0;
TP_DSW_HI_AA_Mean(isnan(TP_DSW_HI_AA_Mean))=0;TP_DSW_HI_AA_SEM(isnan(TP_DSW_HI_AA_SEM))=0;
TP_DLW_HI_AA_Mean(isnan(TP_DLW_HI_AA_Mean))=0;TP_DLW_HI_AA_SEM(isnan(TP_DLW_HI_AA_SEM))=0;
TP_DSW_HI_AB_Mean(isnan(TP_DSW_HI_AB_Mean))=0;TP_DSW_HI_AB_SEM(isnan(TP_DSW_HI_AB_SEM))=0;
TP_DLW_HI_AB_Mean(isnan(TP_DLW_HI_AB_Mean))=0;TP_DLW_HI_AB_SEM(isnan(TP_DLW_HI_AB_SEM))=0;
TP_DSW_NH_UN_N0_Mean(isnan(TP_DSW_NH_UN_N0_Mean))=0;TP_DSW_NH_UN_N0_SEM(isnan(TP_DSW_NH_UN_N0_SEM))=0;
TP_DSW_NH_UN_N60_Mean(isnan(TP_DSW_NH_UN_N60_Mean))=0;TP_DSW_NH_UN_N60_SEM(isnan(TP_DSW_NH_UN_N60_SEM))=0;
TP_DSW_NH_UN_N70_Mean(isnan(TP_DSW_NH_UN_N70_Mean))=0;TP_DSW_NH_UN_N70_SEM(isnan(TP_DSW_NH_UN_N70_SEM))=0;
TP_DLW_NH_UN_N0_Mean(isnan(TP_DLW_NH_UN_N0_Mean))=0;TP_DLW_NH_UN_N0_SEM(isnan(TP_DLW_NH_UN_N0_SEM))=0;
TP_DLW_NH_UN_N60_Mean(isnan(TP_DLW_NH_UN_N60_Mean))=0;TP_DLW_NH_UN_N60_SEM(isnan(TP_DLW_NH_UN_N60_SEM))=0;
TP_DLW_NH_UN_N70_Mean(isnan(TP_DLW_NH_UN_N70_Mean))=0;TP_DLW_NH_UN_N70_SEM(isnan(TP_DLW_NH_UN_N70_SEM))=0;
TP_DSW_NH_AA_N0_Mean(isnan(TP_DSW_NH_AA_N0_Mean))=0;TP_DSW_NH_AA_N0_SEM(isnan(TP_DSW_NH_AA_N0_SEM))=0;
TP_DSW_NH_AA_N60_Mean(isnan(TP_DSW_NH_AA_N60_Mean))=0;TP_DSW_NH_AA_N60_SEM(isnan(TP_DSW_NH_AA_N60_SEM))=0;
TP_DSW_NH_AA_N70_Mean(isnan(TP_DSW_NH_AA_N70_Mean))=0;TP_DSW_NH_AA_N70_SEM(isnan(TP_DSW_NH_AA_N70_SEM))=0;
TP_DLW_NH_AA_N0_Mean(isnan(TP_DLW_NH_AA_N0_Mean))=0;TP_DLW_NH_AA_N0_SEM(isnan(TP_DLW_NH_AA_N0_SEM))=0;
TP_DLW_NH_AA_N60_Mean(isnan(TP_DLW_NH_AA_N60_Mean))=0;TP_DLW_NH_AA_N60_SEM(isnan(TP_DLW_NH_AA_N60_SEM))=0;
TP_DLW_NH_AA_N70_Mean(isnan(TP_DLW_NH_AA_N70_Mean))=0;TP_DLW_NH_AA_N70_SEM(isnan(TP_DLW_NH_AA_N70_SEM))=0;
TP_DSW_NH_AB_N0_Mean(isnan(TP_DSW_NH_AB_N0_Mean))=0;TP_DSW_NH_AB_N0_SEM(isnan(TP_DSW_NH_AB_N0_SEM))=0;
TP_DSW_NH_AB_N60_Mean(isnan(TP_DSW_NH_AB_N60_Mean))=0;TP_DSW_NH_AB_N60_SEM(isnan(TP_DSW_NH_AB_N60_SEM))=0;
TP_DSW_NH_AB_N70_Mean(isnan(TP_DSW_NH_AB_N70_Mean))=0;TP_DSW_NH_AB_N70_SEM(isnan(TP_DSW_NH_AB_N70_SEM))=0;
TP_DLW_NH_AB_N0_Mean(isnan(TP_DLW_NH_AB_N0_Mean))=0;TP_DLW_NH_AB_N0_SEM(isnan(TP_DLW_NH_AB_N0_SEM))=0;
TP_DLW_NH_AB_N60_Mean(isnan(TP_DLW_NH_AB_N60_Mean))=0;TP_DLW_NH_AB_N60_SEM(isnan(TP_DLW_NH_AB_N60_SEM))=0;
TP_DLW_NH_AB_N70_Mean(isnan(TP_DLW_NH_AB_N70_Mean))=0;TP_DLW_NH_AB_N70_SEM(isnan(TP_DLW_NH_AB_N70_SEM))=0;
TP_DSW_HI_UN_N0_Mean(isnan(TP_DSW_HI_UN_N0_Mean))=0;TP_DSW_HI_UN_N0_SEM(isnan(TP_DSW_HI_UN_N0_SEM))=0;
TP_DSW_HI_UN_N60_Mean(isnan(TP_DSW_HI_UN_N60_Mean))=0;TP_DSW_HI_UN_N60_SEM(isnan(TP_DSW_HI_UN_N60_SEM))=0;
TP_DSW_HI_UN_N70_Mean(isnan(TP_DSW_HI_UN_N70_Mean))=0;TP_DSW_HI_UN_N70_SEM(isnan(TP_DSW_HI_UN_N70_SEM))=0;
TP_DLW_HI_UN_N0_Mean(isnan(TP_DLW_HI_UN_N0_Mean))=0;TP_DLW_HI_UN_N0_SEM(isnan(TP_DLW_HI_UN_N0_SEM))=0;
TP_DLW_HI_UN_N60_Mean(isnan(TP_DLW_HI_UN_N60_Mean))=0;TP_DLW_HI_UN_N60_SEM(isnan(TP_DLW_HI_UN_N60_SEM))=0;
TP_DLW_HI_UN_N70_Mean(isnan(TP_DLW_HI_UN_N70_Mean))=0;TP_DLW_HI_UN_N70_SEM(isnan(TP_DLW_HI_UN_N70_SEM))=0;
TP_DSW_HI_AA_N0_Mean(isnan(TP_DSW_HI_AA_N0_Mean))=0;TP_DSW_HI_AA_N0_SEM(isnan(TP_DSW_HI_AA_N0_SEM))=0;
TP_DSW_HI_AA_N60_Mean(isnan(TP_DSW_HI_AA_N60_Mean))=0;TP_DSW_HI_AA_N60_SEM(isnan(TP_DSW_HI_AA_N60_SEM))=0;
TP_DSW_HI_AA_N70_Mean(isnan(TP_DSW_HI_AA_N70_Mean))=0;TP_DSW_HI_AA_N70_SEM(isnan(TP_DSW_HI_AA_N70_SEM))=0;
TP_DLW_HI_AA_N0_Mean(isnan(TP_DLW_HI_AA_N0_Mean))=0;TP_DLW_HI_AA_N0_SEM(isnan(TP_DLW_HI_AA_N0_SEM))=0;
TP_DLW_HI_AA_N60_Mean(isnan(TP_DLW_HI_AA_N60_Mean))=0;TP_DLW_HI_AA_N60_SEM(isnan(TP_DLW_HI_AA_N60_SEM))=0;
TP_DLW_HI_AA_N70_Mean(isnan(TP_DLW_HI_AA_N70_Mean))=0;TP_DLW_HI_AA_N70_SEM(isnan(TP_DLW_HI_AA_N70_SEM))=0;
TP_DSW_HI_AB_N0_Mean(isnan(TP_DSW_HI_AB_N0_Mean))=0;TP_DSW_HI_AB_N0_SEM(isnan(TP_DSW_HI_AB_N0_SEM))=0;
TP_DSW_HI_AB_N60_Mean(isnan(TP_DSW_HI_AB_N60_Mean))=0;TP_DSW_HI_AB_N60_SEM(isnan(TP_DSW_HI_AB_N60_SEM))=0;
TP_DSW_HI_AB_N70_Mean(isnan(TP_DSW_HI_AB_N70_Mean))=0;TP_DSW_HI_AB_N70_SEM(isnan(TP_DSW_HI_AB_N70_SEM))=0;
TP_DLW_HI_AB_N0_Mean(isnan(TP_DLW_HI_AB_N0_Mean))=0;TP_DLW_HI_AB_N0_SEM(isnan(TP_DLW_HI_AB_N0_SEM))=0;
TP_DLW_HI_AB_N60_Mean(isnan(TP_DLW_HI_AB_N60_Mean))=0;TP_DLW_HI_AB_N60_SEM(isnan(TP_DLW_HI_AB_N60_SEM))=0;
TP_DLW_HI_AB_N70_Mean(isnan(TP_DLW_HI_AB_N70_Mean))=0;TP_DLW_HI_AB_N70_SEM(isnan(TP_DLW_HI_AB_N70_SEM))=0;

TP_DGSW_B_Mean(isnan(TP_DGSW_B_Mean))=0;TP_DGSW_B_SEM(isnan(TP_DGSW_B_SEM))=0;
TP_DGLW_B_Mean(isnan(TP_DGLW_B_Mean))=0;TP_DGLW_B_SEM(isnan(TP_DGLW_B_SEM))=0;
TP_DSW_NH_B_Mean(isnan(TP_DSW_NH_B_Mean))=0;TP_DSW_NH_B_SEM(isnan(TP_DSW_NH_B_SEM))=0;
TP_DLW_NH_B_Mean(isnan(TP_DLW_NH_B_Mean))=0;TP_DLW_NH_B_SEM(isnan(TP_DLW_NH_B_SEM))=0;
TP_DSW_N0_B_Mean(isnan(TP_DSW_N0_B_Mean))=0;TP_DSW_N0_B_SEM(isnan(TP_DSW_N0_B_SEM))=0;
TP_DSW_N60_B_Mean(isnan(TP_DSW_N60_B_Mean))=0;TP_DSW_N60_B_SEM(isnan(TP_DSW_N60_B_SEM))=0;
TP_DSW_N70_B_Mean(isnan(TP_DSW_N70_B_Mean))=0;TP_DSW_N70_B_SEM(isnan(TP_DSW_N70_B_SEM))=0;
TP_DLW_N0_B_Mean(isnan(TP_DLW_N0_B_Mean))=0;TP_DLW_N0_B_SEM(isnan(TP_DLW_N0_B_SEM))=0;
TP_DLW_N60_B_Mean(isnan(TP_DLW_N60_B_Mean))=0;TP_DLW_N60_B_SEM(isnan(TP_DLW_N60_B_SEM))=0;
TP_DLW_N70_B_Mean(isnan(TP_DLW_N70_B_Mean))=0;TP_DLW_N70_B_SEM(isnan(TP_DLW_N70_B_SEM))=0;
TP_DSW_NH_N0_B_Mean(isnan(TP_DSW_NH_N0_B_Mean))=0;TP_DSW_NH_N0_B_SEM(isnan(TP_DSW_NH_N0_B_SEM))=0;
TP_DSW_NH_N60_B_Mean(isnan(TP_DSW_NH_N60_B_Mean))=0;TP_DSW_NH_N60_B_SEM(isnan(TP_DSW_NH_N60_B_SEM))=0;
TP_DSW_NH_N70_B_Mean(isnan(TP_DSW_NH_N70_B_Mean))=0;TP_DSW_NH_N70_B_SEM(isnan(TP_DSW_NH_N70_B_SEM))=0;
TP_DLW_NH_N0_B_Mean(isnan(TP_DLW_NH_N0_B_Mean))=0;TP_DLW_NH_N0_B_SEM(isnan(TP_DLW_NH_N0_B_SEM))=0;
TP_DLW_NH_N60_B_Mean(isnan(TP_DLW_NH_N60_B_Mean))=0;TP_DLW_NH_N60_B_SEM(isnan(TP_DLW_NH_N60_B_SEM))=0;
TP_DLW_NH_N70_B_Mean(isnan(TP_DLW_NH_N70_B_Mean))=0;TP_DLW_NH_N70_B_SEM(isnan(TP_DLW_NH_N70_B_SEM))=0;
TP_DSW_HI_B_Mean(isnan(TP_DSW_HI_B_Mean))=0;TP_DSW_HI_B_SEM(isnan(TP_DSW_HI_B_SEM))=0;
TP_DLW_HI_B_Mean(isnan(TP_DLW_HI_B_Mean))=0;TP_DLW_HI_B_SEM(isnan(TP_DLW_HI_B_SEM))=0;
TP_DSW_HI_N0_B_Mean(isnan(TP_DSW_HI_N0_B_Mean))=0;TP_DSW_HI_N0_B_SEM(isnan(TP_DSW_HI_N0_B_SEM))=0;
TP_DSW_HI_N60_B_Mean(isnan(TP_DSW_HI_N60_B_Mean))=0;TP_DSW_HI_N60_B_SEM(isnan(TP_DSW_HI_N60_B_SEM))=0;
TP_DSW_HI_N70_B_Mean(isnan(TP_DSW_HI_N70_B_Mean))=0;TP_DSW_HI_N70_B_SEM(isnan(TP_DSW_HI_N70_B_SEM))=0;
TP_DLW_HI_N0_B_Mean(isnan(TP_DLW_HI_N0_B_Mean))=0;TP_DLW_HI_N0_B_SEM(isnan(TP_DLW_HI_N0_B_SEM))=0;
TP_DLW_HI_N60_B_Mean(isnan(TP_DLW_HI_N60_B_Mean))=0;TP_DLW_HI_N60_B_SEM(isnan(TP_DLW_HI_N60_B_SEM))=0;
TP_DLW_HI_N70_B_Mean(isnan(TP_DLW_HI_N70_B_Mean))=0;TP_DLW_HI_N70_B_SEM(isnan(TP_DLW_HI_N70_B_SEM))=0;
TP_DSW_NH_UN_B_Mean(isnan(TP_DSW_NH_UN_B_Mean))=0;TP_DSW_NH_UN_B_SEM(isnan(TP_DSW_NH_UN_B_SEM))=0;
TP_DLW_NH_UN_B_Mean(isnan(TP_DLW_NH_UN_B_Mean))=0;TP_DLW_NH_UN_B_SEM(isnan(TP_DLW_NH_UN_B_SEM))=0;
TP_DSW_NH_AA_B_Mean(isnan(TP_DSW_NH_AA_B_Mean))=0;TP_DSW_NH_AA_B_SEM(isnan(TP_DSW_NH_AA_B_SEM))=0;
TP_DLW_NH_AA_B_Mean(isnan(TP_DLW_NH_AA_B_Mean))=0;TP_DLW_NH_AA_B_SEM(isnan(TP_DLW_NH_AA_B_SEM))=0;
TP_DSW_NH_AB_B_Mean(isnan(TP_DSW_NH_AB_B_Mean))=0;TP_DSW_NH_AB_B_SEM(isnan(TP_DSW_NH_AB_B_SEM))=0;
TP_DLW_NH_AB_B_Mean(isnan(TP_DLW_NH_AB_B_Mean))=0;TP_DLW_NH_AB_B_SEM(isnan(TP_DLW_NH_AB_B_SEM))=0;
TP_DSW_HI_UN_B_Mean(isnan(TP_DSW_HI_UN_B_Mean))=0;TP_DSW_HI_UN_B_SEM(isnan(TP_DSW_HI_UN_B_SEM))=0;
TP_DLW_HI_UN_B_Mean(isnan(TP_DLW_HI_UN_B_Mean))=0;TP_DLW_HI_UN_B_SEM(isnan(TP_DLW_HI_UN_B_SEM))=0;
TP_DSW_HI_AA_B_Mean(isnan(TP_DSW_HI_AA_B_Mean))=0;TP_DSW_HI_AA_B_SEM(isnan(TP_DSW_HI_AA_B_SEM))=0;
TP_DLW_HI_AA_B_Mean(isnan(TP_DLW_HI_AA_B_Mean))=0;TP_DLW_HI_AA_B_SEM(isnan(TP_DLW_HI_AA_B_SEM))=0;
TP_DSW_HI_AB_B_Mean(isnan(TP_DSW_HI_AB_B_Mean))=0;TP_DSW_HI_AB_B_SEM(isnan(TP_DSW_HI_AB_B_SEM))=0;
TP_DLW_HI_AB_B_Mean(isnan(TP_DLW_HI_AB_B_Mean))=0;TP_DLW_HI_AB_B_SEM(isnan(TP_DLW_HI_AB_B_SEM))=0;
TP_DSW_NH_UN_N0_B_Mean(isnan(TP_DSW_NH_UN_N0_B_Mean))=0;TP_DSW_NH_UN_N0_B_SEM(isnan(TP_DSW_NH_UN_N0_B_SEM))=0;
TP_DSW_NH_UN_N60_B_Mean(isnan(TP_DSW_NH_UN_N60_B_Mean))=0;TP_DSW_NH_UN_N60_B_SEM(isnan(TP_DSW_NH_UN_N60_B_SEM))=0;
TP_DSW_NH_UN_N70_B_Mean(isnan(TP_DSW_NH_UN_N70_B_Mean))=0;TP_DSW_NH_UN_N70_B_SEM(isnan(TP_DSW_NH_UN_N70_B_SEM))=0;
TP_DLW_NH_UN_N0_B_Mean(isnan(TP_DLW_NH_UN_N0_B_Mean))=0;TP_DLW_NH_UN_N0_B_SEM(isnan(TP_DLW_NH_UN_N0_B_SEM))=0;
TP_DLW_NH_UN_N60_B_Mean(isnan(TP_DLW_NH_UN_N60_B_Mean))=0;TP_DLW_NH_UN_N60_B_SEM(isnan(TP_DLW_NH_UN_N60_B_SEM))=0;
TP_DLW_NH_UN_N70_B_Mean(isnan(TP_DLW_NH_UN_N70_B_Mean))=0;TP_DLW_NH_UN_N70_B_SEM(isnan(TP_DLW_NH_UN_N70_B_SEM))=0;
TP_DSW_NH_AA_N0_B_Mean(isnan(TP_DSW_NH_AA_N0_B_Mean))=0;TP_DSW_NH_AA_N0_B_SEM(isnan(TP_DSW_NH_AA_N0_B_SEM))=0;
TP_DSW_NH_AA_N60_B_Mean(isnan(TP_DSW_NH_AA_N60_B_Mean))=0;TP_DSW_NH_AA_N60_B_SEM(isnan(TP_DSW_NH_AA_N60_B_SEM))=0;
TP_DSW_NH_AA_N70_B_Mean(isnan(TP_DSW_NH_AA_N70_B_Mean))=0;TP_DSW_NH_AA_N70_B_SEM(isnan(TP_DSW_NH_AA_N70_B_SEM))=0;
TP_DLW_NH_AA_N0_B_Mean(isnan(TP_DLW_NH_AA_N0_B_Mean))=0;TP_DLW_NH_AA_N0_B_SEM(isnan(TP_DLW_NH_AA_N0_B_SEM))=0;
TP_DLW_NH_AA_N60_B_Mean(isnan(TP_DLW_NH_AA_N60_B_Mean))=0;TP_DLW_NH_AA_N60_B_SEM(isnan(TP_DLW_NH_AA_N60_B_SEM))=0;
TP_DLW_NH_AA_N70_B_Mean(isnan(TP_DLW_NH_AA_N70_B_Mean))=0;TP_DLW_NH_AA_N70_B_SEM(isnan(TP_DLW_NH_AA_N70_B_SEM))=0;
TP_DSW_NH_AB_N0_B_Mean(isnan(TP_DSW_NH_AB_N0_B_Mean))=0;TP_DSW_NH_AB_N0_B_SEM(isnan(TP_DSW_NH_AB_N0_B_SEM))=0;
TP_DSW_NH_AB_N60_B_Mean(isnan(TP_DSW_NH_AB_N60_B_Mean))=0;TP_DSW_NH_AB_N60_B_SEM(isnan(TP_DSW_NH_AB_N60_B_SEM))=0;
TP_DSW_NH_AB_N70_B_Mean(isnan(TP_DSW_NH_AB_N70_B_Mean))=0;TP_DSW_NH_AB_N70_B_SEM(isnan(TP_DSW_NH_AB_N70_B_SEM))=0;
TP_DLW_NH_AB_N0_B_Mean(isnan(TP_DLW_NH_AB_N0_B_Mean))=0;TP_DLW_NH_AB_N0_B_SEM(isnan(TP_DLW_NH_AB_N0_B_SEM))=0;
TP_DLW_NH_AB_N60_B_Mean(isnan(TP_DLW_NH_AB_N60_B_Mean))=0;TP_DLW_NH_AB_N60_B_SEM(isnan(TP_DLW_NH_AB_N60_B_SEM))=0;
TP_DLW_NH_AB_N70_B_Mean(isnan(TP_DLW_NH_AB_N70_B_Mean))=0;TP_DLW_NH_AB_N70_B_SEM(isnan(TP_DLW_NH_AB_N70_B_SEM))=0;
TP_DSW_HI_UN_N0_B_Mean(isnan(TP_DSW_HI_UN_N0_B_Mean))=0;TP_DSW_HI_UN_N0_B_SEM(isnan(TP_DSW_HI_UN_N0_B_SEM))=0;
TP_DSW_HI_UN_N60_B_Mean(isnan(TP_DSW_HI_UN_N60_B_Mean))=0;TP_DSW_HI_UN_N60_B_SEM(isnan(TP_DSW_HI_UN_N60_B_SEM))=0;
TP_DSW_HI_UN_N70_B_Mean(isnan(TP_DSW_HI_UN_N70_B_Mean))=0;TP_DSW_HI_UN_N70_B_SEM(isnan(TP_DSW_HI_UN_N70_B_SEM))=0;
TP_DLW_HI_UN_N0_B_Mean(isnan(TP_DLW_HI_UN_N0_B_Mean))=0;TP_DLW_HI_UN_N0_B_SEM(isnan(TP_DLW_HI_UN_N0_B_SEM))=0;
TP_DLW_HI_UN_N60_B_Mean(isnan(TP_DLW_HI_UN_N60_B_Mean))=0;TP_DLW_HI_UN_N60_B_SEM(isnan(TP_DLW_HI_UN_N60_B_SEM))=0;
TP_DLW_HI_UN_N70_B_Mean(isnan(TP_DLW_HI_UN_N70_B_Mean))=0;TP_DLW_HI_UN_N70_B_SEM(isnan(TP_DLW_HI_UN_N70_B_SEM))=0;
TP_DSW_HI_AA_N0_B_Mean(isnan(TP_DSW_HI_AA_N0_B_Mean))=0;TP_DSW_HI_AA_N0_B_SEM(isnan(TP_DSW_HI_AA_N0_B_SEM))=0;
TP_DSW_HI_AA_N60_B_Mean(isnan(TP_DSW_HI_AA_N60_B_Mean))=0;TP_DSW_HI_AA_N60_B_SEM(isnan(TP_DSW_HI_AA_N60_B_SEM))=0;
TP_DSW_HI_AA_N70_B_Mean(isnan(TP_DSW_HI_AA_N70_B_Mean))=0;TP_DSW_HI_AA_N70_B_SEM(isnan(TP_DSW_HI_AA_N70_B_SEM))=0;
TP_DLW_HI_AA_N0_B_Mean(isnan(TP_DLW_HI_AA_N0_B_Mean))=0;TP_DLW_HI_AA_N0_B_SEM(isnan(TP_DLW_HI_AA_N0_B_SEM))=0;
TP_DLW_HI_AA_N60_B_Mean(isnan(TP_DLW_HI_AA_N60_B_Mean))=0;TP_DLW_HI_AA_N60_B_SEM(isnan(TP_DLW_HI_AA_N60_B_SEM))=0;
TP_DLW_HI_AA_N70_B_Mean(isnan(TP_DLW_HI_AA_N70_B_Mean))=0;TP_DLW_HI_AA_N70_B_SEM(isnan(TP_DLW_HI_AA_N70_B_SEM))=0;
TP_DSW_HI_AB_N0_B_Mean(isnan(TP_DSW_HI_AB_N0_B_Mean))=0;TP_DSW_HI_AB_N0_B_SEM(isnan(TP_DSW_HI_AB_N0_B_SEM))=0;
TP_DSW_HI_AB_N60_B_Mean(isnan(TP_DSW_HI_AB_N60_B_Mean))=0;TP_DSW_HI_AB_N60_B_SEM(isnan(TP_DSW_HI_AB_N60_B_SEM))=0;
TP_DSW_HI_AB_N70_B_Mean(isnan(TP_DSW_HI_AB_N70_B_Mean))=0;TP_DSW_HI_AB_N70_B_SEM(isnan(TP_DSW_HI_AB_N70_B_SEM))=0;
TP_DLW_HI_AB_N0_B_Mean(isnan(TP_DLW_HI_AB_N0_B_Mean))=0;TP_DLW_HI_AB_N0_B_SEM(isnan(TP_DLW_HI_AB_N0_B_SEM))=0;
TP_DLW_HI_AB_N60_B_Mean(isnan(TP_DLW_HI_AB_N60_B_Mean))=0;TP_DLW_HI_AB_N60_B_SEM(isnan(TP_DLW_HI_AB_N60_B_SEM))=0;
TP_DLW_HI_AB_N70_B_Mean(isnan(TP_DLW_HI_AB_N70_B_Mean))=0;TP_DLW_HI_AB_N70_B_SEM(isnan(TP_DLW_HI_AB_N70_B_SEM))=0;

TP_FGSW_Mean(isnan(TP_FGSW_Mean))=0;TP_FGSW_SEM(isnan(TP_FGSW_SEM))=0;
TP_FGLW_Mean(isnan(TP_FGLW_Mean))=0;TP_FGLW_SEM(isnan(TP_FGLW_SEM))=0;
TP_FSW_NH_Mean(isnan(TP_FSW_NH_Mean))=0;TP_FSW_NH_SEM(isnan(TP_FSW_NH_SEM))=0;
TP_FLW_NH_Mean(isnan(TP_FLW_NH_Mean))=0;TP_FLW_NH_SEM(isnan(TP_FLW_NH_SEM))=0;
TP_FSW_N0_Mean(isnan(TP_FSW_N0_Mean))=0;TP_FSW_N0_SEM(isnan(TP_FSW_N0_SEM))=0;
TP_FSW_N60_Mean(isnan(TP_FSW_N60_Mean))=0;TP_FSW_N60_SEM(isnan(TP_FSW_N60_SEM))=0;
TP_FSW_N70_Mean(isnan(TP_FSW_N70_Mean))=0;TP_FSW_N70_SEM(isnan(TP_FSW_N70_SEM))=0;
TP_FLW_N0_Mean(isnan(TP_FLW_N0_Mean))=0;TP_FLW_N0_SEM(isnan(TP_FLW_N0_SEM))=0;
TP_FLW_N60_Mean(isnan(TP_FLW_N60_Mean))=0;TP_FLW_N60_SEM(isnan(TP_FLW_N60_SEM))=0;
TP_FLW_N70_Mean(isnan(TP_FLW_N70_Mean))=0;TP_FLW_N70_SEM(isnan(TP_FLW_N70_SEM))=0;
TP_FSW_NH_N0_Mean(isnan(TP_FSW_NH_N0_Mean))=0;TP_FSW_NH_N0_SEM(isnan(TP_FSW_NH_N0_SEM))=0;
TP_FSW_NH_N60_Mean(isnan(TP_FSW_NH_N60_Mean))=0;TP_FSW_NH_N60_SEM(isnan(TP_FSW_NH_N60_SEM))=0;
TP_FSW_NH_N70_Mean(isnan(TP_FSW_NH_N70_Mean))=0;TP_FSW_NH_N70_SEM(isnan(TP_FSW_NH_N70_SEM))=0;
TP_FLW_NH_N0_Mean(isnan(TP_FLW_NH_N0_Mean))=0;TP_FLW_NH_N0_SEM(isnan(TP_FLW_NH_N0_SEM))=0;
TP_FLW_NH_N60_Mean(isnan(TP_FLW_NH_N60_Mean))=0;TP_FLW_NH_N60_SEM(isnan(TP_FLW_NH_N60_SEM))=0;
TP_FLW_NH_N70_Mean(isnan(TP_FLW_NH_N70_Mean))=0;TP_FLW_NH_N70_SEM(isnan(TP_FLW_NH_N70_SEM))=0;
TP_FSW_HI_Mean(isnan(TP_FSW_HI_Mean))=0;TP_FSW_HI_SEM(isnan(TP_FSW_HI_SEM))=0;
TP_FLW_HI_Mean(isnan(TP_FLW_HI_Mean))=0;TP_FLW_HI_SEM(isnan(TP_FLW_HI_SEM))=0;
TP_FSW_HI_N0_Mean(isnan(TP_FSW_HI_N0_Mean))=0;TP_FSW_HI_N0_SEM(isnan(TP_FSW_HI_N0_SEM))=0;
TP_FSW_HI_N60_Mean(isnan(TP_FSW_HI_N60_Mean))=0;TP_FSW_HI_N60_SEM(isnan(TP_FSW_HI_N60_SEM))=0;
TP_FSW_HI_N70_Mean(isnan(TP_FSW_HI_N70_Mean))=0;TP_FSW_HI_N70_SEM(isnan(TP_FSW_HI_N70_SEM))=0;
TP_FLW_HI_N0_Mean(isnan(TP_FLW_HI_N0_Mean))=0;TP_FLW_HI_N0_SEM(isnan(TP_FLW_HI_N0_SEM))=0;
TP_FLW_HI_N60_Mean(isnan(TP_FLW_HI_N60_Mean))=0;TP_FLW_HI_N60_SEM(isnan(TP_FLW_HI_N60_SEM))=0;
TP_FLW_HI_N70_Mean(isnan(TP_FLW_HI_N70_Mean))=0;TP_FLW_HI_N70_SEM(isnan(TP_FLW_HI_N70_SEM))=0;
TP_FSW_NH_UN_Mean(isnan(TP_FSW_NH_UN_Mean))=0;TP_FSW_NH_UN_SEM(isnan(TP_FSW_NH_UN_SEM))=0;
TP_FLW_NH_UN_Mean(isnan(TP_FLW_NH_UN_Mean))=0;TP_FLW_NH_UN_SEM(isnan(TP_FLW_NH_UN_SEM))=0;
TP_FSW_NH_AA_Mean(isnan(TP_FSW_NH_AA_Mean))=0;TP_FSW_NH_AA_SEM(isnan(TP_FSW_NH_AA_SEM))=0;
TP_FLW_NH_AA_Mean(isnan(TP_FLW_NH_AA_Mean))=0;TP_FLW_NH_AA_SEM(isnan(TP_FLW_NH_AA_SEM))=0;
TP_FSW_NH_AB_Mean(isnan(TP_FSW_NH_AB_Mean))=0;TP_FSW_NH_AB_SEM(isnan(TP_FSW_NH_AB_SEM))=0;
TP_FLW_NH_AB_Mean(isnan(TP_FLW_NH_AB_Mean))=0;TP_FLW_NH_AB_SEM(isnan(TP_FLW_NH_AB_SEM))=0;
TP_FSW_HI_UN_Mean(isnan(TP_FSW_HI_UN_Mean))=0;TP_FSW_HI_UN_SEM(isnan(TP_FSW_HI_UN_SEM))=0;
TP_FLW_HI_UN_Mean(isnan(TP_FLW_HI_UN_Mean))=0;TP_FLW_HI_UN_SEM(isnan(TP_FLW_HI_UN_SEM))=0;
TP_FSW_HI_AA_Mean(isnan(TP_FSW_HI_AA_Mean))=0;TP_FSW_HI_AA_SEM(isnan(TP_FSW_HI_AA_SEM))=0;
TP_FLW_HI_AA_Mean(isnan(TP_FLW_HI_AA_Mean))=0;TP_FLW_HI_AA_SEM(isnan(TP_FLW_HI_AA_SEM))=0;
TP_FSW_HI_AB_Mean(isnan(TP_FSW_HI_AB_Mean))=0;TP_FSW_HI_AB_SEM(isnan(TP_FSW_HI_AB_SEM))=0;
TP_FLW_HI_AB_Mean(isnan(TP_FLW_HI_AB_Mean))=0;TP_FLW_HI_AB_SEM(isnan(TP_FLW_HI_AB_SEM))=0;
TP_FSW_NH_UN_N0_Mean(isnan(TP_FSW_NH_UN_N0_Mean))=0;TP_FSW_NH_UN_N0_SEM(isnan(TP_FSW_NH_UN_N0_SEM))=0;
TP_FSW_NH_UN_N60_Mean(isnan(TP_FSW_NH_UN_N60_Mean))=0;TP_FSW_NH_UN_N60_SEM(isnan(TP_FSW_NH_UN_N60_SEM))=0;
TP_FSW_NH_UN_N70_Mean(isnan(TP_FSW_NH_UN_N70_Mean))=0;TP_FSW_NH_UN_N70_SEM(isnan(TP_FSW_NH_UN_N70_SEM))=0;
TP_FLW_NH_UN_N0_Mean(isnan(TP_FLW_NH_UN_N0_Mean))=0;TP_FLW_NH_UN_N0_SEM(isnan(TP_FLW_NH_UN_N0_SEM))=0;
TP_FLW_NH_UN_N60_Mean(isnan(TP_FLW_NH_UN_N60_Mean))=0;TP_FLW_NH_UN_N60_SEM(isnan(TP_FLW_NH_UN_N60_SEM))=0;
TP_FLW_NH_UN_N70_Mean(isnan(TP_FLW_NH_UN_N70_Mean))=0;TP_FLW_NH_UN_N70_SEM(isnan(TP_FLW_NH_UN_N70_SEM))=0;
TP_FSW_NH_AA_N0_Mean(isnan(TP_FSW_NH_AA_N0_Mean))=0;TP_FSW_NH_AA_N0_SEM(isnan(TP_FSW_NH_AA_N0_SEM))=0;
TP_FSW_NH_AA_N60_Mean(isnan(TP_FSW_NH_AA_N60_Mean))=0;TP_FSW_NH_AA_N60_SEM(isnan(TP_FSW_NH_AA_N60_SEM))=0;
TP_FSW_NH_AA_N70_Mean(isnan(TP_FSW_NH_AA_N70_Mean))=0;TP_FSW_NH_AA_N70_SEM(isnan(TP_FSW_NH_AA_N70_SEM))=0;
TP_FLW_NH_AA_N0_Mean(isnan(TP_FLW_NH_AA_N0_Mean))=0;TP_FLW_NH_AA_N0_SEM(isnan(TP_FLW_NH_AA_N0_SEM))=0;
TP_FLW_NH_AA_N60_Mean(isnan(TP_FLW_NH_AA_N60_Mean))=0;TP_FLW_NH_AA_N60_SEM(isnan(TP_FLW_NH_AA_N60_SEM))=0;
TP_FLW_NH_AA_N70_Mean(isnan(TP_FLW_NH_AA_N70_Mean))=0;TP_FLW_NH_AA_N70_SEM(isnan(TP_FLW_NH_AA_N70_SEM))=0;
TP_FSW_NH_AB_N0_Mean(isnan(TP_FSW_NH_AB_N0_Mean))=0;TP_FSW_NH_AB_N0_SEM(isnan(TP_FSW_NH_AB_N0_SEM))=0;
TP_FSW_NH_AB_N60_Mean(isnan(TP_FSW_NH_AB_N60_Mean))=0;TP_FSW_NH_AB_N60_SEM(isnan(TP_FSW_NH_AB_N60_SEM))=0;
TP_FSW_NH_AB_N70_Mean(isnan(TP_FSW_NH_AB_N70_Mean))=0;TP_FSW_NH_AB_N70_SEM(isnan(TP_FSW_NH_AB_N70_SEM))=0;
TP_FLW_NH_AB_N0_Mean(isnan(TP_FLW_NH_AB_N0_Mean))=0;TP_FLW_NH_AB_N0_SEM(isnan(TP_FLW_NH_AB_N0_SEM))=0;
TP_FLW_NH_AB_N60_Mean(isnan(TP_FLW_NH_AB_N60_Mean))=0;TP_FLW_NH_AB_N60_SEM(isnan(TP_FLW_NH_AB_N60_SEM))=0;
TP_FLW_NH_AB_N70_Mean(isnan(TP_FLW_NH_AB_N70_Mean))=0;TP_FLW_NH_AB_N70_SEM(isnan(TP_FLW_NH_AB_N70_SEM))=0;
TP_FSW_HI_UN_N0_Mean(isnan(TP_FSW_HI_UN_N0_Mean))=0;TP_FSW_HI_UN_N0_SEM(isnan(TP_FSW_HI_UN_N0_SEM))=0;
TP_FSW_HI_UN_N60_Mean(isnan(TP_FSW_HI_UN_N60_Mean))=0;TP_FSW_HI_UN_N60_SEM(isnan(TP_FSW_HI_UN_N60_SEM))=0;
TP_FSW_HI_UN_N70_Mean(isnan(TP_FSW_HI_UN_N70_Mean))=0;TP_FSW_HI_UN_N70_SEM(isnan(TP_FSW_HI_UN_N70_SEM))=0;
TP_FLW_HI_UN_N0_Mean(isnan(TP_FLW_HI_UN_N0_Mean))=0;TP_FLW_HI_UN_N0_SEM(isnan(TP_FLW_HI_UN_N0_SEM))=0;
TP_FLW_HI_UN_N60_Mean(isnan(TP_FLW_HI_UN_N60_Mean))=0;TP_FLW_HI_UN_N60_SEM(isnan(TP_FLW_HI_UN_N60_SEM))=0;
TP_FLW_HI_UN_N70_Mean(isnan(TP_FLW_HI_UN_N70_Mean))=0;TP_FLW_HI_UN_N70_SEM(isnan(TP_FLW_HI_UN_N70_SEM))=0;
TP_FSW_HI_AA_N0_Mean(isnan(TP_FSW_HI_AA_N0_Mean))=0;TP_FSW_HI_AA_N0_SEM(isnan(TP_FSW_HI_AA_N0_SEM))=0;
TP_FSW_HI_AA_N60_Mean(isnan(TP_FSW_HI_AA_N60_Mean))=0;TP_FSW_HI_AA_N60_SEM(isnan(TP_FSW_HI_AA_N60_SEM))=0;
TP_FSW_HI_AA_N70_Mean(isnan(TP_FSW_HI_AA_N70_Mean))=0;TP_FSW_HI_AA_N70_SEM(isnan(TP_FSW_HI_AA_N70_SEM))=0;
TP_FLW_HI_AA_N0_Mean(isnan(TP_FLW_HI_AA_N0_Mean))=0;TP_FLW_HI_AA_N0_SEM(isnan(TP_FLW_HI_AA_N0_SEM))=0;
TP_FLW_HI_AA_N60_Mean(isnan(TP_FLW_HI_AA_N60_Mean))=0;TP_FLW_HI_AA_N60_SEM(isnan(TP_FLW_HI_AA_N60_SEM))=0;
TP_FLW_HI_AA_N70_Mean(isnan(TP_FLW_HI_AA_N70_Mean))=0;TP_FLW_HI_AA_N70_SEM(isnan(TP_FLW_HI_AA_N70_SEM))=0;
TP_FSW_HI_AB_N0_Mean(isnan(TP_FSW_HI_AB_N0_Mean))=0;TP_FSW_HI_AB_N0_SEM(isnan(TP_FSW_HI_AB_N0_SEM))=0;
TP_FSW_HI_AB_N60_Mean(isnan(TP_FSW_HI_AB_N60_Mean))=0;TP_FSW_HI_AB_N60_SEM(isnan(TP_FSW_HI_AB_N60_SEM))=0;
TP_FSW_HI_AB_N70_Mean(isnan(TP_FSW_HI_AB_N70_Mean))=0;TP_FSW_HI_AB_N70_SEM(isnan(TP_FSW_HI_AB_N70_SEM))=0;
TP_FLW_HI_AB_N0_Mean(isnan(TP_FLW_HI_AB_N0_Mean))=0;TP_FLW_HI_AB_N0_SEM(isnan(TP_FLW_HI_AB_N0_SEM))=0;
TP_FLW_HI_AB_N60_Mean(isnan(TP_FLW_HI_AB_N60_Mean))=0;TP_FLW_HI_AB_N60_SEM(isnan(TP_FLW_HI_AB_N60_SEM))=0;
TP_FLW_HI_AB_N70_Mean(isnan(TP_FLW_HI_AB_N70_Mean))=0;TP_FLW_HI_AB_N70_SEM(isnan(TP_FLW_HI_AB_N70_SEM))=0;

% Fill
fill(ax1,[t, flipud(t')'],[(TP_DSW_NH_UN_N0_Mean+TP_DSW_NH_UN_N0_SEM), flipud((TP_DSW_NH_UN_N0_Mean-TP_DSW_NH_UN_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[t, flipud(t')'],[(TP_DSW_NH_UN_N60_Mean+TP_DSW_NH_UN_N60_SEM), flipud((TP_DSW_NH_UN_N60_Mean-TP_DSW_NH_UN_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax1,[t, flipud(t')'],[(TP_DSW_NH_UN_N70_Mean+TP_DSW_NH_UN_N70_SEM), flipud((TP_DSW_NH_UN_N70_Mean-TP_DSW_NH_UN_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[t, flipud(t')'],[(TP_DLW_NH_UN_N0_Mean+TP_DLW_NH_UN_N0_SEM), flipud((TP_DLW_NH_UN_N0_Mean-TP_DLW_NH_UN_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[t, flipud(t')'],[(TP_DLW_NH_UN_N60_Mean+TP_DLW_NH_UN_N60_SEM), flipud((TP_DLW_NH_UN_N60_Mean-TP_DLW_NH_UN_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax2,[t, flipud(t')'],[(TP_DLW_NH_UN_N70_Mean+TP_DLW_NH_UN_N70_SEM), flipud((TP_DLW_NH_UN_N70_Mean-TP_DLW_NH_UN_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[t, flipud(t')'],[(TP_DSW_NH_UN_N0_B_Mean+TP_DSW_NH_UN_N0_B_SEM), flipud((TP_DSW_NH_UN_N0_B_Mean-TP_DSW_NH_UN_N0_B_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[t, flipud(t')'],[(TP_DSW_NH_UN_N60_B_Mean+TP_DSW_NH_UN_N60_B_SEM), flipud((TP_DSW_NH_UN_N60_B_Mean-TP_DSW_NH_UN_N60_B_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax3,[t, flipud(t')'],[(TP_DSW_NH_UN_N70_B_Mean+TP_DSW_NH_UN_N70_B_SEM), flipud((TP_DSW_NH_UN_N70_B_Mean-TP_DSW_NH_UN_N70_B_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[t, flipud(t')'],[(TP_DLW_NH_UN_N0_B_Mean+TP_DLW_NH_UN_N0_B_SEM), flipud((TP_DLW_NH_UN_N0_B_Mean-TP_DLW_NH_UN_N0_B_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[t, flipud(t')'],[(TP_DLW_NH_UN_N60_B_Mean+TP_DLW_NH_UN_N60_B_SEM), flipud((TP_DLW_NH_UN_N60_B_Mean-TP_DLW_NH_UN_N60_B_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax4,[t, flipud(t')'],[(TP_DLW_NH_UN_N70_B_Mean+TP_DLW_NH_UN_N70_B_SEM), flipud((TP_DLW_NH_UN_N70_B_Mean-TP_DLW_NH_UN_N70_B_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax5,[t, flipud(t')'],[(TP_DSW_HI_UN_N0_Mean+TP_DSW_HI_UN_N0_SEM), flipud((TP_DSW_HI_UN_N0_Mean-TP_DSW_HI_UN_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax5,[t, flipud(t')'],[(TP_DSW_HI_UN_N60_Mean+TP_DSW_HI_UN_N60_SEM), flipud((TP_DSW_HI_UN_N60_Mean-TP_DSW_HI_UN_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax5,[t, flipud(t')'],[(TP_DSW_HI_UN_N70_Mean+TP_DSW_HI_UN_N70_SEM), flipud((TP_DSW_HI_UN_N70_Mean-TP_DSW_HI_UN_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax6,[t, flipud(t')'],[(TP_DLW_HI_UN_N0_Mean+TP_DLW_HI_UN_N0_SEM), flipud((TP_DLW_HI_UN_N0_Mean-TP_DLW_HI_UN_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax6,[t, flipud(t')'],[(TP_DLW_HI_UN_N60_Mean+TP_DLW_HI_UN_N60_SEM), flipud((TP_DLW_HI_UN_N60_Mean-TP_DLW_HI_UN_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax6,[t, flipud(t')'],[(TP_DLW_HI_UN_N70_Mean+TP_DLW_HI_UN_N70_SEM), flipud((TP_DLW_HI_UN_N70_Mean-TP_DLW_HI_UN_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax7,[t, flipud(t')'],[(TP_DSW_HI_UN_N0_B_Mean+TP_DSW_HI_UN_N0_B_SEM), flipud((TP_DSW_HI_UN_N0_B_Mean-TP_DSW_HI_UN_N0_B_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax7,[t, flipud(t')'],[(TP_DSW_HI_UN_N60_B_Mean+TP_DSW_HI_UN_N60_B_SEM), flipud((TP_DSW_HI_UN_N60_B_Mean-TP_DSW_HI_UN_N60_B_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax7,[t, flipud(t')'],[(TP_DSW_HI_UN_N70_B_Mean+TP_DSW_HI_UN_N70_B_SEM), flipud((TP_DSW_HI_UN_N70_B_Mean-TP_DSW_HI_UN_N70_B_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax8,[t, flipud(t')'],[(TP_DLW_HI_UN_N0_B_Mean+TP_DLW_HI_UN_N0_B_SEM), flipud((TP_DLW_HI_UN_N0_B_Mean-TP_DLW_HI_UN_N0_B_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax8,[t, flipud(t')'],[(TP_DLW_HI_UN_N60_B_Mean+TP_DLW_HI_UN_N60_B_SEM), flipud((TP_DLW_HI_UN_N60_B_Mean-TP_DLW_HI_UN_N60_B_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax8,[t, flipud(t')'],[(TP_DLW_HI_UN_N70_B_Mean+TP_DLW_HI_UN_N70_B_SEM), flipud((TP_DLW_HI_UN_N70_B_Mean-TP_DLW_HI_UN_N70_B_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')

fill(ax9,[t, flipud(t')'],[(TP_DSW_NH_AA_N0_Mean+TP_DSW_NH_AA_N0_SEM), flipud((TP_DSW_NH_AA_N0_Mean-TP_DSW_NH_AA_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax9,[t, flipud(t')'],[(TP_DSW_NH_AA_N60_Mean+TP_DSW_NH_AA_N60_SEM), flipud((TP_DSW_NH_AA_N60_Mean-TP_DSW_NH_AA_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax9,[t, flipud(t')'],[(TP_DSW_NH_AA_N70_Mean+TP_DSW_NH_AA_N70_SEM), flipud((TP_DSW_NH_AA_N70_Mean-TP_DSW_NH_AA_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax10,[t, flipud(t')'],[(TP_DLW_NH_AA_N0_Mean+TP_DLW_NH_AA_N0_SEM), flipud((TP_DLW_NH_AA_N0_Mean-TP_DLW_NH_AA_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax10,[t, flipud(t')'],[(TP_DLW_NH_AA_N60_Mean+TP_DLW_NH_AA_N60_SEM), flipud((TP_DLW_NH_AA_N60_Mean-TP_DLW_NH_AA_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax10,[t, flipud(t')'],[(TP_DLW_NH_AA_N70_Mean+TP_DLW_NH_AA_N70_SEM), flipud((TP_DLW_NH_AA_N70_Mean-TP_DLW_NH_AA_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax11,[t, flipud(t')'],[(TP_DSW_NH_AA_N0_B_Mean+TP_DSW_NH_AA_N0_B_SEM), flipud((TP_DSW_NH_AA_N0_B_Mean-TP_DSW_NH_AA_N0_B_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax11,[t, flipud(t')'],[(TP_DSW_NH_AA_N60_B_Mean+TP_DSW_NH_AA_N60_B_SEM), flipud((TP_DSW_NH_AA_N60_B_Mean-TP_DSW_NH_AA_N60_B_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax11,[t, flipud(t')'],[(TP_DSW_NH_AA_N70_B_Mean+TP_DSW_NH_AA_N70_B_SEM), flipud((TP_DSW_NH_AA_N70_B_Mean-TP_DSW_NH_AA_N70_B_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax12,[t, flipud(t')'],[(TP_DLW_NH_AA_N0_B_Mean+TP_DLW_NH_AA_N0_B_SEM), flipud((TP_DLW_NH_AA_N0_B_Mean-TP_DLW_NH_AA_N0_B_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax12,[t, flipud(t')'],[(TP_DLW_NH_AA_N60_B_Mean+TP_DLW_NH_AA_N60_B_SEM), flipud((TP_DLW_NH_AA_N60_B_Mean-TP_DLW_NH_AA_N60_B_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax12,[t, flipud(t')'],[(TP_DLW_NH_AA_N70_B_Mean+TP_DLW_NH_AA_N70_B_SEM), flipud((TP_DLW_NH_AA_N70_B_Mean-TP_DLW_NH_AA_N70_B_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax13,[t, flipud(t')'],[(TP_DSW_HI_AA_N0_Mean+TP_DSW_HI_AA_N0_SEM), flipud((TP_DSW_HI_AA_N0_Mean-TP_DSW_HI_AA_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax13,[t, flipud(t')'],[(TP_DSW_HI_AA_N60_Mean+TP_DSW_HI_AA_N60_SEM), flipud((TP_DSW_HI_AA_N60_Mean-TP_DSW_HI_AA_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax13,[t, flipud(t')'],[(TP_DSW_HI_AA_N70_Mean+TP_DSW_HI_AA_N70_SEM), flipud((TP_DSW_HI_AA_N70_Mean-TP_DSW_HI_AA_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax14,[t, flipud(t')'],[(TP_DLW_HI_AA_N0_Mean+TP_DLW_HI_AA_N0_SEM), flipud((TP_DLW_HI_AA_N0_Mean-TP_DLW_HI_AA_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax14,[t, flipud(t')'],[(TP_DLW_HI_AA_N60_Mean+TP_DLW_HI_AA_N60_SEM), flipud((TP_DLW_HI_AA_N60_Mean-TP_DLW_HI_AA_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax14,[t, flipud(t')'],[(TP_DLW_HI_AA_N70_Mean+TP_DLW_HI_AA_N70_SEM), flipud((TP_DLW_HI_AA_N70_Mean-TP_DLW_HI_AA_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax15,[t, flipud(t')'],[(TP_DSW_HI_AA_N0_B_Mean+TP_DSW_HI_AA_N0_B_SEM), flipud((TP_DSW_HI_AA_N0_B_Mean-TP_DSW_HI_AA_N0_B_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax15,[t, flipud(t')'],[(TP_DSW_HI_AA_N60_B_Mean+TP_DSW_HI_AA_N60_B_SEM), flipud((TP_DSW_HI_AA_N60_B_Mean-TP_DSW_HI_AA_N60_B_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax15,[t, flipud(t')'],[(TP_DSW_HI_AA_N70_B_Mean+TP_DSW_HI_AA_N70_B_SEM), flipud((TP_DSW_HI_AA_N70_B_Mean-TP_DSW_HI_AA_N70_B_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax16,[t, flipud(t')'],[(TP_DLW_HI_AA_N0_B_Mean+TP_DLW_HI_AA_N0_B_SEM), flipud((TP_DLW_HI_AA_N0_B_Mean-TP_DLW_HI_AA_N0_B_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax16,[t, flipud(t')'],[(TP_DLW_HI_AA_N60_B_Mean+TP_DLW_HI_AA_N60_B_SEM), flipud((TP_DLW_HI_AA_N60_B_Mean-TP_DLW_HI_AA_N60_B_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax16,[t, flipud(t')'],[(TP_DLW_HI_AA_N70_B_Mean+TP_DLW_HI_AA_N70_B_SEM), flipud((TP_DLW_HI_AA_N70_B_Mean-TP_DLW_HI_AA_N70_B_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')

fill(ax17,[t, flipud(t')'],[(TP_DSW_NH_AB_N0_Mean+TP_DSW_NH_AB_N0_SEM), flipud((TP_DSW_NH_AB_N0_Mean-TP_DSW_NH_AB_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax17,[t, flipud(t')'],[(TP_DSW_NH_AB_N60_Mean+TP_DSW_NH_AB_N60_SEM), flipud((TP_DSW_NH_AB_N60_Mean-TP_DSW_NH_AB_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax17,[t, flipud(t')'],[(TP_DSW_NH_AB_N70_Mean+TP_DSW_NH_AB_N70_SEM), flipud((TP_DSW_NH_AB_N70_Mean-TP_DSW_NH_AB_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax18,[t, flipud(t')'],[(TP_DLW_NH_AB_N0_Mean+TP_DLW_NH_AB_N0_SEM), flipud((TP_DLW_NH_AB_N0_Mean-TP_DLW_NH_AB_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax18,[t, flipud(t')'],[(TP_DLW_NH_AB_N60_Mean+TP_DLW_NH_AB_N60_SEM), flipud((TP_DLW_NH_AB_N60_Mean-TP_DLW_NH_AB_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax18,[t, flipud(t')'],[(TP_DLW_NH_AB_N70_Mean+TP_DLW_NH_AB_N70_SEM), flipud((TP_DLW_NH_AB_N70_Mean-TP_DLW_NH_AB_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax19,[t, flipud(t')'],[(TP_DSW_NH_AB_N0_B_Mean+TP_DSW_NH_AB_N0_B_SEM), flipud((TP_DSW_NH_AB_N0_B_Mean-TP_DSW_NH_AB_N0_B_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax19,[t, flipud(t')'],[(TP_DSW_NH_AB_N60_B_Mean+TP_DSW_NH_AB_N60_B_SEM), flipud((TP_DSW_NH_AB_N60_B_Mean-TP_DSW_NH_AB_N60_B_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax19,[t, flipud(t')'],[(TP_DSW_NH_AB_N70_B_Mean+TP_DSW_NH_AB_N70_B_SEM), flipud((TP_DSW_NH_AB_N70_B_Mean+TP_DSW_NH_AB_N70_B_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax20,[t, flipud(t')'],[(TP_DLW_NH_AB_N0_B_Mean+TP_DLW_NH_AB_N0_B_SEM), flipud((TP_DLW_NH_AB_N0_B_Mean-TP_DLW_NH_AB_N0_B_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax20,[t, flipud(t')'],[(TP_DLW_NH_AB_N60_B_Mean+TP_DLW_NH_AB_N60_B_SEM), flipud((TP_DLW_NH_AB_N60_B_Mean-TP_DLW_NH_AB_N60_B_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax20,[t, flipud(t')'],[(TP_DLW_NH_AB_N70_B_Mean+TP_DLW_NH_AB_N70_B_SEM), flipud((TP_DLW_NH_AB_N70_B_Mean-TP_DLW_NH_AB_N70_B_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax21,[t, flipud(t')'],[(TP_DSW_HI_AB_N0_Mean+TP_DSW_HI_AB_N0_SEM), flipud((TP_DSW_HI_AB_N0_Mean-TP_DSW_HI_AB_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax21,[t, flipud(t')'],[(TP_DSW_HI_AB_N60_Mean+TP_DSW_HI_AB_N60_SEM), flipud((TP_DSW_HI_AB_N60_Mean-TP_DSW_HI_AB_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax21,[t, flipud(t')'],[(TP_DSW_HI_AB_N70_Mean+TP_DSW_HI_AB_N70_SEM), flipud((TP_DSW_HI_AB_N70_Mean-TP_DSW_HI_AB_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax22,[t, flipud(t')'],[(TP_DLW_HI_AB_N0_Mean+TP_DLW_HI_AB_N0_SEM), flipud((TP_DLW_HI_AB_N0_Mean-TP_DLW_HI_AB_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax22,[t, flipud(t')'],[(TP_DLW_HI_AB_N60_Mean+TP_DLW_HI_AB_N60_SEM), flipud((TP_DLW_HI_AB_N60_Mean-TP_DLW_HI_AB_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax22,[t, flipud(t')'],[(TP_DLW_HI_AB_N70_Mean+TP_DLW_HI_AB_N70_SEM), flipud((TP_DLW_HI_AB_N70_Mean-TP_DLW_HI_AB_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax23,[t, flipud(t')'],[(TP_DSW_HI_AB_N0_B_Mean+TP_DSW_HI_AB_N0_B_SEM), flipud((TP_DSW_HI_AB_N0_B_Mean-TP_DSW_HI_AB_N0_B_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax23,[t, flipud(t')'],[(TP_DSW_HI_AB_N60_B_Mean+TP_DSW_HI_AB_N60_B_SEM), flipud((TP_DSW_HI_AB_N60_B_Mean-TP_DSW_HI_AB_N60_B_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax23,[t, flipud(t')'],[(TP_DSW_HI_AB_N70_B_Mean+TP_DSW_HI_AB_N70_B_SEM), flipud((TP_DSW_HI_AB_N70_B_Mean-TP_DSW_HI_AB_N70_B_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax24,[t, flipud(t')'],[(TP_DLW_HI_AB_N0_B_Mean+TP_DLW_HI_AB_N0_B_SEM), flipud((TP_DLW_HI_AB_N0_B_Mean-TP_DLW_HI_AB_N0_B_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax24,[t, flipud(t')'],[(TP_DLW_HI_AB_N60_B_Mean+TP_DLW_HI_AB_N60_B_SEM), flipud((TP_DLW_HI_AB_N60_B_Mean-TP_DLW_HI_AB_N60_B_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax24,[t, flipud(t')'],[(TP_DLW_HI_AB_N70_B_Mean+TP_DLW_HI_AB_N70_B_SEM), flipud((TP_DLW_HI_AB_N70_B_Mean-TP_DLW_HI_AB_N70_B_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')

fill(ax25,[t, flipud(t')'],[(TP_FSW_NH_UN_N0_Mean+TP_FSW_NH_UN_N0_SEM), flipud((TP_FSW_NH_UN_N0_Mean-TP_FSW_NH_UN_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax25,[t, flipud(t')'],[(TP_FSW_NH_UN_N60_Mean+TP_FSW_NH_UN_N60_SEM), flipud((TP_FSW_NH_UN_N60_Mean-TP_FSW_NH_UN_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax25,[t, flipud(t')'],[(TP_FSW_NH_UN_N70_Mean+TP_FSW_NH_UN_N70_SEM), flipud((TP_FSW_NH_UN_N70_Mean-TP_FSW_NH_UN_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax26,[t, flipud(t')'],[(TP_FLW_NH_UN_N0_Mean+TP_FLW_NH_UN_N0_SEM), flipud((TP_FLW_NH_UN_N0_Mean-TP_FLW_NH_UN_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax26,[t, flipud(t')'],[(TP_FLW_NH_UN_N60_Mean+TP_FLW_NH_UN_N60_SEM), flipud((TP_FLW_NH_UN_N60_Mean-TP_FLW_NH_UN_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax26,[t, flipud(t')'],[(TP_FLW_NH_UN_N70_Mean+TP_FLW_NH_UN_N70_SEM), flipud((TP_FLW_NH_UN_N70_Mean-TP_FLW_NH_UN_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax27,[t, flipud(t')'],[(TP_FSW_HI_UN_N0_Mean+TP_FSW_HI_UN_N0_SEM), flipud((TP_FSW_HI_UN_N0_Mean-TP_FSW_HI_UN_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax27,[t, flipud(t')'],[(TP_FSW_HI_UN_N60_Mean+TP_FSW_HI_UN_N60_SEM), flipud((TP_FSW_HI_UN_N60_Mean-TP_FSW_HI_UN_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax27,[t, flipud(t')'],[(TP_FSW_HI_UN_N70_Mean+TP_FSW_HI_UN_N70_SEM), flipud((TP_FSW_HI_UN_N70_Mean-TP_FSW_HI_UN_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax28,[t, flipud(t')'],[(TP_FLW_HI_UN_N0_Mean+TP_FLW_HI_UN_N0_SEM), flipud((TP_FLW_HI_UN_N0_Mean-TP_FLW_HI_UN_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax28,[t, flipud(t')'],[(TP_FLW_HI_UN_N60_Mean+TP_FLW_HI_UN_N60_SEM), flipud((TP_FLW_HI_UN_N60_Mean-TP_FLW_HI_UN_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax28,[t, flipud(t')'],[(TP_FLW_HI_UN_N70_Mean+TP_FLW_HI_UN_N70_SEM), flipud((TP_FLW_HI_UN_N70_Mean-TP_FLW_HI_UN_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')

fill(ax29,[t, flipud(t')'],[(TP_FSW_NH_AA_N0_Mean+TP_FSW_NH_AA_N0_SEM), flipud((TP_FSW_NH_AA_N0_Mean-TP_FSW_NH_AA_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax29,[t, flipud(t')'],[(TP_FSW_NH_AA_N60_Mean+TP_FSW_NH_AA_N60_SEM), flipud((TP_FSW_NH_AA_N60_Mean-TP_FSW_NH_AA_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax29,[t, flipud(t')'],[(TP_FSW_NH_AA_N70_Mean+TP_FSW_NH_AA_N70_SEM), flipud((TP_FSW_NH_AA_N70_Mean-TP_FSW_NH_AA_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax30,[t, flipud(t')'],[(TP_FLW_NH_AA_N0_Mean+TP_FLW_NH_AA_N0_SEM), flipud((TP_FLW_NH_AA_N0_Mean-TP_FLW_NH_AA_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax30,[t, flipud(t')'],[(TP_FLW_NH_AA_N60_Mean+TP_FLW_NH_AA_N60_SEM), flipud((TP_FLW_NH_AA_N60_Mean-TP_FLW_NH_AA_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax30,[t, flipud(t')'],[(TP_FLW_NH_AA_N70_Mean+TP_FLW_NH_AA_N70_SEM), flipud((TP_FLW_NH_AA_N70_Mean-TP_FLW_NH_AA_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax31,[t, flipud(t')'],[(TP_FSW_HI_AA_N0_Mean+TP_FSW_HI_AA_N0_SEM), flipud((TP_FSW_HI_AA_N0_Mean-TP_FSW_HI_AA_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax31,[t, flipud(t')'],[(TP_FSW_HI_AA_N60_Mean+TP_FSW_HI_AA_N60_SEM), flipud((TP_FSW_HI_AA_N60_Mean-TP_FSW_HI_AA_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax31,[t, flipud(t')'],[(TP_FSW_HI_AA_N70_Mean+TP_FSW_HI_AA_N70_SEM), flipud((TP_FSW_HI_AA_N70_Mean-TP_FSW_HI_AA_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax32,[t, flipud(t')'],[(TP_FLW_HI_AA_N0_Mean+TP_FLW_HI_AA_N0_SEM), flipud((TP_FLW_HI_AA_N0_Mean-TP_FLW_HI_AA_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax32,[t, flipud(t')'],[(TP_FLW_HI_AA_N60_Mean+TP_FLW_HI_AA_N60_SEM), flipud((TP_FLW_HI_AA_N60_Mean-TP_FLW_HI_AA_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax32,[t, flipud(t')'],[(TP_FLW_HI_AA_N70_Mean+TP_FLW_HI_AA_N70_SEM), flipud((TP_FLW_HI_AA_N70_Mean-TP_FLW_HI_AA_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')

fill(ax33,[t, flipud(t')'],[(TP_FSW_NH_AB_N0_Mean+TP_FSW_NH_AB_N0_SEM), flipud((TP_FSW_NH_AB_N0_Mean-TP_FSW_NH_AB_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax33,[t, flipud(t')'],[(TP_FSW_NH_AB_N60_Mean+TP_FSW_NH_AB_N60_SEM), flipud((TP_FSW_NH_AB_N60_Mean-TP_FSW_NH_AB_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax33,[t, flipud(t')'],[(TP_FSW_NH_AB_N70_Mean+TP_FSW_NH_AB_N70_SEM), flipud((TP_FSW_NH_AB_N70_Mean-TP_FSW_NH_AB_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax34,[t, flipud(t')'],[(TP_FLW_NH_AB_N0_Mean+TP_FLW_NH_AB_N0_SEM), flipud((TP_FLW_NH_AB_N0_Mean-TP_FLW_NH_AB_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax34,[t, flipud(t')'],[(TP_FLW_NH_AB_N60_Mean+TP_FLW_NH_AB_N60_SEM), flipud((TP_FLW_NH_AB_N60_Mean-TP_FLW_NH_AB_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax34,[t, flipud(t')'],[(TP_FLW_NH_AB_N70_Mean+TP_FLW_NH_AB_N70_SEM), flipud((TP_FLW_NH_AB_N70_Mean-TP_FLW_NH_AB_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax35,[t, flipud(t')'],[(TP_FSW_HI_AB_N0_Mean+TP_FSW_HI_AB_N0_SEM), flipud((TP_FSW_HI_AB_N0_Mean-TP_FSW_HI_AB_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax35,[t, flipud(t')'],[(TP_FSW_HI_AB_N60_Mean+TP_FSW_HI_AB_N60_SEM), flipud((TP_FSW_HI_AB_N60_Mean-TP_FSW_HI_AB_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax35,[t, flipud(t')'],[(TP_FSW_HI_AB_N70_Mean+TP_FSW_HI_AB_N70_SEM), flipud((TP_FSW_HI_AB_N70_Mean-TP_FSW_HI_AB_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax36,[t, flipud(t')'],[(TP_FLW_HI_AB_N0_Mean+TP_FLW_HI_AB_N0_SEM), flipud((TP_FLW_HI_AB_N0_Mean-TP_FLW_HI_AB_N0_SEM)')'],QuietColor,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax36,[t, flipud(t')'],[(TP_FLW_HI_AB_N60_Mean+TP_FLW_HI_AB_N60_SEM), flipud((TP_FLW_HI_AB_N60_Mean-TP_FLW_HI_AB_N60_SEM)')'],N60Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')
fill(ax36,[t, flipud(t')'],[(TP_FLW_HI_AB_N70_Mean+TP_FLW_HI_AB_N70_SEM), flipud((TP_FLW_HI_AB_N70_Mean-TP_FLW_HI_AB_N70_SEM)')'],N70Color,'FaceAlpha',.1,'Edgecolor','none','handlevisibility' ,'off')

xlabel([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10 ax11 ax12 ax13 ax14 ax15 ax16 ax17 ax18 ax19 ax20 ax21 ax22 ax23 ax24 ax25 ax26 ax27 ax28 ax29 ax30 ax31 ax32 ax33 ax34 ax35 ax36],'Time [s]')
ylabel([ax1 ax5 ax9 ax13 ax17 ax21],'Pupil diameter [mm]')
ylabel([ax3 ax7 ax11 ax15 ax19 ax23],'Pupil baseline difference [mm]')
ylabel([ax25 ax27 ax29 ax31 ax33 ax35],'Fixation duration [s]')
xlim([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 ax10 ax11 ax12 ax13 ax14 ax15 ax16 ax17 ax18 ax19 ax20 ax21 ax22 ax23 ax24 ax25 ax26 ax27 ax28 ax29 ax30 ax31 ax32 ax33 ax34 ax35 ax36],[-TimeStartW 3])
% xlim(ax1,[-TimeStartW, min(t_G(cumsum(squeeze(sum(~isnan(DSW_NH),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(DSW_NH),[1 2])))))])
% xlim(ax2,[-TimeStartW, min(t_G(cumsum(squeeze(sum(~isnan(DLW_NH),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(DLW_NH),[1 2])))))])
% xlim(ax3,[-TimeStartW, min(t_G(cumsum(squeeze(sum(~isnan(DSW_HI),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(DSW_HI),[1 2])))))])
% xlim(ax4,[-TimeStartW, min(t_G(cumsum(squeeze(sum(~isnan(DLW_HI),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(DLW_HI),[1 2])))))])
% xlim(ax5,[-TimeStartW, min(t_G(cumsum(squeeze(sum(~isnan(DSW_NH_B),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(DSW_NH_B),[1 2])))))])
% xlim(ax6,[-TimeStartW, min(t_G(cumsum(squeeze(sum(~isnan(DLW_NH_B),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(DLW_NH_B),[1 2])))))])
% xlim(ax7,[-TimeStartW, min(t_G(cumsum(squeeze(sum(~isnan(DSW_HI_B),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(DSW_HI_B),[1 2])))))])
% xlim(ax8,[-TimeStartW, min(t_G(cumsum(squeeze(sum(~isnan(DLW_HI_B),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(DLW_HI_B),[1 2])))))])
% xlim(ax11,[-TimeStartW, min(t_G(cumsum(squeeze(sum(~isnan(FSW_NH),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(FSW_NH),[1 2])))))])
% xlim(ax12,[-TimeStartW, min(t_G(cumsum(squeeze(sum(~isnan(FLW_NH),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(FLW_NH),[1 2])))))])
% xlim(ax13,[-TimeStartW, min(t_G(cumsum(squeeze(sum(~isnan(FSW_HI),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(FSW_HI),[1 2])))))])
% xlim(ax14,[-TimeStartW, min(t_G(cumsum(squeeze(sum(~isnan(FLW_HI),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(FLW_HI),[1 2])))))])
% xlim(ax15,[-TimeStartW, min(t_G_B(cumsum(squeeze(sum(~isnan(FSW_NH_B),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(FSW_NH_B),[1 2])))))])
% xlim(ax16,[-TimeStartW, min(t_G_B(cumsum(squeeze(sum(~isnan(FLW_NH_B),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(FLW_NH_B),[1 2])))))])
% xlim(ax17,[-TimeStartW, min(t_G_B(cumsum(squeeze(sum(~isnan(FSW_HI_B),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(FSW_HI_B),[1 2])))))])
% xlim(ax18,[-TimeStartW, min(t_G_B(cumsum(squeeze(sum(~isnan(FLW_HI_B),[1 2]))) >= VisRatio * sum(squeeze(sum(~isnan(FLW_HI_B),[1 2])))))])
ylim([ax1 ax2],[min([ylim(ax1),ylim(ax2)]), max([ylim(ax1),ylim(ax2)])])
ylim([ax3 ax4],[min([ylim(ax3),ylim(ax4)]), max([ylim(ax3),ylim(ax4)])])
ylim([ax5 ax6],[min([ylim(ax5),ylim(ax6)]), max([ylim(ax5),ylim(ax6)])])
ylim([ax7 ax8],[min([ylim(ax7),ylim(ax8)]), max([ylim(ax7),ylim(ax8)])])
ylim([ax9 ax10],[min([ylim(ax9),ylim(ax10)]), max([ylim(ax9),ylim(ax10)])])
ylim([ax11 ax12],[min([ylim(ax11),ylim(ax12)]), max([ylim(ax11),ylim(ax12)])])
ylim([ax13 ax14],[min([ylim(ax13),ylim(ax14)]), max([ylim(ax13),ylim(ax14)])])
ylim([ax15 ax16],[min([ylim(ax15),ylim(ax16)]), max([ylim(ax15),ylim(ax16)])])
ylim([ax17 ax18],[min([ylim(ax17),ylim(ax18)]), max([ylim(ax17),ylim(ax18)])])
ylim([ax19 ax20],[min([ylim(ax19),ylim(ax20)]), max([ylim(ax19),ylim(ax20)])])
ylim([ax21 ax22],[min([ylim(ax21),ylim(ax22)]), max([ylim(ax21),ylim(ax22)])])
ylim([ax23 ax24],[min([ylim(ax23),ylim(ax24)]), max([ylim(ax23),ylim(ax24)])])
ylim([ax25 ax26],[min([ylim(ax25),ylim(ax26)]), max([ylim(ax25),ylim(ax26)])])
ylim([ax27 ax28],[min([ylim(ax27),ylim(ax28)]), max([ylim(ax27),ylim(ax28)])])
ylim([ax29 ax30],[min([ylim(ax29),ylim(ax30)]), max([ylim(ax29),ylim(ax30)])])
ylim([ax31 ax32],[min([ylim(ax31),ylim(ax32)]), max([ylim(ax31),ylim(ax32)])])
ylim([ax33 ax34],[min([ylim(ax33),ylim(ax34)]), max([ylim(ax33),ylim(ax34)])])
ylim([ax35 ax36],[min([ylim(ax35),ylim(ax36)]), max([ylim(ax35),ylim(ax36)])])

set([ax1 ax3 ax5 ax7 ax9 ax11 ax13 ax15 ax17 ax19 ax21 ax23 ax25 ax27 ax29 ax31 ax33 ax35],'Color',[SColor,0.04])
set([ax2 ax4 ax6 ax8 ax10 ax12 ax14 ax16 ax18 ax20 ax22 ax24 ax26 ax28 ax30 ax32 ax34 ax36],'Color',[LColor,0.04])

% ylmax5_8=cell2mat(ylim([ax5 ax6 ax7 ax8]));
% [max5_8,idx5_8] = max(diff(ylmax5_8,1,2));
% ylim([ax5 ax6 ax7 ax8],ylmax5_8(idx5_8,:)) % Setting ylim to the max out of all
% lgd4=legend(ax4,'N0','N60','N70','Location','southeastoutside');
% lgd4.Title.String = 'Types of windows:';
% lgd8=legend(ax8,'N0','N60','N70','Location','southeastoutside');
% lgd8.Title.String = 'Types of windows:';
% lgd14=legend(ax14,'N0','N60','N70','Location','southeastoutside');
% lgd14.Title.String = 'Types of windows:';
% lgd24=legend(ax24,'N0','N60','N70','Location','southeastoutside');
% lgd24.Title.String = 'Types of windows:';
% lgd28=legend(ax28,'N0','N60','N70','Location','southeastoutside');
% lgd28.Title.String = 'Types of windows:';
% lgd34=legend(ax34,'N0','N60','N70','Location','southeastoutside');
% lgd34.Title.String = 'Types of windows:';
% title([ax1 ax5 ax11 ax21 ax25 ax31],'NH - Speaking')
% title([ax2 ax6 ax12 ax22 ax26 ax32],'NH - Listening')
% title([ax3 ax7 ax13 ax23 ax27 ax33],'HI - Speaking')
% title([ax4 ax8 ax14 ax24 ax28 ax34],'HI - Listening')
% sgtitle(f1,'Global Pupil Diameter - NH','FontWeight','bold')
% sgtitle(f2,'Global Pupil Diameter - HI','FontWeight','bold')
% sgtitle(f3,'Global Pupil Diameter adaptive baselined - NH','FontWeight','bold')
% sgtitle(f4,'Global Pupil Diameter adaptive baselined - HI','FontWeight','bold')
% sgtitle(f5,'Global Fixation duration - NH','FontWeight','bold')
% sgtitle(f6,'Global Fixation duration - HI','FontWeight','bold')
% sgtitle(f9,'Pupil Diameter averaged over TPs - NH','FontWeight','bold')
% sgtitle(f10,'Pupil Diameter averaged over TPs - HI','FontWeight','bold')
% sgtitle(f11,'Pupil Diameter averaged over TPs adaptive baselined - NH','FontWeight','bold')
% sgtitle(f12,'Pupil Diameter averaged over TPs adaptive baselined - HI','FontWeight','bold')
% sgtitle(f13,'Fixation duration averaged over TPs - NH','FontWeight','bold')
% sgtitle(f14,'Fixation duration averaged over TPs - HI','FontWeight','bold')

% flgd=figure;axlgd=gca();hold(axlgd,'on')
% plot(axlgd,nan,color=QuietColor,linewidth=2) % plot nans to show color in legend
% plot(axlgd,nan,color=N60Color,linewidth=2) % plot nans to show color in legend
% plot(axlgd,nan,color=N70Color,linewidth=2) % plot nans to show color in legend
% lgd=legend(axlgd,'N0','N60','N70','Location','southeastoutside');
% lgd.Title.String = 'Noise conditions:';
% 
% 
% flgd2=figure;axlgd2=gca();hold(axlgd2,'on')
% plot(axlgd2,nan,color=[SColor 0.3],linewidth=3) % plot nans to show color in legend
% plot(axlgd2,nan,color=[LColor 0.3],linewidth=3) % plot nans to show color in legend
% lgd2=legend(axlgd2,'Speaking','Listening','Location','southeastoutside');
% lgd2.Title.String = 'Types of windows:';
% 
% flgd3=figure;axlgd3=gca();hold(axlgd3,'on')
% plot(axlgd3,nan,color=UNColor,linewidth=2) % plot nans to show color in legend
% plot(axlgd3,nan,color=AAColor,linewidth=2) % plot nans to show color in legend
% plot(axlgd3,nan,color=ABColor,linewidth=2) % plot nans to show color in legend
% lgd3=legend(axlgd3,'UN','AA','AB','Location','southeastoutside');
% lgd3.Title.String = 'Hearing aid settings:';

% saveLegendToImage(flgd, lgd, 'lgd', 'eps');
% saveLegendToImage(flgd2, lgd2, 'lgd2', 'eps');
% saveLegendToImage(flgd3, lgd3, 'lgd3', 'eps');

% Export figs as pdfs
% tiles = {t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18};
% for i=1:size(tiles,2)
% exportgraphics(tiles{i}, ['output_' num2str(i) '.pdf']);
% end

%% NOTES
% sum(sum(sum(~isnan(DGSW),3)>0,2)); % Used to know how many windows (non-nan rows) are inside each variable (Layers,Rows,Columns)