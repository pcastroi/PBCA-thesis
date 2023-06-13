clear all; clc; close all;

addpath(genpath('C:\git'));

Param.Fs = 50; % Sampling frequency of pupil data
Param.RemoveBeforeAndAfter = [0 0]*1e-3; % Samples within the time range before and after NaNs will set NaNs as well.
Param.MinLengthNaNRepair = 0; % Drop values (i.e., change to NaN) before and after NaNs only for contiguous NaNs of at least __ samples.
MaxVisualAngle = 1.5; % [degrees], critical visual angle for fixation definition
LPWinSize = 0.5; % [s]: Window size of hamming-window for low-pass filtering
FilterWidth = round((LPWinSize*Param.Fs)/2); % [samples]: Width of hamming filter used for fixation duration
LPWindow = hamming(round(LPWinSize*Param.Fs));
LPWindow = LPWindow/sum(LPWindow); % Hamming-window

% Colors
ColRaw = "#39040e";
ColPre = "#c9290e";
ColInt = "#0e0ba9";
ColConv = "#0db93a";

alldata=load('C:\git\PBCA-thesis\data\AMEND_II\Pair03\HI\AAHI_N0.mat');
alldata_mat = cell2mat(alldata.data);
timestamps = vertcat(alldata_mat.timeStamp);
%% PUPIL
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

LDiamPre = LDiamRaw;
RDiamPre = RDiamRaw;

LDiamPre(~LvalOut)=NaN;
RDiamPre(~RvalOut)=NaN;

% Processing - Interpolating NaNs
[LDiamInt,LMetadata] = preprocpupil(LDiamPre,Param);
[RDiamInt,RMetadata] = preprocpupil(RDiamPre,Param);

% Low-Pass Filtering
LDiamConv = conv(LDiamInt,LPWindow,'same'); 
RDiamConv = conv(RDiamInt,LPWindow,'same');

% Remove start/end artifacts (peak/dip) originated from Low-Pass Filtering
LDiamConv(1:round(length(LPWindow)/2-1)) = mean(LDiamInt(1:round(length(LPWindow)/2-1)));
LDiamConv(end-round(length(LPWindow)/2-1):end) = mean(LDiamInt(end-round(length(LPWindow)/2-1):end));
RDiamConv(1:round(length(LPWindow)/2-1)) = mean(RDiamInt(1:round(length(LPWindow)/2-1)));
RDiamConv(end-round(length(LPWindow)/2-1):end) = mean(RDiamInt(end-round(length(LPWindow)/2-1):end));

% Decide 'better' eye results
[Min, idx_decision] = min([sum(LMetadata.Isnan) sum(RMetadata.Isnan)]);
if idx_decision == 1
    DiameterPre = LDiamPre';
    DiameterConv = LDiamConv;
    DiameterInt = LDiamInt;
    DiameterRaw = LDiamRaw';
    eyeChosen = 'Left';
    DiamNaN = sum(LMetadata.Isnan);
elseif idx_decision == 2
    DiameterPre = RDiamPre';
    DiameterConv = RDiamConv;
    DiameterInt = RDiamInt;
    DiameterRaw = RDiamRaw';
    eyeChosen = 'Right';
    DiamNaN = sum(RMetadata.Isnan);
end
t_diam = linspace(0,size(DiameterConv,1)/Param.Fs,size(DiameterConv,1));
%% GAZE
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
FixationRaw = fix_duration([GazeX_fix,GazeY_fix,GazeZ_fix],MaxVisualAngle,Param.Fs);        

% 5) Filtering - Fixation duration with hamming window (F = 3)
Fixation = ndnanfilter(FixationRaw,'hamming',FilterWidth);
%% Plots
f1 = figure;ax1 = gca;
hold(ax1,'on')
% grid(ax1,'on')
plot(ax1, t_diam, DiameterRaw, 'Color', ColRaw, 'LineWidth', 2)
plot(ax1, t_diam, DiameterPre+1, 'Color', ColPre, 'LineWidth', 2)
plot(ax1, t_diam, DiameterInt+2, 'Color', ColInt, 'LineWidth', 2)
plot(ax1, t_diam, DiameterConv+3, 'Color', ColConv, 'LineWidth', 2)
text(ax1,t_diam(5819),DiameterRaw(5819),{'1'},'Color',ColRaw,'EdgeColor',ColRaw,'BackgroundColor',"w",'LineWidth',2)
text(ax1,t_diam(5819),DiameterPre(5819)+1,{'2'},'Color',ColPre,'EdgeColor',ColPre,'BackgroundColor',"w",'LineWidth',2)
text(ax1,t_diam(5819),DiameterInt(5819)+2,{'3'},'Color',ColInt,'EdgeColor',ColInt,'BackgroundColor',"w",'LineWidth',2)
text(ax1,t_diam(5819),DiameterConv(5819)+3,{'4'},'Color',ColConv,'EdgeColor',ColConv,'BackgroundColor',"w",'LineWidth',2)

xlim(ax1,[114,122])
xticklabels(0:8)
set(gca,'ytick',[])
xlabel(ax1,'Time [s]')
ylabel(ax1,{'Pupil size','(separated vertically for clarity)'})
% lgd=legend(ax1,{'Raw data','Pre-processed data','Interpolated data','Low-pass filtered data'},'Location','northeast');

f2 = figure; ax2 = gca;
hold(ax2,'on')
plot(ax2, t_diam, FixationRaw,'Color', ColRaw, 'LineWidth', 2)
plot(ax2, t_diam, Fixation+1,'Color', ColConv, 'LineWidth', 2)
text(ax2,t_diam(8557),FixationRaw(8557),{'1'},'Color',ColRaw,'EdgeColor',ColRaw,'BackgroundColor',"w",'LineWidth',2)
xlim(ax2,[167,175])
xticklabels(0:8)
ylims = get(gca, 'YLim');
height = ylims(2) - ylims(1);
for i = 1:proc_info.number_of_blinks
    rectangle('position', [proc_info.blink_starts_s(i) ylims(1) proc_info.blink_durations(i) height ],...
        'facecolor', [[33,33,33]/255 0.4],...
        'edgecolor', 'none');
end
for i = 1:proc_info.number_of_saccades
        rectangle('position', [proc_info.saccade_starts_s(i) ylims(1) proc_info.saccade_durations(i)  height ],...
        'facecolor', [[4, 117, 111]/255 0.4],...
        'edgecolor', 'none');
end
plot(ax2, t_diam, Fixation+1,'Color', ColConv, 'LineWidth', 2)
text(ax2,t_diam(8557),Fixation(8557)+1,{'2'},'Color',ColConv,'EdgeColor',ColConv,'BackgroundColor',"w",'LineWidth',2)

set(gca,'ytick',[])
xlabel(ax2,'Time [s]')
ylabel(ax2,{'Fixation duration','(separated vertically for clarity)'})