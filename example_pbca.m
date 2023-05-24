clear all; clc; close all;

addpath(genpath('C:\git'));

Param.Fs=50;
alldata=load('C:\git\PBCA-thesis\data\AMEND_II\Pair03\HI\AAHI_N0.mat');
alldata_mat = cell2mat(alldata.data);
timestamps = vertcat(alldata_mat.timeStamp);

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
GazeXRaw = gaze3d(:,1);
GazeYRaw = gaze3d(:,2);
GazeZRaw = gaze3d(:,3);
GazeX = GazeXRaw;
GazeY = GazeYRaw;
GazeZ = GazeZRaw;
%% PUPIL
% Replace blanks '[]' for 'NaN' in fields diameterLeft and diameterRight
[alldata_mat(cellfun(@isempty,{alldata_mat.diameterLeft})).diameterLeft] = deal(NaN);
[alldata_mat(cellfun(@isempty,{alldata_mat.diameterRight})).diameterRight] = deal(NaN);

LDiamRaw = [alldata_mat.diameterLeft];
RDiamRaw = [alldata_mat.diameterRight];

% New artifact-removal method
standardRawSettings = rawDataFilter();
[LvalOut,LspeedFiltData,LdevFiltData] = rawDataFilter(linspace(0,length(LDiamRaw)./Param.Fs,length(LDiamRaw))',LDiamRaw',standardRawSettings);
[RvalOut,RspeedFiltData,RdevFiltData] = rawDataFilter(linspace(0,length(RDiamRaw)./Param.Fs,length(RDiamRaw))',RDiamRaw',standardRawSettings);

LDiamRaw(~LvalOut)=NaN;
RDiamRaw(~RvalOut)=NaN;

% Decide 'better' eye results
[Min, idx_decision] = min([isnan(LDiamRaw) isnan(RDiamRaw)]);
if idx_decision == 1
    Diameter = LDiamRaw';
    eyeChosen = 'Left';
elseif idx_decision == 2
    Diameter = RDiamRaw';
    eyeChosen = 'Right';
end
%% EXAMPLE

t_Gaze = linspace(0,length(alldata_mat)/Param.Fs,length(alldata_mat));

%calculate angular velocity
[az, el] = calcAngular(GazeX,GazeY,GazeZ);
gazeAz_velocity = [0;diff(az)/(1/Param.Fs)];
gazeEl_velocity = [0;diff(el)/(1/Param.Fs)];

gaze_vel=sqrt(gazeAz_velocity.^2.*cosd(el).^2+gazeEl_velocity.^2);

data = [timestamps*Param.Fs GazeX GazeY Diameter gaze_vel];

% This script shows an example of how to call the PUPILS preprocessing
% pipeline and define its options, 

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
options.screen_distance = 500;    % Screen distance in mm
options.dpi = 120;                % pixels/inches

% Inputs:
%  - data     :   dataframe (samples x categories) minimum dataframe should
%                 include [time stamps, x-coordinate, y-coordinate,
%                 pupilsize] in that order. In addition, x and y velocity
%                 values might be included.
%  - options  :   structure defining the options for the event detection
%                 algorithms, missing data interpolation and noise removal. 
%
% Outputs:
%  - data_out :   dataframe that contains the original data with appended
%                 information, in the following order: [blink information,
%                 saccade information, interpolated data, denoised data]
% 
%  - info     :   structure containing metadata regarding events and
%                 quality of the pre-processed data

[proc_data, proc_info] = processPupilData(data, options);

options.screen_distance = 100000;    % Screen distance in mm
options.dpi = 10;                % pixels/inches

[proc_data2, proc_info2] = processPupilData(data, options);

%% Fixation duration
% wsiz=50;for i = 1:size(proc_data,1);ends(i,:)=find(proc_info.fixation_ends_idx>i-wsiz/2 & proc_info.fixation_ends_idx<i+wsiz/2);starts(i,:)=find(proc_info.fixation_starts_idx>i-wsiz/2 & proc_info.fixation_starts_idx<i+wsiz/2);end

%% Data visualization

info_blinks = sprintf('blink loss: %.2f %%', proc_info.percentage_blinks );
info_sacc = sprintf('saccadic movements: %.2f %%',proc_info.percentage_saccades);
info_interp = sprintf('Interpolated data: %.2f %%',proc_info.percentage_interpolated );

cols = [110 87 115;
        212 93 121;
        234 144 133;
        233 226 208;
        112 108 97]./255;

fs = options.fs;

N = length(proc_data);
T = N/fs;
t = 0:(1/fs):T-(1/fs); % Define time vector


figure('Position', [100 100 1000 600])

subplot(3,2, [1,3,5])
sacc_col = size(proc_data, 2) - 3; 
sacc_idx = find(proc_data(:,sacc_col) ~=0);

saccades_x = proc_data(sacc_idx , 2);
saccades_y = proc_data(sacc_idx , 3);

fixations_x =  proc_data(:, 2);
fixations_x(sacc_idx) = [];

fixations_y = proc_data(:, 3);
fixations_y(sacc_idx) = [];

scatter(saccades_x, saccades_y, 'o',...
     'markerfacecolor', cols(3, :),...
     'markeredgecolor', cols(3, :),...
     'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.1); 

 hold on
 scatter(fixations_x, fixations_y, 'o',...
     'markerfacecolor', cols(1, :),...
     'markeredgecolor', cols(1, :),...
     'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.1);


le = legend('saccades', 'fixations');
set(le, 'box', 'on', 'location', 'southwest')


title('Gaze')
xlabel('x coordinates (pixels)')
ylabel('y coordinates (pixels)')


subplot(3,2,2)

plot(t, proc_data(:, 4), 'k')
ylims = get(gca, 'YLim');
axis([t(1) t(end) ylims(1) ylims(2)])

title('Original pupil traze')
xlabel('t(s)')
ylabel('Pupil diameter (\mum)')


subplot(3,2,4)

p = plot(t, proc_data(:, 4), 'k');
ylims = get(gca, 'YLim');

height = ylims(2) - ylims(1);

hold on

for i = 1:proc_info.number_of_saccades
    h(i) = rectangle('position', [proc_info.saccade_starts_s(i) ylims(1) proc_info.saccade_durations(i)  height ],...
        'facecolor', [cols(3, :) 0.4],...
        'edgecolor', 'none');
end

for i = 1:proc_info.number_of_blinks
    h2(i) = rectangle('position', [proc_info.blink_starts_s(i) ylims(1) proc_info.blink_durations(i) height ],...
        'facecolor', [cols(5, :) 1],...
        'edgecolor', 'none');
end

axis([t(1) t(end) ylims(1) ylims(2)])

hg = hggroup;
% set(h2,'Parent',hg) 
set(hg,'Displayname','Blinks')
hg2 = hggroup();
% set(h,'Parent',hg) 
set(hg2,'Displayname','Saccades')

axP = get(gca,'Position');
le = legend([hg hg2]);
set(le,'location', 'eastoutside', 'box', 'on');
set(gca, 'Position', axP)

text(5, 600, info_blinks)
text(5, 1200, info_sacc)

title('Events')
xlabel('t(s)')
ylabel('Pupil diameter (\mum)')

subplot(3,2,6)

plot(t, proc_data(:, size(proc_data, 2)), 'color', cols(2, :), 'linewidth', 1)
ylims = get(gca, 'YLim');
axis([t(1) t(end) ylims(1) ylims(2)])

text(5, 1500, info_interp)

title('Processed pupil traze')

xlabel('time (s)');
ylabel('Pupil Diameter (\mum)');