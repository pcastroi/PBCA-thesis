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

% % Preprocessing gaze3d - Setting outliers as NaNs (remove artifacts)
% XThresh = [mean(GazeX,'omitnan')-std(GazeX,'omitnan'),mean(GazeX,'omitnan')+std(GazeX,'omitnan')];
% YThresh = [mean(GazeY,'omitnan')-std(GazeY,'omitnan'),mean(GazeY,'omitnan')+std(GazeY,'omitnan')];
% ZThresh = [mean(GazeZ,'omitnan')-std(GazeZ,'omitnan'),mean(GazeZ,'omitnan')+std(GazeZ,'omitnan')];
% F_NOutl = 1; % number of outliers per file
% 
% for s=1:length(GazeX)
%     if GazeX(s) < XThresh(1) || GazeX(s) > XThresh(2)
%         GazeX(s)=NaN;
%         F_NOutl = F_NOutl + 1;
%     end
%     if GazeY(s) < YThresh(1) || GazeY(s) > YThresh(2)
%         GazeY(s)=NaN;
%         F_NOutl = F_NOutl + 1;
%     end
%     if GazeZ(s) < ZThresh(1) || GazeZ(s) > ZThresh(2)
%         GazeZ(s)=NaN;
%         F_NOutl = F_NOutl + 1;
%     end
% end

MaxVisualAngle = 2; % [degrees], critical visual angle for fixation definition
LPWinSize = 0.5; % [s]: Window size of hamming-window for low-pass filtering
FilterWidth = round((LPWinSize*Param.Fs)/2); % [samples]: Width of hamming filter used for fixation duration

% Fixation duration
Fixation = fix_duration([GazeX,GazeY,GazeZ],MaxVisualAngle,Param.Fs);

% Filtering - Fixation duration with hamming window (F = 3)
Fixation = ndnanfilter(Fixation,'hamming',FilterWidth);

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
[Min, idx_decision] = min([sum(isnan(LDiamRaw)) sum(isnan(RDiamRaw))]);
if idx_decision == 1
    Diameter = LDiamRaw';
    eyeChosen = 'Left';
elseif idx_decision == 2
    Diameter = RDiamRaw';
    eyeChosen = 'Right';
end
%% EXAMPLE
%calculate angular velocity
[az, el] = calcAngular(GazeX,GazeY,GazeZ);
gazeAz_velocity = [0;diff(az)/(1/Param.Fs)];
gazeEl_velocity = [0;diff(el)/(1/Param.Fs)];
gaze_vel=sqrt(gazeAz_velocity.^2.*cosd(el).^2+gazeEl_velocity.^2);

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

% Inputs:
%  - data     :   dataframe (samples x categories) minimum dataframe should
%                 include [time stamps, x-coordinate, y-coordinate,
%                 pupilsize, ] in that order. In addition, x and y velocity
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

data = [timestamps*Param.Fs GazeX GazeY Diameter gaze_vel];
[proc_data, proc_info] = processPupilData(data, options);

%% Fixation duration
FixWDur=2;
FixWN=FixWDur*Param.Fs;
FixData = [zeros(ceil(FixWN/2),1); proc_data(:,8) ; zeros(ceil(FixWN/2),1)];
fix_hel = zeros(length(proc_data(:,8)),1);
j = ceil(FixWN/2):length(FixData)-ceil(FixWN/2);
for i = 1:size(proc_data,1)
    fix_hel(i,1)=nnz(FixData(j(i)-ceil(FixWN/2)+1:j(i)+ceil(FixWN/2)-1))/Param.Fs; % Count non-zero elements
end

GazeX_fix = GazeX; GazeY_fix = GazeY; GazeZ_fix = GazeZ;
GazeX_fix(proc_data(:,8)==0)=NaN; GazeY_fix(proc_data(:,8)==0)=NaN; GazeZ_fix(proc_data(:,8)==0)=NaN;
Combined_Fixation = fix_duration([GazeX_fix,GazeY_fix,GazeZ_fix],MaxVisualAngle,Param.Fs);
Combined_Fixation = ndnanfilter(Combined_Fixation,'hamming',FilterWidth);
% plot(fix_hel/max(fix_hel))
figure;grid on;hold on;plot(Fixation/max(Fixation));plot(Combined_Fixation/max(Combined_Fixation));title(['Difference between Fix and Combined fix = ',sprintf('%2f',mean(Fixation-Combined_Fixation,'omitnan'))])
figure;grid on;hold on;plot(Diameter);plot(proc_data(:,10))
%% Data visualization

info_blinks = sprintf('blink loss: %.2f %%', proc_info.percentage_blinks );
info_sacc = sprintf('saccadic movements: %.2f %%',proc_info.percentage_saccades);
info_interp = sprintf('Interpolated data: %.2f %%',proc_info.percentage_interpolated );

cols = [53 155 67;
        212 93 121;
        234 144 133;
        233 226 208;
        112 108 97;
        30, 30, 31;
        111 185 191;
        242 183 5]./255;

fs = options.fs;

N = length(proc_data);
T = N/fs;
t = 0:(1/fs):T-(1/fs); % Define time vector


figure('Position', [100 100 1000 600])

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


figure
% blinks
subplot(3,1,1)
hold on
plot(t, proc_data(:, 4), 'k','LineWidth',1,'HandleVisibility','off')
% yline(mean(proc_data(:, 4),'omitnan'), '--k')
axis([t(1) t(end) nanmean(proc_data(1:500, 4))-3*nanstd(proc_data(1:500, 4))-0.1 max(proc_data(1:Param.Fs*10, 4))+0.1])
xlim([0 8])
ylims1 = get(gca, 'YLim');
height = ylims1(2) - ylims1(1);
for i = 1:proc_info.number_of_blinks
    rectangle('position', [proc_info.blink_starts_s(i) ylims1(1) proc_info.blink_durations(i) height ],...
        'facecolor', [cols(6, :) 0.4],...
        'edgecolor', 'none');
end
ylim(ylims1)
% title('Blink detection')
ylabel('Pupil diameter [mm]')
grid on
line(NaN,NaN,'linewidth',5,'Color',[cols(6, :) 0.4])
yline((nanmean(proc_data(1:500, 4))-3*nanstd(proc_data(1:500, 4))),'--') %,'Threshold','LabelHorizontalAlignment', 'center'
legend('Blinks','Threshold')

% saccades
subplot(3,1,2)
hold on
plot(t,gaze_vel,'k','LineWidth',1,'HandleVisibility','off');
axis([t(1) t(end) min(gaze_vel(1:Param.Fs*10))-10 max(gaze_vel(1:Param.Fs*10))+10])
xlim([0 8])
ylims = get(gca, 'YLim');
height = ylims(2) - ylims(1);
for i = 1:proc_info.number_of_saccades
        rectangle('position', [proc_info.saccade_starts_s(i) ylims(1) proc_info.saccade_durations(i)  height ],...
        'facecolor', [cols(7, :) 0.4],...
        'edgecolor', 'none');
end
% title('Saccade detection')
ylabel('Angular velocity [°/s]')
grid on
line(NaN,NaN,'linewidth',5,'Color',[cols(7, :) 0.4])
yline(options.vel_threshold,'--') %,'Threshold','LabelHorizontalAlignment', 'center'
legend('Saccades','Threshold')
% legend('$\textsf{Angular velocity: } \: v_{gaze}$','$\textsf{Saccade events}$','Interpreter','latex')

% Combined events
subplot(3,1,3)
hold on
plot(t, proc_data(:, size(proc_data, 2)),'Color', cols(2, :),'LineWidth',1,'HandleVisibility','off')
axis([t(1) t(end) min(proc_data(1:Param.Fs*10, size(proc_data, 2)))-0.1 max(proc_data(1:Param.Fs*10, size(proc_data, 2)))+0.1])
xlim([0 8])
ylims3 = ylims1;
height = ylims3(2) - ylims3(1);
% for i = 1:proc_info.number_of_blinks
%     rectangle('position', [proc_info.blink_starts_s(i) ylims(1) proc_info.blink_durations(i) height ],...
%         'facecolor', [cols(5, :) 0.4],...
%         'edgecolor', 'none');
% end
% for i = 1:proc_info.number_of_saccades
%         rectangle('position', [proc_info.saccade_starts_s(i) ylims(1) proc_info.saccade_durations(i)  height ],...
%         'facecolor', [cols(3, :) 0.4],...
%         'edgecolor', 'none');
% end
for i = 1:proc_info.number_of_fixations
    rectangle('position', [proc_info.fixation_starts_s(i) ylims3(1) proc_info.fixation_durations(i) height ],...
        'facecolor', [cols(8, :) 0.2],...
        'edgecolor', 'none');
end
ylim(ylims3)
grid on
% line(NaN,NaN,'linewidth',5,'Color',[cols(5, :) 0.4])
% line(NaN,NaN,'linewidth',5,'Color',[cols(3, :) 0.4])
line(NaN,NaN,'linewidth',5,'Color',[cols(8, :) 0.4])
legend('Fixations')
% title('Fixation detection')
xlabel('Time [s]')
ylabel('Pupil diameter [mm]')
