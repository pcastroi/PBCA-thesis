clear all; clc; close all;

addpath(genpath('C:\git'));

Param.Fs=50;
alldata=load('C:\git\PBCA-thesis\data\AMEND_II\Pair03\HI\AAHI_N0.mat');
alldata_mat = cell2mat(alldata.data);

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
    Diameter = LDiamRaw;
    eyeChosen = 'Left';
elseif idx_decision == 2
    Diameter = RDiamRaw;
    eyeChosen = 'Right';
end
%% EXAMPLE

t_Gaze = linspace(0,length(alldata_mat)/Param.Fs,length(alldata_mat));

% common parameters
nTrials = 120;
sample_window = 0.02; %2ms i.e., 500Hz tracking
nEdgeSamples = 50; %number of samples to discard either side of blink
tVelocity = 30; %threshold velocity for saccades (deg/s)
tAcceleration = 8000; %threshold acceleration for saccades (deg/s^2)
tDist = 0.5; %threshold (Euclidean) distance for saccades (deg)

WExt = 10;

%%%%%%% (0) %%%%%%%%%%

% find start and end times
t0 = t_Gaze(WExt); %start time 
tEnd = t_Gaze(end-WExt); %end time

% ...and sample indices closest t0/tEnd
t0_time = abs(t_Gaze - t0); %t_Gaze is common to all trials
tEnd_time = abs(t_Gaze - tEnd);
[~, t0_ind] = min(t0_time);
[~, tEnd_ind] = min(tEnd_time);
t0_ind = t0_ind - 4;
tEnd_ind = tEnd_ind + 4;
% ...extend window by 4 samples either side to ensure velocity
% and acceleration calculation below covers entire stimulation
% period (vel/acc sliding windows chop off 2 samples each
% at either end of calculation)

%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% (1) %%%%%%%%%%

% find baseline eye position (must remove NaNs from xPos
% samples later for calculation of median, but not for
% subtraction from and extraction of relevant samples below)
baselineX = GazeX(1:(t0_ind - 1));
baselineY = GazeY(1:(t0_ind - 1));
baselineX_NaN = isnan(baselineX);
baselineY_NaN = isnan(baselineY);
baselineX(baselineX_NaN) = [];
baselineY(baselineY_NaN) = [];

% get relevant eye position data aligned to baseline, and
% also relevant pupil samples
xPos = GazeX(t0_ind:tEnd_ind) - median(baselineX);
yPos = GazeY(t0_ind:tEnd_ind) - median(baselineY);
pupil = Diameter(t0_ind:tEnd_ind);

%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% (2) %%%%%%%%%%

% blink detection/removal
noPupil = isnan(pupil);
eyedat.blink = 0;

% if trial has blink
if max(noPupil) == 1

    % note instance of blink
    eyedat.blink = 1;

    % find blink boundaries
    blinkStart = find(diff(noPupil) == 1) + 1; %(index+1 was first NaN, due to differentiation
    blinkEnd = find(diff(noPupil) == -1);

    % blink that started before first interval will not show up
    % in blinkStart boundary search above
    if isempty(blinkStart) == 1
        blinkStart = 1;
    end
    % blink that ended after second interval will not show up
    % in blinkEnd boundary search above
    if isempty(blinkEnd) == 1
        blinkEnd = length(noPupil);
    end

    % include samples just before/after blink,
    % accounting for timepoints that move outside analysis
    % window
    blinkStart = blinkStart - nEdgeSamples;
    blinkEnd = blinkEnd + nEdgeSamples;
    for blink = 1:length(blinkStart)
        if blinkStart(blink) < 1
            blinkStart(blink) = 1;
        end
    end
    for blink = 1:length(blinkEnd)
        if blinkEnd(blink) > length(noPupil)
            blinkEnd(blink) = length(noPupil);
        end
    end

    % set samples within the boundaries defined above to NaN;
    % code below should handle both single and multi-blink trials, and
    % trials in which the first or last blink extends outside
    % boundaries
    if length(blinkStart) == length(blinkEnd)

        for blink = 1:length(blinkStart)
            noPupil(blinkStart(blink):blinkEnd(blink)) = 1;
        end

    elseif length(blinkStart) > length(blinkEnd)

        for blink = 1:(length(blinkStart) - 1) % all but final blink in trial
            noPupil(blinkStart(blink):blinkEnd(blink)) = 1;
        end
        for blink = length(blinkStart) % final blink
            noPupil(blinkStart(blink):length(noPupil)) = 1;
        end

    elseif length(blinkStart) < length(blinkEnd)

        for blink = 1 % first blink in trial
            noPupil(1:blinkEnd(blink)) = 1;
        end
        for blink = 2:(length(blinkEnd)) % any additional blinks
            noPupil(blinkStart(blink - 1):blinkEnd(blink)) = 1;
        end

    end

    %...and blink samples are removed below (after velocity
    %and acceleration calculation). Removing prior to this
    %would create possibility of blink window edges being
    %accidentally marked as saccades.

end

%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% (3) %%%%%%%%%%

% calculate velocity/acceleration using 5-sample window. See
% Engbert and Kliegl, 2003. Denominator accounts for the
% six sample 'differences' used in numerator (i.e., n-2 to
% n+2 = 4 samples, n-1 to n+1 = 2 samples).
xVel = zeros(size(xPos)); yVel = zeros(size(yPos));
for ii = 3:(size(xPos, 1) - 2) % 2 additional samples chopped off either end (see ~line 230 above)
    xVel(ii) = (xPos(ii + 2) + xPos(ii + 1) - xPos(ii - 1) - xPos(ii - 2))/(6*sample_window);
    yVel(ii) = (yPos(ii + 2) + yPos(ii + 1) - yPos(ii - 1) - yPos(ii - 2))/(6*sample_window);
end
euclidVel = sqrt((xVel.*xVel) + (yVel.*yVel));

xAcc = zeros(size(xPos)); yAcc = zeros(size(yPos));
for ii = 3:(size(xVel, 1) - 2) % 2 additional samples chopped off either end (see ~line 230 above)
    xAcc(ii) = (xVel(ii + 2) + xVel(ii + 1) - xVel(ii - 1) - xVel(ii - 2))/(6*sample_window);
    yAcc(ii) = (yVel(ii + 2) + yVel(ii + 1) - yVel(ii - 1) - yVel(ii - 2))/(6*sample_window);
end
euclidAcc = sqrt((xAcc.*xAcc) + (yAcc.*yAcc));

% remove blink samples from pos/velocity/acceleration
% variables, so that blink edges are not accidentally
% marked as saccades below
xPos(noPupil == 1) = [];
yPos(noPupil == 1) = [];
xVel(noPupil == 1) = [];
yVel(noPupil == 1) = [];
xAcc(noPupil == 1) = [];
yAcc(noPupil == 1) = [];
euclidVel(noPupil == 1) = [];
euclidAcc(noPupil == 1) = [];

% save to post-processed eye data structure (pupil
% variables may now contain more samples than the others)
eyedat.pupil = pupil;
eyedat.xPos = xPos;
eyedat.yPos = yPos;
eyedat.xVel = xVel;
eyedat.yVel = yVel;
eyedat.xAcc = xAcc;
eyedat.yAcc = yAcc;
eyedat.euclidAcc = euclidAcc;
eyedat.euclidVel = euclidVel;
eyedat.noPupil = noPupil;

%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% (4) %%%%%%%%%%

% saccade detection
eyedat.saccades = [];
candidates = find(euclidVel > tVelocity);
if candidates

    % check for multiple candidate saccades in single
    % trial, using threshold parameters defined at top
    % (see Engbert & Kliegl papers, and Eyelink manual)
    saccades = [];
    diffCandidates = diff(candidates);
    breaks = [0 find(diffCandidates > 1) size(candidates, 2)];
    for jj = 1:(size(breaks, 2) - 1)

        % find individual candidate saccades
        saccade = [candidates(breaks(jj) + 1) candidates(breaks(jj + 1))];

        % exceeds acceleration threshold?
        peakAcceleration = max(euclidAcc(saccade(1):saccade(2)));
        if peakAcceleration > tAcceleration

            % exceeds amplitude threshold?
            xDist = xPos(saccade(2)) - xPos(saccade(1));
            yDist = yPos(saccade(2)) - yPos(saccade(1));
            euclidDist = sqrt((xDist*xDist) + (yDist*yDist));
            if euclidDist > tDist

                % store saccade info
                peakVelocity = max(euclidVel(saccade(1):saccade(2)));
                saccades = [saccades; saccade xDist yDist euclidDist peakVelocity];
            end

        end

    end

    % save saccade data
    eyedat.saccades = saccades;
end

%% Plot
figure(1)

subplot(3,1,1)
plot(eyedat.euclidVel, 'b', 'LineWidth', 2);
ylabel('Velocity (deg/s)')
ylim([0 150]);
line([0 700], [30 30], 'LineStyle', '--', 'Color', 'r')

subplot(3,1,2)
plot(eyedat.euclidAcc, 'm', 'LineWidth', 2);
ylabel('Acceleration (deg/s^2)')
ylim([0 15000]);
line([0 700], [8000 8000], 'LineStyle', '--', 'Color', 'r')

subplot(3,1,3)
plot(eyedat.xPos, 'k', 'LineWidth', 2);
xlabel('Sample')
ylabel('Position (deg)')
ylim([-2 2]);