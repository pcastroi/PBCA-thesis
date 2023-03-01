function [y] = eyeSaccade(x,pitch,param)


% find fixation and saccade
y.vector = x > param.threshold;


% set first and last 2 values to 0
y.vector(1:2) = 0;
y.vector(end-1:end) = 0;

% differentiate the saccade vector
dy = diff(y.vector);

% find start and end saccade point
y.saccadeStart = find([dy; 0]==1)+1;
y.saccadeEnd = find([dy; 0]==-1);

% add offset toi
y.saccadeStart = y.saccadeStart;
y.saccadeEnd   = y.saccadeEnd;


% preallocation for saccade amplitude
sA = zeros(size(y.saccadeStart,1),1);

time= [1:1:length(pitch)]./param.fsr;
%
if param.debug == 1
    figure(3);
    clf;
    p1 = plot(time,pitch,'b','LineWidth',2);
    hold on; grid on;
    p2 = plot(repelem(y.saccadeStart,1,2)'./param.fsr+param.startT,...
        repelem([-20 20],size(y.saccadeStart,1),1)',...
        '.--k','LineWidth',1);
    p3 = plot(repelem(y.saccadeEnd,1,2)'./param.fsr+param.startT,...
        repelem([-20 20],size(y.saccadeEnd,1),1)',...
        '.--r','LineWidth',1);
    xlabel('Time [sec]','Interpreter','latex');
    ylabel('pitch [deg]','Interpreter','latex');
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    ax.XLim = [0 260];
    ax.YLim = [-60 60];
    ax.FontSize = 14;
    ax.LineWidth = 2;
    legend([p1 p2(1) p3(1)],'pitch','Saccade Start','Saccade End','Interpreter','latex','Location','north','Orientation','horizontal');

    %     print('-f3','C:\Programming\Matlab\adrl-eeg-saeb\rawdata\code\results\figures\fig_eye_saccade','-dpng');
end

% loop over pitch angles (sum of the angles)
for n = 1:size(y.saccadeStart,1)
    sA(n) = abs(max(pitch(y.saccadeStart(n):y.saccadeEnd(n)))-min(pitch(y.saccadeStart(n):y.saccadeEnd(n))));
end

% boolean vector to remove micro & outlier saccades
microSaccade = sA > param.angleThreshold(1);
outlierSaccade = sA < param.angleThreshold(2);
keepSaccade = microSaccade & outlierSaccade;

% remove micro saccades
y.sA = sA(keepSaccade);

%
if isempty(y.sA)
    y.sA = NaN;
end

% calculate average saccade amplitude
y.asa = mean(y.sA,1);

% remove start and end values for micro saccades
y.saccadeStart = y.saccadeStart(keepSaccade);
y.saccadeEnd = y.saccadeEnd(keepSaccade);

% calculate saccade duration in samples
y.saccadeSamples = y.saccadeEnd-y.saccadeStart;

% calculate saccade duration in time
y.saccadeDur = y.saccadeSamples./param.fsr;

% get number of saccades
y.saccadeNum = size(y.saccadeDur,1);

% preallocate velocity amplitude
vA = zeros(size(y.saccadeStart,1),1);

% loop over saccades and find velocity amplitude
for n = 1:size(y.saccadeStart,1)
    vA(n) = max(x(y.saccadeStart(n):y.saccadeEnd(n)));
end

% save velocity amplitude
y.vA = vA;

%
if isempty(y.vA)
    y.vA = 0;
end

% calculate average velocity amplitude
y.ava = mean(y.vA);

end

% eof