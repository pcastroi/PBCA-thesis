function [y] = eyeFixation(x,param)

% differentiate the saccade vector
dy = diff(x);

% find start and end fixation point
y.fixationStart = [1; find([0; dy]==-1)];
y.fixationEnd = [find([dy; 0]==1); size(x,1)];

% calculate fixation duration in samples
y.fixationSamples = y.fixationEnd-y.fixationStart;

% calculate fixation duration in time
y.fixationDur = y.fixationSamples./param.fsr;

% remove fixation duration below fixation threshold
shortFixations = y.fixationDur > param.fixThreshold(1);
longFixations = y.fixationDur < param.fixThreshold(2);
keepFixation = shortFixations & longFixations;

% filter out fixation values
y.fixationStart = y.fixationStart(keepFixation);
y.fixationEnd = y.fixationEnd(keepFixation);
y.fixationSamples = y.fixationSamples(keepFixation);
y.fixationDur = y.fixationDur(keepFixation);

% get number of fixations
y.fixationNum = size(y.fixationDur,1);

% calculate the average fixation duration
y.afd = mean(y.fixationDur,1);

end

% eof