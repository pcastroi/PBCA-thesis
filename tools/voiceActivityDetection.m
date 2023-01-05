function [actArr, t,buff] = voiceActivityDetection(sig, fs,varargin)

% VOICEACTIVITYDETECTION
%    Detects voice activity in the waveform "sig" by computing the
%    squared RMS level (power) in overlapping windows of 5 ms duration with 
%    1 ms overlap, and determining speech dominant windows based on a 
%    threshold.
%    
%    Inputs:
%    - sig:      waveforms of talker                          [1x* double] 
%    - fs:       sampling rate                                [1x1 double]
%    - varargin: if true, the squared RMS level is plotted    [true/false]
%                together with the array indicating speech
%                dominant or silence domaninant windows.
%    Outputs:
%    - actArr:   binary array indicating speech dominant (1)
%                or silence domaninant (0) windows.           [2x* logical]
%    - t:        time array with mid-points of each of the 
%                windows                                      [1x* double]
%% INPUT CHECK

% Default: don't plot figures
if nargin <= 2
  plotON = false;
else 
  plotON = varargin{1};
end

if nargin < 2
  error('activity:argCheck','Too few input arguments');
else if nargin > 3
  error('activity:argCheck','Too many input arguments');
end

%% COMPUTE RMS VALUES OF THE TWO SIGNALS
% Define analysis window
tLength = 5e-3;                                                             % Window size in [s]
window  = round(fs * tLength);                                              % Window size in samples
overlap = round(fs / 1e3);                                                  % Overlap (1 ms) in samples
  
% Compute RMS^2 = power
[buff,~] = buffer(sig,window,overlap);
power    = rms(buff,1).^2;

timeVec = 0 : 1/fs : length(sig)/fs-1/fs;                                   % Time vector
[buffTime,~] =  buffer(timeVec,window,overlap);
t = mean(buffTime,1);

%% ACTIVITY DETECTION
% Thresholds
% actThresh   = 1.5e-5;
% Previusly was : actThresh   = (1/2)*0.1* mean(power)

% power(find(power>3*std(power)))= 0;
% actThresh   = 0.1* mean(power(find(power)))
actThresh   = 0.1* mean(power)
% actThresh= thr;
burstThresh = 90e-3/tLength;                                                % Threshold for small bursts that are unlikely to be speech
silThresh   = 180e-3/tLength;                                               % Threshold for defining silent gaps

% %%%%% Select a theresholf for power%%%%%%%%%%%%%%%%%%%%%

% figure(888)
% plot(power); axis tight;
% title(['click to zoom y axis'])
% [xzoom,yzoom]=ginput(1);
% figure(888)
% plot(power); axis tight; ylim([0 yzoom])
% title(['click to input the threshold for onset detection from Speech'])
% 
% [THF,THF2]=ginput(1); % MAUNALY INPUT THE THRESHOLD FOR MOVEMENT DETECTION FROM THE EMG BY CLICKING ON THE PLOT (FIGURE 888)

% actThresh= THF2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Activity
actArr = power > actThresh;
% Defining activity in channel (binary)

% Set first sample to zero to detect changes in case it's active from the
% beginning:
actArr(1) = 0;

% Bridge pauses with dur < 180 ms, remove bursts of speech with dur < 90 ms 
% and bursts of speech with RMS value smaller than the noise.
  
% BRIDGE GAPS:
% 
[on, off] = detectChanges(actArr);
%
% Define gap duration
gapDuration = abs(off(1:end-1)-on(2:end));   
%
% Remove pauses < 180 ms
for j = 1:length(gapDuration)
  if gapDuration(j) < silThresh     
    actArr(off(j):on(j+1)) = 1;
  end
end


% REMOVE IRRELEVANT BURSTS OF SPEECH:
%
% Detect changes after 180 ms bursts within a stream are bridged
%
[on, off] = detectChanges(actArr);
%
% Define burst duration
burst = off-on;
%
% Set bursts < 90 ms to zero
for i = 1:length(burst)
  if burst(i) < burstThresh 
    actArr(on(i):off(i)) = 0;

  % If burst is below a mean RMS level of 
    elseif mean(power(on(i):off(i))) < actThresh*2
      actArr(on(i):off(i)) = 0;
  end
end



%% ACTIVITY ARRAYS AND RMS WAVEFORMS STACKED
if plotON
  
  figure;
  plot(t,power); hold on; 
  plot(t,actArr.*max(power)*0.8);
  xlabel('Time [s]');

end

end


