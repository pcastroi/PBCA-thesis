function [y,param] = adrl_eye_resample(x,param)

% create steep acausal lowpass with 30Hz cut-off and order of 226
fc      = 30;
n       = 226;
causal  = 0;
debug   = 0;
type    = 'low';
y       = matfiltfirzerophase(x,param.fs,n,fc,causal,type,debug);

% include new fs
[y,param.time] = preproc_resample(y,param.fsr,param.fs,param.time(1));

param.time = param.time';

end

function [xr,tn]=preproc_resample(xx,fsr,fs,ti)

% normalize to zero mean
mu      = mean(xx);
xx      = bsxfun(@minus,xx,mu);

% perform resampling
xr      = resample(xx,fsr,fs);

% get new time vector
tn      = deftime(xr,fsr,ti);

% add the mean value
xr      = bsxfun(@plus,xr,mu);

end

function t = deftime(xx,fs,ti)

% the more straightforward way:
% t = [0:1/fs:size(xx,1)/fs-1/fs]'+round(ti*fs)/fs;

t  = (0:1/fs:size(xx,1)/fs-1/fs)';

if ti<t(1) || ti>t(end) 
    error('Please consider using another resampling approach than <preproc_resample>')
end

[~,idb] = min(abs(t-ti));
t       = t + t(idb);

end

% eof