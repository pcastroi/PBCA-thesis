function [param] = eyeBandpass(x,param)

% Filter data in range of [1 15] Hz

% apply zerosphase highpass filter
fc      = param.hp.fc;
n       = param.hp.n;
causal  = 0;
debug   = 0;
type    = 'high';
yh      = matfiltfirzerophase(x,param.fsr,n,fc,causal,type,debug);

% apply zerosphase lowpass filter
fc      = param.lp.fc;
n       = param.lp.n;
causal  = 0;
debug   = 0;
type    = 'low';
yl      = matfiltfirzerophase(yh,param.fsr,n,fc,causal,type,debug);

% save bp filtered pitch
param.pitchBP = yl;

end

% eof