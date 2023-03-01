function [param] = eyeFilter(x,param)

% create steep acausal lowpass with 5Hz cut-off and order of 113
fc                  = 9;
n                   = 42;
causal              = 0;
debug               = 0;
type                = 'low';
param.velocityFilt  = matfiltfirzerophase(x,param.fsr,n,fc,causal,type,debug);

end

% eof