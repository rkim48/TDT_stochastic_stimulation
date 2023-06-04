stimRate = 100;
minISI = 0.001;
multiChStimProb = [0.22 0.26 0.3 0.22];
[stim_ts,stim_ch] = generateStimTimes(stimRate,minISI,0.5,multiChStimProb);
%%



%%
function [stim_ts,stim_ch] = generateStimTimes(stimRate,minISI,trialLength,multiChStimProb)
% stimRate is set to spike rate of homogenous Poisson process
% 1/stimRate or mu is mean of exponential distribution underlying ISIs
% min ISI (s) is the minimum ISI that can occur (arises from refractory period
% of neurons) 
% trialLength (s) defines the max stimulation timestamp such that stim_ts occurs
% in interval [0, trialLength) 
mu = 1/stimRate;
ISIs = max(minISI,exprnd(mu,1,500)); % min ISI 
stim_ts = cumsum(ISIs);
stim_ts = stim_ts(stim_ts < trialLength); % get stim before trialLength

nSTIM = numel(stim_ts);
s = RandStream('mlfg6331_64');

stim_ch = nan(nSTIM,4);
for i = 1:nSTIM
    nStimCh = datasample(s,1:4,1,'Weights',multiChStimProb);
    stim_ch(i,1:nStimCh) = randsample(s,1:32,nStimCh);
end


end