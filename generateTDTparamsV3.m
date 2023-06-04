stimRate = 100;
minISI = 0.001;
multiChStimProb = [0.22 0.26 0.3 0.22];
[stim_ts,stim_ch] = generateStimTimes(stimRate,minISI,0.5,multiChStimProb);

%%
load('columnNamesAll.mat','columnNames');

nROWS = size(stim_ch,1);
period = 0.4; % ms, want period to be longer than pulseDuration 
current = 5;
pulseDuration = period;
count = 1;
delay = 0;

periodColumns = nan(nROWS,4);
countColumns = nan(nROWS,4);
durColumns = nan(nROWS,4);
ampColumns = nan(nROWS,4);
delayColumns = nan(nROWS,4);
chColumns = nan(nROWS,4);

for i = 1:nROWS
    insertCh = stim_ch{i};
    numInsertCh = numel(insertCh);
    
    periodColumns(i,1:numInsertCh) = period;
    ampColumns(i,1:numInsertCh) = current;
    countColumns(i,1:numInsertCh) = count;
    durColumns(i,1:numInsertCh) = pulseDuration;
    delayColumns(i,1:numInsertCh) = delay;
    chColumns(i,1:numInsertCh) = insertCh(1:numInsertCh);
end

% Create params array to store parameters for single trial 
paramsArr = nan(nROWS,16);
periodColIdx = find(contains(columnNames,'Period')==1);
countColIdx = find(contains(columnNames,'Count')==1);
ampColIdx = find(contains(columnNames,'Amp')==1);
durColIdx = find(contains(columnNames,'Dur')==1);
delayColIdx = find(contains(columnNames,'Delay')==1);
chColIdx = find(contains(columnNames,'Chan')==1);

paramsArr(:,periodColIdx) = periodColumns;
paramsArr(:,ampColIdx) = ampColumns;
paramsArr(:,countColIdx) = countColumns;
paramsArr(:,durColIdx) = durColumns;
paramsArr(:,delayColIdx) = delayColumns;
paramsArr(:,chColIdx) = chColumns;

paramTable = array2table(paramsArr,'RowNames',string(1:nROWS));
paramTable.Properties.VariableNames = string(columnNames);
writetable(paramTable,'paramArray.par.csv','WriteRowNames',true);

%% Generate sequence file
timestamp_convert = [stim_ts(1) ; diff(stim_ts)]*1000;
columnNames = {'Seq-1','Time-1'};
seqArr = [(1:nROWS)' timestamp_convert];
seqTable = array2table(seqArr,'RowNames',string(1:nROWS));
seqTable.Properties.VariableNames = string(columnNames);
writetable(seqTable,'paramArray.seq.csv','WriteRowNames',true);

%% Plot trial 

plotTrial(stim_ts,stim_ch)

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
stim_ts = stim_ts(stim_ts < trialLength)'; % get stim before trialLength

nSTIM = numel(stim_ts);
s = RandStream('mlfg6331_64');

stim_ch = cell(nSTIM,1);
for i = 1:nSTIM
    nStimCh = datasample(s,1:4,1,'Weights',multiChStimProb);
    stim_ch{i} = randsample(s,1:32,nStimCh);
end

end

function plotTrial(stim_ts,stim_ch)
    nSTIM = size(stim_ch,1);
    for i = 1:nSTIM
        stim_time = stim_ts(i);
        ch_stim = stim_ch{i};
        nCH = numel(ch_stim);
        scatter(repmat(stim_time,1,nCH),ch_stim,'k','Marker','|','LineWidth',1);
        hold on;
    end
    ylim([0.5 32.5])
    set(gca, 'YDir','reverse')
    set(gca().YAxis,'TickLength',[0 0])
    set(gca().XAxis,'TickLength',[0 0])
    ylabel('Stim channel')
    xlabel('Time (s)')

end