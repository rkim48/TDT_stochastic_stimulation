timestamps = [1000 2000 3000]; % ms
ch = cell(3,1);
ch{1} = 1; 
ch{2} = 2;
ch{3} = [1 2 3];
nROWS = length(timestamps);
%% Define experiment variables
nTRIALS = 100;
nCH = 32;
currents = 1:5; % uA
pulseDurations = [0.4,0.5,0.6,0.7]; % ms (includes phases and IPI)
stimRates = 5:10; % Hz
trialLength = 0.5; % s
minISI = 0.001; % s, minimum inter-stim interval
%% Randomize 

T = combvec(currents,pulseDurations,stimRates)';
T = T(randperm(size(T,1)), :); 
T = T(1:nTRIALS,:);

% Stim timestamps for a single trial
bigStimTsArr = cell(nTRIALS,1);
for i = 1:nTRIALS
    row = T(i,:);
    stimRate = row(3);
    stimTsArr = cell(32,1);
    for j = 1:nCH
        stimTsArr{j} = generateStimTimes(stimRate,minISI,0.5);
    end
    bigStimTsArr{i} = stimTsArr;
end

[chVec,tsVec] = transformStimTsArr(stimTsArr);

[chColumns,tsColumns,nROWS,rowTime] = groupChAndTs(chVec,tsVec*1000);

%% Generate paramater file
load('columnNamesAll.mat','columnNames');

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
    insertCh = ch{i};
    numInsertCh = numel(insertCh);
    
    periodColumns(i,1:numInsertCh) = period;
    ampColumns(i,1:numInsertCh) = current;
    countColumns(i,1:numInsertCh) = count;
    durColumns(i,1:numInsertCh) = pulseDuration;
    delayColumns(i,1:numInsertCh) = delay;
    chColumns(i,1:numInsertCh) = insertCh;
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
timestamp_convert = [timestamps(1) diff(timestamps)];
columnNames = {'Seq-1','Time-1'};
seqArr = [(1:nROWS)' timestamp_convert'];
seqTable = array2table(seqArr,'RowNames',string(1:nROWS));
seqTable.Properties.VariableNames = string(columnNames);
writetable(seqTable,'paramArray.seq.csv','WriteRowNames',true);

%% Helper functions

function stim_ts = generateStimTimes(stimRate,minISI,trialLength)
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
end

function [chVec,tsVec] = transformStimTsArr(stimTsArr)

    chVec = [];
    tsVec = [];
    for i = 1:size(stimTsArr,1)
        stimTs = stimTsArr{i};
        tsVec = [tsVec stimTs];
        chVec = [chVec i * ones(1,numel(stimTs))];
    end
    [~,b] = sort(tsVec);
    tsVec = tsVec(b);
    chVec = chVec(b);
end

function [chColumns,tsColumns,nROWS,rowTime] = groupChAndTs(chVec,tsVec)

    nTS = numel(tsVec);
    nROWS = ceil(nTS/4);
    nPAD = nROWS * 4 - nTS;
    tsVec = [tsVec nan(1,nPAD)];
    chVec = [chVec nan(1,nPAD)];
    chColumns = reshape(chVec,4,[])';
    tsColumns = reshape(tsVec,4,[])';

    lastTs = [0 ; tsColumns(1:end-1,4)];
    tsColumns = tsColumns - lastTs;
    rowTime = [0 ; tsColumns(1:end-1,4)];
end