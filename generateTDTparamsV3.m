stimRate = 50;
stimRate = 25;
% stimRate = 10;
% scaleParam = 0.4;
minISI = 0.0015;
multiChStimProb = [0.22 0.26 0.3 0.22];
nTRIALS = 10; % number of stim trials
goodChannels = 1:32;
[stim_ts,stim_ch] = generateStimTimes(stimRate,minISI,nTRIALS,multiChStimProb,goodChannels);
% [stim_ts,stim_ch] = generateStimTimesGamma(stimRate,scaleParam,minISI,window,multiChStimProb);

% subplot(1,2,1)
plotTrial(stim_ts,stim_ch)
xline([0:0.5:10])
% subplot(1,2,2)
% histogram(gamrnd(stimRate,scaleParam,1,1000))

% high sparsity: mu = 10
% medium sparsity: mu = 25
% low sparsity: mu = 50

%% Experiment 1 - single trial
good_channels = 1:32;
stimRate = 100;
nTRIALS = 1;
[stim_ts_arr,stim_ch_arr,trial_id_arr] = generateStimTimes(stimRate,minISI,nTRIALS,multiChStimProb,good_channels);
plotTrial(stim_ts_arr,stim_ch_arr)

groupedArr = groupStimByTrial(stim_ch_arr,stim_ts_arr,trial_id_arr)
%% Experiment 2 - 30 stim trials (half min duration)
stimRate = 50;
nTRIALS = 30;
[stim_ts_arr,stim_ch_arr,trial_id_arr] = generateStimTimes(stimRate,minISI,nTRIALS,multiChStimProb,good_channels);
% plotTrial(stim_ts_arr,stim_ch_arr)
% xline(0:0.5:nTRIALS);
groupedArr = groupStimByTrial(stim_ch_arr,stim_ts_arr,trial_id_arr);

%% Experiment 3 - 60 stim trials (1 min duration) with different "sparsitites" parametrized by exp mu
stimRate1 = 50;
stimRate2 = 25;
stimRate3 = 10;
nTRIALS = 10;
[stim_ts_arr1,stim_ch_arr1,trial_id_arr1] = generateStimTimes(stimRate1,minISI,nTRIALS,multiChStimProb,good_channels);
[stim_ts_arr2,stim_ch_arr2,trial_id_arr2] = generateStimTimes(stimRate2,minISI,nTRIALS,multiChStimProb,good_channels);
[stim_ts_arr3,stim_ch_arr3,trial_id_arr3] = generateStimTimes(stimRate3,minISI,nTRIALS,multiChStimProb,good_channels);
groupedArr1 = groupStimByTrial(stim_ts_arr1,stim_ch_arr1,trial_id_arr1);
groupedArr2 = groupStimByTrial(stim_ts_arr2,stim_ch_arr2,trial_id_arr2);
groupedArr3 = groupStimByTrial(stim_ts_arr3,stim_ch_arr3,trial_id_arr3);

allGroupedArr = [groupedArr1; groupedArr2; groupedArr3];
[m,n] = size(allGroupedArr);
idx = randperm(m);
allGroupedArr = allGroupedArr(idx,:);

for i = 1:3
subplot(3,1,i)
if i == 1
    plotTrial(stim_ts_arr1,stim_ch_arr1)
elseif i == 2
    plotTrial(stim_ts_arr2,stim_ch_arr2)
elseif i == 3
    plotTrial(stim_ts_arr3,stim_ch_arr3)
end
xline(0:0.5:nTRIALS);
end
%%
% shuffled
stim_ch_arr = {allGroupedArr{:,2}};
stim_ch_arr = vertcat(stim_ch_arr{:});
stim_ts_arr = cell2mat(allGroupedArr(:,1));
plotTrial(stim_ts_arr,stim_ch_arr)
xline(0:0.5:30);

%%
load('columnNamesAll.mat','columnNames');

nROWS = size(stim_ch_arr,1);
period = 0.4; % ms, want period to be longer than pulseDuration 
current = -10;
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
    insertCh = stim_ch_arr{i};
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
timestamp_convert = [diff(stim_ts_arr)*1000; 0];
columnNames = {'Seq-1','Time-1'};
seqArr = [(1:nROWS)' timestamp_convert];
seqTable = array2table(seqArr,'RowNames',string(1:nROWS));
seqTable.Properties.VariableNames = string(columnNames);
writetable(seqTable,'paramArray.seq.csv','WriteRowNames',true);

%% Plot trial 

plotTrial(stim_ts_arr,stim_ch_arr)

%%
function [stim_ts_arr,stim_ch_arr,trial_id_arr] = generateStimTimes(stimRate,minISI,nTRIALS,multiChStimProb,goodChannels)
% stimRate is set to spike rate of homogenous Poisson process
% 1/stimRate or mu is mean of exponential distribution underlying ISIs
% min ISI (s) is the minimum ISI that can occur (arises from refractory period
% of neurons) 
% trialLength (s) defines the max stimulation timestamp such that stim_ts occurs
% in interval [0, trialLength) 
mu = 1/stimRate;
trialLength = 0.5;

stim_ts_arr = [];
stim_ch_arr = {};
trial_id_arr = [];
s = RandStream('mlfg6331_64');
for trial = 1:nTRIALS
    ISIs = max(minISI,exprnd(mu,1,1000)); % min ISI 
    stim_ts = cumsum(ISIs);
    stim_ts = stim_ts(stim_ts < trialLength)'; % get stim before trialLength
    stim_ts = stim_ts + trial - 1; 
    nSTIM = numel(stim_ts);

    trial_id = trial * ones(nSTIM,1);

    stim_ch = cell(nSTIM,1);
    for i = 1:nSTIM
        nStimCh = datasample(s,1:4,1,'Weights',multiChStimProb);
        stim_ch{i} = randsample(s,goodChannels,nStimCh);
    end
    
    stim_ch_arr = [stim_ch_arr; stim_ch];
    stim_ts_arr = [stim_ts_arr; stim_ts];
    trial_id_arr = [trial_id_arr; trial_id];
end
end

function groupedArr = groupStimByTrial(stim_ch_arr,stim_ts_arr,trial_id_arr)

    uniqueTrialID = unique(trial_id_arr);
    nTRIALS = numel(uniqueTrialID);
    groupedArr = cell(nTRIALS,2);
    for i = 1:nTRIALS
        trialID = uniqueTrialID(i);
        idx = find(trial_id_arr == trialID);
        groupedArr{i,1} = stim_ch_arr(idx);
        groupedArr{i,2} = stim_ts_arr(idx);
    end
end

function [stim_ts,stim_ch] = generateStimTimesGamma(stimRate,scaleParam,minISI,trialLength,multiChStimProb)
% stimRate is set to spike rate of homogenous Poisson process
% 1/stimRate or mu is mean of exponential distribution underlying ISIs
% min ISI (s) is the minimum ISI that can occur (arises from refractory period
% of neurons) 
% trialLength (s) defines the max stimulation timestamp such that stim_ts occurs
% in interval [0, trialLength) 
a = 1/stimRate; % shape parameter
b = scaleParam; % scale parameter
ISIs = max(minISI,gamrnd(a,b,1,trialLength*1000)); % min ISI 
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
    % xlim([0 0.5])
    set(gca, 'YDir','reverse')
    set(gca().YAxis,'TickLength',[0 0])
    set(gca().XAxis,'TickLength',[0 0])
    ylabel('Stim channel')
    xlabel('Time (s)')

end