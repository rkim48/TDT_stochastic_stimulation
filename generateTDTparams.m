load('columnNames.mat','columnNames');

%% Define experiment variables
nTRIALS = 100;
nCH = 32;
currents = 1:5; % uA
pulseDurations = [0.4,0.5,0.6,0.7]; % ms (includes phases and IPI)
stimRates = 5:10; % Hz
trialLength = 0.5; % s
minISI = 0.001; % s, minimum inter-stim interval
%% Randomize 

T = combinations(currents,pulseDurations,stimRates);
T = T(randperm(size(T,1)), :); 
T = T(1:nTRIALS,:);

% Stim timestamps for a single trial
bigStimTsArr = cell(nTRIALS,1);
for i = 1:nTRIALS
    row = T(i,:);
    stimRate = row.stimRates;
    stimTsArr = cell(32,1);
    for j = 1:nCH
        stimTsArr{j} = generateStimTimes(stimRate,minISI,0.5);
    end
    bigStimTsArr{i} = stimTsArr;
end

subplot(1,2,1);
plotStimTimes(stimTsArr,stimRate)
subplot(1,2,2);
plotISI(stimTsArr);
%%
[chVec,tsVec] = transformStimTsArr(stimTsArr);

%%
[chColumns,tsColumns,nROWS,rowTime] = groupChAndTs(chVec,tsVec*1000);

%% Trial parameters
trial_idx = 1;
row = T(trial_idx,:);
stimRate = row.stimRates;
pulseDuration = row.pulseDurations;
current = row.currents;
period = 1; % ms

% Create params array to store parameters for single trial 
paramsArr = nan(nROWS,24);
periodColIdx = find(contains(columnNames,'Period')==1);
countColIdx = find(contains(columnNames,'Count')==1);
ampColIdx = find(contains(columnNames,'Amp')==1);
durColIdx = find(contains(columnNames,'Dur')==1);
delayColIdx = find(contains(columnNames,'Delay')==1);
chColIdx = find(contains(columnNames,'Chan')==1);

paramsArr(:,periodColIdx) = period * ones(nROWS,4);
paramsArr(:,countColIdx) = 1 * ones(nROWS,4);
paramsArr(:,ampColIdx) = current * ones(nROWS,4);
paramsArr(:,durColIdx) = pulseDuration * ones(nROWS,4);
paramsArr(:,delayColIdx) = tsColumns;
paramsArr(:,chColIdx) = chColumns;

%% Write .par file

paramTable = array2table(paramsArr,'RowNames',string(1:nROWS));
paramTable.Properties.VariableNames = string(columnNames);
writetable(paramTable,'paramArray.par.csv','WriteRowNames',true);

%% Write .seq file

columnNames = {'Seq-1','Time-1'};
seqArr = [(1:nROWS)' rowTime];
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

function plotStimTimes(arr,stimRate)

for i = 1:size(arr,1)
    ch_stim = arr{i};
    scatter(ch_stim,i*ones(1,numel(ch_stim)),'k','Marker','|');
    hold on;
end
ylim([0.5 32.5])
set(gca, 'YDir','reverse')
set(gca().YAxis,'TickLength',[0 0])
set(gca().XAxis,'TickLength',[0 0])
ylabel('Stim channel')
xlabel('Time (s)')
title(sprintf('Example stimulation raster plot for lambda = %d',stimRate))
end

function plotISI(arr)
    ISIarr = [];
    for i = 1:size(arr,1)
        chISI = diff(arr{i});
        ISIarr = [ISIarr chISI];
    end
    histogram(ISIarr,'BinWidth',0.01);
    title('ISI distribution of generated stimulation timestamps');
    xlabel('Time (s)')
    ylabel('Frequency')
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