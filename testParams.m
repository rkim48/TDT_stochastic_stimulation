%% Test latencies 

load('columnNames.mat','columnNames');
% Create params array to store parameters for single trial 
nROWS = 9;
current = 1;
period = 0.4;
pulseDurations = 80:40:400;
tsColumns = [0 200 400 600]; % delays
chColumns = [1 2 3 4];

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
paramsArr(:,durColIdx) = repmat(pulseDurations',1,4);
paramsArr(:,delayColIdx) = repmat(tsColumns,nROWS,1);
paramsArr(:,chColIdx) = repmat(chColumns,nROWS,1);

%% Write par file 

paramTable = array2table(paramsArr,'RowNames',string(1:nROWS));
paramTable.Properties.VariableNames = string(columnNames);
writetable(paramTable,'test.par.csv','WriteRowNames',true);

%% Write seq file
seq = 1:9;
nSEQ = numel(seq);
time = 801.6 * ones(nSEQ,1);
seqArr = [seq' time];
columnNames = ["Seq-1" "Time-1"];
paramTable = array2table(seqArr,'RowNames',string(1:nSEQ));
paramTable.Properties.VariableNames = string(columnNames);
writetable(paramTable,'test.seq.csv','WriteRowNames',true);

