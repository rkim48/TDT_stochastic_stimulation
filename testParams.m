%% Test latencies 

load('columnNames.mat','columnNames');
% Create params array to store parameters for single trial 
nROWS = 2;
current = 1;
period = 400;
pulseDuration = 400;
tsColumns = [0 0.2 0.4 0.6; 0 0.2 0.4 0.6]; % delays
chColumns = [1 2 3 4; 1 2 3 4];

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

%% Write par file 

paramTable = array2table(paramsArr,'RowNames',string(1:nROWS));
paramTable.Properties.VariableNames = string(columnNames);
writetable(paramTable,'test.par.csv','WriteRowNames',true);

%% Write seq file
nSEQ = 4;
seq = [1 2 1 2]';
time = [1000 2000 3000 4000]';
seqArr = [seq time];
columnNames = ["Seq-1" "Time-1"];
paramTable = array2table(seqArr,'RowNames',string(1:nSEQ));
paramTable.Properties.VariableNames = string(columnNames);
writetable(paramTable,'test.seq.csv','WriteRowNames',true);

