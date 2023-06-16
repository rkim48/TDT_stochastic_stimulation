addpath(genpath('C:\TDT\TDTMatlabSDK'));
tank = 'D:\Synapse\Tanks\TEST-230506-205111\';
d = dir(tank);
d = d([d(:).isdir]);
d = d(~ismember({d(:).name},{'.','..'}));
session = d(end).name;
data = TDTbin2mat(fullfile(tank,session));
streamNames = data.streams
%% Param out
% Confirm stimulation matches input

eS1p = data.scalars.eS1p;
timestamps = eS1p.ts;
params = eS1p.data';
load('TDT_stochastic_stimulation\columnNamesAll.mat','columnNames');

%% Create table from param output
nROWS = size(params,1);
for i = 1:nROWS
    chColIdx = find(contains(columnNames,'Chan')==1);
    paramRow = params(i,:);
    channels = paramRow(chColIdx);
    stim_idx = find(channels > 0);
    startNanIdx = chColIdx(stim_idx(end))+1;
    paramRow(startNanIdx:end) = nan;
    params(i,:) = paramRow;
end
paramTable = array2table(params,"VariableNames",columnNames,'RowNames',string(1:nROWS));
% writetable(paramTable,'outputParamArray.par.csv','WriteRowNames',true);


stim_ch_out = table2array(paramTable(:,{'ChanA','ChanB','ChanC','ChanD'}));

%%
ts = [1000*diff(timestamps) 1];
columnNames = {'Seq-1','Time-1'};
seqArr = [(1:nROWS)' ts'];
seqTable = array2table(seqArr,'RowNames',string(1:nROWS));
seqTable.Properties.VariableNames = string(columnNames);
% writetable(seqTable,'outputParamArray.seq.csv','WriteRowNames',true);

%%

stim_ts_out = timestamps';
%%
downsampleN = 100;
N = 9000;
subplot(1,2,1);
h1 = plotTrial(stim_ts(1:N),stim_ch(1:N),downsampleN);
hold on;
h2=plotTrial(stim_ts_out(1:N),stim_ch_out(1:N,:),downsampleN);
title('With delay')
subplot(1,2,2);
h1 = plotTrial(stim_ts(1:N),stim_ch(1:N),downsampleN);
hold on;
h2=plotTrial(stim_ts_out(1:N)-stim_ts_out(1)+0.001,stim_ch_out(1:N,:),downsampleN);
legend([h1(1),h2(1)],'Requested timestamps','Actual timestamps')
title('Delay subtracted')
%%
function h = plotTrial(stim_ts,stim_ch,downsampleN)
    stim_ch = stim_ch(1:downsampleN:end,:);
    stim_ts = stim_ts(1:downsampleN:end);
    nSTIM = size(stim_ch,1); 
    for i = 1:nSTIM
        stim_time = stim_ts(i);
        if iscell(stim_ch)
            ch_stim = stim_ch{i};
            nCH = numel(ch_stim);
            h = scatter(repmat(stim_time,1,nCH),ch_stim,'k','Marker','|','LineWidth',1);
        else
            ch_stim = stim_ch(i,:);
            ch_stim = ch_stim(~isnan(ch_stim));
            nCH = numel(ch_stim);
            h = scatter(repmat(stim_time,1,nCH),ch_stim,'r','Marker','|','LineWidth',1);
        end
        hold on;
    end
    ylim([0.5 32.5])
    set(gca, 'YDir','reverse')
    set(gca().YAxis,'TickLength',[0 0])
    set(gca().XAxis,'TickLength',[0 0])
    ylabel('Stim channel')
    xlabel('Time (s)')
    hold on;
end