function [spikeT]= getBlockSpikes;

cells =1;
[fname, pname] = uigetfile('*.mat','cluster data');
load(fullfile(pname,fname));
for i =1:length(Block_Name);
    sprintf('%d : %s ',i,Block_Name{i})
end
b = input('which block to analyze ? ');
Block_Name = Block_Name{b}

[afname, apname] = uigetfile('*.mat','analysis data');
if afname~=0
    afile = fullfile(apname,afname);
    load(afile);
    use_afile=1;
else
    use_afile=0;
    %%% select which units to analyze (channel, cluster number)
end

event_times_all = event_times_all-(b-1)*10^5;


for c = 1:length(cells);

    ch= cells(c,1);
    cl = cells(c,2);
    channel_times =squeeze(event_times_all(ch,:));
    times = channel_times(channel_times>0 & channel_times<10^5);
    idx = squeeze(idx_all(ch,:));
    idx = idx(channel_times>0 & channel_times<10^5);
    spikes = times(idx == cl);
    spikeT{c} = spikes(spikes>0);
end    
