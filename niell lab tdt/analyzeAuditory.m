pname = uigetdir('C:\data\tdt tanks','block data')
delims = strfind(pname,'\');
selected_path = pname(1 :delims(length(delims))-1)
Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
Block_Name = pname(delims(length(delims))+1 :length(pname))

nChan=input('# channels : ');
tic
flags =struct('lfpTseries',0,'lfpSpectra',0,'MUspike',1,'mouseOn',0,'newCluster',1,'analog',1);
tdtData= getTDTdata(Tank_Name, Block_Name, 1:nChan, flags);
toc

%%% sync is first voltage channel
sync = tdtData.analogV{1};

%%% find edges in sync
edge = find(diff(sync)>0.1);
startTimes = tdtData.analogT(edge);

%%% which epoc corresponds to wich stim

duration =0.1;

load('C:\data\tuning-curve-tones090512.mat')

for p = 2:length(stimuli);
    stimlist(p-1) = stimuli(p).param.amplitude/10+1;
end


for ch = length(tdtData.spikeT);  %%% loop through 8 tetrodes
    
    ch
    
    %t = tdtData.MUspikeT{(ch-1)*4+1};
    t = tdtData.spikeT{ch};
    %%% collect spikes for each epoc
    for s = 1:length(edge);
        sp = find(t>startTimes(s) & t<startTimes(s)+duration);
        if ~isempty(sp)
        R(s,ch) = length(sp)/duration;
        spikes{s,ch} = t(sp)-startTimes(s);
        else
            R(s,ch) = 0;
        spikes{s,ch} = []; 
        end
      
    end
    
    %%% collect all epocs that correspond to same stim parameter
    for stim=1:max(stimlist);
        eps = find(stimlist==stim);
        tuning(stim,ch) = mean(R(eps,ch));
    end
    figure
    plot(tuning(:,ch))
end

figure
plot(tuning)

stimlist = [1 2 3 1 2 3 1 2 3 1 2 3 1 2 3];
