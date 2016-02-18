function spikes = getSpikes(clustfile,afile, block)
%%% read in single unit spike times for a given block
%%% this is mostly just a matter of sorting the spikes from one block 
%%% but nice to do it just in one line!


load(clustfile,'Block_Name','Tank_Name');
blocknum = find(strcmp(block,Block_Name));
Block_Name = Block_Name{blocknum}

load(afile,'spikes');

if ~exist('spikes','var') | isempty(spikes(blocknum))
    try
        load(afile,'spikeT')
        for c = 1:length(spikeT)
            sp = spikeT{c};
            sp = sp-(blocknum-1)*10^5;
            sp = sp(sp>0 & sp<10^5);
            spikes(blocknum).sp{c} = sp;
        end       
    catch
        spikes(blocknum).sp = NaN;
    end
    save(afiles,'spikes','-append')
end

spikes = spikes(blocknum);
