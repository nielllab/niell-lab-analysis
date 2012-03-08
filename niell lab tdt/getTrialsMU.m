function [Spike_Timing index trial_no] = getTrialsMU(allEpocs,Spike_TS,cond,duration);
    
            condEpocs=find(allEpocs(1,:)==cond);
            trial_no = length(condEpocs);
            Epocs_TS = allEpocs(2,condEpocs);
            
            
            % Now for all events, we find out when the xTrig was:
            index = zeros(1,N);
            TS_xTrg=index;
            for i = 1:size(Epocs_TS,2)-1;
                epochSpikes = find(Spike_TS>=Epocs_TS(i) & Spike_TS<Epocs_TS(i+1));
                index(epochSpikes)=i;
                TS_xTrg(epochSpikes)=Epocs_TS(i);
            end
            
            %%% subtract off trial start time, and only keep spikes for this
            %%% orientation i.e. index>0
            Spike_Timing=Spike_TS(index>0)-TS_xTrg(index>0);
            index = index(index>0);