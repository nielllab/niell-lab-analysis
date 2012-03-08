function [t index Ntrial] = getTrials(allEpocs,times,cond,duration);
        Epocs=find(allEpocs(1,:)==cond);
        Epocs_TS = allEpocs(2,Epocs);  
        % Now for all events, we find out when the xTrig was:
        index = zeros(1,length(times));
        TS_xTrg=index;
        Ntrial = length(Epocs_TS);

        for i = 1:Ntrial;
            if i<Ntrial
                epochSpikes = find(times>=Epocs_TS(i) & times<Epocs_TS(i)+duration);
            else
                epochSpikes = find(times>=Epocs_TS(i) & times<Epocs_TS(i)+duration);
            end
            index(epochSpikes)=i;
            TS_xTrg(epochSpikes)=Epocs_TS(i);
        end

       
        t=times(index>0)-TS_xTrg(index>0);
        index = index(index>0);
        