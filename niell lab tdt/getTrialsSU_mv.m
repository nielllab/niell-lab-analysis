function [t index Ntrial Epocs_TS] = getTrialsSU_mv(allEpocs,times,cond,duration,tsamp,vsmooth,rep,thresh);
Epocs=find(allEpocs(1,:)==cond);
Epocs_TS = allEpocs(2,Epocs);
% Now for all events, we find out when the xTrig was:
index = zeros(1,length(times));
TS_xTrg=index;
Ntrial = length(Epocs_TS);

for i = 1:Ntrial
    sp(i) = mean(vsmooth(tsamp>=Epocs_TS(i) & tsamp<Epocs_TS(i)+duration));
end

if rep ==1
    useTrials = find(sp<thresh);
    if isempty(useTrials);
        [m useTrials] = min(sp);
    end
elseif rep==2
    useTrials = find(sp>thresh);
    if isempty(useTrials);
        [m useTrials] = max(sp);
    end
end
useTrials
for i = 1:length(useTrials);
    epochSpikes = find(times>=Epocs_TS(useTrials(i)) & times<Epocs_TS(useTrials(i))+duration);
    
    index(epochSpikes)=i;
    TS_xTrg(epochSpikes)=Epocs_TS(useTrials(i));
end

Ntrial = length(useTrials);

size( times(index>0));
size(TS_xTrg(index>0));
t=times(index>0)-TS_xTrg(index>0);
index = index(index>0);
