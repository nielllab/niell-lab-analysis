function [onset_hist resp err_est onset_bins] = getFlashHist(epSpikes,eps,onset)
use_eps = eps(:,eps(1,:)<length(onset));
useEpSpikes = epSpikes(eps(1,:)<length(onset));

inds = find(onset(use_eps(1,:)));
n= length(inds);
ontime = use_eps(2,inds);

dt=0.05;
histbins = 0:dt:0.8;
stimbins = histbins>=0.2 & histbins<=0.5;
h=0;
clear nsp

allT=[];
for i = 1:n

t = useEpSpikes{inds(i)};
allT = [allT t];
    nsp(i) = sum(t>=0.3 & t<=0.6);

end

    if ~isempty(allT)
        h =histc(allT, histbins);
    end
    
if length(h)>1
    
    %                     hold on
    %                     plot(histbins(1:end-1)+dt/2,h(1:end-1)/(n*dt),c)
    err_est = std(nsp)/(0.3*sqrt(n));
    
    onset_hist = h(1:end-1)/(n*dt);
    
else
    err_est=0;
    onset_hist = histc(0,histbins); onset_hist=onset_hist(1:end-1);
end
   onset_bins=histbins;
   resp = prctileMean(nsp',10)/0.3;
    