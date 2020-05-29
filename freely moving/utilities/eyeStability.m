%%% calculate stability of eyes when animal is totally still (no locomotion or head rotation)

sp_thresh = 1;  %%% speed threshold, cm/sec
dHead_thresh = 1;  %%% rotation threshold, deg/frame

close all
clear n
for i = 1:length(appEpoch)
    vid = useData(i);
    
    %%% get head rotaiton and speed
    acc_dth = sqrt(accelChannels{vid}(1:end-1,4).^2 + accelChannels{vid}(1:end-1,5).^2+ accelChannels{vid}(1:end-1,6).^2);
    sp = mouseSp{vid}(1:end-1);
    
    figure
    plot(acc_dth,sp,'.'); xlabel('rotation'); ylabel('speed')

    
    %%% calculate stdev of eyes during stationary
    stationary = acc_dth<dHead_thresh & sp'<sp_thresh;
    nanstd(dLtheta{vid}(stationary));
    sum(stationary);
    length(stationary);
    
    stability(i,1) = nanstd(dLtheta{vid}(stationary));
    stability(i,2) = nanstd(dRtheta{vid}(stationary));
    n(i) = sum(stationary);
end

%%% how many stationary timepoints needed to include video in average?
nthresh = 30;

sprintf('data from %d videos with >%d frames',sum(n>nthresh),nthresh)
sprintf('left = %0.2f +/- %0.2f',mean(stability(n>nthresh,1),1),std(stability(n>nthresh,1),[],1))
sprintf('right = %0.2f +/- %0.2f',mean(stability(n>nthresh,2),1),std(stability(n>nthresh,2),[],1))

