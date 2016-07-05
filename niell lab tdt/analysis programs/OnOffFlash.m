
clear onset_hist err_est baseline_hist baseline_resp flash_resp baseline_est onset err base
nx =128; ny= 64;

%%% find spikes
for i=1:length(eps)
    epSpikes{i} = times(times>=eps(2,i) & times<eps(2,i)+0.8)-eps(2,i);
    nEpSpikes(i)=sum(epSpikes{i}>=0.3 & epSpikes{i}<=0.6);
end

%%% get timecourse at each location
tic
parfor xi=1:nx;
   xi
x =xi*2;
   [onset{xi} resp{xi} err{xi} base{xi} baseresp{xi} onset_bins] = doOnOffFlash(squeeze(moviedata(x,:,:)),squeeze(sz_mov(x,:,:)),eps,epSpikes,ny);
   
   
end
toc

%%% recompile parallel variables
for xi=1:nx
    onset_hist(xi,:,:,:,:) = onset{xi};
    err_est(xi,:,:,:) = err{xi};
    baseline_hist(xi,:,:,:,:) = base{xi};
    flash_resp(xi,:,:,:) = resp{xi};
    baseline_resp(xi,:)= baseresp{xi};
    
end
clear onset err base


%%% show all traces

stimbins = onset_bins>=0.2 & onset_bins<=0.7;

%%% get min/max
clear resps
for rep = 1:2
    resps(:,:,rep,:) = squeeze(mean(onset_hist(:,:,rep,2:5,stimbins),4))-baseline_hist(:,:,stimbins);
end


sfactor=8;
rf_fig = figure;
for x = 1:nx/sfactor
    for y = 1:ny/sfactor
        xi= x*sfactor; yi = y*sfactor;
        subplot(nx/sfactor,ny/sfactor,(x-1)*ny/sfactor + y)
        plot(squeeze(mean(onset_hist(xi,yi,1,2:5,stimbins),4))-squeeze(baseline_hist(xi,yi,stimbins)));
        hold on
        plot(squeeze(mean(onset_hist(xi,yi,2,2:5,stimbins),4))-squeeze(baseline_hist(xi,yi,stimbins)),'r');
        ylim([min(resps(:)) max(resps(:))]); xlim([1 sum(stimbins)])
        axis off
        set(gca,'LooseInset',get(gca,'TightInset'));
    end
end

pix_trace = squeeze(mean(onset_hist(:,:,:,2:5,stimbins),4));
pix_baseline = baseline_hist(:,:,stimbins);

%%% get on/off responses, with baseline subtracted off
stimbins = onset_bins>=0.3 & onset_bins<=0.55;
clear resps
for rep = 1:2
    resps(:,:,rep,:) = squeeze(mean(onset_hist(:,:,rep,2:5,stimbins),4))-baseline_hist(:,:,stimbins);
end


datafig=figure;

%flash_resp = squeeze(mean(onset_hist(:,:,:,:,histbins>0.25 & histbins<=0.5),5) - mean(onset_hist(:,:,:,:,histbins>=0.15 & histbins<=0.25),5));
% for rep = 1:2
%     for sz = 1:6
%         flash_resp(:,:,rep,sz) = squeeze(mean(onset_hist(:,:,rep,sz,histbins>=0.35 & histbins<=0.6),5) - mean(baseline_hist(:,:,histbins>=0.35 & histbins<=0.6),3));
%     end
% end

%%% amplitude of response
for rep = 1:2
    for sz = 1:6
        flash_resp(:,:,rep,sz) = flash_resp(:,:,rep,sz) - baseline_resp(:,:);
    end
end


lims = max(abs([min(flash_resp(:)) max(flash_resp(:)) 10]));
lims = [-lims lims];

for rep = 1:2
    for sz = 1:6
        figure(datafig)
        subplot(2,6,(rep-1)*6+sz);
        imagesc(squeeze(flash_resp(:,:,rep,sz)),lims)
        axis off; axis equal; set(gca,'LooseInset',get(gca,'TightInset'));
        
%         figure(errfig)
%                 subplot(2,6,(rep-1)*6+sz);
%         imagesc(squeeze(flash_resp(:,:,rep,sz)./err_est(:,:,rep,sz)),[-2 2])
%         axis off; axis equal
    end
end

%%% on and off RFs, averaged over size 2:5
figure
for rep = 1:2
    subplot(2,2,rep)
    rf(:,:,rep) =medfilt2( mean(flash_resp(:,:,rep,2:5),4));
    imagesc(rf(:,:,rep)',lims/2);
    if rep==1; title('off'); else title('on'); end
    axis off; axis equal; set(gca,'LooseInset',get(gca,'TightInset'));
   
end
 err= mean(err_est(:,:,:,1:4),4)/sqrt(4);
subplot(2,2,3);
diffRF = rf(:,:,2)'- rf(:,:,1)';
imagesc(diffRF,lims/2)
title('diff')
axis off ;axis equal; set(gca,'LooseInset',get(gca,'TightInset'));
subplot(2,2,4)
errormap = mean(mean(err_est,4),3)/sqrt(4);
zscore = diffRF./errormap'; diffRF(abs(zscore)<6)=0;
imagesc(diffRF,lims/2); title('zscore >10')
axis off; axis equal; set(gca,'LooseInset',get(gca,'TightInset'));

cmap = cbrewer('div','RdBu',64);
cmap = flipud(cmap);


%%% get on/off bias over space
figure
rfz = rf./err;
rfzAbs = abs(rfz(:,:,1))+ abs(rfz(:,:,2));
amp = rfzAbs/10; amp(amp>1)=1; amp(amp<0.5)=0.;
amp = repmat(amp,[1 1 3]);
rfpos = rf; rfpos(rf<0.001)=0.001;
onoffbias = (rfpos(:,:,2)-rfpos(:,:,1))./(rfpos(:,:,1)+rfpos(:,:,2));
onoffIm = mat2im(onoffbias,cmap,[-1 1]);

subplot(1,2,1)
imshow(imresize(onoffIm .*amp,10));
title('on off')

%%% get sustain ratio over space 

sustain = mean(resps(:,:,:,1:end),4)./max(resps(:,:,:,1:end),[],4);
sustain(:,:,1)=medfilt2(sustain(:,:,1)); sustain(:,:,2)=medfilt2(sustain(:,:,2)); 
sustain = (rfpos(:,:,1).*sustain(:,:,1) + rfpos(:,:,2).*sustain(:,:,2) )./(rfpos(:,:,1) + rfpos(:,:,2));
sustainIm = mat2im(sustain,cmap,[0 1]);

subplot(1,2,2)
imshow(imresize(sustainIm.*amp,10));
title('sustain')
