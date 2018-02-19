load('developmentB3Nwt_all111015.mat')

PostnatalAge = [17 17 16 18 18 22 16 17 17 17 22 22 60 60 16 16 18 18 19 19 60 14 14 14 19 19 20 20 18 23 25 25 17 17 18 18 19 19 60 60 60 60 60 60 60 60 25 22 24 24 27 28  29 25 26 26 27 29 60 ];

for i = 1:length(wx)
    genotype(i) = geno(site(i));  %%% 0 = ko, 1 = het, 2 = wt
    age(i) = PostnatalAge(site(i));
end

figure
hist(geno)


age(age>28)=28;
figure
hist(age)
xlabel('age')

figure
hist(age(genotype==2),14:30);
xlabel('age')

% ageBins = [14 16; 17 21; 22 25; 26 29; 29.5 30.5]';
% ageBins = [14 16; 17 21; 22 25; 26 30]';


ageBins = [14 16; 17 21; 22 24; 27.5 28.5]';
%ageBins = [14 16; 17 18; 19 21; 22 24; 27.5 28.5]';
%ageBins = [14 16; 17 19; 20 22; 23 25; 27.5 28.5]';

load('c:\wn016alpha1_10hzLg60Hz.mat');

spectrum_thresh=0.2; %%%to avoid bringing up high freqs
%spectrum_thresh=0.05 %%% to correct all freqs, but not rounding error
n_frames = size(moviedata,3);
tic
spectrum = zeros(size(moviedata(1:76,:,:)),'single');
for f= 1:n_frames
    spectrum(:,:,f) = fft2(single(moviedata(1:76,:,f))-128);
end
mean_spectrum = mean(abs(spectrum(:,:,:)),3);
mean_spectrum = mean_spectrum / max(max(mean_spectrum));




% 
% %%% example rfs
% for a = [1 4]
%     a
%     figure
%    set(gcf,'Name',sprintf('age %d',ageBins(1,a)));
%    use = find(genotype==2 & age>=ageBins(1,a) & age<=ageBins(2,a));
%    length(use)
%     if length(use)>100
%         use = use(1:100);
%     end
%     for n= 1:length(use);
%         subplot(10,10,n);
%         if ~isempty(wn_all(use(n)).svd_xy)
%             %%%  imagesc(wn_all(use(n)).sta_final,[-0.1 0.1]); axis equal; colormap jet
%             sta = squeeze(wn_all(use(n)).svd_xy(1,:,:));
%             sta_fix = ifft2(fft2(sta).*mean_spectrum);
%             
%             imagesc(sta_fix,[-0.025 0.025]); axis equal; colormap gray
%         else          
%             imagesc(0,[-0.1 0.1]), colormap gray;
%         end
%         axis off
%         set(gca,'LooseInset',get(gca,'TightInset'))
%         title(sprintf('%d',use(n)))
%     end
%     drawnow
% end
% 
% 
% 

%%% extract flash values
for i = 1:length(fl_all);
 
    err = fl_all(i).rf ./ fl_all(i).RFzcore;
    err = nanmedian(err(:));
    fl_all(i).zscore = fl_all(i).rf ./err;
    fl_all(i).zabs = max(abs(fl_all(i).zscore(:,:,1)'), abs(fl_all(i).zscore(:,:,2)'));
    
    usepts = find(fl_all(i).zabs'>5); %%% was 7.5
    npts(i) = length(usepts);
    fl_all(i).onoffoverlap = nanmean(fl_all(i).onoffbias(usepts));
    fl_all(i).sustVariation = nanstd(fl_all(i).sustainBias(usepts));
    onoffoverlap(i) = fl_all(i).onoffoverlap;
    sustVariation(i) =  fl_all(i).sustVariation ;
    meanSust(i) = nanmean(fl_all(i).sustainBias(usepts));
end

sum(npts==0)

cmap = cbrewer('div','RdBu',64);
cmap = flipud(cmap);

%%% show examples - can select STA, onoff, or sustain

adult = [184 187];
mid = [14 60 65];
early = [39 214];
all = [39 214 14 60 65 184 187];  %%% sta examples
figure
for i= 1:length(all)
    subplot(3,3,i);
    imagesc(squeeze(wn_all(all(i)).svd_xy(1,:,:)),[-0.1 0.1]); axis equal; colormap jet
    title(sprintf('%f',wx(all(i))))
    axis off
    set(gca,'LooseInset',get(gca,'TightInset'))
end

adult = [184 187];
mid = [14 60 65];
early = [256 361 39 214];
%all = [256 361 184]

all = [early mid adult];
figure
for i= 1:length(all)
    subplot(3,3,i);
           sta = squeeze(wn_all(all(i)).svd_xy(1,:,:));
            sta_fix = ifft2(fft2(sta).*mean_spectrum);
            
            imagesc(sta_fix,[-0.025 0.025]); axis equal; colormap gray
    title(sprintf('%d %f',all(i),wx(all(i))))
    axis off
    set(gca,'LooseInset',get(gca,'TightInset'))
end



adult = [184 187];
mid = [14 60 65];
early = [39 214];
all = [39 214 14 60 65 184 187];  %%% sta examples
all = [101 212 159 185]
figure
for i= 1:length(all)
    subplot(3,3,i);
    onoffIm = mat2im(fl_all(all(i)).onoffbias',cmap,[-1.5 1.5]);
    amp = fl_all(all(i)).zabs/15;  amp(amp>1)=1; amp = repmat(amp,[1 1 3]);
    imshow(onoffIm.*amp);  title(sprintf('%0.2f',fl_all(all(i)).onoffoverlap));
    axis off
    set(gca,'LooseInset',get(gca,'TightInset'))
end

adult = [184 187];
mid = [14 60 65];
early = [39 214];
all = [39 214 14 60 65 184 187];  %%% sta examples
all = [101 212 159 185]
figure
for i= 1:length(all)
    subplot(3,3,i);
    amp = fl_all(all(i)).zabs/15;  amp(amp>1)=1; amp = repmat(amp,[1 1 3]);
    sustIm = mat2im(fl_all(all(i)).sustainBias',cmap,[0.2 0.8]);
    imshow(sustIm.*amp); title(sprintf('%0.2f',fl_all(all(i)).sustVariation));
    axis off
    set(gca,'LooseInset',get(gca,'TightInset'))
end

% %%% show all on/off vs age
% for a = 1:4
%     a
%     figure
%        set(gcf,'Name',sprintf('age %d sust',ageBins(1,a)));
%     use = find(genotype==2 & age>=ageBins(1,a) & age<=ageBins(2,a));
%     if length(use)>100
%         use = use(1:100);
%     end
%     for n= 1:length(use);
%         subplot(10,10,n);
%         
%         rf = squeeze(fl_all(use(n)).rf(:,:,2) - fl_all(use(n)).rf(:,:,1))';
%         m = max(max(abs(rf)));
%         onoffIm = mat2im(fl_all(use(n)).onoffbias',cmap,[-1.5 1.5]);
%         amp = fl_all(use(n)).zabs/15;  amp(amp>1)=1; amp = repmat(amp,[1 1 3]);
%         imshow(onoffIm.*amp);  title(sprintf('%0.2f %d',fl_all(use(n)).onoffoverlap,use(n)));
%         sustIm = mat2im(fl_all(use(n)).sustainBias',cmap,[0 0.75]);
%         imshow(sustIm.*amp); title(sprintf('%0.2f %d',fl_all(use(n)).sustVariation,use(n)));
%         % imagesc(rf,[-10 10]);
%         % imagesc(squeeze(wn_all(use(n),1).svd_xy(1,:,:)),[-0.1 0.1]);
%         % imagesc(abs(fl_all(use(n)).zscore(:,:,1)') + abs(fl_all(use(n)).zscore(:,:,2)')   ,[2 20]);
%         % imagesc(fl_all(use(n)).zabs,[2 20])
%         
%         %imagesc(fl_all(use(n)).onoffbias',[-1 1]);
%         axis equal; colormap jet;  axis off
%         set(gca,'LooseInset',get(gca,'TightInset'))
%   
%     end
%     drawnow
% end

plotDevGenoData(abs(onoffoverlap'),age,agelist, ageBins,~isnan(onoffoverlap), genotype,'onoff overlap')
scatterDevGenoData(abs(onoffoverlap'),age,agelist, ageBins,~isnan(onoffoverlap), genotype,'onoff overlap')

plotDevGenoData(abs(sustVariation'),age,agelist, ageBins,~isnan(sustVariation), genotype,'sustain variation')
scatterDevGenoData(abs(sustVariation'),age,agelist, ageBins,~isnan(sustVariation), genotype,'sustain variation')

plotDevGenoData(meanSust',age,agelist, ageBins,~isnan(sustVariation), genotype,'mean sustain')
scatterDevGenoData(meanSust',age,agelist, ageBins,~isnan(sustVariation), genotype,'mean sustain')



%%% bursting



plotDevGenoData(burst_fraction(:,1),age,agelist, ageBins,~isnan(burst_fraction(:,1)'), genotype,'burst stationary')
ylim([0 0.25])
scatterDevGenoData(burst_fraction(:,1),age,agelist, ageBins,~isnan(burst_fraction(:,1)'), genotype,'burst stationary')
ylim([0 0.75])

plotDevGenoData(burst_fraction(:,2),age,agelist, ageBins,~isnan(burst_fraction(:,2)'), genotype,'burst moving')
ylim([0 0.25])
scatterDevGenoData(burst_fraction(:,2),age,agelist, ageBins,~isnan(burst_fraction(:,2)'), genotype,'burst moving')
ylim([0 0.75])


%%% size suppression

plotDevGenoData(1-fl_supp',age,agelist, ageBins,fl_amp>fl_thresh, genotype,'size suppression')
ylim([0 0.75])
scatterDevGenoData(1-fl_supp',age,agelist, ageBins,fl_amp>fl_thresh, genotype,'size suppression')

plotDevwtData(1-fl_supp',age,agelist, ageBins,fl_amp>fl_thresh, genotype,'size suppression')
ylim([0 0.75])

%%% preferred spot size

plotDevGenoData(fl_sz',age,agelist, ageBins,fl_amp>fl_thresh, genotype,'pref spot size')
set(gca,'Ytick',1:6); ylim([ 3 6])
set(gca,'Yticklabel',{'2','4','8','16','32','full'});

plotDevwtData(fl_sz',age,agelist, ageBins,fl_amp>fl_thresh, genotype,'pref spot radius')
set(gca,'Ytick',1:6); ylim([ 3 6])
set(gca,'Yticklabel',{'1','2','4','8','16','full'});

scatterDevGenoData(fl_sz'+rand(size(fl_sz')),age,agelist, ageBins,fl_amp>fl_thresh, genotype,'pref spot size')
set(gca,'Ytick',1:6); ylim([ 0 6.5])
set(gca,'Yticklabel',{'2','4','8','16','32','full'});



%%% RF width


plotDevGenoData(wx,age,agelist, ageBins,1, genotype,'wx from wn (sp/sec)')

scatterDevGenoData(wx,age,agelist, ageBins,1, genotype,'wx from wn (sp/sec)')

scatterDevwtData(wx,age,agelist, ageBins,genotype==2 , genotype,'RF radius')
xlim([13.5 25.5])

plotDevwtData(wx,age,agelist, ageBins,genotype==2, genotype,'RF radius')
ylim([0 16])

plotDevwtData(~isnan(wx),age,agelist, ageBins,genotype==2 , genotype,'fraction center surround')
ylim([0 0.75])

%%% RF area
plotDevGenoData(wx.*wy,age,agelist, ageBins,1, genotype,'area from wn (sp/sec)')
scatterDevGenoData(wx.*wy,age,agelist, ageBins,genotype==2, genotype,'area from wn (sp/sec)')

%%% grating response amp
plotDevwtData(sf_amp',age,agelist, ageBins,1, genotype,'grating sp/sec')
ylim([0 5])
scatterDevGenoData(sf_amp',age,agelist, ageBins,1, genotype,'grating sp/sec')

plotDevwtData(sf_amp',age,agelist, ageBins,1, genotype,'grating sp/sec')
ylim([0 10])
plotDevwtData(fl_amp',age,agelist, ageBins,1, genotype,'grating sp/sec',1)
legend('grating','spots')

%%% grrating responsive fraction

plotDevwtData(fl_amp'>2,age,agelist, ageBins,1, genotype,'fraction responsive',0)
ylim([0 1])
plotDevwtData(sf_amp'>2,age,agelist, ageBins,1, genotype,'fraction responsive',1)
ylim([0 1])
legend('spots','gratings')



%%% peak SF
%%% need to fix
for i = 1:length(drift_all)
    d= drift_all(i);
    sftune = squeeze(nanmean(d.sf_tune(1:2,1,:),1));  %%% sftune(1:2 tfs, 1:3 f0 f1 f2, sf); (for some reason sftune(3,:,:)=0 not sure why this is not just 1:2
   sftune(1) = sftune(1)-dr_spont(i); % if using f0, need to subtract
   %%% spont
    [sf_amp(i) peak_sf(i)] = max(sftune);
       % [sf_amp(i) peak_sf(i)] = max(squeeze(mean(nanmean(d.sf_tune(:,1:2,:),2),1)));    
       df1f0(i) = max(d.orient_tune(peak_tf(i),2,:)) / max(d.orient_tune(peak_tf(i),1,:));
end


figure
hist(df1f0(drift_amp>2 & df1f0>0 & df1f0<2))

figure
hist(peak_sf(genotype==2),1:7)

sfs = [0 0.01 0.02 0.04 0.08 0.16 0.32 nan];
sf_inds = peak_sf;
sf_inds(peak_sf ==0) = 7;
sf_inds(isnan(peak_sf))=8;

plotDevGenoData(sfs(sf_inds)',age,agelist, ageBins,sf_amp>2, genotype,'peak sf (cpd)')
scatterDevGenoData(sfs(sf_inds)',age,agelist, ageBins,sf_amp>2, genotype,'peak sf (cpd)')


plotDevGenoData(driftF1F0',age,agelist, ageBins,sf_amp>2, genotype,'F1F0')
scatterDevGenoData(driftF1F0',age,agelist, ageBins,sf_amp>2, genotype,'F1F0')



plotDevGenoData(sf_amp',age,agelist, ageBins,~isnan(sf_amp), genotype,'sf amp')
scatterDevGenoData(sf_amp',age,agelist, ageBins,~isnan(sf_amp), genotype,'sf amp')


plotDevwtData(sfs(sf_inds)',age,agelist, ageBins,sf_amp>2 & genotype==2, genotype,'peak sf (cpd)')

%%% fraction lowpass (greatest response to 0cpd)

plotDevGenoData((peak_sf==1)',age,agelist, ageBins,sf_amp>1, genotype,'fraction lowpass')
ylim([0 1])


%%% fraction OS

plotDevwtData(driftOSI(:,1)>0.25,age,agelist, ageBins,sf_amp>2, genotype,'fraction')
ylim([ 0 1])
%%% fraction DS
plotDevwtData(driftDSI(:,1)>0.25,age,agelist, ageBins,sf_amp>2, genotype,'fraction',1)
ylim([0 1])
legend('OS','DS')


plotDevwtData(driftOSI(:,1),age-0.1,agelist-0.1, ageBins-0.1,sf_amp>2, genotype,'mean OSI')
ylim([ 0 0.5])
%%% fraction DS
plotDevwtData(driftDSI(:,1),age+0.1,agelist+0.1, ageBins+0.1,sf_amp>2, genotype,'mean selectivity',1)
ylim([0 0.5])
legend('OSI','DSI')


scatterDevGenoData(driftOSI(:,1),age,agelist, ageBins,sf_amp>2, genotype,'fraction orientation selective >0.2')
ylim([ 0 1])

%%% fraction DS
plotDevGenoData(driftDSI(:,1)>0.25,age,agelist, ageBins,sf_amp>2, genotype,'fraction direction selective >0.2')
ylim([0 1])

%%% fraction DS
plotDevGenoData(driftDSI(:,1),age,agelist, ageBins,sf_amp>2, genotype,'mean DSI')
ylim([0 0.5])

%%% fraction SBC
sbc= (wn_cr_dom<-1*10^3);
plotDevGenoData(sbc',age,agelist, ageBins,1, genotype,'fraction sbc')
ylim([0 0.5])




%%% figure 1
figure
all = [256 361 184];
all = [256 361 39  186 184  530]
for i= 1:length(all)
    subplot(1,6,i);
           sta = squeeze(wn_all(all(i)).svd_xy(1,:,:));
            sta_fix = ifft2(fft2(sta).*mean_spectrum);
            if i==5 | i==4, range = [-0.05 0.05], else range = [-0.025 0.025], end
            if i==6, range = [-0.04 0.04]; end
            imagesc(sta_fix,range); axis equal; colormap gray
            %if i==3, title('adult'), else title('P16'), end
            hold on
           if i ==1
               plot([ 68 81],[68 68],'w','LineWidth',4) %%% 10 deg scale bar
           end
    axis off
    axis([17 17+75 1 76])
    if i==6, axis([1 76 1 76]); end
    set(gca,'LooseInset',get(gca,'TightInset'))
end


figure

all = [256 361 184];
for i= 1:length(all)
    subplot(2,3,i);
           sta = squeeze(wn_all(all(i)).svd_xy(1,:,:));
            sta_fix = ifft2(fft2(sta).*mean_spectrum);
            if i==3, range = [-0.05 0.05], elseif i ==2, range = [-0.035 0.035],else range = [-0.025 0.025], end
            imagesc(sta_fix,range); axis equal; colormap gray
            if i==3, title('adult'), else title('P16'), end
            hold on
            plot([ 68 81],[68 68],'w','LineWidth',4) %%% 10 deg scale bar
    axis off
    axis([17 17+75 1 76])
    set(gca,'LooseInset',get(gca,'TightInset'))
end


subplot(2,3,4)
[resp err] = plotDevGenoData(~isnan(wx),age,agelist, ageBins,genotype==2 , genotype,'fraction STA',1)
ylim([0 0.75])
set(gca,'Xtick',[14 18 22 28]);
set(gca,'Xticklabel',{'14','18','22','adult'})


subplot(2,3,5);
scatterDevGenoData(wx,age,agelist, ageBins,genotype==2, genotype,'RF radius (deg)',1)
%plotDevGenoData(wy,age,agelist, ageBins,genotype==2, genotype,'RF radius',1)

subplot(2,3,6);
[resp err] = plotDevGenoData(wx,age,agelist, ageBins,genotype==2, genotype,'RF radius (deg)',1)
ylim([0 18])


[data group] = labelbyage(double(~isnan(wx)),age,ageBins,genotype==2);
[p tbl stats] = kruskalwallis(data(group>0),group(group>0))
multcompare(stats,'Ctype','dunn-sidak')


[data group] = labelbyage(wx,age,ageBins,genotype==2 & ~isnan(wx'));
[p tbl stats] = kruskalwallis(data,group)
multcompare(stats,'Ctype','dunn-sidak')

%%% figure 2
figure
subplot(2,3,1);
use = age>=ageBins(1,1) & age<=ageBins(2,1) & genotype==2;
errorbar(mean(fl_sztune(use & fl_amp>fl_thresh,:)),std(fl_sztune(use & fl_amp>fl_thresh,:))/sqrt(sum(use & fl_amp>fl_thresh)),'k');
xlim([0.5 6.5]); ylim([0 5.5])
ylabel('sp/sec'); xlabel('spot radius (deg)')
title('p14-16');
set(gca,'Xtick',1:6);
set(gca,'Xticklabel',{'1','2','4','8','16','full'});
box off
% figure
% kruskalwallis(fl_sztune(use & fl_amp>fl_thresh,:))

  
subplot(2,3,2);
use = age>=ageBins(1,end) & age<=ageBins(2,end) & genotype==2;
errorbar(mean(fl_sztune(use & fl_amp>fl_thresh,:)),std(fl_sztune(use & fl_amp>fl_thresh,:))/sqrt(sum(use & fl_amp>fl_thresh)),'k');
axis([0.5 6.5 0 5.5])
ylabel('sp/sec'); xlabel('spot radius (deg)')
title('adult')
set(gca,'Xtick',1:6);
set(gca,'Xticklabel',{'1','2','4','8','16','full'});
box off
% figure
% kruskalwallis(fl_sztune(use & fl_amp>fl_thresh,:))

% figure
% subplot(2,3,2);
% use = age>=ageBins(1,1) & age<=ageBins(2,1) & genotype==2;
% errorbar((1:6)-0.1,mean(fl_sztune(use & fl_amp>fl_thresh,:)),std(fl_sztune(use & fl_amp>fl_thresh,:))/sqrt(sum(use & fl_amp>fl_thresh)),'k');
% hold on
% use = age>=ageBins(1,end) & age<=ageBins(2,end) & genotype==2;
% errorbar((1:6)+0.1,mean(fl_sztune(use & fl_amp>fl_thresh,:)),std(fl_sztune(use & fl_amp>fl_thresh,:))/sqrt(sum(use & fl_amp>fl_thresh)),'k');
% xlim([0.5 6.5]); ylim([0 5.5])
% ylabel('sp/sec'); xlabel('spot radius (deg)')
% legend('P14-16','adult')
% set(gca,'Xtick',1:6);
% set(gca,'Xticklabel',{'1','2','4','8','16','full'});
% box off
% figure
% kruskalwallis(fl_sztune(use & fl_amp>fl_thresh,:))


% subplot(2,3,3)
% use = age>=ageBins(1,2) & age<=ageBins(2,2) & genotype==2;
% m = max(mean(fl_sztune(use & fl_amp>fl_thresh,:)));
% errorbar(mean(fl_sztune(use & fl_amp>fl_thresh,:))/m,std(fl_sztune(use & fl_amp>fl_thresh,:))/(m*sqrt(sum(use & fl_amp>fl_thresh))));
% hold on
% use = age>=ageBins(1,1) & age<=ageBins(2,1) & genotype==2;
% m = max(mean(fl_sztune(use & fl_amp>fl_thresh,:)));
% errorbar(mean(fl_sztune(use & fl_amp>fl_thresh,:))/m,std(fl_sztune(use & fl_amp>fl_thresh,:))/(m*sqrt(sum(use & fl_amp>fl_thresh))));
% xlim([0.5 6.5]); ylim([0 1.2])
% ylabel('sp/sec'); xlabel('spot radius')
% title('p14-16');
% set(gca,'Xtick',1:6);
% set(gca,'Xticklabel',{'1','2','4','8','16','full'});
% box off

subplot(2,3,3)
[resp err] = plotDevwtData(fl_sz',age,agelist, ageBins,fl_amp>fl_thresh, genotype,'pref spot radius',1)
set(gca,'Ytick',1:6); ylim([ 3 6])
set(gca,'Yticklabel',{'1','2','4','8','16','full'});

subplot(2,3,4)
[resp err] = plotDevwtData(1-fl_supp',age,agelist, ageBins,fl_amp>fl_thresh, genotype,'size suppression',1)
ylim([0 0.75])


subplot(2,3,5)
[resp err] = plotDevwtData(sfs(sf_inds)',age,agelist, ageBins,sf_amp>2 & genotype==2, genotype,'peak sf (cpd)',1)


% %%%DS OS mean
subplot(2,3,6)
[resp err] = plotDevwtData(driftOSI(:,1),age-0.1,agelist-0.1, ageBins-0.1,sf_amp>2, genotype,'mean OSI',1)
[resp err] =plotDevwtData(driftDSI(:,1),age+0.1,agelist+0.1, ageBins+0.1,sf_amp>2, genotype,'mean selectivity',1)
ylim([0 0.5])
legend('OSI','DSI')

figure
subplot(2,3,6)
plotDevwtData(fl_amp'>2,age,agelist, ageBins,1, genotype,'fraction responsive',1)
ylim([0 1])
plotDevwtData(sf_amp'>2,age,agelist, ageBins,1, genotype,'fraction responsive',1)
ylim([0 1])
legend('spots','gratings')
axis square

%%% DS/OS fraction
% figure
% subplot(2,3,6)
% plotDevwtData(driftOSI(:,1)>0.25,age,agelist, ageBins,sf_amp>2, genotype,'fraction',1)
% plotDevwtData(driftDSI(:,1)>0.25,age,agelist, ageBins,sf_amp>2, genotype,'fraction',1)
% ylim([0 1])
% legend('OS','DS')
% axis square



[data group] = labelbyage(fl_sz',age,ageBins,genotype==2 & fl_amp>fl_thresh);
[p tbl stats] = kruskalwallis(data(group>0),group(group>0))
multcompare(stats,'Ctype','dunn-sidak')


[data group] = labelbyage(1-fl_supp',age,ageBins,genotype==2 & fl_amp>fl_thresh);
[p tbl stats] = kruskalwallis(data(group>0),group(group>0))
multcompare(stats,'Ctype','dunn-sidak')

[data group] = labelbyage(sfs(sf_inds)',age,ageBins,genotype==2 & sf_amp>2);
[p tbl stats] = kruskalwallis(data(group>0),group(group>0))
multcompare(stats,'Ctype','dunn-sidak')


[data group] = labelbyage(driftOSI(:,1),age,ageBins,genotype==2 & sf_amp>2);
[p tbl stats] = kruskalwallis(data(group>0),group(group>0))
multcompare(stats,'Ctype','dunn-sidak')

[data group] = labelbyage(driftDSI(:,1),age,ageBins,genotype==2 & sf_amp>2);
[p tbl stats] = kruskalwallis(data(group>0),group(group>0))
multcompare(stats,'Ctype','dunn-sidak')

[data group] = labelbyage(fl_amp>2,age,ageBins,genotype==2 );
[p tbl stats] = kruskalwallis(data(group>0),group(group>0))
multcompare(stats,'Ctype','dunn-sidak')

[data group] = labelbyage(sf_amp>2,age,ageBins,genotype==2 );
[p tbl stats] = kruskalwallis(data(group>0),group(group>0))
multcompare(stats,'Ctype','dunn-sidak')


%%% figure 3
figure
all = [212 159 185]
all = [ 212 101  338 185 199]
for i= 1:4
    subplot(2,4,i);
    onoffIm = mat2im(fl_all(all(i)).onoffbias',cmap,[-1.35 1.35]);
      amp = fl_all(all(i)).zabs/8;  if i ==5 amp = amp-0.5; else amp = amp-0.3;end;
      amp(amp>1)=1; amp(amp<0)=0; amp = repmat(amp,[1 1 3]); 

    imshow(imresize(onoffIm.*amp,2));  %title(sprintf('%0.2f',fl_all(all(i)).onoffoverlap));
    if i ==1
        axis([21 143 1 123]); title('P16')
           elseif i ==5
          axis([123 244 1 123]);
    elseif i ==2
        axis([61 183 1 123])
    else       
        axis([ 81 203 1 123]); title('adult')
    end
      if i==2, title('P22'), end
    axis off
    set(gca,'LooseInset',get(gca,'TightInset'))
end


all = [ 212 159 185]
all = [ 212 101 338 185  199]
for i= 1:4
    subplot(2,4,4+i);
    amp = fl_all(all(i)).zabs/8;  if i ==5 amp = amp-0.5; else amp = amp-0.3;end;
    amp(amp>1)=1; amp(amp<0)=0; amp = repmat(amp,[1 1 3]); 
    sustIm = mat2im(fl_all(all(i)).sustainBias',cmap,[0 .75]);
    imshow(imresize(sustIm.*amp,2)); % title(sprintf('%0.2f',fl_all(all(i)).sustVariation));
      if i ==1
        axis([21 143 1 123]); title('P16')
      elseif i ==5
          axis([123 244 1 123]);
              elseif i ==2
        axis([61 183 1 123])
    else       
        axis([ 81 203 1 123]); title('adult')
      end
    if i==2, title('P22'), end
    axis off  
    set(gca,'LooseInset',get(gca,'TightInset'))
end

figure
subplot(1,3,1)
[resp err]=plotDevGenoData(~isnan(onoffoverlap'),age,agelist, ageBins, genotype==2, genotype,'fraction mapped',1)
ylim([ 0 1])

subplot(1,3,2)
[resp err]=plotDevGenoData(abs(onoffoverlap'),age,agelist, ageBins,~isnan(onoffoverlap) & genotype==2, genotype,'on/off uniformity',1)
ylim([ 0 1])

subplot(1,3,3)
[resp err]= plotDevGenoData(1-abs(sustVariation'),age,agelist, ageBins,~isnan(sustVariation)& genotype==2, genotype,'sustain uniformity',1)
ylim([ 0 1])

[data group] = labelbyage(abs(onoffoverlap'),age,ageBins,genotype==2 &~isnan(onoffoverlap));
[p tbl, stats] = kruskalwallis(data(group>0),group(group>0))
multcompare(stats,'Ctype','dunn-sidak')

[data group] = labelbyage(1-abs(sustVariation'),age,ageBins,genotype==2 & ~isnan(sustVariation));
[p tbl stats] = kruskalwallis(data(group>0),group(group>0))
multcompare(stats,'Ctype','dunn-sidak')


figure
colormap(cmap);
colorbar
