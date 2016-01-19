PostnatalAge = [17 17 16 18 18 22 16 17 17 17 22 22 60 60 16 16 18 18 19 19 60 14 14 14 19 19 20 20 18 23 25 25 17 17 18 18 19 19 60 60 60 60 60 60 60 60 25 22 24 24 27 28  29 25 26 26 27 29 60 ];

for i = 1:length(wx)
    genotype(i) = geno(site(i));
    age(i) = PostnatalAge(site(i));
end

figure
hist(geno)


age(age>30)=30;
figure
hist(age)

ageBins = [14 16; 17 21; 22 25; 26 29; 29.5 30.5]';
ageBins = [14 16; 17 21; 22 25; 26 30]';


ageBins = [14 16; 17 21; 22 29; 29.5 30.5]';

for a = 1:length(ageBins)
    a
    figure
    use = find(genotype==2 & age>=ageBins(1,a) & age<=ageBins(2,a));
    if length(use)>49 
        use = use(1:49);
    end
    for n= 1:length(use);
   subplot(7,7,n);
   if ~isempty(wn_all(use(n)).svd_xy)
 %%%  imagesc(wn_all(use(n)).sta_final,[-0.1 0.1]); axis equal; colormap jet
    imagesc(squeeze(wn_all(use(n)).svd_xy(1,:,:)),[-0.1 0.1]); axis equal; colormap jet
   else
       
       imagesc(0,[-0.1 0.1]), colormap gray;
   end
   axis off
   set(gca,'LooseInset',get(gca,'TightInset'))
   title(sprintf('%d',use(n)))
    end
    drawnow
end

    



for i = 1:length(fl_all);
 i
err = fl_all(i).rf ./ fl_all(i).RFzcore;  
 err = nanmedian(err(:));
 fl_all(i).zscore = fl_all(i).rf ./err;
 fl_all(i).zabs = max(abs(fl_all(i).zscore(:,:,1)'), abs(fl_all(i).zscore(:,:,2)'));
 
 usepts = find(fl_all(i).zabs'>7.5);
 fl_all(i).onoffoverlap = nanmean(fl_all(i).onoffbias(usepts));
 fl_all(i).sustVariation = nanstd(fl_all(i).sustainBias(usepts));
 onoffoverlap(i) = fl_all(i).onoffoverlap;
 sustVariation(i) =  fl_all(i).sustVariation ;
 meanSust(i) = nanmean(fl_all(i).sustainBias(usepts));
end

cmap = cbrewer('div','RdBu',64);
cmap = flipud(cmap);


adult = [184 187];
mid = [14 60 65];
early = [39 214];
all = [39 214 14 60 65 184 187];  %%% sta examples
all = [101 212 159 185]
figure
for i= 1:length(all)
    subplot(3,3,i);
    %imagesc(squeeze(wn_all(all(i)).svd_xy(1,:,:)),[-0.1 0.1]); axis equal; colormap jet
    title(sprintf('%f',wx(all(i))))
       onoffIm = mat2im(fl_all(all(i)).onoffbias',cmap,[-1.5 1.5]);
        amp = fl_all(all(i)).zabs/15;  amp(amp>1)=1; amp = repmat(amp,[1 1 3]);
        imshow(onoffIm.*amp);  title(sprintf('%0.2f',fl_all(all(i)).onoffoverlap));
         sustIm = mat2im(fl_all(all(i)).sustainBias',cmap,[0.2 0.8]);
        imshow(sustIm.*amp); title(sprintf('%0.2f',fl_all(use(n)).sustVariation));
      axis off
   set(gca,'LooseInset',get(gca,'TightInset'))  
end


for a = 1:4
    a
    figure
    use = find(genotype==2 & age>=ageBins(1,a) & age<=ageBins(2,a));
    if length(use)>36
        use = use(1:36);
    end
    for n= 1:length(use);
        subplot(6,6,n);
     
        rf = squeeze(fl_all(use(n)).rf(:,:,2) - fl_all(use(n)).rf(:,:,1))';        
        m = max(max(abs(rf)));
       onoffIm = mat2im(fl_all(use(n)).onoffbias',cmap,[-1.5 1.5]);
        amp = fl_all(use(n)).zabs/15;  amp(amp>1)=1; amp = repmat(amp,[1 1 3]);
        imshow(onoffIm.*amp);  title(sprintf('%0.2f',fl_all(use(n)).onoffoverlap));
         sustIm = mat2im(fl_all(use(n)).sustainBias',cmap,[0.25 0.75]);
        %imshow(sustIm.*amp); title(sprintf('%0.2f',fl_all(use(n)).sustVariation));
        % imagesc(rf,[-10 10]); 
            % imagesc(squeeze(wn_all(use(n),1).svd_xy(1,:,:)),[-0.1 0.1]); 
          % imagesc(abs(fl_all(use(n)).zscore(:,:,1)') + abs(fl_all(use(n)).zscore(:,:,2)')   ,[2 20]);
        % imagesc(fl_all(use(n)).zabs,[2 20])
         
          %imagesc(fl_all(use(n)).onoffbias',[-1 1]);
        axis equal; colormap jet;  axis off
        set(gca,'LooseInset',get(gca,'TightInset'))
        title(sprintf('%d',use(n)))
    end
    drawnow
end

plotDevGenoData(abs(onoffoverlap'),age,agelist, ageBins,~isnan(onoffoverlap), genotype,'onoff overlap')
scatterDevGenoData(abs(onoffoverlap'),age,agelist, ageBins,~isnan(onoffoverlap), genotype,'onoff overlap')
plotDevGenoData(abs(sustVariation'),age,agelist, ageBins,~isnan(sustVariation), genotype,'sustain variation')
    scatterDevGenoData(abs(sustVariation'),age,agelist, ageBins,~isnan(sustVariation), genotype,'sustain variation')
plotDevGenoData(meanSust',age,agelist, ageBins,~isnan(sustVariation), genotype,'sustain variation')
    scatterDevGenoData(meanSust',age,agelist, ageBins,~isnan(sustVariation), genotype,'sustain variation')
    
    
    
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

plotDevwtData(~isnan(wx),age,agelist, ageBins,genotype==2 & ~isnan(histox), genotype,'fraction center surround')
ylim([0 0.75])

%%% RF area 
plotDevGenoData(wx.*wy,age,agelist, ageBins,1, genotype,'wx from wn (sp/sec)')
scatterDevGenoData(wx.*wy,age,agelist, ageBins,1, genotype,'wx from wn (sp/sec)')

%%% grating response amp
plotDevGenoData(sf_amp',age,agelist, ageBins,1, genotype,'grating sp/sec')
scatterDevGenoData(sf_amp',age,agelist, ageBins,1, genotype,'grating sp/sec')

%%% grrating responsive fraction

plotDevGenoData(sf_amp'>1,age,agelist, ageBins,1, genotype,'% responsive gratings (>1sp/sec)')



%%% peak SF
%%% need to fix
for i = 1:length(drift_all)
    d= drift_all(i);
    sftune = squeeze(mean(nanmean(d.sf_tune(:,1:2,:),2),1));
sftune(1) = sftune(1)-dr_spont(i);
[sf_amp(i) peak_sf(i)] = max(squeeze(mean(nanmean(d.sf_tune(:,1:2,:),2),1)));
end


sfs = [0 0.01 0.02 0.04 0.08 0.16 0.32 nan];
sf_inds = peak_sf;
sf_inds(peak_sf ==0) = 7;
sf_inds(isnan(peak_sf))=8;

plotDevGenoData(sfs(sf_inds)',age,agelist, ageBins,sf_amp>1, genotype,'peak sf (cpd)')
scatterDevGenoData(sfs(sf_inds)',age,agelist, ageBins,sf_amp>1, genotype,'peak sf (cpd)')

plotDevwtData(sfs(sf_inds)',age,agelist, ageBins,sf_amp>1 & genotype==2, genotype,'peak sf (cpd)')

%%% fraction lowpass (greatest response to 0cpd)

plotDevGenoData((peak_sf==1)',age,agelist, ageBins,sf_amp>2, genotype,'fraction lowpass')
ylim([0 1])


%%% fraction OS

plotDevGenoData(driftOSI(:,1)>0.2,age,agelist, ageBins,sf_amp>2, genotype,'fraction orientation selective >0.2')
ylim([ 0 0.5])


%%% fraction DS
plotDevGenoData(driftDSI(:,1)>0.2,age,agelist, ageBins,sf_amp>2, genotype,'fraction direction selective >0.2')
ylim([0 0.5])

%%% fraction SBC
sbc= (wn_cr_dom<-1*10^3);
plotDevGenoData(sbc',age,agelist, ageBins,1, genotype,'fraction sbc')
ylim([0 0.5])

