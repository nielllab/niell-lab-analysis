clear all
load('rgcdata04282015.mat')
age(age>30)=31;
% figure
% hist(age,14:30)
% xlabel('age')

genotype = zeros(size(age))+2;


% ageBins = [14 16; 17 21; 22 25; 26 29; 29.5 30.5]';
% ageBins = [14 16; 17 21; 22 25; 26 30]';


ageBins = [16 20; 21 30; 31 33]';


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





%%% example rfs
for a = 1:length(ageBins)
    a
    figure
   set(gcf,'Name',sprintf('age %d',ageBins(1,a)));
   use = find(genotype==2 & age>=ageBins(1,a) & age<=ageBins(2,a));
   length(use)
    if length(use)>36
        use = use(1:36);
    end
    for n= 1:length(use);
        subplot(6,6,n);
        if ~isempty(wn_all(use(n)).svd_xy)
            %%%  imagesc(wn_all(use(n)).sta_final,[-0.1 0.1]); axis equal; colormap jet
            sta = squeeze(wn_all(use(n)).svd_xy(1,:,:));
            sta_fix = ifft2(fft2(sta).*mean_spectrum);
            
            imagesc(sta_fix,[-0.05 0.05]); axis equal; colormap gray
        else          
            imagesc(0,[-0.1 0.1]), colormap gray;
        end
        axis off
        set(gca,'LooseInset',get(gca,'TightInset'))
        title(sprintf('%d',use(n)))
    end
    drawnow
end




%%% figure 1

%ageBins = [15.5 20.5; 20.5 30.5; 27.5 28.5]';


%
ageBins = [16 18; 20 24;25 30; 30.5 31.5]';
%ageBins = [16 20; 21 30; 31 33]';
figure

all = [45 32 2 ]
for i= 1:length(all)
    subplot(2,3,i);
           sta = squeeze(wn_all(all(i)).svd_xy(1,:,:));
            sta_fix = ifft2(fft2(sta).*mean_spectrum);
            if i==3, range = [-0.04 0.04], elseif i ==2, range = [-0.04 0.04],else range = [-0.04 0.04], end
            imagesc(sta_fix,range); axis equal; colormap gray
            if i==3, title('adult'), elseif i==1, title('P16'), else title('P18'), end
            if i ==2, axis xy, end
            hold on
            if i ==1, plot([ 25 38],[68 68],'w','LineWidth',4), end %%% 10 deg scale bar
    axis off
    axis([17 17+75 1 76])
    set(gca,'LooseInset',get(gca,'TightInset'))
end

subplot(2,3,4)
plotDevwtData(~isnan(wx),age,agelist, ageBins,genotype==2 , genotype,'fraction center surround',1)
ylim([0 1.05])
xlim([ 13.5 33])
set(gca,'Xtick',[14 18 22 26 31]);
set(gca,'Xticklabel',{'14','18','22','26','adult'})

subplot(2,3,5);
scatterDevGenoData(0.5*(wx + wy),age,agelist, ageBins,genotype==2, genotype,'RF radius (deg)',1)
%plotDevGenoData(wy,age,agelist, ageBins,genotype==2, genotype,'RF radius',1)
xlim([ 13.5 33])
ylim([0 30])
set(gca,'Xtick',[14 18 22 26 31]);
set(gca,'Xticklabel',{'14','18','22','26','adult'})

subplot(2,3,6);
[resp err]= plotDevGenoData(0.5*(wx+wy),age,agelist, ageBins,genotype==2, genotype,'RF radius (deg)',1)
ylim([0 18])
xlim([ 13.5 33])
set(gca,'Xtick',[14 18 22 26 31]);
set(gca,'Xticklabel',{'14','18','22','26','adult'})

[data group] = labelbyage(double(~isnan(wx)),age,ageBins,genotype==2 );
[p tbl stats] = kruskalwallis(data(group>0),group(group>0))
multcompare(stats,'Ctype','dunn-sidak')

[data group] = labelbyage(0.5*(wx+wy),age,ageBins,genotype==2 & ~isnan(wx'));
[p tbl stats] = kruskalwallis(data(group>0),group(group>0))
multcompare(stats,'Ctype','dunn-sidak')

