close all; clear all
pathname = {'D:\gdrive\cohort5_100919\Cohort5\timestamp info\'};

fileList = [];
for i = 1:length(pathname)
    fileList = [fileList ; dir([pathname{i} '*DLC*.csv'])];
end

for f = 1:length(fileList)
    csvFile = fileList(f).name;
    syncFile= strrep(csvFile,'DLC_resnet50_preycapSep30shuffle2_250000filtered.csv','_sync.mat');
    thispath = fileList(i).folder;
    
    trackData = csvread(fullfile(thispath,csvFile),3,0);
    load(fullfile(thispath,syncFile),'ttlSync');
    ttlSync(end+1) = ttlSync(end); %%% trimmed one timepoint off
    
       crickP = trackData(:,[4 7 10]);

    clear crick
    for i = 1:2
        crickPos = trackData(:,[1 4 7]+i);
    crickPos(crickP<0.95) = NaN;
    crick(:,i) = nanmean(crickPos,2);
    end
    
    nose = trackData(:,11:12);
    earL = trackData(:,14:15);
    earR = trackData(:,17:18);
    
    headP = trackData(:,[13 16 19]);
    
    
    figure
    subplot(3,1,1);
    plot(crickP); ylim([0 1]); title('cricket P'); xlim([1 length(ttlSync)])
    subplot(3,1,2);
    plot(headP); ylim([0 1]); title('head P');xlim([1 length(ttlSync)])
    subplot(3,1,3);
    plot(ttlSync); ylim([-0.1 1.1]); title('ttl'); xlim([1 length(ttlSync)])
    
    head = 0.5*(earL + earR);
    
    range = sqrt((crick(:,1)-head(:,1)).^2 + (crick(:,2) - head(:,2)).^2);
    
    rangeVec = crick-head;
    rangeTheta = atan2d(rangeVec(:,1),rangeVec(:,2));
    
    earVec = earR-earL;
    earTheta = atan2d(earVec(:,1), earVec(:,2));
    headTheta = earTheta+90;
    
    az = rangeTheta-headTheta;
    
    az = mod(az+180,360) - 180;
    
    
figure
subplot(3,1,1);
plot(range); title('range');xlim([1 length(ttlSync)])
subplot(3,1,2);
plot(az); title('azimuth'); xlim([1 length(ttlSync)])
subplot(3,1,3);
plot(ttlSync); ylim([-0.1 1.1]); title(csvFile(1:30)); ylabel('TTL'); xlim([1 length(ttlSync)])
    
rangebins = 50:100:1500;
azbins = -175:10:175;

rangeHist(:,1,f) = hist(range(ttlSync==0),rangebins)/sum(ttlSync==0);
rangeHist(:,2,f) = hist(range(ttlSync==1),rangebins)/sum(ttlSync==1);

azHist(:,1,f) = hist(az(ttlSync==0 & range>100),azbins)/sum(ttlSync==0 & range>100);
azHist(:,2,f) = hist(az(ttlSync==1 & range>100),azbins)/sum(ttlSync==1 & range>100);

figure
subplot(1,2,1);
plot(rangebins, rangeHist(:,:,f)); xlabel('range'); title(csvFile(1:30));
subplot(1,2,2);
plot(azbins,azHist(:,:,f)); xlabel('az');

    
end

azAll = nanmean(azHist,3);
rangeAll = nanmean(rangeHist,3);

figure
subplot(1,2,1);
plot(rangebins,rangeAll); xlabel('range'); title('mean across sessions')
subplot(1,2,2);
plot(azbins,azAll); xlabel('az');

