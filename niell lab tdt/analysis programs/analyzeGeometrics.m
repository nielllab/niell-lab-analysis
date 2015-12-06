%%%analyzeGeometrics
%%%use to get spiking relative to geometric stimulus movie
%%% P.R.L. Parker, Cris Niell Lab, 11/12/2015

%%%ordered legend info: column 1 = number, 2 = shape/size, 3 = color, 4 =
%%%position/direction

%use getMovieData to pull spike times, mouse velocity/time, frame info
[spikes mouseT mouseV framenums frametimes] = getMovieData;

%load movie legend
load C:\Users\nlab\Desktop\Ephys\Stimuli\GeomStim
fpsec = 60; %movie frames per second
binSize = 0.05 %in seconds

ordleg = repmat(ordleg,5,1); %movie plays five times sequentially
%stimulus parameters: black, white, shape/size, location/direction
blStim = find(ordleg(:,3)==1);
whStim = find(ordleg(:,3)==2);
triBigStim = find(ordleg(:,2)==1);
triSmlStim = find(ordleg(:,2)==2);
sqrBigStim = find(ordleg(:,2)==3);
sqrSmlStim = find(ordleg(:,2)==4);
crcBigStim = find(ordleg(:,2)==5);
crcSmlStim = find(ordleg(:,2)==6);
starBigStim = find(ordleg(:,2)==7);
starSmlStim = find(ordleg(:,2)==8);
centStatStim = find(ordleg(:,4)==1);
topStatStim = find(ordleg(:,4)==2);
rightStatStim = find(ordleg(:,4)==3);
botStatStim = find(ordleg(:,4)==4);
leftStatStim = find(ordleg(:,4)==5);
bot2topMovStim = find(ordleg(:,4)==6);
left2rightMovStim = find(ordleg(:,4)==7);
top2botMovStim = find(ordleg(:,4)==8);
right2leftMovStim = find(ordleg(:,4)==9);

firstStim = ceil(framenums(1)/90) + 1; %90 frames per stimulus, start on the first full stimulus
firstFrame = firstStim*90; % first usable frame number
firstFrameind = find(framenums==firstFrame); firstFrameind = firstFrameind(1); %index of it
%% 

%bin spikes by desired binSize
spikebins = frametimes(firstFrameind)-0.5:binSize:frametimes(firstFrameind)+1080-binSize; %create binning edges w/half second baseline before first stimulus
spikeHist = zeros(size(spikes,2),size(spikebins,2));
for n = 1:size(spikes,2)
    spikeHist(n,:) = histc(spikes{n},spikebins);
end

%bin spikes by pre, stim, post (0.5 sec each)
epochbins = frametimes(firstFrameind)-0.5:0.5:frametimes(firstFrameind)+1080-0.5; %create binning edges w/half second baseline before first stimulus
epochHist = zeros(size(spikes,2),size(epochbins,2));
for n = 1:size(spikes,2)
    epochHist(n,:) = histc(spikes{n},epochbins);
end

%rearrange spikes to baseline + stim + post for each full stimulus (5
%repeats); rows are stimulus number, columns are spike rate, 3rd dim is
%cell number)
sortedSpikes = zeros(size(ordleg,1),1500/(binSize*1000),size(spikes,2));
for unit = 1:size(spikes,2)
    count = 1;
    for stimNum = 1:size(ordleg,1)
        sortedSpikes(stimNum,:,unit) = spikeHist(unit,count:count+(1500/(binSize*1000))-1);
        count = count + (1500/(binSize*1000));
    end
end
sortedSpikes = sortedSpikes/binSize; %get rate

%generate averaged histograms for each type of stimulus
for n = 1:size(spikes,2)
    bl(n,:) = mean(sortedSpikes(blStim,:,n));
    maxFR(1) = max(bl(n,:));
    wh(n,:) = mean(sortedSpikes(whStim,:,n));
    maxFR(2) = max(wh(n,:));
    triBig(n,:) = mean(sortedSpikes(triBigStim,:,n));
    maxFR(3) = max(triBig(n,:));
    triSml(n,:) = mean(sortedSpikes(triSmlStim,:,n));
    maxFR(4) = max(triSml(n,:));
    sqrBig(n,:) = mean(sortedSpikes(sqrBigStim,:,n));
    maxFR(5) = max(sqrBig(n,:));
    sqrSml(n,:) = mean(sortedSpikes(sqrSmlStim,:,n));
    maxFR(6) = max(sqrSml(n,:));
    crcBig(n,:) = mean(sortedSpikes(crcBigStim,:,n));
    maxFR(7) = max(crcBig(n,:));
    crcSml(n,:) = mean(sortedSpikes(crcSmlStim,:,n));
    maxFR(8) = max(crcSml(n,:));
    starBig(n,:) = mean(sortedSpikes(starBigStim,:,n));
    maxFR(9) = max(starBig(n,:));
    starSml(n,:) = mean(sortedSpikes(starSmlStim,:,n));
    maxFR(10) = max(starSml(n,:));
    centStat(n,:) = mean(sortedSpikes(centStatStim,:,n));
    maxFR(11) = max(centStat(n,:));
    topStat(n,:) = mean(sortedSpikes(topStatStim,:,n));
    maxFR(12) = max(topStat(n,:));
    rightStat(n,:) = mean(sortedSpikes(rightStatStim,:,n));
    maxFR(13) = max(rightStat(n,:));
    botStat(n,:) = mean(sortedSpikes(botStatStim,:,n));
    maxFR(14) = max(botStat(n,:));
    leftStat(n,:) = mean(sortedSpikes(leftStatStim,:,n));
    maxFR(15) = max(leftStat(n,:));
    bot2topMov(n,:) = mean(sortedSpikes(bot2topMovStim,:,n));
    maxFR(16) = max(bot2topMov(n,:));
    left2rightMov(n,:) = mean(sortedSpikes(left2rightMovStim,:,n));
    maxFR(17) = max(left2rightMov(n,:));
    top2botMov(n,:) = mean(sortedSpikes(top2botMovStim,:,n));
    maxFR(18) = max(top2botMov(n,:));
    right2leftMov(n,:) = mean(sortedSpikes(right2leftMovStim,:,n));
    maxFR(19) = max(right2leftMov(n,:));
    allstim(n,:) = mean(sortedSpikes(:,:,n));
    maxFR(20) = max(allstim(n,:));
    ymax(n) = max(maxFR);
end

%%
%same rearranging except for firing rate per entire stimulus period (pre stim post)
sortedEpochs = zeros(size(ordleg,1),1500/(0.5*1000),size(spikes,2));
for unit = 1:size(spikes,2)
    count = 1;
    for stimNum = 1:size(ordleg,1)
        sortedEpochs(stimNum,:,unit) = epochHist(unit,count:count+(1500/(0.5*1000))-1);
        count = count + (1500/(0.5*1000));
    end
end
sortedEpochs = sortedEpochs/0.5; %get rate

% get change in firing rate stim/baseline
for n = 1:size(spikes,2)
    blFRmean(n,:) = mean(sortedEpochs(blStim,:,n));
    blFRse(n,:) = std(sortedEpochs(blStim,:,n))/sqrt(size(blStim,1));
    whFRmean(n,:) = mean(sortedEpochs(whStim,:,n));
    whFRse(n,:) = std(sortedEpochs(whStim,:,n))/sqrt(size(whStim,1));
    triBigFRmean(n,:) = mean(sortedEpochs(triBigStim,:,n));
    triBigFRstd(n,:) = std(sortedEpochs(triBigStim,:,n))/sqrt(size(triBigStim,1));
    triSmlFRmean(n,:) = mean(sortedEpochs(triSmlStim,:,n));
    triSmlFRstd(n,:) = std(sortedEpochs(triSmlStim,:,n))/sqrt(size(triSmlStim,1));
    sqrBigFRmean(n,:) = mean(sortedEpochs(sqrBigStim,:,n));
    sqrBigFRstd(n,:) = std(sortedEpochs(sqrBigStim,:,n))/sqrt(size(sqrBigStim,1));
    sqrSmlFRmean(n,:) = mean(sortedEpochs(sqrSmlStim,:,n));
    sqrSmlFRstd(n,:) = std(sortedEpochs(sqrSmlStim,:,n))/sqrt(size(sqrSmlStim,1));
    crcBigFRmean(n,:) = mean(sortedEpochs(crcBigStim,:,n));
    crcBigFRstd(n,:) = std(sortedEpochs(crcBigStim,:,n))/sqrt(size(crcBigStim,1));
    crcSmlFRmean(n,:) = mean(sortedEpochs(crcSmlStim,:,n));
    crcBigFRstd(n,:) = std(sortedEpochs(crcSmlStim,:,n))/sqrt(size(crcSmlStim,1));
    starBigFRmean(n,:) = mean(sortedEpochs(starBigStim,:,n));
    starBigFRstd(n,:) = std(sortedEpochs(starBigStim,:,n))/sqrt(size(starBigStim,1));
    starSmlFRmean(n,:) = mean(sortedEpochs(starSmlStim,:,n));
    starSmlFRstd(n,:) = std(sortedEpochs(starSmlStim,:,n))/sqrt(size(starSmlStim,1));
    centStatFRmean(n,:) = mean(sortedEpochs(centStatStim,:,n));
    centStatFRstd(n,:) = std(sortedEpochs(centStatStim,:,n))/sqrt(size(centStatStim,1));
    topStatFRmean(n,:) = mean(sortedEpochs(topStatStim,:,n));
    topStatFRstd(n,:) = std(sortedEpochs(topStatStim,:,n))/sqrt(size(topStatStim,1));
    rightStatFRmean(n,:) = mean(sortedEpochs(rightStatStim,:,n));
    rightStatFRstd(n,:) = std(sortedEpochs(rightStatStim,:,n))/sqrt(size(rightStatStim,1));
    botStatFRmean(n,:) = mean(sortedEpochs(botStatStim,:,n));
    botStatFRstd(n,:) = std(sortedEpochs(botStatStim,:,n))/sqrt(size(botStatStim,1));
    leftStatFRmean(n,:) = mean(sortedEpochs(leftStatStim,:,n));
    leftStatFRstd(n,:) = std(sortedEpochs(leftStatStim,:,n))/sqrt(size(leftStatStim,1));
    bot2topMovFRmean(n,:) = mean(sortedEpochs(bot2topMovStim,:,n));
    bot2topMovFRstd(n,:) = std(sortedEpochs(bot2topMovStim,:,n))/sqrt(size(bot2topMovStim,1));
    left2rightMovFRmean(n,:) = mean(sortedEpochs(left2rightMovStim,:,n));
    left2rightMovFRstd(n,:) = std(sortedEpochs(left2rightMovStim,:,n))/sqrt(size(left2rightMovStim,1));
    top2botMovFRmean(n,:) = mean(sortedEpochs(top2botMovStim,:,n));
    top2botMovFRstd(n,:) = std(sortedEpochs(top2botMovStim,:,n))/sqrt(size(top2botMovStim,1));
    right2leftMovFRmean(n,:) = mean(sortedEpochs(right2leftMovStim,:,n));
    right2leftMovFRstd(n,:) = std(sortedEpochs(right2leftMovStim,:,n))/sqrt(size(right2leftMovStim,1));
    allFRmean(n,:) = mean(sortedEpochs(:,:,n));
    allFRse(n,:) = std(sortedEpochs(:,:,n))/sqrt(size(ordleg,1));
end

%create ratio measures for stimulus-specific preferences, normalized to max
for n = 1:size(spikes,2)
    colorPref(n,:) = [blFRmean(n,2)-blFRmean(n,1) whFRmean(n,2)-whFRmean(n,1)];
%     colorPref(n,:) = colorPref(n,:)/max(colorPref(n,:));
    shapePref(n,:) = [triBigFRmean(n,2)-triBigFRmean(n,1) triSmlFRmean(n,2)-triSmlFRmean(n,1) sqrBigFRmean(n,2)-sqrBigFRmean(n,1) sqrSmlFRmean(n,2)-sqrSmlFRmean(n,1) crcBigFRmean(n,2)-crcBigFRmean(n,1) crcSmlFRmean(n,2)-crcSmlFRmean(n,1) starBigFRmean(n,2)-starBigFRmean(n,1) starSmlFRmean(n,2)-starSmlFRmean(n,1)];
%     shapePref(n,:) = shapePref(n,:)/max(shapePref(n,:));
    sizePref(n,:) = [mean([triBigFRmean(n,2)-triBigFRmean(n,1) sqrBigFRmean(n,2)-sqrBigFRmean(n,1) crcBigFRmean(n,2)-crcBigFRmean(n,1) starBigFRmean(n,2)-starBigFRmean(n,1)]) mean([triSmlFRmean(n,2)-triSmlFRmean(n,1) sqrSmlFRmean(n,2)-sqrSmlFRmean(n,1) crcSmlFRmean(n,2)-crcSmlFRmean(n,1) starSmlFRmean(n,2)-starSmlFRmean(n,1)])];
%     sizePref(n,:) = sizePref(n,:)/max(sizePref(n,:));
    posPref(n,:) = [centStatFRmean(n,2)-centStatFRmean(n,1) topStatFRmean(n,2)-topStatFRmean(n,1) rightStatFRmean(n,2)-rightStatFRmean(n,1) botStatFRmean(n,2)-botStatFRmean(n,1) leftStatFRmean(n,2)-leftStatFRmean(n,1)];
%     posPref(n,:) = posPref(n,:)/max(posPref(n,:));
    movPref(n,:) = [bot2topMovFRmean(n,2)-bot2topMovFRmean(n,1) left2rightMovFRmean(n,2)-left2rightMovFRmean(n,1) top2botMovFRmean(n,2)-top2botMovFRmean(n,1) right2leftMovFRmean(n,2)-right2leftMovFRmean(n,1)];
%     movPref(n,:) = movPref(n,:)/max(movPref(n,:));
    movstatPref(n,:) = [mean([bot2topMovFRmean(n,2)-bot2topMovFRmean(n,1) left2rightMovFRmean(n,2)-left2rightMovFRmean(n,1) top2botMovFRmean(n,2)-top2botMovFRmean(n,1) right2leftMovFRmean(n,2)-right2leftMovFRmean(n,1)]) mean([centStatFRmean(n,2)-centStatFRmean(n,1) topStatFRmean(n,2)-topStatFRmean(n,1) rightStatFRmean(n,2)-rightStatFRmean(n,1) botStatFRmean(n,2)-botStatFRmean(n,1) leftStatFRmean(n,2)-leftStatFRmean(n,1)])];
%     movstatPref(n,:) = movstatPref(n,:)/max(movstatPref(n,:));
end

[fname pname] =uiputfile('*.ps'); psfilename=fullfile(pname,fname);
if exist(psfilename,'file')==2;delete(psfilename);end

%make plots for every cell, each subplot is different stimulus, see key
%file to make life easier

xvals = [-0.5:binSize:1-binSize];
xrange = [-0.5 1];
yrange = [0 40];

for n = 1:size(spikes,2)
    figure('Name',sprintf('unit %d',n))
    subplot(4,5,1)
    plot(xvals,bl(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,2)
    plot(xvals,wh(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,3)
    plot(xvals,triBig(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,4)
    plot(xvals,triSml(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,5)
    plot(xvals,sqrBig(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,6)
    plot(xvals,sqrSml(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,7)
    plot(xvals,crcBig(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,8)
    plot(xvals,crcSml(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,9)
    plot(xvals,starBig(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,10)
    plot(xvals,starSml(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,11)
    plot(xvals,centStat(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,12)
    plot(xvals,topStat(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,13)
    plot(xvals,rightStat(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,14)
    plot(xvals,botStat(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,15)
    plot(xvals,leftStat(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,16)
    plot(xvals,bot2topMov(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,17)
    plot(xvals,left2rightMov(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,18)
    plot(xvals,top2botMov(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,19)
    plot(xvals,right2leftMov(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,20)
    plot(xvals,allstim(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
    
    figure('Name',sprintf('unit %d',n))
    subplot(3,3,1)
    bar(1:2,colorPref(n,:))
    xlim([0 3])
%     axis([0 3 0 1])
    set(gca,'Xtick',1:2,'XTickLabel',{'black', 'white'})
    subplot(3,3,2)
    bar(1:8,shapePref(n,:))
    xlim([0 9])
%     axis([0 9 0 1])
    set(gca,'Xtick',1:8,'XTickLabel',{'triBig', 'triSml', 'sqrBig', 'sqrSml', 'crcBig', 'crcSml', 'starBig', 'starSml'})
    subplot(3,3,3)
    bar(1:2,sizePref(n,:))
    xlim([0 3])
%     axis([0 3 0 1])
    set(gca,'Xtick',1:2,'XTickLabel',{'big', 'small'})
    subplot(3,3,4)
    bar(1:5,posPref(n,:))
    xlim([0 6])
%     axis([0 6 0 1])
    set(gca,'Xtick',1:5,'XTickLabel',{'centStat', 'topStat', 'rightStat', 'botStat', 'leftStat'})    
    subplot(3,3,5)
    bar(1:4,movPref(n,:))
    xlim([0 5])
%     axis([0 5 0 1])
    set(gca,'Xtick',1:4,'XTickLabel',{'bot2topMov', 'left2rightMov', 'top2botMov', 'right2leftMov'})
    subplot(3,3,6)
    bar(1:2,movstatPref(n,:))
    xlim([0 3])
%     axis([0 3 0 1])
    set(gca,'Xtick',1:2,'XTickLabel',{'moving', 'static'})
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
end

ps2pdf('psfile', psfilename, 'pdffile', [psfilename(1:(end-2)) 'pdf']);
delete(psfilename);

save(fullfile(pname,fname(1:end-3)),'fpsec','sortedSpikes','sortedEpochs','colorPref','shapePref','sizePref','posPref','movPref','movstatPref');
