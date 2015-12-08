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

%get relevant measures for each stimulus: histograms and FR changes, mean
%and SEM for each
ep1 = 1:0.5/binSize; %indices for baseline period
ep2 = ep1(end)+1:1/binSize; %stimulus period
ep3 = ep2(end)+1:1.5/binSize; %post-stimulus period
for n = 1:size(spikes,2)
    bl(:,:,n) = sortedSpikes(blStim,:,n); %hist w/all trials for this stim
    blHistMean(n,:) = mean(bl(:,:,n)); %mean hist for this stim
    blHistSem(n,:) = nansem(bl(:,:,n)); %SEM hist for this stim
    blMeanDF(n,:) = mean(mean(bl(:,ep2,n),2)-mean(bl(:,ep1,n),2)); %mean firing rate change
    blSemDF(n,:) = nansem(mean(bl(:,ep2,n),2)-mean(bl(:,ep1,n),2)); %SEM firing rate change
    maxFR(1) = max(blHistMean(n,:)); %using to calc max stim-evoked firing for plotting
    wh(:,:,n) = sortedSpikes(whStim,:,n);
    whHistMean(n,:) = mean(wh(:,:,n));
    whHistSem(n,:) = nansem(wh(:,:,n));
    whMeanDF(n,:) = mean(mean(wh(:,ep2,n),2)-mean(wh(:,ep1,n),2));
    whSemDF(n,:) = nansem(mean(wh(:,ep2,n),2)-mean(wh(:,ep1,n),2));
    maxFR(2) = max(whHistMean(n,:));
    triBig(:,:,n) = sortedSpikes(triBigStim,:,n);
    triBigHistMean(n,:) = mean(triBig(:,:,n));
    triBigHistSem(n,:) = nansem(triBig(:,:,n));
    triBigMeanDF(n,:) = mean(mean(triBig(:,ep2,n),2)-mean(triBig(:,ep1,n),2));
    triBigSemDF(n,:) = nansem(mean(triBig(:,ep2,n),2)-mean(triBig(:,ep1,n),2));
    maxFR(3) = max(triBigHistMean(n,:));
    triSml(:,:,n) = sortedSpikes(triSmlStim,:,n);
    triSmlHistMean(n,:) = mean(triSml(:,:,n));
    triSmlHistSem(n,:) = nansem(triSml(:,:,n));
    triSmlMeanDF(n,:) = mean(mean(triSml(:,ep2,n),2)-mean(triSml(:,ep1,n),2));
    triSmlSemDF(n,:) = nansem(mean(triSml(:,ep2,n),2)-mean(triSml(:,ep1,n),2));
    maxFR(4) = max(triSmlHistMean(n,:));
    sqrBig(:,:,n) = sortedSpikes(sqrBigStim,:,n);
    sqrBigHistMean(n,:) = mean(sqrBig(:,:,n));
    sqrBigHistSem(n,:) = nansem(sqrBig(:,:,n));
    sqrBigMeanDF(n,:) = mean(mean(sqrBig(:,ep2,n),2)-mean(sqrBig(:,ep1,n),2));
    sqrBigSemDF(n,:) = nansem(mean(sqrBig(:,ep2,n),2)-mean(sqrBig(:,ep1,n),2));
    maxFR(5) = max(sqrBigHistMean(n,:));
    sqrSml(:,:,n) = sortedSpikes(sqrSmlStim,:,n);
    sqrSmlHistMean(n,:) = mean(sqrSml(:,:,n));
    sqrSmlHistSem(n,:) = nansem(sqrSml(:,:,n));
    sqrSmlMeanDF(n,:) = mean(mean(sqrSml(:,ep2,n),2)-mean(sqrSml(:,ep1,n),2));
    sqrSmlSemDF(n,:) = nansem(mean(sqrSml(:,ep2,n),2)-mean(sqrSml(:,ep1,n),2));
    maxFR(6) = max(sqrSmlHistMean(n,:));
    crcBig(:,:,n) = sortedSpikes(crcBigStim,:,n);
    crcBigHistMean(n,:) = mean(crcBig(:,:,n));
    crcBigHistSem(n,:) = nansem(crcBig(:,:,n));
    crcBigMeanDF(n,:) = mean(mean(crcBig(:,ep2,n),2)-mean(crcBig(:,ep1,n),2));
    crcBigSemDF(n,:) = nansem(mean(crcBig(:,ep2,n),2)-mean(crcBig(:,ep1,n),2));
    maxFR(7) = max(crcBigHistMean(n,:));
    crcSml(:,:,n) = sortedSpikes(crcSmlStim,:,n);
    crcSmlHistMean(n,:) = mean(crcSml(:,:,n));
    crcSmlHistSem(n,:) = nansem(crcSml(:,:,n));
    crcSmlMeanDF(n,:) = mean(mean(crcSml(:,ep2,n),2)-mean(crcSml(:,ep1,n),2));
    crcSmlSemDF(n,:) = nansem(mean(crcSml(:,ep2,n),2)-mean(crcSml(:,ep1,n),2));
    maxFR(8) = max(crcSmlHistMean(n,:));
    starBig(:,:,n) = sortedSpikes(starBigStim,:,n);
    starBigHistMean(n,:) = mean(starBig(:,:,n));
    starBigHistSem(n,:) = nansem(starBig(:,:,n));
    starBigMeanDF(n,:) = mean(mean(starBig(:,ep2,n),2)-mean(starBig(:,ep1,n),2));
    starBigSemDF(n,:) = nansem(mean(starBig(:,ep2,n),2)-mean(starBig(:,ep1,n),2));
    maxFR(9) = max(starBigHistMean(n,:));
    starSml(:,:,n) = sortedSpikes(starSmlStim,:,n);
    starSmlHistMean(n,:) = mean(starSml(:,:,n));
    starSmlHistSem(n,:) = nansem(starSml(:,:,n));
    starSmlMeanDF(n,:) = mean(mean(starSml(:,ep2,n),2)-mean(starSml(:,ep1,n),2));
    starSmlSemDF(n,:) = nansem(mean(starSml(:,ep2,n),2)-mean(starSml(:,ep1,n),2));
    maxFR(10) = max(starSmlHistMean(n,:));
    centStat(:,:,n) = sortedSpikes(centStatStim,:,n);
    centStatHistMean(n,:) = mean(centStat(:,:,n));
    centStatHistSem(n,:) = nansem(centStat(:,:,n));
    centStatMeanDF(n,:) = mean(mean(centStat(:,ep2,n),2)-mean(centStat(:,ep1,n),2));
    centStatSemDF(n,:) = nansem(mean(centStat(:,ep2,n),2)-mean(centStat(:,ep1,n),2));
    maxFR(11) = max(centStatHistMean(n,:));
    topStat(:,:,n) = sortedSpikes(topStatStim,:,n);
    topStatHistMean(n,:) = mean(topStat(:,:,n));
    topStatHistSem(n,:) = nansem(topStat(:,:,n));
    topStatMeanDF(n,:) = mean(mean(topStat(:,ep2,n),2)-mean(topStat(:,ep1,n),2));
    topStatSemDF(n,:) = nansem(mean(topStat(:,ep2,n),2)-mean(topStat(:,ep1,n),2));
    maxFR(12) = max(topStatHistMean(n,:));
    rightStat(:,:,n) = sortedSpikes(rightStatStim,:,n);
    rightStatHistMean(n,:) = mean(rightStat(:,:,n));
    rightStatHistSem(n,:) = nansem(rightStat(:,:,n));
    rightStatMeanDF(n,:) = mean(mean(rightStat(:,ep2,n),2)-mean(rightStat(:,ep1,n),2));
    rightStatSemDF(n,:) = nansem(mean(rightStat(:,ep2,n),2)-mean(rightStat(:,ep1,n),2));
    maxFR(13) = max(rightStatHistMean(n,:));
    botStat(:,:,n) = sortedSpikes(botStatStim,:,n);
    botStatHistMean(n,:) = mean(botStat(:,:,n));
    botStatHistSem(n,:) = nansem(botStat(:,:,n));
    botStatMeanDF(n,:) = mean(mean(botStat(:,ep2,n),2)-mean(botStat(:,ep1,n),2));
    botStatSemDF(n,:) = nansem(mean(botStat(:,ep2,n),2)-mean(botStat(:,ep1,n),2));
    maxFR(14) = max(botStatHistMean(n,:));
    leftStat(:,:,n) = sortedSpikes(leftStatStim,:,n);
    leftStatHistMean(n,:) = mean(leftStat(:,:,n));
    leftStatHistSem(n,:) = nansem(leftStat(:,:,n));
    leftStatMeanDF(n,:) = mean(mean(leftStat(:,ep2,n),2)-mean(leftStat(:,ep1,n),2));
    leftStatSemDF(n,:) = nansem(mean(leftStat(:,ep2,n),2)-mean(leftStat(:,ep1,n),2));
    maxFR(15) = max(leftStatHistMean(n,:));
    bot2topMov(:,:,n) = sortedSpikes(bot2topMovStim,:,n);
    bot2topMovHistMean(n,:) = mean(bot2topMov(:,:,n));
    bot2topMovHistSem(n,:) = nansem(bot2topMov(:,:,n));
    bot2topMovMeanDF(n,:) = mean(mean(bot2topMov(:,ep2,n),2)-mean(bot2topMov(:,ep1,n),2));
    bot2topMovSemDF(n,:) = nansem(mean(bot2topMov(:,ep2,n),2)-mean(bot2topMov(:,ep1,n),2));
    maxFR(16) = max(bot2topMovHistMean(n,:));
    left2rightMov(:,:,n) = sortedSpikes(left2rightMovStim,:,n);
    left2rightMovHistMean(n,:) = mean(left2rightMov(:,:,n));
    left2rightMovHistSem(n,:) = nansem(left2rightMov(:,:,n));
    left2rightMovMeanDF(n,:) = mean(mean(left2rightMov(:,ep2,n),2)-mean(left2rightMov(:,ep1,n),2));
    left2rightMovSemDF(n,:) = nansem(mean(left2rightMov(:,ep2,n),2)-mean(left2rightMov(:,ep1,n),2));
    maxFR(17) = max(left2rightMovHistMean(n,:));
    top2botMov(:,:,n) = sortedSpikes(top2botMovStim,:,n);
    top2botMovHistMean(n,:) = mean(top2botMov(:,:,n));
    top2botMovHistSem(n,:) = nansem(top2botMov(:,:,n));
    top2botMovMeanDF(n,:) = mean(mean(top2botMov(:,ep2,n),2)-mean(top2botMov(:,ep1,n),2));
    top2botMovSemDF(n,:) = nansem(mean(top2botMov(:,ep2,n),2)-mean(top2botMov(:,ep1,n),2));
    maxFR(18) = max(top2botMovHistMean(n,:));
    right2leftMov(:,:,n) = sortedSpikes(right2leftMovStim,:,n);
    right2leftMovHistMean(n,:) = mean(right2leftMov(:,:,n));
    right2leftMovHistSem(n,:) = nansem(right2leftMov(:,:,n));
    right2leftMovMeanDF(n,:) = mean(mean(right2leftMov(:,ep2,n),2)-mean(right2leftMov(:,ep1,n),2));
    right2leftMovSemDF(n,:) = nansem(mean(right2leftMov(:,ep2,n),2)-mean(right2leftMov(:,ep1,n),2));
    maxFR(19) = max(right2leftMovHistMean(n,:));
    allTrials(:,:,n) = sortedSpikes(:,:,n);
    allTrialsHistMean(n,:) = mean(allTrials(:,:,n));
    allTrialsHistSem(n,:) = nansem(allTrials(:,:,n));
    allTrialsMeanDF(n,:) = mean(mean(allTrials(:,ep2,n),2)-mean(allTrials(:,ep1,n),2));
    allTrialsSemDF(n,:) = nansem(mean(allTrials(:,ep2,n),2)-mean(allTrials(:,ep1,n),2));
    maxFR(20) = max(allTrialsHistMean(n,:));
    
    ymax(n) = max(maxFR); %each cell's max FR for plotting
end

% %%
% %create ratio measures for stimulus-specific preferences, normalized to max
% for n = 1:size(spikes,2)
%     colorPrefmean(n,:) = [blMeanDF whMean(n,2)-whMean(n,1)];
% %     colorPref(n,:) = colorPref(n,:)/max(colorPref(n,:));
%     shapePref(n,:) = [triBigFRmean(n,2)-triBigFRmean(n,1) triSmlFRmean(n,2)-triSmlFRmean(n,1) sqrBigFRmean(n,2)-sqrBigFRmean(n,1) sqrSmlFRmean(n,2)-sqrSmlFRmean(n,1) crcBigFRmean(n,2)-crcBigFRmean(n,1) crcSmlFRmean(n,2)-crcSmlFRmean(n,1) starBigFRmean(n,2)-starBigFRmean(n,1) starSmlFRmean(n,2)-starSmlFRmean(n,1)];
% %     shapePref(n,:) = shapePref(n,:)/max(shapePref(n,:));
%     sizePref(n,:) = [mean([triBigFRmean(n,2)-triBigFRmean(n,1) sqrBigFRmean(n,2)-sqrBigFRmean(n,1) crcBigFRmean(n,2)-crcBigFRmean(n,1) starBigFRmean(n,2)-starBigFRmean(n,1)]) mean([triSmlFRmean(n,2)-triSmlFRmean(n,1) sqrSmlFRmean(n,2)-sqrSmlFRmean(n,1) crcSmlFRmean(n,2)-crcSmlFRmean(n,1) starSmlFRmean(n,2)-starSmlFRmean(n,1)])];
% %     sizePref(n,:) = sizePref(n,:)/max(sizePref(n,:));
%     posPref(n,:) = [centStatFRmean(n,2)-centStatFRmean(n,1) topStatFRmean(n,2)-topStatFRmean(n,1) rightStatFRmean(n,2)-rightStatFRmean(n,1) botStatFRmean(n,2)-botStatFRmean(n,1) leftStatFRmean(n,2)-leftStatFRmean(n,1)];
% %     posPref(n,:) = posPref(n,:)/max(posPref(n,:));
%     movPref(n,:) = [bot2topMovFRmean(n,2)-bot2topMovFRmean(n,1) left2rightMovFRmean(n,2)-left2rightMovFRmean(n,1) top2botMovFRmean(n,2)-top2botMovFRmean(n,1) right2leftMovFRmean(n,2)-right2leftMovFRmean(n,1)];
% %     movPref(n,:) = movPref(n,:)/max(movPref(n,:));
%     movstatPref(n,:) = [mean([bot2topMovFRmean(n,2)-bot2topMovFRmean(n,1) left2rightMovFRmean(n,2)-left2rightMovFRmean(n,1) top2botMovFRmean(n,2)-top2botMovFRmean(n,1) right2leftMovFRmean(n,2)-right2leftMovFRmean(n,1)]) mean([centStatFRmean(n,2)-centStatFRmean(n,1) topStatFRmean(n,2)-topStatFRmean(n,1) rightStatFRmean(n,2)-rightStatFRmean(n,1) botStatFRmean(n,2)-botStatFRmean(n,1) leftStatFRmean(n,2)-leftStatFRmean(n,1)])];
% %     movstatPref(n,:) = movstatPref(n,:)/max(movstatPref(n,:));
% end

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
    shadedErrorBar(xvals,blHistMean(n,:),blHistSem(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,2)
    shadedErrorBar(xvals,whHistMean(n,:),whHistSem(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,3)
    shadedErrorBar(xvals,triBigHistMean(n,:),triBigHistSem(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,4)
    shadedErrorBar(xvals,triSmlHistMean(n,:),triSmlHistSem(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,5)
    shadedErrorBar(xvals,sqrBigHistMean(n,:),sqrBigHistSem(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,6)
    shadedErrorBar(xvals,sqrSmlHistMean(n,:),sqrSmlHistSem(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,7)
    shadedErrorBar(xvals,crcBigHistMean(n,:),crcBigHistSem(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,8)
    shadedErrorBar(xvals,crcSmlHistMean(n,:),crcBigHistSem(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,9)
    shadedErrorBar(xvals,starBigHistMean(n,:),starBigHistSem(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,10)
    shadedErrorBar(xvals,starSmlHistMean(n,:),starSmlHistSem(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,11)
    shadedErrorBar(xvals,centStatHistMean(n,:),centStatHistSem(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,12)
    shadedErrorBar(xvals,topStatHistMean(n,:),topStatHistSem(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,13)
    shadedErrorBar(xvals,rightStatHistMean(n,:),rightStatHistSem(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,14)
    shadedErrorBar(xvals,botStatHistMean(n,:),botStatHistSem(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,15)
    shadedErrorBar(xvals,leftStatHistMean(n,:),leftStatHistSem(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,16)
    shadedErrorBar(xvals,bot2topMovHistMean(n,:),bot2topMovHistSem(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,17)
    shadedErrorBar(xvals,left2rightMovHistMean(n,:),left2rightMovHistSem(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,18)
    shadedErrorBar(xvals,top2botMovHistMean(n,:),top2botMovHistSem(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,19)
    shadedErrorBar(xvals,right2leftMovHistMean(n,:),right2leftMovHistSem(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    subplot(4,5,20)
    shadedErrorBar(xvals,allTrialsHistMean(n,:),allTrialsHistSem(n,:))
    axis([xrange(1) xrange(2) yrange(1) ymax(n)]);
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
    
    figure('Name',sprintf('unit %d',n))
    subplot(3,3,1)
    errorbar(1:2,[blMeanDF(n) whMeanDF(n)],[blSemDF(n) whSemDF(n)])
    axis([0 3 -10 10])
    set(gca,'Xtick',1:2,'XTickLabel',{'black', 'white'})
    subplot(3,3,2)
    errorbar(1:8,[triBigMeanDF(n) triSmlMeanDF(n) sqrBigMeanDF(n) sqrSmlMeanDF(n) crcBigMeanDF(n) crcSmlMeanDF(n) starBigMeanDF(n) starSmlMeanDF(n)],...
        [triBigSemDF(n) triSmlSemDF(n) sqrBigSemDF(n) sqrSmlSemDF(n) crcBigSemDF(n) crcSmlSemDF(n) starBigSemDF(n) starSmlSemDF(n)])
    axis([0 9 -10 10])
    set(gca,'Xtick',1:8,'XTickLabel',{'TRI', 'tri', 'SQR', 'sqr', 'CRC', 'crc', 'STAR', 'star'})
    subplot(3,3,3)
    errorbar(1:2,[mean([triBigMeanDF(n) sqrBigMeanDF(n) crcBigMeanDF(n) starBigMeanDF(n)]) mean([triSmlMeanDF(n) sqrSmlMeanDF(n) crcSmlMeanDF(n) starSmlMeanDF(n)])],...
        [nansem([triBigMeanDF(n) sqrBigMeanDF(n) crcBigMeanDF(n) starBigMeanDF(n)]) nansem([triSmlMeanDF(n) sqrSmlMeanDF(n) crcSmlMeanDF(n) starSmlMeanDF(n)])])
    axis([0 3 -10 10])
    set(gca,'Xtick',1:2,'XTickLabel',{'big', 'small'})
    subplot(3,3,4)
    errorbar(1:5,[centStatMeanDF(n) topStatMeanDF(n) rightStatMeanDF(n) botStatMeanDF(n) leftStatMeanDF(n)],...
        [centStatSemDF(n) topStatSemDF(n) rightStatSemDF(n) botStatSemDF(n) leftStatSemDF(n)])
    axis([0 6 -10 10])
    set(gca,'Xtick',1:5,'XTickLabel',{'center', 'top', 'right', 'bottom', 'left'})    
    subplot(3,3,5)
    errorbar(1:4,[bot2topMovMeanDF(n) left2rightMovMeanDF(n) top2botMovMeanDF(n) right2leftMovMeanDF(n)],...
        [bot2topMovSemDF(n) left2rightMovSemDF(n) top2botMovSemDF(n) right2leftMovSemDF(n)])
    axis([0 5 -10 10])
    set(gca,'Xtick',1:4,'XTickLabel',{'B2T', 'L2R', 'T2B', 'R2L'})
    subplot(3,3,6)
    errorbar(1:2,[mean([bot2topMovMeanDF(n) left2rightMovMeanDF(n) top2botMovMeanDF(n) right2leftMovMeanDF(n)]) mean([centStatMeanDF(n) topStatMeanDF(n) rightStatMeanDF(n) botStatMeanDF(n) leftStatMeanDF(n)])],...
        [nansem([bot2topMovMeanDF(n) left2rightMovMeanDF(n) top2botMovMeanDF(n) right2leftMovMeanDF(n)]) mean([centStatMeanDF(n) topStatMeanDF(n) rightStatMeanDF(n) botStatMeanDF(n) leftStatMeanDF(n)])])
    axis([0 3 -10 10])
    set(gca,'Xtick',1:2,'XTickLabel',{'mov', 'stat'})
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
end

ps2pdf('psfile', psfilename, 'pdffile', [psfilename(1:(end-2)) 'pdf']);
delete(psfilename);

save(fullfile(pname,fname(1:end-3)),'fpsec','binSize','sortedSpikes');
