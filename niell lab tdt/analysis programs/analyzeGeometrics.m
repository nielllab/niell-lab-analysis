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

%bin spikes
bins = frametimes(firstFrameind)-0.5:binSize:frametimes(firstFrameind)+1080-binSize; %create binning edges w/half second baseline before first stimulus
spikehist = zeros(size(spikes,2),size(bins,2));
for n = 1:size(spikes,2)
    spikehist(n,:) = histc(spikes{n},bins);
end

%rearrange spikes to baseline + stim + post for each full stimulus (5
%repeats); rows are stimulus number, columns are spike rate, 3rd dim is
%cell number)
sortedSpikes = zeros(5*size(ordleg,1),1500/(binSize*1000),size(spikes,2));
for unit = 1:size(spikes,2)
    count = 1;
    for stimNum = 1:5*size(ordleg,1)
        sortedSpikes(stimNum,:,unit) = spikehist(unit,count:count+(1500/(binSize*1000))-1);
        count = count + (1500/(binSize*1000));
    end
end
sortedSpikes = sortedSpikes/binSize; %get rate

%get indices for different stimuli
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

% get change in firing rate stim/baseline
for n = 1:size(spikes,2)
    blFR(n) = mean(bl(n,1+500/(binSize*1000):1000/(binSize*1000)))/...
        mean(bl(n,1:500/(binSize*1000)));
    whFR(n) = mean(wh(n,1+500/(binSize*1000):1000/(binSize*1000)))/...
        mean(wh(n,1:500/(binSize*1000)));
    triBigFR(n) = mean(triBig(n,1+500/(binSize*1000):1000/(binSize*1000)))/...
        mean(triBig(n,1:500/(binSize*1000)));
    triSmlFR(n) = mean(triSml(n,1+500/(binSize*1000):1000/(binSize*1000)))/...
        mean(triSml(n,1:500/(binSize*1000)));
    sqrBigFR(n) = mean(sqrBig(n,1+500/(binSize*1000):1000/(binSize*1000)))/...
        mean(sqrBig(n,1:500/(binSize*1000)));
    sqrSmlFR(n) = mean(sqrSml(n,1+500/(binSize*1000):1000/(binSize*1000)))/...
        mean(sqrSml(n,1:500/(binSize*1000)));
    crcBigFR(n) = mean(crcBig(n,1+500/(binSize*1000):1000/(binSize*1000)))/...
        mean(crcBig(n,1:500/(binSize*1000)));
    crcSmlFR(n) = mean(crcSml(n,1+500/(binSize*1000):1000/(binSize*1000)))/...
        mean(crcSml(n,1:500/(binSize*1000)));
    starBigFR(n) = mean(starBig(n,1+500/(binSize*1000):1000/(binSize*1000)))/...
        mean(starBig(n,1:500/(binSize*1000)));
    starSmlFR(n) = mean(starSml(n,1+500/(binSize*1000):1000/(binSize*1000)))/...
        mean(starSml(n,1:500/(binSize*1000)));
    centStatFR(n) = mean(centStat(n,1+500/(binSize*1000):1000/(binSize*1000)))/...
        mean(centStat(n,1:500/(binSize*1000)));
    topStatFR(n) = mean(topStat(n,1+500/(binSize*1000):1000/(binSize*1000)))/...
        mean(topStat(n,1:500/(binSize*1000)));
    rightStatFR(n) = mean(rightStat(n,1+500/(binSize*1000):1000/(binSize*1000)))/...
        mean(rightStat(n,1:500/(binSize*1000)));
    botStatFR(n) = mean(botStat(n,1+500/(binSize*1000):1000/(binSize*1000)))/...
        mean(botStat(n,1:500/(binSize*1000)));
    leftStatFR(n) = mean(leftStat(n,1+500/(binSize*1000):1000/(binSize*1000)))/...
        mean(leftStat(n,1:500/(binSize*1000)));
    bot2topMovFR(n) = mean(bot2topMov(n,1+500/(binSize*1000):1000/(binSize*1000)))/...
        mean(bot2topMov(n,1:500/(binSize*1000)));
    left2rightMovFR(n) = mean(left2rightMov(n,1+500/(binSize*1000):1000/(binSize*1000)))/...
        mean(left2rightMov(n,1:500/(binSize*1000)));
    top2botMovFR(n) = mean(top2botMov(n,1+500/(binSize*1000):1000/(binSize*1000)))/...
        mean(top2botMov(n,1:500/(binSize*1000)));
    right2leftMovFR(n) = mean(right2leftMov(n,1+500/(binSize*1000):1000/(binSize*1000)))/...
        mean(right2leftMov(n,1:500/(binSize*1000)));
    allstimFR(n) = mean(allstim(n,1+500/(binSize*1000):1000/(binSize*1000)))/...
        mean(allstim(n,1:500/(binSize*1000)));
end

%create ratio measures for stimulus-specific preferences, normalized to max
for n = 1:size(spikes,2)
    colorPref(n,:) = [blFR(n) whFR(n)];
    colorPref(n,:) = colorPref(n,:)/max(colorPref(n,:));
    shapePref(n,:) = [triBigFR(n) triSmlFR(n) sqrBigFR(n) sqrSmlFR(n) crcBigFR(n) crcSmlFR(n) starBigFR(n) starSmlFR(n)];
    shapePref(n,:) = shapePref(n,:)/max(shapePref(n,:));
    sizePref(n,:) = [mean([triBigFR(n) sqrBigFR(n) crcBigFR(n) starBigFR(n)]) mean([triSmlFR(n) sqrSmlFR(n) crcSmlFR(n) starSmlFR(n)])];
    sizePref(n,:) = sizePref(n,:)/max(sizePref(n,:));
    posPref(n,:) = [centStatFR(n) topStatFR(n) rightStatFR(n) botStatFR(n) leftStatFR(n)];
    posPref(n,:) = posPref(n,:)/max(posPref(n,:));
    movPref(n,:) = [bot2topMovFR(n) left2rightMovFR(n) top2botMovFR(n) right2leftMovFR(n)];
    movPref(n,:) = movPref(n,:)/max(movPref(n,:));
    movstatPref(n,:) = [mean([bot2topMovFR(n) left2rightMovFR(n) top2botMovFR(n) right2leftMovFR(n)]) mean([centStatFR(n) topStatFR(n) rightStatFR(n) botStatFR(n) leftStatFR(n)])];
    movstatPref(n,:) = movstatPref(n,:)/max(movstatPref(n,:));
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
end

ps2pdf('psfile', psfilename, 'pdffile', [psfilename(1:(end-2)) 'pdf']);
delete(psfilename);

save(fullfile(pname,fname(1:end-3)),'fpsec','sortedSpikes');
