%%%analyzeGeometrics
%%%use to get spiking relative to geometric stimulus movie
%%% P.R.L. Parker, Cris Niell Lab, 11/12/2015

%%%ordered legend info: column 1 = number, 2 = shape/size, 3 = color, 4 =
%%%position/direction

%use getMovieData to pull spike times, mouse velocity/time, frame info
[spikes mouseT mouseV framenums frametimes] = getMovieData;

%load movie legend
load GeomStim
fpsec = 60;

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

%bin spikes by movie frame
spikehist = zeros(size(spikes,2),size(frametimes,2));
for n = 1:size(spikes,2)
    spikehist(n,:) = histc(spikes{n},frametimes);
end

%rearrange spikes to baseline + stim + post for each full stimulus (5
%repeats); rows are stimulus number, columns are spike rate, 3rd dim is
%cell number)
sortedSpikes = zeros(5*size(ordleg,1),90,size(spikes,2));
for a = 1:size(spikes,2)
    count = firstFrameind;
    for b = 1:5*size(ordleg,1)
        sortedSpikes(b,:,a) = spikehist(a,count-30:count+59);
        count = count + 90;
    end
end
sortedSpikes = sortedSpikes*fpsec; %get rate

%get indices for different stimuli
for n = 1:size(spikes,2)
    bl(n,:) = mean(sortedSpikes(blStim,:,n));
    wh(n,:) = mean(sortedSpikes(whStim,:,n));
    triBig(n,:) = mean(sortedSpikes(triBigStim,:,n));
    triSml(n,:) = mean(sortedSpikes(triSmlStim,:,n));
    sqrBig(n,:) = mean(sortedSpikes(sqrBigStim,:,n));
    sqrSml(n,:) = mean(sortedSpikes(sqrSmlStim,:,n));
    crcBig(n,:) = mean(sortedSpikes(crcBigStim,:,n));
    crcSml(n,:) = mean(sortedSpikes(crcSmlStim,:,n));
    starBig(n,:) = mean(sortedSpikes(starBigStim,:,n));
    starSml(n,:) = mean(sortedSpikes(starSmlStim,:,n));
    centStat(n,:) = mean(sortedSpikes(centStatStim,:,n));
    topStat(n,:) = mean(sortedSpikes(topStatStim,:,n));
    rightStat(n,:) = mean(sortedSpikes(rightStatStim,:,n));
    botStat(n,:) = mean(sortedSpikes(botStatStim,:,n));
    leftStat(n,:) = mean(sortedSpikes(leftStatStim,:,n));
    bot2topMov(n,:) = mean(sortedSpikes(bot2topMovStim,:,n));
    left2rightMov(n,:) = mean(sortedSpikes(left2rightMovStim,:,n));
    top2botMov(n,:) = mean(sortedSpikes(top2botMovStim,:,n));
    right2leftMov(n,:) = mean(sortedSpikes(right2leftMovStim,:,n));
    all(n,:) = mean(sortedSpikes(:,:,n));
end

[fname pname] =uiputfile('*.ps'); psfilename=fullfile(pname,fname);
if exist(psfilename,'file')==2;delete(psfilename);end

%make plots for every cell, each subplot is different stimulus, see key
%file to make life easier
xrange = [-0.5 1];
for n = 1:size(spikes,2)
    figure('Name',sprintf('unit %d',n))
    subplot(4,5,1)
    plot([-30:59]*1/fpsec,bl(n,:))
    xlim(xrange)
    subplot(4,5,2)
    plot([-30:59]*1/fpsec,wh(n,:))
    xlim(xrange)
    subplot(4,5,3)
    plot([-30:59]*1/fpsec,triBig(n,:))
    xlim(xrange)
    subplot(4,5,4)
    plot([-30:59]*1/fpsec,triSml(n,:))
    xlim(xrange)
    subplot(4,5,5)
    plot([-30:59]*1/fpsec,sqrBig(n,:))
    xlim(xrange)
    subplot(4,5,6)
    plot([-30:59]*1/fpsec,sqrSml(n,:))
    xlim(xrange)
    subplot(4,5,7)
    plot([-30:59]*1/fpsec,crcBig(n,:))
    xlim(xrange)
    subplot(4,5,8)
    plot([-30:59]*1/fpsec,crcSml(n,:))
    xlim(xrange)
    subplot(4,5,9)
    plot([-30:59]*1/fpsec,starBig(n,:))
    xlim(xrange)
    subplot(4,5,10)
    plot([-30:59]*1/fpsec,starSml(n,:))
    xlim(xrange)
    subplot(4,5,11)
    plot([-30:59]*1/fpsec,centStat(n,:))
    xlim(xrange)
    subplot(4,5,12)
    plot([-30:59]*1/fpsec,topStat(n,:))
    xlim(xrange)
    subplot(4,5,13)
    plot([-30:59]*1/fpsec,rightStat(n,:))
    xlim(xrange)
    subplot(4,5,14)
    plot([-30:59]*1/fpsec,botStat(n,:))
    xlim(xrange)
    subplot(4,5,15)
    plot([-30:59]*1/fpsec,leftStat(n,:))
    xlim(xrange)
    subplot(4,5,16)
    plot([-30:59]*1/fpsec,bot2topMov(n,:))
    xlim(xrange)
    subplot(4,5,17)
    plot([-30:59]*1/fpsec,left2rightMov(n,:))
    xlim(xrange)
    subplot(4,5,18)
    plot([-30:59]*1/fpsec,top2botMov(n,:))
    xlim(xrange)
    subplot(4,5,19)
    plot([-30:59]*1/fpsec,right2leftMov(n,:))
    xlim(xrange)
    subplot(4,5,20)
    plot([-30:59]*1/fpsec,all(n,:))
    xlim(xrange)
    
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
end

ps2pdf('psfile', psfilename, 'pdffile', [psfilename(1:(end-2)) 'pdf']);
delete(psfilename);

save(fullfile(pname,fname(1:end-3)),'fpsec','sortedSpikes');
