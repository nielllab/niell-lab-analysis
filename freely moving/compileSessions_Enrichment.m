clear all; close all
load('J470a_control.mat');
savePDF=1;
if savePDF
    psfilename = 'C:\analysisPSEnrichment.ps';
    if exist(psfilename,'file')==2;delete(psfilename);end
end

mouse_xy=[];cricket_xy=[]; az=[];

for i=1:length(Data)
    mouse_xy{i,1,1}=Data(i).mouse_xy;
    mouseV{i,:}= Data(i).mouseV;
    cricket_xy{i,1}=Data(i).cricketxy;
    cricketV{i,:}= Data(i).cricketV;
    theta{i,:}= Data(i).theta;
    dTheta{i,:}= diff(Data(i).theta);
    range{i,:}= Data(i).range;
    az = mod(Data(i).az,2*pi); az(az>pi) = az(az>pi)-2*pi;
    azT{i,:,:}= az;
end

nVid=1:length(Data);
useData=3:length(Data);
%%
figure
for vid = 1:length(useData)
    subplot(9,4,vid)  ;
    bar([mean(isnan(mouse_xy{vid,1}(1,:))) mean(isnan(cricket_xy{vid,1}(1,:)))])
    ylabel('% error'); xlim([0.5 2.5]); ylim([0 1]);axis square
    set(gca,'XTick',[1 2])
    set(gca,'XTickLabel',{'mouse','crick'})
end
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-dpsc',psfilename,'-append'); close(gcf); end

%%
figure
for vid=1:length(useData)
    clear dur
    subplot(9,4,vid)
    plot(mouse_xy{vid,1}(1,:),mouse_xy{vid,1}(2,:)); hold on;
    plot(cricket_xy{vid,1}(1,:),cricket_xy{vid,1}(2,:))
    dur=num2str(length(mouse_xy{vid,1}(1,:))/30,'%.2f');
    title([dur,'sec']); axis square
end

if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-dpsc',psfilename,'-append'); close(gcf); end

%%
clear time2cap
figure;% fr2cap =[];time2cap=[]
for vid=1:length(useData)
    clear vidEnd lastQ dur
    dist=range{vid}
    vidEnd=length(range{vid});
    lastQ=round(vidEnd/4);
    app=vidEnd-lastQ:vidEnd;
    plot(dist); hold on; axis square
    plot([app(1),app(1)],[app(1),0],'--'); hold on

 if sum(~isnan(dist(app)))>20
        
        fr2cap(vid,:) = find(dist==min(dist(app)));
        time2cap(vid,:) = (find(dist==min(dist(app))))/30;
        dur=num2str(time2cap(vid,:),'%.2f');
        subplot(9,4,vid)
%         plot([fr2cap(vid,:),fr2cap(vid,:)],[fr2cap(vid,:),0],'go');
        plot([fr2cap(vid,:),fr2cap(vid,:)],[fr2cap(vid,:),0],'g-','linewidth',.75);
        title([dur,'sec']); xlabel('frame num');ylabel('dist to cricket (pix)')
   
 else
fr2cap(vid,:)=NaN;
 end
end
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-dpsc',psfilename,'-append'); close(gcf); end


%%
figure;
mnCapT=nanmean(time2cap);
h=histogram(time2cap,20);
hold on; plot([mnCapT,mnCapT],[mnCapT,0],'g-','linewidth',3);
ylim([0 length(nVid)])
title('time to capure'); xlabel('seconds');ylabel('num of videos')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-dpsc',psfilename,'-append'); close(gcf); end

%%

if savePDF
    pSname='T:\PreyCaptureAnalysis\Data\ControlAnalysis';
    filen=sprintf('%s',ani,'AnalyzedControl','.pdf')
    pdfilename=fullfile(pSname,filen)
    dos(['ps2pdf ' psfilename ' ' pdfilename]);
    delete(psfilename);
end


afilename=sprintf('%s',ani,'Analyzed','.mat')
save(fullfile(pSname, afilename))

