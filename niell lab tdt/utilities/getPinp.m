function pinp = getPinp(clustfile,afile, block,trainPeriod,freq,showImg,makePDF,redo)
%%% read in single unit spike times for a given block
%%% this is mostly just a matter of sorting the spikes from one block
%%% but nice to do it just in one line!

load(clustfile,'Block_Name','Tank_Name','ch');
blocknum = find(strcmp(block,Block_Name));
Block_Name = Block_Name{blocknum}

load(afile,'pinp');

if makePDF
    psfname = 'C:\temp.ps';
    if exist(psfname,'file')==2;delete(psfname);end
end




if ~exist('pinp','var') | length(pinp)<blocknum  | isempty(pinp(blocknum).psth) | redo
    
    spikes = getSpikes(clustfile,afile, block,redo);
    flags = struct('laserOn',1);
    tdtData= getTDTdata(Tank_Name, Block_Name, [], flags);
    laserT = tdtData.laserT;
    dV = diff(tdtData.laserTTL);
    
    edges = laserT(find(dV>2.5)+1);  %%% find rising edge
    edgeDiff=zeros(size(edges));
    edgeDiff(2:end)=diff(edges);
    edgeDiff(1)=100;
    
    trainEdges=edges(edgeDiff>1); %%% find edges that have greater than 1sec interval between - these will be beginning of a train
    
    dt=0.001;
    histbins = -0.05:dt:0.05;
    longdt=2;
    longBins=-trainPeriod:longdt:2*trainPeriod;
    psth=zeros(length(spikes.sp),length(hist(0,histbins)));
    longPsth = zeros(length(spikes.sp),length(hist(0,longBins)));
    
    mainfig=figure;
    nfig=0;
    for cell_n = 1:length(spikes.sp);
        times = spikes.sp{cell_n};
        
        histfig(cell_n)=figure;
        set(gcf,'position',[200 200 600 700]);
        subplot(3,2,1);
        plot(laserT,tdtData.laserTTL,'g');
        hold on
        plot(0:1:max(times),hist(times,0:1:max(times)));
        
        title(sprintf('unit %d',cell_n))
        %title(sprintf('ch %d ',cell_n))
        
        subplot(3,2,2)
        hold on
        for t = 1:length(edges);
            tdiff = times-edges(t);
            tdiff = tdiff(abs(tdiff)<max(histbins));
            plot(tdiff*1000,t*ones(size(tdiff)),'ks','MarkerSize',1)
            psth(cell_n,:) = psth(cell_n,:) + hist(tdiff,histbins);
        end
        xlim([min(histbins)*1000 max(histbins)*1000])
        
        
        psth(cell_n,:)=psth(cell_n,:)/(dt*length(edges));
        subplot(3,2,4);
        plot(histbins*1000,psth(cell_n,:));
        hold on
        plot([0 0],[0 max(psth(cell_n,:))],'g');
        xlabel('msec')
        ylabel('sp/sec')
        
        
        subplot(3,2,3)
        hold on
        for t = 1:length(trainEdges);
            tdiff = times-trainEdges(t);
            tdiff = tdiff(tdiff<max(longBins) & tdiff>min(longBins));
            plot([tdiff;tdiff],[t*ones(size(tdiff));t*ones(size(tdiff))-1],'k')
            longPsth(cell_n,:) = longPsth(cell_n,:) + hist(tdiff,longBins);
        end
        longPsth(cell_n,:)=longPsth(cell_n,:)/(longdt*length(trainEdges));
        xlim([min(longBins) max(longBins)])
        subplot(3,2,5);
        plot(longBins,longPsth(cell_n,:));
        hold on
        plot([0 0],[0 max(longPsth(cell_n,:))],'g');
        xlabel('sec')
        ylabel('sp/sec')
        xlim([min(longBins) max(longBins)])
        % ylim([0 7]); plot([0 17],[6 6],'g','Linewidth',4);
        
        if makePDF
            set(gcf, 'PaperPositionMode', 'auto');
            print('-dpsc',psfname,'-append');
        end
        
        if ~showImg
            close(gcf)
        end
        
        figure(mainfig);
        plot(histbins*1000,psth(:,:));hold on
        plot([0 0],[0 (max(psth(:)))],'g');
        xlabel('msec')
        ylabel('sp/sec')
        
        
        psth(isnan(psth))=0;
        pinp(blocknum).psth(cell_n,:) = psth;
        drawnow
        

    end
    
    
    
    %%% calculate baseline, evoked, and artifact by choosing windows
    baseline = mean(psth(:,5:45),2);
    baseStd = std(psth(:,5:45),[],2);
    artifact = psth(:,50:52);
    ev = max(psth(:,53:55),[],2);
    evoked = ev- baseline;
    zscore =evoked./baseStd;
    zscore(zscore>20)=20;
    junk = (psth(:,52)>50);
    
    pinp(blocknum).pinped = (zscore>10& evoked>20 & ~junk);
    
    save(afile,'pinp','-append')
    
    if makePDF
        ps2pdf('psfile', psfname, 'pdffile', [clustfile(1:end-11) 'PINP.pdf']);
    end
    
end
pinp = pinp(blocknum);
