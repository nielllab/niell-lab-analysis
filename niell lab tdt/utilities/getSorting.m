function getSorting(clustfile,afile,label);
display('loading data')
load(clustfile,'wave_all','mean_wvform','idx_all','event_times_all','Block_Name');
load(afile,'cells')
Block_Name

for i = 1:length(cells)
    tet = ceil(cells(i,1)/4);
    tet_ch = cells(i,1);
    c = cells(i,2);
    figure
%     subplot(2,2,1);
%     plot(mean_wvform(:,tet_ch:tet_ch+3,c));
%    legend('1','2','3','4')
    
    subplot(2,2,2);
    t1= squeeze(event_times_all{tet_ch}(find(idx_all{tet_ch} == c)));
    length(t1);
    dt = diff(t1);
    dt =dt(dt<.02);
    hist(dt,.0005:0.001:.02);

    for b = 1:length(Block_Name);
        allSp = event_times_all{tet_ch};
        theseSp = event_times_all{tet_ch}(idx_all{tet_ch}==c);
        blockMax = max( allSp(allSp> (b-1)*10^5 & allSp < b*10^5) - (b-1)*10^5);
        theseSp = theseSp(theseSp> (b-1)*10^5 & theseSp < b*10^5) - (b-1)*10^5;
        rate(b) = length(theseSp)/blockMax;
        nm=Block_Name{b};
        nm(nm=='_') = ' ';
        Block_Name{b}=nm;
    end
    
    
    subplot(2,2,1);
    bar(rate); xlim([0.5 b+0.5]);
    set(gca,'Xtick',1:b); set(gca,'XtickLabel',Block_Name); set(gca,'Xticklabelrotation',90)
    title(sprintf('%s cell %d ch %d cl %d',label,i,tet_ch,c));   
    yl = get(gca,'Ylim'); ylim([0 max(yl(2),2)]);
    
    if exist('wave_all','var')
        wvall = wave_all{tet};
        wvclust = wvall(find(idx_all{tet_ch}==c),:,:);
        
        dt = diff(t1);
        breaks = find(dt>5*10^4);
        breaks(end+1)=length(dt);
        dt(dt>5*10^4)=1;
        tmerge = cumsum(dt);
        binwidth=60;
        amps =squeeze(min(wvclust(:,5:10,:),[],2));
        clear ampmean
        for t =1:floor(max(tmerge)/binwidth);
            ampmean(t,:) = median(amps(tmerge>(t-1)*binwidth & tmerge<t*binwidth,:),1);
        end
        
        subplot(2,2,3:4)
        if length(amps)~=0  %%% happens from merge
            plot([0 tmerge],amps,'.','MarkerSize',2 ); hold on
            plot((binwidth:binwidth:max(tmerge))-binwidth/2,ampmean,'LineWidth',2);
            for b = 1:length(breaks)
                plot([tmerge(breaks(b)) tmerge(breaks(b))],[-10^-4 2*10^-5],'Linewidth',4,'color','b')
                bl = Block_Name{b}; bl(bl=='_')=' ';
                text(tmerge(breaks(b))-120,2*10^-5,bl,'Rotation',90,'HorizontalAlignment','right','FontSize',14,'Color','b')
            end
        end
        
    end
end

