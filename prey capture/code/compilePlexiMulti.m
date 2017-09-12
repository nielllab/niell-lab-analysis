
close all
clear all
%contrastBatch;% change to size for size data 
SizeBatch
%MuscimolBatch
%plexiBatch
%V1MuscimolBatch % X:\Preycapture\Muscimol V1 Saved Tracks

for i = 1:length(files)

    [latency{i} tracks{i} Time1stApp(i) IAppI(i) mouseTouch{i} cricketEdge{i} mouseTouchEd{i} MouseTouchL{i} Rstarts{i} hits{i} range{i} tdhist{i} cricketTouch{i}] = analyzePlexi([pathname files(i).trackpts],files(i).fps,28,1);

end

group = [files.group]; lighting = [files.lighting];%chnge to size for size
%contrast
%grouplabels = {'100','0','50','25','12.5','6.25'};
%size
grouplabels = {'1x','1xdark','2x','3x','.5x','.25x','.125x','n1x','n2x','n3x','n.5x','n.25x','n.125x'};

%lightlabels = {'dark','light'};

close all
distbins = -25:3.5:25;
%distbinR=4:3:36;
accuracyFig = figure;
col = 'bkrgcmybrgcmy';
clear accuracy mouseAll cricketAll
pctAccurate=NaN(length(mouseTouch),length(grouplabels)); %preallocate size of percent correct error

figure; hold on

    for cond = 1:7 %expand to fit num conditions you have
        if cond==1
            trials=find(group==1 & lighting==1); %size is 1x
%         if cond ==1
%             trials = find(group==1 & lighting==1);  %%%  scaled conditions
% %         elseif cond==2
% %             trials=find(group==1 & lighting==0);
%         elseif cond ==2
%             trials = find(group ==2 );
%         elseif cond ==4
%             trials = find(group ==3 );
        elseif cond ==2
            trials = find(group ==2 & lighting==1);%2x size
        elseif cond ==3
            trials = find(group ==3 & lighting==1 );%4x size
        elseif cond ==4
            trials = find(group ==4 & lighting==1 );%.5x
        elseif cond ==5
            trials = find(group ==5 & lighting==1 );%.25x
        elseif cond ==6
            trials = find(group ==7 & lighting==1 ); %naive 1x
        elseif cond ==7
            trials = find(group ==10 & lighting==1 );%naive .5x
%         elseif cond ==11
%             trials = find(group ==10 );
%         elseif cond ==12
%             trials = find(group ==11 );
%         elseif cond ==13
%             trials = find(group ==12 );
        end  
          
       
        %trials = find(lighting == cond-1);
        clear diff Left mousePts cricketPts diff_trans

        for j=1:length(trials)
        mousePts(j) = mouseTouch{trials(j)}(1); cricketPts(j) = cricketTouch{trials(j)}(1);Left(j)=MouseTouchL{trials(j)}(1);
        end
        
        diff=mousePts-cricketPts;
        
        for i=1:length(diff)
            if Left(i)==1
            diff_trans(i)=-diff(i)
            else
            diff_trans(i)=diff(i)
            end
        end
        
        
        accuracy(cond,1:length(distbins))=hist(diff_trans,distbins)/length(diff);   
        mouseAll(cond,:) = hist(mousePts,1.25:2.5:40);
        cricketAll(cond,:) = hist(cricketPts,1.25:2.5:40);
        numApp(cond)=length(trials);
        
         for i=1:length(trials)
             mouseErr(i)=abs(mouseTouch{trials(i)}(1)-cricketTouch{trials(i)}(1));
             accurate=find(mouseErr(i)<3);
             pctAccurate(i,cond)=length(accurate)/length(mouseErr(i));      
         end
         
%         figure(accuracyFig);
%         subplot(4,2,2 + 7-cond); plot(distbins,squeeze(accuracy(cond,:))); axis([-30 30 0 1]);     
%         
       % figure; hold on
        for tr = 1:length(trials);
            for n = 1%1:length(mouseTouch{trials(tr)});
                   plot(tracks{trials(tr)}{n}(:,1),tracks{trials(tr)}{n}(:,2),col(cond));
%                    if cricketTouch{trials(tr)}(n)==1
%                        plot(tracks{trials(tr)}{n}(:,1),tracks{trials(tr)}{n}(:,2),'r','LineWidth',2);
%                    end
                   
            end
        end
        axis([0 40 -40 40]); title(sprintf('%s %s',grouplabels{cond}));         
        
    end

%col = 'kbrgcmy';
% overlay of approachError by light (blue) vs. dark (black)
%% virtual versus real and virtual exp vs. virtual naive

%plot accuracy as size increase from standard

figure
        plot(distbins,squeeze(accuracy(1,:)),'LineWidth',2,'Color',[0 0 0]+0.05*10); hold on
        plot(distbins,squeeze(accuracy(3,:)),'y','LineWidth',2); hold on
        plot(distbins,squeeze(accuracy(2,:)),'r','LineWidth',2); hold on
%plot accuracy as size decrease from standard
figure
       
        plot(distbins,squeeze(accuracy(1,:)),'LineWidth',2,'Color',[0 0 0]+0.05*10); hold on
        plot(distbins,squeeze(accuracy(4,:)),'g','LineWidth',2); hold on
         plot(distbins,squeeze(accuracy(5,:)),'c','LineWidth',2); hold on
        title ''
%plot accuracy as size increase experience vs. naive
figure

        plot(distbins,squeeze(accuracy(1,:)),'LineWidth',2,'Color',[0 0 0]+0.05*10); hold on
         plot(distbins,squeeze(accuracy(4,:)),'g','LineWidth',2); hold on
        title ''
        
 figure

        plot(distbins,squeeze(accuracy(1,:)),'LineWidth',2,'Color',[0 0 0]+0.05*10); hold on
         plot(distbins,squeeze(accuracy(4,:)),'g','LineWidth',2); hold on
        title ''
        

        figure

        plot(distbins,squeeze(accuracy(6,:)),'LineWidth',2,'Color',[0 0 0]+0.05*10); hold on
         plot(distbins,squeeze(accuracy(7,:)),'g','LineWidth',2); hold on
        title ''
        
figure
        plot(distbins,squeeze(accuracy(3,:)),'b','LineWidth',2); hold on
        plot(distbins,squeeze(accuracy(9,:)),'LineWidth',2,'Color',[0 0 0]+0.05*10); hold on
        title '2x virtual Exp vs. Virtual Naive'
        
figure
        plot(distbins,squeeze(accuracy(4,:)),'b','LineWidth',2); hold on
        plot(distbins,squeeze(accuracy(10,:)),'LineWidth',2,'Color',[0 0 0]+0.05*10); hold on
        title '4x virtual Exp vs. Virtual Naive'
        
figure
        plot(distbins,squeeze(accuracy(5,:)),'b','LineWidth',2); hold on
        plot(distbins,squeeze(accuracy(11,:)),'LineWidth',2,'Color',[0 0 0]+0.05*10); hold on
        title '.5x virtual Exp vs. Virtual Naive'
figure      
        plot(distbins,squeeze(accuracy(1,:)),'b','LineWidth',2); hold on
        plot(distbins,squeeze(accuracy(8,:)),'LineWidth',2,'Color',[0 0 0]+0.05*10); hold on
        title '1x'
figure      
        plot(distbins,squeeze(accuracy(3,:)),'r','LineWidth',2); hold on
        plot(distbins,squeeze(accuracy(9,:)),'LineWidth',2,'Color',[0 0 0]+0.05*10); hold on
        title '2x'
        
figure      
        plot(distbins,squeeze(accuracy(4,:)),'g','LineWidth',2); hold on
        plot(distbins,squeeze(accuracy(10,:)),'LineWidth',2,'Color',[0 0 0]+0.05*10); hold on
        title '4x'
        
figure      
        plot(distbins,squeeze(accuracy(5,:)),'c','LineWidth',2); hold on
        plot(distbins,squeeze(accuracy(11,:)),'LineWidth',2,'Color',[0 0 0]+0.05*10); hold on
        title '.5x'
        
 figure      
        plot(distbins,squeeze(accuracy(6,:)),'m','LineWidth',2); hold on
        plot(distbins,squeeze(accuracy(12,:)),'LineWidth',2,'Color',[0 0 0]+0.05*10); hold on
        title '.25x'
figure
        plot(distbins,squeeze(accuracy(7,:)),'y','LineWidth',2);hold on
        plot(distbins,squeeze(accuracy(13,:)),'LineWidth',2,'Color',[0 0 0]+0.05*10); hold on
        title '.125x'

        
%
%% plot shaded error bar plot for approachError paths light vs. dark
  
% pctAccurate=fliplr(pctAccurate)
% data_mu(1,1) = nanmean(pctAccurate(:,1));
% err(1,1) = nanstd(pctAccurate(:,1))/sqrt(sum(~isnan(pctAccurate(:,1))));
% data_mu(1,2) = nanmean(pctAccurate(:,2));
% err(1,2) = nanstd(pctAccurate(:,2))/sqrt(sum(~isnan(pctAccurate(:,1))));
nhit=sum(~isnan(pctAccurate));%n number of trials
phit=sum(pctAccurate==1)./nhit;%prob you get an accurate approach
[m v]=binostat(nhit,phit);
mHit=m./nhit; vHit=v./nhit;
stdDev=sqrt(vHit);
SEMhit=stdDev./sqrt(nhit);

[tbl, chi2stat,pval]=crosstab(pctAccurate(:,1),pctAccurate(:,2));

figure
bar(mHit); hold on; errorbar(1:7,mHit,SEMhit,'o')
ylabel('Percent of approaches that were accurate')
set(gca,'Xtick',1:6); set(gca,'XtickLabel',{'100','0','50','25','12.5','6.25'});

hold on
for ii=1:6
    tmp=(pctAccurate(:,ii)+(rand(size(pctAccurate(:,ii)))-0.5)*0.1); %temporarily store data in variable "tmp"
    x = repmat(ii,1,length(tmp)); %the x axis location
    x = (x+(rand(size(x))-0.5)*0.1-.2); %add a little random "jitter" to aid visibility

    plot(x,tmp,'ok')
end
hold off

%% light trials
for cond=1:4
    
    clear tr
    if cond ==1
        trials = find(group==1 & lighting==1);  %%%  scaled conditions
%         elseif cond==2
%             trials=find(group==1 & lighting==0);
        elseif cond ==2
            trials = find(group ==2 );
        elseif cond ==3
            trials = find(group ==3 );
        elseif cond ==4
            trials = find(group ==4 );
%         elseif cond ==6
%             trials = find(group ==5 );
%         elseif cond ==7
%             trials = find(group ==6 );
%         elseif cond ==8
%             trials = find(group ==7 );
%         elseif cond ==9
%             trials = find(group ==8 );
%         elseif cond ==10
%             trials = find(group ==9 );
%         elseif cond ==11
%             trials = find(group ==10 );
%         elseif cond ==12
%             trials = find(group ==11 );
%         elseif cond ==13
%             trials = find(group ==12 );
        end  

figure; hold on; 
%axis ([0 40 -40 40]); 

for tr = 1:length(trials);
    if ~isempty(mouseTouch{trials(tr)})
            for n = 1:length(mouseTouch{trials(tr)});                 
                       plot(tracks{trials(tr)}{n}(:,1),tracks{trials(tr)}{n}(:,2),'b','LineWidth',1); 
                        
            end
    end 
end

clear tracksG goodApp tracksG AppErr
goodApp=0;

for tr = 1:length(trials);
    %if ~isempty(mouseTouch{trials(tr)})
            for i = 1:length(mouseTouch{trials(tr)});
                %if cricketEdge{trials(tr)}(i)==0;
                       goodApp=goodApp+1;    
                       tracksG{goodApp}=tracks{trials(tr)}{i};    
                       AppErr(goodApp)=abs(tracks{trials(tr)}{i}(end,2));
                  % end       
            end
    %end
end

ErrApproach{cond}=AppErr;
ErrApproachMu(cond)=nanmean(AppErr);
ErrApproachStdEM(cond)=nanstd(AppErr)/sqrt(length(AppErr));
% make arrays all the same size so that distance and error can be averaged
binD=1:2:34;
binE=-35:5:35; 

clear errHist avg_y

for i=1:length(tracksG)
 [errHist(:,:,i),avg_y(i,:)]=myHist2Avg(tracksG{i}(:,1),tracksG{i}(:,2),binD,binE);
end

avg_y_sym{cond}=abs(avg_y);
stdD{cond}=nanstd(avg_y_sym{cond},[],1);
stdErr{cond}=stdD{cond}/sqrt(length(avg_y_sym{cond}));

% figure;
% imagesc(mean(errHist,3));

%plot tracks of approaches
figure
col=[0,0,1,0.6]; %4th entry makes line transparent


for i=1:length(tracksG)
    
    plot(tracksG{i}(:,1),tracksG{i}(:,2),'color',col,'LineWidth',1.5);hold on
    %plot(tracksG{end}(:,1),tracksG{end}(:,2),'r');
end
set(gca,'xdir','reverse','XLim',[1 40],'YLim',[-35 35])
    
end
      
%% 
figure
bar(ErrApproachMu(1,1:4)); hold on; errorbar(1:4,ErrApproachMu(1,1:4),ErrApproachStdEM(1,1:4),'o')
ylabel('Pre Contact prey azimuth')
set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'light','dark'});

hold on
for ii=1:4
    tmp=ErrApproach{ii}; %temporarily store data in variable "tmp"
    x = repmat(ii,1,length(tmp)); %the x axis location
    x = (x+(rand(size(x))-0.5)*0.1-.1); %add a little random "jitter" to aid visibility

    plot(x,tmp,'ok')
end
hold off

% Plot mean error with std Error mean  for the two conditions

figure
      shadedErrorBar(binD,nanmean(avg_y_sym{4}),stdErr{4},'b'); hold on  
      plot(binD,nanmean(avg_y_sym{4}),'b','LineWidth',2); hold on
      
      shadedErrorBar(binD,nanmean(avg_y_sym{2}),stdErr{2},'k'); hold on
      plot(binD,nanmean(avg_y_sym{2}),'k','LineWidth',2); hold on
      
      shadedErrorBar(binD,nanmean(avg_y_sym{3}),stdErr{3},'r'); hold on
      plot(binD,nanmean(avg_y_sym{3}),'k','LineWidth',2); hold on
        
      shadedErrorBar(binD,nanmean(avg_y_sym{4}),stdErr{4},'r'); hold on
      plot(binD,nanmean(avg_y_sym{4}),'k','LineWidth',2); hold on

      shadedErrorBar(binD,nanmean(avg_y_sym{5}),stdErr{5},'g'); hold on
      plot(binD,nanmean(avg_y_sym{5}),'k','LineWidth',2); hold on
      
      shadedErrorBar(binD,nanmean(avg_y_sym{6}),stdErr{6},'m'); hold on
      plot(binD,nanmean(avg_y_sym{6}),'k','LineWidth',2); hold on
      
      shadedErrorBar(binD,nanmean(avg_y_sym{7}),stdErr{7},'y'); hold on
      plot(binD,nanmean(avg_y_sym{7}),'k','LineWidth',2); hold on
     
      title 'Mean Error Symetric' 
      set(gca,'xdir','reverse','XLim',[4 25],'YLim',[0 18])
      

 
      %% virtual experienced versus naive mice
      
      figure
      shadedErrorBar(binD,nanmean(avg_y_sym{1}),stdErr{1},'b'); hold on  
      plot(binD,nanmean(avg_y_sym{1}),'b','LineWidth',2); hold on
     
      shadedErrorBar(binD,nanmean(avg_y_sym{8}),stdErr{8},'k'); hold on  
      plot(binD,nanmean(avg_y_sym{8}),'k','LineWidth',2); hold on
      
      title '1x Mean Error Symetric Experienced v Naive' 
      set(gca,'xdir','reverse','XLim',[4 35],'YLim',[0 18])
      
      figure
      shadedErrorBar(binD,nanmean(avg_y_sym{9}),stdErr{9},'k'); hold on  
      plot(binD,nanmean(avg_y_sym{9}),'k','LineWidth',2); hold on
      
      shadedErrorBar(binD,nanmean(avg_y_sym{3}),stdErr{3},'b'); hold on  
      plot(binD,nanmean(avg_y_sym{3}),'b','LineWidth',2); hold on
      
      title '2x Mean Error Symetric Experienced v Naive' 
      set(gca,'xdir','reverse','XLim',[4 35],'YLim',[0 18])
      
      figure
      shadedErrorBar(binD,nanmean(avg_y_sym{10}),stdErr{10},'k'); hold on  
      plot(binD,nanmean(avg_y_sym{10}),'k','LineWidth',2); hold on
      
      shadedErrorBar(binD,nanmean(avg_y_sym{4}),stdErr{4},'b'); hold on  
      plot(binD,nanmean(avg_y_sym{4}),'b','LineWidth',2); hold on
      
      title '4x Mean Error Symetric Experienced v Naive' 
      set(gca,'xdir','reverse','XLim',[4 35],'YLim',[0 18])
      
      figure
      shadedErrorBar(binD,nanmean(avg_y_sym{11}),stdErr{11},'k'); hold on  
      plot(binD,nanmean(avg_y_sym{11}),'k','LineWidth',2); hold on
      
      shadedErrorBar(binD,nanmean(avg_y_sym{5}),stdErr{5},'b'); hold on  
      plot(binD,nanmean(avg_y_sym{5}),'b','LineWidth',2); hold on
      
      title '.5x Mean Error Symetric Experienced v Naive' 
      set(gca,'xdir','reverse','XLim',[4 35],'YLim',[0 18])
      
       figure
      shadedErrorBar(binD,nanmean(avg_y_sym{12}),stdErr{12},'k'); hold on  
      plot(binD,nanmean(avg_y_sym{12}),'k','LineWidth',2); hold on
      
      shadedErrorBar(binD,nanmean(avg_y_sym{6}),stdErr{6},'b'); hold on  
      plot(binD,nanmean(avg_y_sym{6}),'b','LineWidth',2); hold on
      
      title '.25x Mean Error Symetric Experienced v Naive' 
      set(gca,'xdir','reverse','XLim',[4 35],'YLim',[0 18])