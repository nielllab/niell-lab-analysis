close all
clear all
%plexiBatch
%contrastBatch
%SizeBatch

for i = 1:length(files)
%     if files(i).Tnum==1
    [latency{i} tracks{i} mouseTouch{i} cricketEdge{i} mouseTouchEd{i} MouseTouchL{i} Rstarts{i} hits{i} range{i} tdhist{i} cricketTouch{i}] = analyzePlexi([pathname files(i).trackpts],files(i).fps,28);
%     end
end

group = [files.group]; lighting = [files.lighting];earPlug=[files.EP];
grouplabels = {'','','','','','',''};
lightlabels = {'dark','light'};

close all
distbins = -25:3.5:25;
%distbinR=4:3:36;
accuracyFig = figure;
col = 'rgbcmk';
clear accuracy mouseAll cricketAll
pctAccurate=NaN(length(mouseTouch),2); %preallocate size of percent correct error
% for light on/off only
    for light = 1:2
        trials = find(lighting == light-1);
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
        
        
        accuracy(light,1:length(distbins))=hist(diff_trans,distbins)/length(diff);   
        mouseAll(light,:) = hist(mousePts,1.25:2.5:40);
        cricketAll(light,:) = hist(cricketPts,1.25:2.5:40);
        
        
         for i=1:length(trials)
             mouseErr(i)=abs(mouseTouch{trials(i)}(1)-cricketTouch{trials(i)}(1));
             accurate=find(mouseErr(i)<3);
             pctAccurate(i,light)=length(accurate)/length(mouseErr(i));      
         end
         
        figure(accuracyFig);
        subplot(4,2,2 + 3-light); plot(distbins,squeeze(accuracy(light,:))); axis([-30 30 0 1]);     
        
        figure; hold on
        for tr = 1:length(trials);
            for n = 1%1:length(mouseTouch{trials(tr)});
                   plot(tracks{trials(tr)}{n}(:,1),tracks{trials(tr)}{n}(:,2),'b');
%                    if cricketTouch{trials(tr)}(n)==1
%                        plot(tracks{trials(tr)}{n}(:,1),tracks{trials(tr)}{n}(:,2),'r','LineWidth',2);
%                    end
                   
            end
        end
        axis([0 40 -40 40]); title(sprintf('%s %s',lightlabels{light}));         
        
    end


% overlay of approachError by light (blue) vs. dark (black)
figure
        plot(distbins,squeeze(accuracy(1,:)),'k','LineWidth',2); hold on
        plot(distbins,squeeze(accuracy(2,:)),'b','LineWidth',2);    
%
%% plot shaded error bar plot for approachError paths light vs. dark
  
% pctAccurate=fliplr(pctAccurate)
% data_mu(1,1) = nanmean(pctAccurate(:,1));
% err(1,1) = nanstd(pctAccurate(:,1))/sqrt(sum(~isnan(pctAccurate(:,1))));
% data_mu(1,2) = nanmean(pctAccurate(:,2));
% err(1,2) = nanstd(pctAccurate(:,2))/sqrt(sum(~isnan(pctAccurate(:,1))));
nhit=sum(~isnan(pctAccurate));%n number of trials
phit=sum(pctAccurate==1)./sum(~isnan(pctAccurate));%prob you get an accurate approach
[m v]=binostat(nhit,phit);
mHit=m./nhit; vHit=v./nhit;
stdDev=sqrt(vHit);
SEMhit=stdDev./sqrt(nhit);

[tbl, chi2stat,pval]=crosstab(pctAccurate(:,1),pctAccurate(:,2));

figure
bar(fliplr(mHit)); hold on; errorbar(1:2,fliplr(mHit),fliplr(SEMhit),'o')
ylabel('Percent of approaches that were accurate')
set(gca,'Xtick',1:2); set(gca,'XtickLabel',{'light','dark'});

hold on
for ii=1:2
    tmp=(pctAccurate(:,ii)+(rand(size(pctAccurate(:,ii)))-0.5)*0.1); %temporarily store data in variable "tmp"
    x = repmat(ii,1,length(tmp)); %the x axis location
    x = (x+(rand(size(x))-0.5)*0.1-.2); %add a little random "jitter" to aid visibility

    plot(x,tmp,'ok')
end
hold off

%% light trials

    
% for light on w/wout Ear plug
EPlabels = {'w/EP','w/out EP'};
accuracyFig = figure;
col = 'rgbcmk';
clear accuracy mouseAll cricketAll
pctAccurate=NaN(length(mouseTouch),2); %preallocate size of percent correct error

    for cond = 1:2
        if cond==1     
        trials = find(lighting ==1 & earPlug==1);
        
        clear diff Left mousePts cricketPts diff_trans
        elseif cond==2
        trials=find(lighting ==1 & earPlug==0);
        end
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
        
        
         for i=1:length(trials)
             mouseErr(i)=abs(mouseTouch{trials(i)}(1)-cricketTouch{trials(i)}(1));
             accurate=find(mouseErr(i)<3);
             pctAccurate(i,cond)=length(accurate)/length(mouseErr(i));      
         end
         
        figure(accuracyFig);
        subplot(4,2,2 + 3-cond); plot(distbins,squeeze(accuracy(cond,:))); axis([-30 30 0 1]);     
        
        figure; hold on
        for tr = 1:length(trials);
            for n = 1%1:length(mouseTouch{trials(tr)});
                   plot(tracks{trials(tr)}{n}(:,1),tracks{trials(tr)}{n}(:,2),'b');
%                    if cricketTouch{trials(tr)}(n)==1
%                        plot(tracks{trials(tr)}{n}(:,1),tracks{trials(tr)}{n}(:,2),'r','LineWidth',2);
%                    end
                   
            end
        end
        axis([0 40 -40 40]); title(sprintf('%s %s',EPlabels{cond}));         
        
    end
    
    
   % overlay of approachError by light (blue) vs. dark (black)
figure
        plot(distbins,squeeze(accuracy(1,:)),'k','LineWidth',2); hold on
        plot(distbins,squeeze(accuracy(2,:)),'b','LineWidth',2);     

nhit=sum(~isnan(pctAccurate));%n number of trials
phit=sum(pctAccurate==1)./sum(~isnan(pctAccurate));%prob you get an accurate approach
[m v]=binostat(nhit,phit);
mHit=m./nhit; vHit=v./nhit;
stdDev=sqrt(vHit);
SEMhit=stdDev./sqrt(nhit);

[tbl, chi2stat,pval]=crosstab(pctAccurate(:,1),pctAccurate(:,2));

figure
bar(fliplr(mHit)); hold on; errorbar(1:2,fliplr(mHit),fliplr(SEMhit),'o')
ylabel('Percent of approaches that were accurate')
set(gca,'Xtick',1:2); set(gca,'XtickLabel',{'EP light','NoEP'});

hold on
for ii=1:2
    tmp=(pctAccurate(:,ii)+(rand(size(pctAccurate(:,ii)))-0.5)*0.1); %temporarily store data in variable "tmp"
    x = repmat(ii,1,length(tmp)); %the x axis location
    x = (x+(rand(size(x))-0.5)*0.1-.2); %add a little random "jitter" to aid visibility

    plot(x,tmp,'ok')
end
hold off

% overlay of approachError by light (blue) vs. dark (black)
figure
        plot(distbins,squeeze(accuracy(1,:)),'k','LineWidth',2); hold on
        plot(distbins,squeeze(accuracy(2,:)),'b','LineWidth',2);    
%% 
for cond=1:4
    
    clear tr
    if cond==1
trials= find(lighting == 1);
    elseif cond==2
trials=  find(lighting == 0);
    elseif cond==3
        trials=find(lighting==1 & earPlug==1);
    elseif cond==4
        trials=find(lighting==1 & earPlug==0);
    end
figure; hold on; 
axis ([0 40 -40 40]); 

for tr = 1:length(trials);
    if ~isempty(mouseTouch{trials(tr)})
            for n = 1%:length(mouseTouch{trials(tr)});                 
                       plot(tracks{trials(tr)}{n}(:,1),tracks{trials(tr)}{n}(:,2),'b','LineWidth',1); 
                        
            end
    end 
end

clear tracksG goodApp tracksG AppErr
goodApp=0;

for tr = 1:length(trials);
    %if ~isempty(mouseTouch{trials(tr)})
            for i = 1:1%length(mouseTouch{trials(tr)});
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
set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'light','dark','light w/EP', 'light w/out EP'});

hold on
for ii=1:1
    tmp=ErrApproach{ii}; %temporarily store data in variable "tmp"
    x = repmat(ii,1,length(tmp)); %the x axis location
    x = (x+(rand(size(x))-0.5)*0.1-.1); %add a little random "jitter" to aid visibility

    plot(x,tmp,'ok')
end
hold off

%% Plot mean error with std Error mean  for the two conditions

figure
      shadedErrorBar(binD,nanmean(avg_y_sym{1}),stdErr{1},'b'); hold on  
      plot(binD,nanmean(avg_y_sym{1}),'b','LineWidth',2); hold on
      
      shadedErrorBar(binD,nanmean(avg_y_sym{2}),stdErr{2},'k'); hold on
      plot(binD,nanmean(avg_y_sym{2}),'k','LineWidth',2); hold on
      
      shadedErrorBar(binD,nanmean(avg_y_sym{3}),stdErr{3},'g'); hold on
      plot(binD,nanmean(avg_y_sym{3}),'g','LineWidth',2); hold on
      
      

      title 'Mean Error Symetric' 
      set(gca,'xdir','reverse','XLim',[0 25],'YLim',[0 20])
      
      
 %% plot bar graph of small error likelihood
 
 