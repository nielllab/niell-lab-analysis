
% close all
% clear all
% 
 MuscimolBatch

for i =1:length(files)

   [latency{i} tracks{i} Time1stApp(i) IAppI(i) mouseTouch{i} cricketEdge{i} mouseTouchEd{i} MouseTouchL{i} Rstarts{i} hits{i} range{i} tdhist{i} cricketTouch{i}] = analyzePlexi([pathname files(i).trackpts],files(i).fps,28,1);
%     
%     fn = [pathname files(i).trackpts]
%     if ~isempty(files(i).trackpts)
%     [ location(:,:,i) targ{i} body{i} r{i} thetaR{i} thetaSac{i} vhead{i} vbody{i} abody{i} vtarg{i} atarg{i} head{i} sub{i}] = analyzeCapture(fn,files(i).fps,files(i).scale,files(i).subj);
%      dist(i) = mean(r{i});   
%     end
end

group = [files.group]; lighting = [files.lighting];%chnge to size for size

grouplabels = {'preControl/muscimol','postControl','postMuscimol','iDREADD','Dctl','V1PostM'};

%grouplabels = {'1xLight','2x','4x','.5x','.25x'};

%lightlabels = {'dark','light'};

close all
distbins = -25:3.5:25;
%distbinR=4:3:36;
accuracyFig = figure;
col = 'bkrgcmybrgcmy';
clear accuracy mouseAll cricketAll
pctAccurate=NaN(length(mouseTouch),length(grouplabels)); %preallocate size of percent correct error
avgStart=NaN(length(mouseTouch),length(grouplabels));
figure; hold on

    for cond = 1:6 %expand to fit num conditions you have
        if cond==1
            %trials=find(group==1 | group==2);
         trials=find(group==1 & lighting==1);

        elseif cond ==2
            trials = find(group ==2 );
        elseif cond ==3
            trials = find(group ==3 );
        elseif cond ==4
             trials = find(group ==4 );     
        elseif cond ==5
             trials = find(group ==5 );
             elseif cond ==6
             trials = find(group ==7 );
        end  
          
       
        %trials = find(lighting == cond-1);
        clear diff Left mousePts cricketPts diff_trans mouseErr a accurate StartingD
        
        mousePts=[];cricketPts=[];Left=[];
        for j=1:length(trials)
        x{j} = [mouseTouch{trials(j)}]; y{j}=[cricketTouch{trials(j)}];z{j}=[MouseTouchL{trials(j)}];DistStart{j}=[Rstarts{trials(j)}];
        mousePts=[x{:}];
        cricketPts=[y{:}];
        Left=[z{:}];
        end
       
        % - value equals cricket is on the mouse's left if the cricket is
        % on the left
       
        %histogram differences between mouse head and cricket when mouse
        %touches the glass
        
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
        numApp(cond)=length(mousePts);
       
        
         for i=1:length(trials)
             mouseErr{i}=abs(mouseTouch{trials(i)}-cricketTouch{trials(i)});
             a=[mouseErr{i}(:)];
             accurate=a<3;
             pctAccurate(i,cond)=sum(accurate)/length(a);      
         end
         
         for i=1:length(trials)
             s=Rstarts{trials(i)};   
             avgStart(i,cond)=nanmean(s);      
         end
         
        figure(accuracyFig);
        subplot(6,1,cond); plot(distbins,squeeze(accuracy(cond,:)));% axis([-30 30 0 1]);     
        
        figure; hold on
        for tr = 1:length(trials);
            for n = 1:length(mouseTouch{trials(tr)});
                   plot(tracks{trials(tr)}{n}(:,1),tracks{trials(tr)}{n}(:,2),col(cond));
%                    if cricketTouch{trials(tr)}(n)==1
%                        plot(tracks{trials(tr)}{n}(:,1),tracks{trials(tr)}{n}(:,2),'r','LineWidth',2);
%                    end
                   
            end
        end
        %axis([0 35 -30 30]); title(sprintf('%s %s',grouplabels{cond}));         
        
    end

%col = 'kbrgcmy';
% overlay of approachError by light (blue) vs. dark (black)
%% virtual versus real and virtual exp vs. virtual naive
figure
        plot(distbins,squeeze(accuracy(2,:)),'b','LineWidth',2); hold on
        plot(distbins,squeeze(accuracy(6,:)),'LineWidth',2,'Color',[0 0 0]+0.05*10); hold on
        title 'preControl vs. postV1Muscimol'
        
figure
        plot(distbins,squeeze(accuracy(1,:)),'b','LineWidth',2); hold on
        plot(distbins,squeeze(accuracy(2,:)),'m','LineWidth',2); hold on
        plot(distbins,squeeze(accuracy(3,:)),'c','LineWidth',2); hold on
        plot(distbins,squeeze(accuracy(4,:)),'g','LineWidth',2); hold on
        plot(distbins,squeeze(accuracy(5,:)),'y','LineWidth',2); hold on
        title 'postDiI vs. PostV1'
        

figure
        plot(distbins,squeeze(accuracy(5,:)),'b','LineWidth',2); hold on
        plot(distbins,squeeze(accuracy(4,:)),'LineWidth',2,'Color',[0 0 0]+0.05*10); hold on
        title 'ctl vs. iDREADD'
       
%
%% plot shaded error bar plot for approachError paths light vs. dark
  
%this is the average start distance of the tacks, but we want avg of
%inflection points based on head angle


% data_mu(1,1) = nanmean(avgStart(:,1));
% err(1,1) = nanstd(avgStart(:,1)/sqrt(sum(~isnan(avgStart(:,1)))));
% data_mu(1,2) = nanmean(avgStart(:,2));
% err(1,2) = nanstd(avgStart(:,2)/sqrt(sum(~isnan(avgStart(:,2)))));
% data_mu(1,3) = nanmean(avgStart(:,3));
% err(1,3) = nanstd(avgStart(:,3)/sqrt(sum(~isnan(avgStart(:,3)))));
% data_mu(1,4) = nanmean(avgStart(:,4));
% err(1,4) = nanstd(avgStart(:,4)/sqrt(sum(~isnan(avgStart(:,4)))));
% data_mu(1,5) = nanmean(avgStart(:,5));
% err(1,5) = nanstd(avgStart(:,5)/sqrt(sum(~isnan(avgStart(:,5)))));

nhit=sum(~isnan(pctAccurate));%n number of trials
phit=sum(pctAccurate==1)./nhit;%prob you get an accurate approach
[m v]=binostat(nhit,phit);
mHit=m./nhit; vHit=v./nhit;
stdDev=sqrt(vHit);
SEMhit=stdDev./sqrt(nhit);

[tbl, chi2stat,pval]=crosstab(pctAccurate(:,1),pctAccurate(:,2));

figure
bar(mHit(1,[2,6])); hold on; errorbar(1:2,mHit(1,[2,6]),SEMhit(1,[2,6]),'o')
ylabel('Percent of approaches that were accurate')
set(gca,'Xtick',1:2); set(gca,'XtickLabel',{'PostControl','PostMusc'});

hold on
for ii=1:2
    tmp=(pctAccurate(:,ii+1)+(rand(size(pctAccurate(:,ii+1)))+0.01)*0.1); %temporarily store data in variable "tmp"
    x = repmat(ii,1,length(tmp)); %the x axis location
    x = (x+(rand(size(x))-0.5)*0.1-.2); %add a little random "jitter" to aid visibility

    plot(x,tmp,'ok')
end
hold off


for cond=1:6
    
    clear tr
    if cond ==1
        trials = find(group==1 & lighting ==1);  %%%  scaled conditions
        elseif cond==2
            trials=find(group==2);
        elseif cond ==3
            trials = find(group ==3 );
            elseif cond ==4
            trials = find(group ==4 );
       elseif cond ==5
            trials = find(group ==5 );
            elseif cond ==6
            trials = find(group ==7 );
    end  

figure; hold on; 
axis ([0 40 -40 40]); 

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
   
%group data for time to 1st approach and Inter-approach Interval

FirstAppT(cond,1)=nanmean(Time1stApp(trials));
FirstAppT(cond,2)=nanstd(Time1stApp(trials))/sqrt(length(Time1stApp(trials)));

IAI(cond,1)=nanmean(IAppI(trials));
IAI(cond,2)=nanstd(IAppI(trials))/sqrt(length(IAppI(trials)));


end

%% plot latency to first approach and inter-approach Interval
figure
bar(FirstAppT([4,6],1)); hold on; errorbar(1:2,FirstAppT([4,6],1),FirstAppT([4,6],2),'o')
ylabel('Time to First Approach')
set(gca,'Xtick',1:2); set(gca,'XtickLabel',{'PostCtl','PostMusc'});

figure
bar(FirstAppT(4:5,1)); hold on; errorbar(1:2,FirstAppT(4:5,1),FirstAppT(4:5,2),'o')
ylabel('Time to First Approach')
set(gca,'Xtick',1:2); set(gca,'XtickLabel',{'iDRDp+CNO','iDRDn+CNO'});

figure
bar(IAI([4,6],1)); hold on; errorbar(1:2,IAI([4,6],1),IAI([4,6],2),'o')
ylabel('Inter-Approach interval')
set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'PostCtl','PostMusc'});

figure
bar(IAI([2,6],1)); hold on; errorbar(1:2,IAI([2,6],1),IAI([2,6],2),'o')
ylabel('Inter-Approach interval')
set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'iDRDp+CNO','iDRDn+CNO'});

%%

figure
bar(ErrApproachMu(1,2:6)); hold on; errorbar(1:5,ErrApproachMu(1,2:6),ErrApproachStdEM(1,2:6),'o')
ylabel('Pre Contact prey azimuth')
set(gca,'Xtick',1:6); set(gca,'XtickLabel',{'PostCtl','PostMusc','iDRDp+CNO','iDRDn+CNO','V1Muscimol'});

hold on
for ii=1:5
    tmp=ErrApproach{ii+1}; %temporarily store data in variable "tmp"
    x = repmat(ii,1,length(tmp)); %the x axis location
    x = (x+(rand(size(x))-0.5)*0.1-.1); %add a little random "jitter" to aid visibility

    plot(x,tmp,'ok')
end
hold off

% Plot mean error with std Error mean  for the two conditions

figure

      shadedErrorBar(binD,nanmean(avg_y_sym{2}),stdErr{2},'k'); hold on
      plot(binD,nanmean(avg_y_sym{2}),'k','LineWidth',2); hold on 
      title 'Mean Error Symetric' 
      set(gca,'xdir','reverse','XLim',[2 30],'YLim',[0 20])

      shadedErrorBar(binD,nanmean(avg_y_sym{6}),stdErr{6},'r'); hold on
      plot(binD,nanmean(avg_y_sym{6}),'r','LineWidth',2); hold on
      title 'Mean Error Symetric' 
      set(gca,'xdir','reverse','XLim',[2 30],'YLim',[0 20])
      
figure
      shadedErrorBar(binD,nanmean(avg_y_sym{4}),stdErr{4},'g'); hold on
      plot(binD,nanmean(avg_y_sym{4}),'g','LineWidth',2); hold on 
      title 'Mean Error Symetric' 
      set(gca,'xdir','reverse','XLim',[2 30],'YLim',[0 20])
      
      shadedErrorBar(binD,nanmean(avg_y_sym{6}),stdErr{6},'c'); hold on
      plot(binD,nanmean(avg_y_sym{6}),'c','LineWidth',2); hold on 
%       shadedErrorBar(binD,nanmean(avg_y_sym{1}),stdErr{1},'b'); hold on  
%       plot(binD,nanmean(avg_y_sym{1}),'b','LineWidth',2); hold on
                 
      title 'Mean Error Symetric' 
      set(gca,'xdir','reverse','XLim',[3 20],'YLim',[0 18])
      
figure

      shadedErrorBar(binD,nanmean(avg_y_sym{2}),stdErr{2},'k'); hold on
      plot(binD,nanmean(avg_y_sym{2}),'k','LineWidth',2); hold on
      
      shadedErrorBar(binD,nanmean(avg_y_sym{3}),stdErr{3},'r'); hold on
      plot(binD,nanmean(avg_y_sym{3}),'r','LineWidth',2); hold on  
      
      title 'Dyd injected vs. Muscimol injectedt' 
      set(gca,'xdir','reverse','XLim',[4 25],'YLim',[0 18])
      
      
      figure

      shadedErrorBar(binD,nanmean(avg_y_sym{5}),stdErr{5},'k'); hold on
      plot(binD,nanmean(avg_y_sym{5}),'k','LineWidth',2); hold on
      
      shadedErrorBar(binD,nanmean(avg_y_sym{4}),stdErr{4},'g'); hold on
      plot(binD,nanmean(avg_y_sym{4}),'g','LineWidth',2); hold on  
      
      title 'Mean Error Symetric CNO+/DREADD- vs. CNO+/DREADD+' 
      set(gca,'xdir','reverse','XLim',[4 25],'YLim',[0 18])
 
      
      %%
   
    
      
