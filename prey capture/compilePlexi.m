close all
clear all
plexiBatch

for i = 1:length(files)
%     if files(i).Tnum==1
    [latency{i} tracks{i} mouseTouch{i} cricketEdge{i} mouseTouchEd{i} MouseTouchL{i} Rstarts{i} hits{i} range{i} tdhist{i} cricketTouch{i}] = analyzePlexi([pathname files(i).trackpts],files(i).fps,28);
%     end
end

group = [files.group]; lighting = [files.lighting];
grouplabels = {'1 side','2 side','','','','',''};
lightlabels = {'dark','light'};

distbins = -35:5:35;
%distbinR=4:3:36;
accuracyFig = figure;
col = 'rgbcmk';
clear accuracy mouseAll cricketAll
pctAccurate=NaN(length(mouseTouch),2); %preallocate size of percent correct error
%for m = 1:1 %  one side or 2-side
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
        
        %get std deviation of percent correct from each mouse
        
%          for i=1:length(trials)
%              mouseErr{i}=abs(mouseTouch{trials(i)}-cricketTouch{trials(i)});
%              accurate=find(mouseErr{i}<3);
%              pctAccurate(i,light)=length(accurate)/length(mouseErr{i});      
%          end
        
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
%end

% overlay of approachError by light (blue) vs. dark (black)
figure
        plot(distbins,squeeze(accuracy(1,:)),'k','LineWidth',2); hold on
        plot(distbins,squeeze(accuracy(2,:)),'b','LineWidth',2);    
% 
% %plot approachError by range 
%     % 2D-histogram
% for t=1:length(tdhist)
%     if ~isempty(tdhist{t})
% tdAll(:,:,t)=mean(tdhist{t},3);
%     end
%    
% end
% 
%  for light=1:2
%  trials = find(lighting == light-1);
%  figure
%  h(:,:,light)=imagesc((mean(tdAll(:,:,trials),3)));
%  end
 
%% plot shaded error bar plot for approachError paths light vs. dark
  
data(1,2) = nanmean(pctAccurate(:,1));
err(1,2) = nanstd(pctAccurate(:,1))/sqrt(sum(~isnan(pctAccurate(:,1))));
data(1,1) = nanmean(pctAccurate(:,2));
err(1,1) = nanstd(pctAccurate(:,2))/sqrt(sum(~isnan(pctAccurate(:,1))));

figure
bar(data(1,:)); hold on; errorbar(1:2,data(1,:),err(1,:),'o')
ylabel('Percent of approaches that were accurate')
set(gca,'Xtick',1:2); set(gca,'XtickLabel',{'light','dark'});
title('figure 4E')
%% light trials
trials= find(lighting == 1);
figure; hold on; 
axis ([0 40 -40 40]); 

for tr = 1:length(trials);
    if ~isempty(mouseTouch{trials(tr)})
            for n = 1:length(mouseTouch{trials(tr)});                 
                       plot(tracks{trials(tr)}{n}(:,1),tracks{trials(tr)}{n}(:,2),'b','LineWidth',1); 
                        
            end
    end 
end

clear tracksG goodApp tracksG
goodApp=0;

for tr = 1:length(trials);
    %if ~isempty(mouseTouch{trials(tr)})
            for i = 1:1%length(mouseTouch{trials(tr)});
                %if cricketEdge{trials(tr)}(i)==0;
                       goodApp=goodApp+1;    
                       tracksG{goodApp}=tracks{trials(tr)}{i};               
                  % end       
            end
    %end
end

% make arrays all the same size so that distance and error can be averaged
binD=1:2:34;
binE=-35:5:35; 

clear errHist avg_y

for i=1:length(tracksG)
 [errHist(:,:,i),avg_y(i,:)]=myHist2Avg(tracksG{i}(:,1),tracksG{i}(:,2),binD,binE);
end

avg_y_sym_L=abs(avg_y);
stdD_L=nanstd(avg_y_sym_L,[],1);
stdErr_L=stdD_L/sqrt(length(avg_y_sym_L));

% figure;
% imagesc(mean(errHist,3));

%plot tracks of approaches
figure
col=[0,0,1,0.6] %4th entry makes line transparent


for i=1:length(tracksG)
    
    plot(tracksG{i}(:,1),tracksG{i}(:,2),'color',col,'LineWidth',1.5);hold on
    %plot(tracksG{end}(:,1),tracksG{end}(:,2),'r');
end
set(gca,'xdir','reverse','XLim',[1 40],'YLim',[-35 35])

    
      
%% dark trials

trials= find(lighting == 0);

figure; hold on; axis ([0 40 -40 40]); 

for tr = 1:length(trials);
    %if ~isempty(mouseTouch{trials(tr)})
            for n = 1:length(mouseTouch{trials(tr)});
                
                  % if mouseTouchEd{trials(tr)}(n)==0
                       plot(tracks{trials(tr)}{n}(:,1),tracks{trials(tr)}{n}(:,2),'k','LineWidth',1);
                   %end                 
            end
   % end
end

goodApp=0; clear tracksG

for tr = 1:length(trials);
   % if ~isempty(mouseTouch{trials(tr)})
            for i = 1%:length(mouseTouch{trials(tr)});
                
                  % if mouseTouchEd{trials(tr)}(i)==0;
                      
                       goodApp=goodApp+1;
                       tracksG{goodApp}=tracks{trials(tr)}{i};               
                  % end
                   
            end
    %end
end

% binD=3:3:30;
% binE=-24:3:24; 
clear errHist avg_y

for i=1:length(tracksG)
 [errHist(:,:,i),avg_y(i,:)]=myHist2Avg(tracksG{i}(:,1),tracksG{i}(:,2),binD,binE);
end

avg_y_sym_D=abs(avg_y);
stdD_D=nanstd(avg_y_sym_D,[],1);
stdErr_D=stdD_D/sqrt(length(avg_y_sym_D));

figure;hold on
col=[0,0,0,0.6] %4th entry makes line transparent
figure;hold on
plot(tracksG{25}(:,1),tracksG{25}(:,2),'r','LineWidth',1.5);
for i=1:length(tracksG)
    
    plot(tracksG{i}(:,1),tracksG{i}(:,2),'color',col,'LineWidth',1.5);
       % plot(tracksG{end}(:,1),tracksG{end}(:,2),'r');

end

set(gca,'xdir','reverse','XLim',[1 40],'YLim',[-35 35])
title 'approachError in Dark'

% Plot mean error with std Error mean  for the two conditions

figure
      shadedErrorBar(binD,nanmean(avg_y_sym_L),stdErr_L,'b'); hold on  
      plot(binD,nanmean(avg_y_sym_L),'b','LineWidth',2); hold on
      
      shadedErrorBar(binD,nanmean(avg_y_sym_D),stdErr_D,'k'); hold on
      plot(binD,nanmean(avg_y_sym_D),'k','LineWidth',2); hold on

      title 'Mean Error Symetric' 
      set(gca,'xdir','reverse','XLim',[0 35],'YLim',[0 20])
      
      
 %% plot bar graph of small error likelihood
 
 