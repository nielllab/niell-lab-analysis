close all;

savePDF = 1;
if savePDF
    psfilename = 'D:\tempPS.ps';
    if exist(psfilename,'file')==2;delete(psfilename);end 
end

dataAll = goodAll==1 & hasDrift==1;
useN =dataAll;
for c= 1:length(useN)
   % for prepost =1:2
  
        [respmax oind] = max(mean(drift_orient(c ,:,1,:),4)); %stationary
       % [g h]=(max(respmax(:)));
       prefOri(c)=oind;
   % end
end

useN =dataAll;
for c= 1:length(useN)
   % for prepost =1:2
        [respmax oind] = max(mean(drift_sf(c ,:,1,:),4)); %stationary
       % [g h]=(max(respmax(:)));
        prefSF(c)=oind;
    %end
end


%%% drift_cond_tcourse(cells,move,prepost,orient,sf,time);
for c=1:length(useN)
    tcourse_pref(c,:,:,:) = squeeze(nanmean(drift_cond_tcourse(c,:,:,prefOri(c),:,:),5));
    %tcourse_pref(c,:,:,:) = squeeze(mean(nanmean(drift_cond_tcourse(c,:,:,mod(prefOri(c)+([0 6])-1,12)+1,:,:),5),4));
    %tcourse_pref(c,:,:,:) = squeeze(mean(nanmean(drift_cond_tcourse(c,:,:,mod(prefOri(c)+([-1:1 5:7])-1,12)+1,:,:),5),4));
    
    tcourse_orth(c,:,:,:) = squeeze(mean(nanmean(drift_cond_tcourse(c,:,:,mod(prefOri(c)+3-1,12)+1,:,:),5),4));
    %tcourse_orth(c,:,:,:) = squeeze(mean(nanmean(drift_cond_tcourse(c,:,:,mod(prefOri(c)+([-3 3])-1,12)+1,:,:),5),4));
   % tcourse_orth(c,:,:,:) = squeeze(mean(nanmean(drift_cond_tcourse(c,:,:,mod(prefOri(c)+([-4:2 2:4])-1,12)+1,:,:),5),4));
    
    tcourse_all(c,:,:,:) = squeeze(mean(nanmean(drift_cond_tcourse(c,:,:,:,:,:),5),4));
end


% for c=1:length(useN)
%     tcourse_pref(c,:,:,:) = squeeze(nanmean(drift_cond_tcourse(c,:,:,:,prefSF(c),:),4));
%  %   tcourse_all(c,:,:,:) = squeeze(mean(nanmean(drift_cond_tcourse(c,:,:,:,:,:),5),4));
% end
%%
peakresp = squeeze(max(drift_orient,[],2))%-squeeze(drift_spont);%subtract driftspont
useResp = peakresp(:,1,1)>1 | peakresp(:,1,2)>1;
use = goodAll & ~inhAll & hasDrift==1 & useResp' ;%which cells to include
useInh = goodAll & inhAll==1 & hasDrift==1 & useResp';
layers = {'Layer 2/3','Layer 4','Layer 5','Inh. Units'};

dt=.1
titles={'Saline pref', 'DOI pref'};
for t=1:2
   figure('units','normalized','outerposition',[.5 .5 .75 .75])

       if t==1, set(gcf,'Name','Saline tcoursepreferred'); else set(gcf,'Name','DOI tcourse preferred')
    end
for i = 1:4
    clear mnpre mnpost prenorm postnorm normrange prerange
    if i==1
        data = find(use & layerAll==5 & treatment==t);
     %   set(gcf,'Name','grating lyr5 stat');
    elseif i==2
        clear data
      %  set(gcf,'Name','grating lyr 2/3 stat');
        data = find(use & (layerAll==2 |layerAll==3) & treatment==t);
   elseif i==3
        clear data
       % set(gcf,'Name','grating lyr 4 stat');
        data = find(use & (layerAll==4) & treatment==t);    
    else i==4
        clear data
      %  set(gcf,'Name','grating inh stat');
        data = find(useInh  & treatment==t);
end
         for c = 1:length(data)
            %normalize response of all cells used, then average
            prerange = squeeze(tcourse_pref(data(c),1,1,:));
            normrange(c,:)= max(prerange) %- min(tcourse_pref(data(c),1,1,:),[],1); %set range w/pre data to normalize to
            prenorm(c,:) =(tcourse_pref(data(c),1,1,:)- min(tcourse_pref(data(c),1,1,:)))...
                ./normrange(c,1);
            postnorm(c,:) = (tcourse_pref(data(c),1,2,:)- min(tcourse_pref(data(c),1,2,:)))...
                ./normrange(c,1);
         %  prenorm(isnan(prenorm))=0; postnorm(isnan(postnorm)) = 0 ;
            prenorm(isinf(prenorm))=nan; postnorm(isinf(postnorm)) = nan ;
            
           prenormds = downsamplebin(prenorm,2,2,1); postnormds = downsamplebin(postnorm,2,2,1);
            mnpre= nanmean(prenormds,1); mnpost = nanmean(postnormds,1);
           % mnpre= nanmean(prenorm,1); mnpost = nanmean(postnorm,1);

            sempre = nanstd(prenormds)/sqrt(length(data));sempost = nanstd(postnormds)/sqrt(length(data));
            %sempre = nanstd(prenorm)/sqrt(length(data));sempost = nanstd(postnorm)/sqrt(length(data));

            
           mnpre = mnpre - repmat(mnpre(:,1),[1 25]); mnpost = mnpost - repmat(mnpost(:,1),[1 25]);
           sempre = sempre - repmat(sempre(:,1),[1 25]); sempost = sempost - repmat(sempost(:,1),[1 25])
          %  mnpre = mnpre - repmat(mnpre(:,1),[1 50]); mnpost = mnpost - repmat(mnpost(:,1),[1 50]);
          % sempre = sempre - repmat(sempre(:,1),[1 50]); sempost = sempost - repmat(sempost(:,1),[1 50]) 
           
            mnpre = circshift(mnpre,[0 5]);
            mnpost = circshift(mnpost,[0 5]);
            sempre = circshift(sempre,[0 5]);
            sempost = circshift(sempost,[0 5]);
            
%             mnpre = circshift(mnpre,[0 10]);
%             mnpost = circshift(mnpost,[0 10]);
%             sempre = circshift(sempre,[0 10]);
%             sempost = circshift(sempost,[0 10])
            ncells = length(data)
        
         end
    x = ((1:length(mnpre)-2.5)*dt -dt/2); x2 = ((1:length(mnpost)-2.5)*dt -dt/2);
    subplot(1,4,i)
%     preerr=shadedErrorBar(x,mnpre(1:45),sempre(1:45)','b',1); hold on; %(3:24)
%     posterr=shadedErrorBar(x2,mnpost(1:45),sempost(1:45)','r',1);
    
     preerr=shadedErrorBar(x,mnpre(1:22),sempre(1:22)','b',1); hold on; %(3:24)
    posterr=shadedErrorBar(x2,mnpost(1:22),sempost(1:22)','r',1);
    axis square; set(gca,'fontsize', 18); hold on; %xlim([0 45])
   % title(titles{t}); 
    title({(titles{t});(layers{i})});
    xlabel('time (s)'); ylabel('spikes/sec');
    text(2 ,.5,max(max((mnpre))+.75), ['n = ' num2str(ncells)],'FontSize',20); 
    xlim([0 2.25]);
end
if savePDF
     set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
end
end
%%
titles={'Saline orth', 'DOI orth'};
for t=1:2
   figure('units','normalized','outerposition',[.5 .5 .75 .75])
       if t==1, set(gcf,'Name','Saline tcourse orth'); else set(gcf,'Name','DOI tcourse orth')
    end
for i = 1:4
    clear mnpre mnpost prenorm postnorm normrange
    if i==1
        data = find(use & layerAll==5 & treatment==t);
     %   set(gcf,'Name','grating lyr5 stat');
    elseif i==2
        clear data
      %  set(gcf,'Name','grating lyr 2/3 stat');
        data = find(use & (layerAll==2 |layerAll==3) & treatment==t);
   elseif i==3
        clear data
       % set(gcf,'Name','grating lyr 4 stat');
        data = find(use & (layerAll==4) & treatment==t);    
    else i==4
        clear data
      %  set(gcf,'Name','grating inh stat');
        data = find(useInh  & treatment==t);
end
         for c = 1:length(data)
            %normalize response of all cells used, then average
            prerange = squeeze(tcourse_orth(data(c),2,1,:))
            normrange(c,:)= max(prerange) %- min(tcourse_pref(data(c),1,1,:),[],1); %set range w/pre data to normalize to
            prenorm(c,:) =(tcourse_orth(data(c),2,1,:)- min(tcourse_orth(data(c),2,1,:)))...
                ./normrange(c,1);
            postnorm(c,:) = (tcourse_orth(data(c),2,2,:)- min(tcourse_orth(data(c),2,2,:)))...
                ./normrange(c,1);
         %   prenorm(isnan(prenorm))=0; postnorm(isnan(postnorm)) = 0 ;
            prenorm(isinf(prenorm))=nan; postnorm(isinf(postnorm)) = nan ;

            prenormds = downsamplebin(prenorm,2,2,1); postnormds = downsamplebin(postnorm,2,2,1);
            mnpre= nanmean(prenormds,1); mnpost = nanmean(postnormds,1);
            sempre = nanstd(prenormds)/sqrt(length(data));sempost = nanstd(postnormds)/sqrt(length(data));
            
%             mnpre= nanmean(prenorm,1);mnpost = nanmean(postnorm,1);
%             sempre = nanstd(prenorm)/sqrt(length(data)); sempost = nanstd(postnorm)/sqrt(length(data));
            
            mnpre = mnpre - repmat(mnpre(:,1),[1 25]); mnpost = mnpost - repmat(mnpost(:,1),[1 25]);
            sempre = sempre - repmat(sempre(:,1),[1 25]); sempost = sempost - repmat(sempost(:,1),[1 25])

            mnpre = circshift(mnpre,[0 5]); mnpost = circshift(mnpost,[0 5]);
            sempre = circshift(sempre,[0 5]);
            sempost = circshift(sempost,[0 5]);
            ncells = length(data)
         end
    x = ((1:length(mnpre)-2.5)*dt -dt/2); x2 = ((1:length(mnpost)-2.5)*dt -dt/2);
    subplot(1,4,i)
    preerr=shadedErrorBar(x,mnpre(3:24),sempre(3:24)','b',1); hold on;
    posterr=shadedErrorBar(x2,mnpost(3:24),sempost(3:24)','r',1);
    axis square; set(gca,'fontsize', 18); hold on; xlim([0 2.5])
   % title(titles{t}); 
    title({(titles{t});(layers{i})});
    xlabel('time (s)'); ylabel('spikes/sec');
    text(2 ,.5,max(max((mnpre))+.75), ['n = ' num2str(ncells)],'FontSize',20); 
end
if savePDF
     set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
end
end

%%

%use = goodAll & ~inhAll & hasDrift==1 & useResp' ;%which cells to include
%useInh = goodAll & inhAll==1 & hasDrift==1 & useResp';
layers = {'Layer 2/3','Layer 4','Layer 5','Inh. Units'};
titles = {'saline all','doi all'};
for t=1:2
   figure('units','normalized','outerposition',[.5 .5 .75 .75])
    if t==1, set(gcf,'Name','Saline tcourse norm all'); else set(gcf,'Name','DOI tcourse norm all')
    end
    for i = 1:4
        % figure
        clear mnpre mnpost prenorm postnorm normrange ncells
        if i==1
            clear data
            % set(gcf,'Name','grating lyr 2/3 stat');
            data = find(use & (layerAll==2 |layerAll==3) & treatment==t);
       %     title('layer 2/3 ')
        elseif i==2
            data = find(use & layerAll==4 & treatment==t);
            %set(gcf,'Name','grating lyr4 stat');
           % title('layer 4')           
        elseif i==3
            
            data = find(use & layerAll==5 & treatment==t);
            %set(gcf,'Name','grating lyr5 stat');
          %  title('layer 5')
        else i==4
            clear data
            % set(gcf,'Name','grating inh stat');
            data = find(useInh  & treatment==t);
        %    title('inh units')     
        end
        for c = 1:length(data)
            %normalize response of all cells used, then average
            prerange = squeeze(tcourse_all(data(c),2,1,:))
            normrange(c,:)= max(prerange) %- min(tcourse_all(data(c),1,:),[],1); %set range w/pre data to normalize to
            prenorm(c,:) =(tcourse_all(data(c),2,1,:)- min(tcourse_all(data(c),2,1,:)))...
                ./normrange(c,1);
            postnorm(c,:) = (tcourse_all(data(c),2,2,:)- min(tcourse_all(data(c),2,2,:)))...
                ./normrange(c,1);
            prenorm(isinf(prenorm))=nan; postnorm(isinf(postnorm)) = nan ;
          
            
            prenormds = downsamplebin(prenorm,2,2,1); postnormds = downsamplebin(postnorm,2,2,1);
            mnpre= nanmean(prenormds,1); mnpost = nanmean(postnormds,1);
            sempre = nanstd(prenormds)/sqrt(length(data));sempost = nanstd(postnormds)/sqrt(length(data));
            %mnpre= nanmean(prenorm,1); mnpost = nanmean(postnorm,1); 
            %sempre = nanstd(prenorm)/sqrt(length(data)); sempost = nanstd(postnorm)/sqrt(length(data)); 
          
            mnpre = mnpre - repmat(mnpre(:,1),[1 25]); mnpost = mnpost - repmat(mnpost(:,1),[1 25]);
            sempre = sempre - repmat(sempre(:,1),[1 25]); sempost = sempost - repmat(sempost(:,1),[1 25])
            
            
            mnpre = circshift(mnpre,[0 5]);mnpost = circshift(mnpost,[0 5]); 
            sempre = circshift(sempre,[0 5]);sempost = circshift(sempost,[0 5]);
  
            ncells = length(data)
        end
        
        x = ((1:length(mnpre)-2.5)*dt -dt/2); x2 = ((1:length(mnpost)-2.5)*dt -dt/2);
        subplot(1,4,i)
        preerr=shadedErrorBar(x,mnpre(3:24),sempre(3:24)','b',1); hold on;
        posterr=shadedErrorBar(x2,mnpost(3:24),sempost(3:24)','r',1);%tightfig;
        axis square; set(gca,'fontsize', 18); hold on; %ylim([0 1]);
        xlim([0 2.5]);
        xlabel('time (s)'); %ylabel('spikes/sec');
        text(1.5 ,.85,max(max((mnpre))+.75), ['n = ' num2str(ncells)],'FontSize',20);
        title({(titles{t});(layers{i})})         
    end
     if savePDF
     set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
  end
end
%%
% figure;
% for c = 1:length(normrange)
%  subplot(6,5,c)
% plot(prenorm(c,:)); hold on; plot(postnorm(c,:));axis square
% end

%%
clear useResp 
peakresp = squeeze(max(drift_orient,[],2))%-squeeze(drift_spont);%subtract driftspont
%peakresp = squeeze(mean(drift_orient,2))-squeeze(drift_spont);%subtract driftspont

%useResp = peakresp(:,1,1)>(min(peakresp(:,1,1)+4)) | peakresp(:,1,2)>(min(peakresp(:,1,2)+4));
useResp = peakresp(:,2,1)>2 | peakresp(:,2,2)>2;
%useOsi_pre = drift_osi(:,1,1)>.3; useOsi_post = drift_osi(:,1,2)>.3; useOsi=(useOsi_pre & useOsi_post)'

use = goodAll &  inhAll==0 & hasDrift==1& useResp' %&useOsi==1;
useInh = goodAll &  inhAll==1 & hasDrift==1& useResp';


% ampmean = squeeze(mean(peakresp,2))
% useResp = ampmean(:,1)>4|ampmean(:,2)>4;
titles=({'saline mn norm','doi mn norm'});
% %dataAll =goodAll & useResp' & hasDrift==1 ;
layers = {'Layer 2/3','Layer 4','Layer 5','Inh. Units'};
for t=1:2
   figure('units','normalized','outerposition',[.5 .5 .75 .75])
    if t==1, set(gcf,'Name','Saline mn stat mv norm'); else set(gcf,'Name','DOI mn stat mv norm')
    end
    for i = 1:4 
        % figure
        clear mnpre mnpost prenorm postnorm normrange ncells
        if i==1
            clear data
            % set(gcf,'Name','grating lyr 2/3 stat');
            data = find(use & (layerAll==2 |layerAll==3) & treatment==t);
       %     title('layer 2/3 ')
        elseif i==2
            data = find(use & layerAll==4 & treatment==t);
            %set(gcf,'Name','grating lyr4 stat');
           % title('layer 4')           
        elseif i==3
            data = find(use & layerAll==5 & treatment==t);
            %set(gcf,'Name','grating lyr5 stat');
          %  title('layer 5')
        else i==4
            clear data
            % set(gcf,'Name','grating inh stat');
            data = find(useInh  & treatment==t);
        %    title('inh units')     
        end
        for c = 1:length(data)
            %normalize response of all cells used, then average
            prerange = squeeze(mnPsth(data(c),1,:))
            normrange(c,:)= max(prerange) %- min(mnPsth(data(c),1,:),[],1); %set range w/pre data to normalize to
            prenorm(c,:) =(mnPsth(data(c),1,:)- min(mnPsth(data(c),1,:)))...
                ./normrange(c,1);
            postnorm(c,:) = (mnPsth(data(c),2,:)- min(mnPsth(data(c),2,:)))...
                ./normrange(c,1);
            prenorm(isinf(prenorm))=nan; postnorm(isinf(postnorm)) = nan ;
            
             prenormds = downsamplebin(prenorm,2,2,1); postnormds = downsamplebin(postnorm,2,2,1);
            mnpre= nanmean(prenormds,1); mnpost = nanmean(postnormds,1);
            sempre = nanstd(prenormds)/sqrt(length(data));sempost = nanstd(postnormds)/sqrt(length(data));
%             mnpre= nanmean(prenorm,1);mnpost = nanmean(postnorm,1); 
%             sempre = nanstd(prenorm)/sqrt(length(data));sempost = nanstd(postnorm)/sqrt(length(data));
            
             mnpre = mnpre - repmat(mnpre(:,1),[1 25]); mnpost = mnpost - repmat(mnpost(:,1),[1 25]);
             sempre = sempre - repmat(sempre(:,1),[1 25]); sempost = sempost - repmat(sempost(:,1),[1 25])
            
            
            mnpre = circshift(mnpre,[0 5]); mnpost = circshift(mnpost,[0 5]); 
            sempre = circshift(sempre,[0 5]); sempost = circshift(sempost,[0 5]);
   
            ncells = length(data)
        end
        

        x = ((1:length(mnpre)-2.5)*dt -dt/2); x2 = ((1:length(mnpost)-2.5)*dt -dt/2);
        subplot(1,4,i)
        preerr=shadedErrorBar(x,mnpre(3:24),sempre(3:24)','b',1); hold on;
        posterr=shadedErrorBar(x2,mnpost(3:24),sempost(3:24)','r',1);%tightfig;
        axis square; set(gca,'fontsize', 18); hold on; %ylim([0 1]);
        xlim([0 2.5]);
        xlabel('time (s)'); %ylabel('spikes/sec');
        text(1.5 ,.85,max(max((mnpre))+.75), ['n = ' num2str(ncells)],'FontSize',20);
        title({(titles{t});(layers{i})});

    end
    if savePDF
     set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
  end
end
%% stats test for temporal components
% layers = {'Layer 2/3','Layer 4','Layer 5','Inh. Units'};
% 
% for t=1:2
%     figure
% for i=1:4
%     if t==1, set(gcf,'Name','Saline'); else set(gcf,'Name','DOI')
%     end
%         % figure
%         clear mnpre mnpost prenorm postnorm normrange ncells pre post
%     clear preTrans postTrans preSus postSus preSpont postSpont 
%        if i==1
%             clear data
%             % set(gcf,'Name','grating lyr 2/3 stat');
%             data = find(use & (layerAll==2 |layerAll==3) & treatment==t);
%        %     title('layer 2/3 ')
%         elseif i==2
%             data = find(use & layerAll==4 & treatment==t);
%             %set(gcf,'Name','grating lyr4 stat');
%            % title('layer 4')        
%         elseif i==3
%             data = find(use & layerAll==5 & treatment==t);
%             %set(gcf,'Name','grating lyr5 stat');
%           %  title('layer 5')
%         else i==4
%             clear data
%             % set(gcf,'Name','grating inh stat');
%             data = find(useInh  & treatment==t);
%         %    title('inh units')     
%         end
%         for c = 1:length(data)
%             %normalize response of all cells used, then average
%             prerange = squeeze(mnPsth(data(c),1,:))
%             normrange(c,:)= mean(prerange)% - min(mnPsth(data(c),1,:),[],1); %set range w/pre data to normalize to
%             
%             prenorm(c,:) =(mnPsth(data(c),1,:)- min(mnPsth(data(c),1,:)))...
%                 ./normrange(c,1);
%             postnorm(c,:) = (mnPsth(data(c),2,:)- min(mnPsth(data(c),2,:)))...
%                 ./normrange(c,1);
%             prenorm(isinf(prenorm))=nan; postnorm(isinf(postnorm)) = nan ;
%             pre= circshift(prenorm,[0 10]); post = circshift(postnorm,[0 10]);        
%             %              preTrans= pre(c,5:15); postTrans= post(c,5:15);
%             %              preSus= pre(c,16:31); postSus = post(c,16:31);
%             preSpont(c,:) = pre(c,1:10);postSpont(c,:)=post(c,1:10);
%             preTrans(c,:) = pre(c,11:14); postTrans(c,:)= post(c,11:14);
%             preSus(c,:)= pre(c,15:39); postSus(c,:) = post(c,15:39);
%             [hSpont p] = ttest(preSpont,postSpont);
%             [hTrans p] = ttest(preTrans,postTrans);
%             [hSus p] = ttest(preSus,postSus);
%             %%0 = reject null if no difference
%                ncells = length(data);
% %            corrTrans = corrcoef(preTrans,postTrans,'rows','pairwise'); corrTrans = corrTrans(2,1);
% %            corrSus = corrcoef(preSus,postSus,'rows','pairwise'); corrSus = corrSus(2,1);
% 
%         end
% 
%         subplot(3,4,i)
%         spont = hist(hSpont,0:1:1); bSp=bar(0:1:1,spont/length(hSpont));
%         ylim([0 1]); xlim([-.50 1.5]);axis square
%         title(layers{i});%ylabel('proportion of cells');
%         set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%         
%         subplot(3,4,i+4)
%         
%         trans = hist(hTrans,0:1:1); bT=bar(0:1:1,trans/length(hTrans));
%         ylim([0 1]); xlim([-.50 1.5]);axis square
%         title(layers{i});%ylabel('proportion of cells');
%         set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%         
%         subplot(3,4,i+8)
%         sus = hist(hSus,0:1:1);
%         bS=bar(0:1:1,sus/length(hSus)); ylim([0 1]);xlim([-.5 1.5]);
%         axis square;set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%         title(layers{i});%ylabel('proportion of cells');
% 
%  end
% end


%% pref psth, not normalized

peakresp = squeeze(max(drift_orient,[],2))%-squeeze(drift_spont);%subtract driftspont
useResp = peakresp(:,1,1)>2 | peakresp(:,1,2)>2;
useOsi_pre = drift_osi(:,2,1)>.3;
useOsi_post = drift_osi(:,2,2)>.3; useOsi=(useOsi_pre & useOsi_post)'

use = goodAll &  inhAll==0 & hasDrift==1& useResp'% &useOsi==1;
useInh = goodAll &  inhAll==1 & hasDrift==1& useResp'% &useOsi==1;
layers = {'Layer 2/3','Layer 4','Layer 5','Inh. Units'};

titles= {'saline pref','doi pref'};
dt = 0.1;
clear used sempre
for t=1:2
   figure('units','normalized','outerposition',[.5 .5 .75 .75])
    if t==1, set(gcf,'Name','Saline not-norm pref'); else set(gcf,'Name','DOI not-norm pref')
    end
    for i = 1:4
         %figure
        clear mnpre mnpost prenorm postnorm normrange ncells preS postS 
        clear preSpont postSpont preTrans postTrans preSus postSus preOff postOff
        if i==1
            clear data
            % set(gcf,'Name','grating lyr 2/3 stat');
            data = find(use & (layerAll==2 |layerAll==3) & treatment==t);
            %     title('layer 2/3 ')
        elseif i==2
            data = find(use & layerAll==4 & treatment==t);
            %set(gcf,'Name','grating lyr4 stat');
            % title('layer 4')
            
        elseif i==3
            data = find(use & layerAll==5 & treatment==t);
            %set(gcf,'Name','grating lyr5 stat');
            %  title('layer 5')
        else i==4
            clear data
            % set(gcf,'Name','grating inh stat');
            data = find(useInh  & treatment==t);
            %    title('inh units')
        end
        for c = 1:length(data)
            pre(c,:) =(tcourse_pref(data(c),1,1,:));post(c,:) =(tcourse_pref(data(c),1,2,:));
            
            preds = downsamplebin(pre,2,2,1); postds = downsamplebin(post,2,2,1);
            mnpre = squeeze(nanmean(preds,1));
            mnpost = squeeze(nanmean(postds,1));
          
            sempre = nanstd(preds)/sqrt(length(data)); %sempre = downsamplebin(sempre,2,2,1);
            sempost = nanstd(postds)/sqrt(length(data));% sempost = downsamplebin(sempost,2,2,1);

            ncells = length(data);
            preS = circshift(pre,[0 5]); postS= circshift(post,[0 5]);
            
%             preSpont(c,:) = preS(c,2:10);postSpont(c,:)=postS(c,2:10);
%             preTrans(c,:) = preS(c,11:14); postTrans(c,:)= postS(c,11:14);
%             preSus(c,:)= preS(c,15:38); postSus(c,:) = postS(c,15:38);
%             preOff(c,:)= preS(c,39:46); postOff(c,:) = postS(c,39:46);
% 
%             [hSpont p] = ttest(preSpont,postSpont);
%             [hTrans p] = ttest(preTrans,postTrans);
%             [hSus p] = ttest(preSus,postSus);
%             [hOff p] = ttest(preOff,postOff);

        end
        
    %   mnpre = mnpre - repmat(mnpre(:,1),[1 25]); mnpost = mnpost - repmat(mnpost(:,1),[1 25]);
     %  sempre = sempre - repmat(sempre(:,1),[1 25]); sempost = sempost - repmat(sempost(:,1),[1 25])
            
    
        mnpre = circshift(mnpre,[0 5]);
        mnpost = circshift(mnpost,[0 5]);
        sempre = circshift(sempre,[0 5]);
        sempost = circshift(sempost,[0 5]);
       
        subplot(1,4,i)
        x = ((1:length(mnpre)-2.5)*dt -dt/2); x2 = ((1:length(mnpost)-2.5)*dt -dt/2);
      %  plot(mnpre); hold on; plot(mnpost,'r');
        preerr=shadedErrorBar(x,mnpre(1:22),sempre(1:22)','b',1); hold on;
        posterr=shadedErrorBar(x2,mnpost(1:22),sempost(1:22)','r',1);
        axis square; set(gca,'fontsize', 18); hold on; %xlim([0 2.5]);
        title({(titles{t});(layers{i})});
        xlabel('time (s)'); ylabel('spikes/sec');
        text(1.75,max(max((mnpre))+.75), ['n = ' num2str(ncells)],'FontSize',20);
         xlim([ 0 2.25]);
        %         subplot(5,4,i+4)
        %         spont = hist(hSpont,0:1:1); bSp=bar(0:1:1,spont/length(hSpont));
        %         ylim([0 1]); xlim([-.50 1.5]);axis square
        %         title(layers{i});%ylabel('proportion of cells');
        %         set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%         axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
% 
%      subplot(5,4,i+8) 
% 
%      trans = hist(hTrans,0:1:1); bT=bar(0:1:1,trans/length(hTrans)); 
%      ylim([0 1]); xlim([-.50 1.5]);axis square
%      title(layers{i});%ylabel('proportion of cells');
%      set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%      axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
% 
%      subplot(5,4,i+12) 
%      sus = hist(hSus,0:1:1);
%      bS=bar(0:1:1,sus/length(hSus)); ylim([0 1]);xlim([-.5 1.5]);
%      axis square;set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%      axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
% 
%      subplot(5,4,i+16) 
%      off = hist(hOff,0:1:1);
%      bS=bar(0:1:1,off/length(hOff)); ylim([0 1]);xlim([-.5 1.5]);
%      axis square;set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%      axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
     
    title(layers{i});%ylabel('proportion of cells');
    end
if savePDF
     set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
end
end

%%
useResp = peakresp(:,1,1)>2 | peakresp(:,1,2)>2;
use = goodAll &  inhAll==0 & hasDrift==1& useResp';
for t=1:2
    data = find(use & (layerAll==2|layerAll==3) & treatment==t); 
    figure
 if t==1, set(gcf,'Name','Saline not-norm all'); else set(gcf,'Name','DOI not-norm all')
    end
for c=1:length(data)
   % set (gcf,'Name','Lyr 2/3 psth/unit') % set up in loop
    subplot(10,5,c)
    plot(squeeze(tcourse_pref(data(c),1,1,:)));hold on;
    plot(squeeze(tcourse_pref(data(c),1,2,:)),'r'); axis square
end
end


%% orthogonal psth

peakresp = squeeze(max(drift_orient,[],2))%-squeeze(drift_spont);%subtract driftspont
useResp = peakresp(:,1,1)>2 | peakresp(:,1,2)>2;
useOsi_pre = drift_osi(:,2,1)>.3;
useOsi_post = drift_osi(:,2,2)>.3; useOsi=(useOsi_pre & useOsi_post)'

use = goodAll &  inhAll==0 & hasDrift==1& useResp'% &useOsi==1;
useInh = goodAll &  inhAll==1 & hasDrift==1& useResp'% &useOsi==1;
layers = {'Layer 2/3','Layer 4','Layer 5','Inh. Units'};
titles = {'saline orth','doi orth'};
dt = 0.1;
clear used sempre
for t=1:2
   figure('units','normalized','outerposition',[.5 .5 .75 .75])
    if t==1, set(gcf,'Name','Saline not-norm orth'); else set(gcf,'Name','DOI not-norm orth')
    end
    for i = 1:4
         %figure
        clear mnpre mnpost prenorm postnorm normrange ncells preS postS  pre post sempre sempost
        clear preSpont postSpont preTrans postTrans preSus postSus preOff postOff
        if i==1
            clear data
            % set(gcf,'Name','grating lyr 2/3 stat');
            data = find(use & (layerAll==2 |layerAll==3) & treatment==t);
            %     title('layer 2/3 ')
        elseif i==2
            data = find(use & layerAll==4 & treatment==t);
            %set(gcf,'Name','grating lyr4 stat');
            % title('layer 4')
            
        elseif i==3
            data = find(use & layerAll==5 & treatment==t);
            %set(gcf,'Name','grating lyr5 stat');
            %  title('layer 5')
        else i==4
            clear data
            % set(gcf,'Name','grating inh stat');
            data = find(useInh  & treatment==t);
            %    title('inh units')
        end
        for c = 1:length(data)
            pre(c,:) =(tcourse_orth(data(c),2,1,:));post(c,:) =(tcourse_orth(data(c),2,2,:));
            preds = downsamplebin(pre,2,2,1);postds = downsamplebin(post,2,2,1);

            mnpre = squeeze(nanmean(preds,1));
            mnpost = squeeze(nanmean(postds,1));
            
       %   sempre = downsamplebin(sempre,2,2,1);
            
            sempre = nanstd(preds)/sqrt(length(data));
            sempost = nanstd(postds)/sqrt(length(data));
            ncells = length(data);
            preS = circshift(pre,[0 5]); postS = circshift(post,[0 5]);
            
            preSpont(c,:) = preS(c,2:10);postSpont(c,:)=postS(c,2:10);
            preTrans(c,:) = preS(c,11:14); postTrans(c,:)= postS(c,11:14);
            preSus(c,:)= preS(c,15:38); postSus(c,:) = postS(c,15:38);
            preOff(c,:)= preS(c,39:46); postOff(c,:) = postS(c,39:46);

            [hSpont p] = ttest(preSpont,postSpont);
            [hTrans p] = ttest(preTrans,postTrans);
            [hSus p] = ttest(preSus,postSus);
            [hOff p] = ttest(preOff,postOff);

        end
        
         mnpre = mnpre - repmat(mnpre(:,1),[1 25]); mnpost = mnpost - repmat(mnpost(:,1),[1 25]);
         sempre = sempre - repmat(sempre(:,1),[1 25]); sempost = sempost - repmat(sempost(:,1),[1 25])
          
        
        mnpre = circshift(mnpre,[0 5]);
        mnpost = circshift(mnpost,[0 5]);
        sempre = circshift(sempre,[0 5]);
        sempost = circshift(sempost,[0 5]);
       
        subplot(1,4,i)
        x = ((1:length(mnpre)-2.5)*dt -dt/2); x2 = ((1:length(mnpost)-2.5)*dt -dt/2);
        preerr=shadedErrorBar(x,mnpre(3:24),sempre(3:24)','b',1); hold on;
        posterr=shadedErrorBar(x2,mnpost(3:24),sempost(3:24)','r',1);
        axis square; set(gca,'fontsize', 18); hold on; %xlim([0 2.5]);
           title({(titles{t});(layers{i})});
        xlabel('time (s)'); ylabel('spikes/sec');
        text(1.75,max(max((mnpre))+.75), ['n = ' num2str(ncells)],'FontSize',20);
     %   xlim([ 0 2.5]);
%         subplot(5,4,i+4)
%         spont = hist(hSpont,0:1:1); bSp=bar(0:1:1,spont/length(hSpont));
%         ylim([0 1]); xlim([-.50 1.5]);axis square
%         title(layers{i});%ylabel('proportion of cells');
%         set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%         axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
% 
%      subplot(5,4,i+8) 
% 
%      trans = hist(hTrans,0:1:1); bT=bar(0:1:1,trans/length(hTrans)); 
%      ylim([0 1]); xlim([-.50 1.5]);axis square
%      title(layers{i});%ylabel('proportion of cells');
%      set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%      axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
% 
%      subplot(5,4,i+12) 
%      sus = hist(hSus,0:1:1);
%      bS=bar(0:1:1,sus/length(hSus)); ylim([0 1]);xlim([-.5 1.5]);
%      axis square;set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%      axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
% 
%      subplot(5,4,i+16) 
%      off = hist(hOff,0:1:1);
%      bS=bar(0:1:1,off/length(hOff)); ylim([0 1]);xlim([-.5 1.5]);
%      axis square;set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%      axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
     
   %ylabel('proportion of cells');
    end
if savePDF
     set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
end
end
%%
peakresp = squeeze(max(drift_orient,[],2))%-squeeze(drift_spont);%subtract driftspont
useResp = peakresp(:,2,1)>2 & peakresp(:,2,2)>2;
%useResp =( peakresp(:,1,1)>2 & peakresp(:,1,1)<15)& (peakresp(:,1,2)>2 & peakresp(:,1,2)<15) ;

useOsi_pre = drift_osi(:,2,1)>.2;
useOsi_post = drift_osi(:,2,2)>.2; useOsi=(useOsi_pre & useOsi_post)'

use = goodAll &  inhAll==0 & hasDrift==1& useResp' %&useOsi==1;
useInh = goodAll &  inhAll==1 & hasDrift==1& useResp' %&useOsi==1;
layers = {'Layer 2/3','Layer 4','Layer 5','Inh. Units'};
titles= {'saline all','doi all'};
   
    
dt = 0.05;
clear used sempre
for t=1:2
figure('units','normalized','outerposition',[.5 .5 .75 .75]); 
    if t==1, set(gcf,'Name','Saline not-norm all'); else set(gcf,'Name','DOI not-norm all')
    end
    for i = 1:4
         %figure
        clear mnpre mnpost prenorm postnorm normrange ncells preS postS pre post sempre sempost
        clear preSpont postSpont preTrans postTrans preSus postSus preOff postOff
        if i==1
            clear data
            % set(gcf,'Name','grating lyr 2/3 stat');
            data = find(use & (layerAll==2 |layerAll==3) & treatment==t);
            %     title('layer 2/3 ')
        elseif i==2
            data = find(use & layerAll==4 & treatment==t);
            %set(gcf,'Name','grating lyr4 stat');
            % title('layer 4')
        elseif i==3
            data = find(use & layerAll==5 & treatment==t);
            %set(gcf,'Name','grating lyr5 stat');
            %  title('layer 5')
        else i==4
            clear data
            % set(gcf,'Name','grating inh stat');
            data = find(useInh  & treatment==t);
            %    title('inh units')
        end
        for c = 1:length(data)
            pre(c,:) =(tcourse_all(data(c),2,1,:));post(c,:) =(tcourse_all(data(c),2,2,:));
          %  preds = downsamplebin(pre,2,2,1);postds = downsamplebin(post,2,2,1);

%            mnpre = squeeze(nanmean(preds,1));
%            mnpost = squeeze(nanmean(postds,1));
%             
            mnpre = squeeze(nanmean(pre,1));
            mnpost = squeeze(nanmean(post,1));
            
%             sempre = nanstd(preds)/sqrt(length(data));
%             sempost = nanstd(postds)/sqrt(length(data));
            sempre = nanstd(pre)/sqrt(length(data));
            sempost = nanstd(post)/sqrt(length(data));
            ncells = length(data);
            preS = circshift(pre,[0 10]); postS= circshift(post,[0 10]);
            
            preSpont(c,:) = preS(c,1:9);postSpont(c,:)=postS(c,1:9);
            preTrans(c,:) = preS(c,11:15); postTrans(c,:)= postS(c,11:15);
            preSus(c,:)= preS(c,17:37); postSus(c,:) = postS(c,17:37);
            preOff(c,:)= preS(c,40:end); postOff(c,:) = postS(c,40:end);

            [hSpont p] = ttest(preSpont,postSpont);
            [hTrans p] = ttest(preTrans,postTrans);
            [hSus p] = ttest(preSus,postSus);
            [hOff p] = ttest(preOff,postOff);

        end
%         
%         mnpre = mnpre - repmat(mnpre(:,1),[1 25]); mnpost = mnpost - repmat(mnpost(:,1),[1 25]);
%         sempre = sempre - repmat(sempre(:,1),[1 25]); sempost = sempost - repmat(sempost(:,1),[1 25]);
%
%         mnpre = mnpre - repmat(mnpre(:,1),[1 50]); mnpost = mnpost - repmat(mnpost(:,1),[1 50]);
         % sempre = sempre - repmat(sempre(:,1),[1 50]); sempost = sempost - repmat(sempost(:,1),[1 50])

%         mnpre = circshift(mnpre,[0 5]);
%         mnpost = circshift(mnpost,[0 5]);
%         sempre = circshift(sempre,[0 5]);
%         sempost = circshift(sempost,[0 5]);

            mnpre = circshift(mnpre,[0 10]);
            mnpost = circshift(mnpost,[0 10]);
            sempre = circshift(sempre,[0 10]);
            sempost = circshift(sempost,[0 10])
       

        subplot(1,4,i)
        x = ((1:length(mnpre)-5)*dt -dt/2); x2 = ((1:length(mnpost)-5)*dt -dt/2); 
        preerr=shadedErrorBar(x,mnpre(1:45),sempre(1:45)','b',1); hold on;%(3:24)
        posterr=shadedErrorBar(x2,mnpost(1:45),sempost(1:45)','r',1);
        axis square; set(gca,'fontsize', 24); hold on; 
        title({(titles{t});(layers{i})});
        xlabel('time (s)'); ylabel('spikes/sec');
        text(1.75,max(max((mnpre))+.75), ['n = ' num2str(ncells)],'FontSize',20);
       xlim([0 2.25]); 

%         subplot(4,3,i)
%         spont = hist(hSpont,0:1:1); bSp=bar(0:1:1,spont/length(hSpont));
%         ylim([0 1]); xlim([-.50 1.5]);axis square
%         title(layers{i});%ylabel('proportion of cells');
%         set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%         axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
% 
%      subplot(4,3,i+3) 
% 
%      trans = hist(hTrans,0:1:1); bT=bar(0:1:1,trans/length(hTrans)); 
%      ylim([0 1]); xlim([-.50 1.5]);axis square
%      title(layers{i});%ylabel('proportion of cells');
%      set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%      axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
% 
%      subplot(4,3,i+6) 
%      sus = hist(hSus,0:1:1);
%      bS=bar(0:1:1,sus/length(hSus)); ylim([0 1]);xlim([-.5 1.5]);
%      axis square;set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%      axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
% 
%      subplot(4,3,i+9) 
%      off = hist(hOff,0:1:1);
%      bS=bar(0:1:1,off/length(hOff)); ylim([0 1]);xlim([-.5 1.5]);
%      axis square;set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%      axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
%      
   % title(layers{i});%ylabel('proportion of cells');
      end
if savePDF
     set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
end
end





%%
clear useResp 
peakresp = squeeze(max(drift_orient,[],2))%-squeeze(drift_spont);%subtract driftspont
%peakresp = squeeze(mean(drift_orient,2))-squeeze(drift_spont);%subtract driftspont

%useResp = peakresp(:,1,1)>(min(peakresp(:,1,1)+4)) | peakresp(:,1,2)>(min(peakresp(:,1,2)+4));
% useResp = peakresp(:,2,1)>2 | peakresp(:,2,2)>2;
% useOsi_pre = drift_osi(:,1,1)>.3; useOsi_post = drift_osi(:,1,2)>.3; useOsi=(useOsi_pre & useOsi_post)'

ampmean = squeeze(mean(peakresp,2))
useResp = ampmean(:,1)>3|ampmean(:,2)>3;


%useResp = peakresp(:,1,1)>0| peakresp(:,1,2)>0;
%useOsi_pre = drift_osi(:,1,1)>0; useOsi_post = drift_osi(:,1,2)>0; useOsi=(useOsi_pre & useOsi_post)';
use = goodAll &  inhAll==0 & hasDrift==1& useResp'; %&useOsi==1;
useInh = goodAll &  inhAll==1 & hasDrift==1& useResp';

titles= {'saline mn','doi mn'};
layers = {'Layer 2/3','Layer 4','Layer 5','Inh. Units'};
dt = .1;
clear used sempre
for t=1:2
   figure('units','normalized','outerposition',[.5 .5 .75 .75])
    if t==1, set(gcf,'Name','Saline mn not-norm'); else set(gcf,'Name','DOI mean not-norm'), end
    for i = 1
        % figure
        clear mnpre mnpost prenorm postnorm normrange ncells preS postS pre post
        clear preSpont postSpont preTrans postTrans preSus postSus
        if i==1
            clear data
            % set(gcf,'Name','grating lyr 2/3 stat');
            data = find(use & (layerAll==2 |layerAll==3) & treatment==t);
            %     title('layer 2/3 ')
        elseif i==2
            data = find(use & layerAll==4 & treatment==t);
            %set(gcf,'Name','grating lyr4 stat');
            % title('layer 4')  
        elseif i==3
            data = find(use & layerAll==5 & treatment==t);
            %set(gcf,'Name','grating lyr5 stat');
            %  title('layer 5')
        else i==4
            clear data
            % set(gcf,'Name','grating inh stat');
            data = find(useInh  & treatment==t);
            %    title('inh units')
        end
        
        for c = 1:length(data)
            pre(c,:) =(mnPsth(data(c),1,:));post(c,:) =(mnPsth(data(c),2,:));
            preds = downsamplebin(pre,2,2,1);postds = downsamplebin(post,2,2,1);

            mnpre = squeeze(nanmean(preds,1));
            mnpost = squeeze(nanmean(postds,1));
            
            sempre = nanstd(preds)/sqrt(length(data));
            sempost = nanstd(postds)/sqrt(length(data));
            ncells = length(data);
%             preS = circshift(pre,[0 10]); postS= circshift(post,[0 10]);
%             
%             preSpont(c,:) = preS(c,1:10);postSpont(c,:)=postS(c,1:10);
%             preTrans(c,:) = preS(c,11:14); postTrans(c,:)= postS(c,11:14);
%             preSus(c,:)= preS(c,15:38); postSus(c,:) = postS(c,15:38);
%             preOff(c,:)= preS(c,39:46); postOff(c,:) = postS(c,39:46);
% 
%             [hSpont p] = ttest(preSpont,postSpont);
%             [hTrans p] = ttest(preTrans,postTrans);
%             [hSus p] = ttest(preSus,postSus);
%             [hOff p] = ttest(preOff,postOff);

        end
        mnpre = mnpre - repmat(mnpre(:,1),[1 25]); mnpost = mnpost - repmat(mnpost(:,1),[1 25]);
       sempre = sempre - repmat(sempre(:,1),[1 25]); sempost = sempost - repmat(sempost(:,1),[1 25]);
        
        mnpre = circshift(mnpre,[0 5]);
        mnpost = circshift(mnpost,[0 5]);
        sempre = circshift(sempre,[0 5]);
        sempost = circshift(sempost,[0 5]);
        
        subplot(1,4,i)
        x = ((1:length(mnpre)-2.5)*dt -dt/2); x2 = ((1:length(mnpost)-2.5)*dt -dt/2);
        preerr=shadedErrorBar(x,mnpre(3:24),sempre(3:24)','b',1); hold on;
        posterr=shadedErrorBar(x2,mnpost(3:24),sempost(3:24)','r',1);
        axis square; set(gca,'fontsize', 18); hold on; %xlim([0 2.5]);
        title({(titles{t});(layers{i})});
        xlabel('time (s)'); ylabel('spikes/sec');
        text(1.75,max(max((mnpre))+.75), ['n = ' num2str(ncells)],'FontSize',20);
        xlim([0 1.2]);
        
%         subplot(5,4,i+4)
%         spont = hist(hSpont,0:1:1); bSp=bar(0:1:1,spont/length(hSpont));
%         ylim([0 1]); xlim([-.50 1.5]);axis square
%         title(layers{i});%ylabel('proportion of cells');
%         set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%         axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
%         
%         subplot(5,4,i+8)  
%         trans = hist(hTrans,0:1:1); bT=bar(0:1:1,trans/length(hTrans));
%         ylim([0 1]); xlim([-.50 1.5]);axis square
%         title(layers{i});%ylabel('proportion of cells');
%         set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%         axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
%         
%         subplot(5,4,i+12)
%         sus = hist(hSus,0:1:1);
%         bS=bar(0:1:1,sus/length(hSus)); ylim([0 1]);xlim([-.5 1.5]);
%         axis square;set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%         axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
%         
%       subplot(5,4,i+16) 
%      off = hist(hOff,0:1:1);
%      bS=bar(0:1:1,off/length(hOff)); ylim([0 1]);xlim([-.5 1.5]);
%      axis square;set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%      axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
              
        title(layers{i});%ylabel('proportion of cells');
    end
if savePDF
     set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
end
end

%%

% layers = {'Layer 2/3','Layer 4','Layer 5','Inh. Units'};
% for t=1:2
%     figure 
%     if t==1, set(gcf,'Name','Saline mean not-norm'); else set(gcf,'Name','DOI mean not-norm')
%     end
%     for i = 1:4
%         clear mnpre mnpost prenorm postnorm normrange ncells pre post
%         clear preTrans postTrans preSus postSus preSpont postSpont
%         if i==1
%             clear data
%             set(gcf,'Name','grating lyr 2/3 stat');
%             data = find(use & (layerAll==2 |layerAll==3) & treatment==t);
%                 title('layer 2/3 ')
%         elseif i==2
%             data = find(use & layerAll==4 & treatment==t);
%             set(gcf,'Name','grating lyr4 stat');
%             title('layer 4')         
%         elseif i==3
%             data = find(use & layerAll==5 & treatment==t);
%             set(gcf,'Name','grating lyr5 stat');
%              title('layer 5')
%         else i==4
%             clear data
%             set(gcf,'Name','grating inh stat');
%             data = find(useInh  & treatment==t);
%                title('inh units')
%         end   
%         for c = 1:length(data)
%             pre(c,:) =(mnPsth(data(c),1,:)); post(c,:) =(mnPsth(data(c),2,:));
%             
%             mnpre = squeeze(nanmean(pre,1));
%             mnpost = squeeze(nanmean(post,1));
%             
%             sempre = nanstd(pre)/sqrt(length(data));
%             sempost = nanstd(post)/sqrt(length(data));
%             ncells = length(data);
%              pre = pre - repmat(pre(1,:),[50 1]);%mnpre = circshift(mnpre,[0 10]);  
%             
%            mnpre = mnpre - repmat(mnpre(:,1),[1 50]);
%             mnpre = circshift(mnpre,[0 10]);
%            mnpost = mnpost - repmat(mnpost(:,1),[1 50]);
%             mnpost = circshift(mnpost,[0 10]);
%            sempre = sempre - repmat(sempre(:,1),[1 50]);
%             sempre = circshift(sempre,[0 10]);
%            sempost = sempost - repmat(sempost(:,1),[1 50]);
%             sempost = circshift(sempost,[0 10]);
%             subplot(1,4,i) 
%             pre = pre - repmat(pre(c,:),[50 1]);
%             pre = circshift(pre,[0 10]); post= circshift(post,[0 10]);    
%             %              preTrans= pre(c,5:15); postTrans= post(c,5:15);
%             %              preSus= pre(c,16:31); postSus = post(c,16:31);
%             preSpont(c,:) = pre(c,2:10);postSpont(c,:)=post(c,2:10);
%             preTrans(c,:) = pre(c,11:14); postTrans(c,:)= post(c,11:14);
%             preSus(c,:)= pre(c,15:38); postSus(c,:) = post(c,15:38);
%             preOff(c,:)= preS(c,39:46); postOff(c,:) = postS(c,39:46);
% 
%             [hSpont p] = ttest(preSpont,postSpont);
%             [hTrans p] = ttest(preTrans,postTrans);
%             [hSus p] = ttest(preSus,postSus);
%             [hOff p] = ttest(preOff,postOff);
%                         %%0 = reject null if no difference
%             ncells = length(data);
%                        corrTrans = corrcoef(preTrans,postTrans,'rows','pairwise'); corrTrans = corrTrans(2,1);
%                        corrSus = corrcoef(preSus,postSus,'rows','pairwise'); corrSus = corrSus(2,1);
%         end
%         
%         subplot(4,4,i)
%         spont = hist(hSpont,0:1:1); bSp=bar(0:1:1,spont/length(hSpont));
%         ylim([0 1]); xlim([-.50 1.5]);axis square
%         title(layers{i});%ylabel('proportion of cells');
%         set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%         axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
%         
%         subplot(4,4,i+4)
%         
%         trans = hist(hTrans,0:1:1); bT=bar(0:1:1,trans/length(hTrans));
%         ylim([0 1]); xlim([-.50 1.5]);axis square
%         title(layers{i});%ylabel('proportion of cells');
%         set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%         axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
%         
%         subplot(4,4,i+8)
%         sus = hist(hSus,0:1:1);
%         bS=bar(0:1:1,sus/length(hSus)); ylim([0 1]);xlim([-.5 1.5]);
%         axis square;set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%         axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
%         
%         title(layers{i});%ylabel('proportion of cells');
%         
%     subplot(4,4,i+12) 
%      off = hist(hOff,0:1:1);
%      bS=bar(0:1:1,off/length(hOff)); ylim([0 1]);xlim([-.5 1.5]);
%      axis square;set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
%      axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
%      
%     title(layers{i});%ylabel('proportion of cells');
% 
%     end
% end




%%
peakresp = squeeze(max(drift_orient,[],2))%-squeeze(drift_spont);%subtract driftspont
% ampmean = squeeze(mean(peakresp,2))
%useResp =( peakresp(:,1,1)>2 & peakresp(:,1,1)<15)| (peakresp(:,1,2)>2 & peakresp(:,1,2)<15) ;
useResp = peakresp(:,1,1)>2 & peakresp(:,1,2)>.5;

% mnOsi = squeeze(nanmean(drift_osi,2))
% useOsi_pre = mnOsi(:,1)>.3;
% useOsi_post =  mnOsi(:,2)>.3; useOsi=(useOsi_pre & useOsi_post)'
%useResp = ampmean(:,1)>2 &ampmean(:,2)>2;
 % &useOsi; ;%which cells to include

%useN= find(useResp' & goodAll==1 & hasDrift & treatment==doi &layerAll==2)

use = goodAll &  inhAll==0 & hasDrift==1& useResp';
for t=1:2
    data = find(use & (layerAll==4) & treatment==t); 
    figure
 if t==1, set(gcf,'Name','Saline not-norm all'); else set(gcf,'Name','DOI not-norm all')
    end
for c=1:length(data)
   % set (gcf,'Name','Lyr 2/3 psth/unit') % set up in loop
    subplot(10,5,c)
    plot(squeeze(tcourse_all(data(c),1,1,:)));hold on;
    plot(squeeze(tcourse_all(data(c),1,2,:)),'r'); axis square
end
end

% use = goodAll &  inhAll==0 & hasDrift==1& useResp';% &useOsi==1;
%%
useResp = ampmean(:,1)>5&ampmean(:,2)>5;
use = goodAll &   hasDrift==1& useResp'; %&useOsi==1;inhAll==0 
useInh = goodAll &  inhAll==1 & hasDrift==1& useResp';
mnPsthDs = downsamplebin(mnPsth,3,2,1)

titles = {'Saline','DOI'};
for t=1:2
   figure%('units','normalized','outerposition',[.8 0 1 1])
    colormap jet
    clear mnpre mnpost prenorm postnorm normrange ncells preS postS pre post
    clear preSpont postSpont preTrans postTrans preSus postSus
    data = find(use & treatment==t);
     for c = 1:length(data)
            pre =(mnPsthDs(data(c),1,:)); post =(mnPsthDs(data(c),2,:)); %
            mnpre = squeeze(nanmean(pre,1));
            mnpost = squeeze(nanmean(post,1));  
            ncells = length(data);
            preS(c,:) = circshift(pre,[0 5]); postS(c,:)= circshift(post,[0 5]);
            
            preSpont(c,:) = preS(c,1:5);postSpont(c,:)=postS(c,1:5);
            preTrans(c,:) = preS(c,6:7); postTrans(c,:)= postS(c,6:7);
            preSus(c,:)= preS(c,8:16); postSus(c,:) = postS(c,8:16);
            preOff(c,:)= preS(c,17:23); postOff(c,:) = postS(c,17:23);
     end
%             subplot(3,2,1)
%             imagesc(corrcoef(preSpont'));axis square; 
%             subplot(3,2,2)
%             imagesc(corrcoef(postSpont'));axis square
%             ccSpont = corrcoef(preSpont,postSpont); ccSpont=ccSpont(2,1);
%             title({(titles{t});(ccSpont)});
% 
%              
%             subplot(3,2,3)
%             imagesc(corrcoef(preTrans'));axis square
%             subplot(3,2,4)
%             imagesc(corrcoef(postTrans'));axis square
%             ccTrans = corrcoef(preTrans,postTrans); ccTrans=ccTrans(2,1);
%             title({(titles{t});(ccTrans)});
% 
% 
%             
%             subplot(3,2,5)
%             imagesc(corrcoef(preSus')); axis square
%             subplot(3,2,6)
%             imagesc(corrcoef(postSus'));axis square
%             ccSus = corrcoef(preSus,postSus); ccSus=ccSus(2,1);
%             title({(titles{t});(ccSus)});
            


preSpontmn = squeeze(nanmean(preSpont,2)); postSpontmn = squeeze(nanmean(postSpont,2));
subplot(1,4,1);hold on
plot([1 2],[preSpontmn postSpontmn],'k-')
errorbar([1 2],[mean(preSpontmn) mean(postSpontmn)],[std(preSpontmn)/length(preSpont) std(postSpontmn)/length(postSpont)],'r','linewidth',2)
%axis([0 3 -2 2]); 
axis square
set(gca,'xtick',[1 2],'xticklabel',{'pre','post'})

preTransmn = squeeze(nanmean(preTrans,2)); postTransmn = squeeze(nanmean(postTrans,2));
subplot(1,4,2);hold on
plot([1 2],[preTransmn postTransmn],'k-')
errorbar([1 2],[mean(preTransmn) mean(postTransmn)],[std(preTransmn)/length(preTrans) std(postTransmn)/length(postTrans)],'r','linewidth',2)
%axis([0 3 -2 2]); 
axis square
set(gca,'xtick',[1 2],'xticklabel',{'pre','post'})

preSusmn = squeeze(nanmean(preSus,2)); postSusmn = squeeze(nanmean(postSus,2));
subplot(1,4,3); hold on
plot([1 2],[preSusmn postSusmn],'k-')
errorbar([1 2],[mean(preSusmn) mean(postSusmn)],[std(preSusmn)/length(preSus) std(postSusmn)/length(postSus)],'r','linewidth',2)
%axis([0 3 -2 2]); 
axis square
set(gca,'xtick',[1 2],'xticklabel',{'pre','post'})

preOffmn = squeeze(nanmean(preOff,2)); postOffmn = squeeze(nanmean(postOff,2));
subplot(1,4,4);hold on
plot([1 2],[preOffmn postOffmn],'k-')
errorbar([1 2],[mean(preOffmn) mean(postOffmn)],[std(preOffmn)/length(preOff) std(postOffmn)/length(postOff)],'r','linewidth',2)
%axis([0 3 -2 2]); 
axis square
set(gca,'xtick',[1 2],'xticklabel',{'pre','post'});





if savePDF
     set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
            end


end


%%

if savePDF
 [f p] = uiputfile('*.pdf','pdf name');
   ps2pdf('psfile', psfilename, 'pdffile', fullfile(p,f));
   delete(psfilename);
end




