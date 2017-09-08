useResp = peakresp(:,1,1)>2 | peakresp(:,1,2)>2;
use = goodAll & ~inhAll & hasDrift==1 & useResp' ;%which cells to include
useInh = goodAll & inhAll==1 & hasDrift==1 & useResp';

t=2

for i = 1:4

    figure
    clear mnpre mnpost prenorm postnorm normrange
    if i==1
        data = find(use & layerAll==5 & treatment==t);
        set(gcf,'Name','grating lyr5 stat');
   
    elseif i==2
        clear data
        set(gcf,'Name','grating lyr 2/3 stat');
        data = find(use & (layerAll==2 |layerAll==3) & treatment==t);
  
       
    else i==3
        clear data
        set(gcf,'Name','grating inh stat');
        data = find(useInh  & treatment==t);
 
end
         for c = 1:length(data)
            %normalize response of all cells used, then average
            prerange = squeeze(tcourse_pref(data(c),1,1,:))
            normrange(c,:)= max(prerange) %- min(tcourse_pref(data(c),1,1,:),[],1); %set range w/pre data to normalize to
            prenorm(c,:) =(tcourse_pref(data(c),1,1,:)- min(tcourse_pref(data(c),1,1,:)))...
                ./normrange(c,1);
            postnorm(c,:) = (tcourse_pref(data(c),1,2,:)- min(tcourse_pref(data(c),1,2,:)))...
                ./normrange(c,1);
         %   prenorm(isnan(prenorm))=0; postnorm(isnan(postnorm)) = 0 ;
            prenorm(isinf(prenorm))=nan; postnorm(isinf(postnorm)) = nan ;

            mnpre= nanmean(prenorm,1);mnpre = circshift(mnpre,[0 10]);
            mnpost = nanmean(postnorm,1);mnpost = circshift(mnpost,[0 10]);
            
            sempre = nanstd(prenorm)/sqrt(length(data));sempre = circshift(sempre,[0 10]);
            sempost = nanstd(postnorm)/sqrt(length(data));sempost = circshift(sempost,[0 10]);
            ncells = length(data)
         end
       
    x = ((1:length(mnpre)-5)*dt -dt/2); x2 = ((1:length(mnpost)-5)*dt -dt/2);
    preerr=shadedErrorBar(x,mnpre(1:45),sempre(1:45)','b',1); hold on;
    posterr=shadedErrorBar(x2,mnpost(1:45),sempost(1:45)','r',1);
    axis square; set(gca,'fontsize', 18); hold on; %xlim([0 45])
    title(titles{t}); xlabel('time (s)'); ylabel('spikes/sec');
    text(2 ,.5,max(max((mnpre))+.75), ['n = ' num2str(ncells)],'FontSize',20); 
end

%%
figure;
for c = 1:length(normrange)

 subplot(6,5,c)
plot(prenorm(c,:)); hold on; plot(postnorm(c,:));axis square
end

%%
clear useResp 
peakresp = squeeze(max(drift_orient,[],2))-squeeze(drift_spont);%subtract driftspont
%peakresp = squeeze(mean(drift_orient,2))-squeeze(drift_spont);%subtract driftspont

%useResp = peakresp(:,1,1)>(min(peakresp(:,1,1)+4)) | peakresp(:,1,2)>(min(peakresp(:,1,2)+4));
%useResp = peakresp(:,1,1)>0 | peakresp(:,1,2)>0;
%use = goodAll & ~inhAll & hasDrift==1 & useResp' ;%which cells to include
%useInh = goodAll & inhAll==1 & hasDrift==1 & useResp';
% 
ampmean = squeeze(mean(peakresp,2))
useResp = ampmean(:,1)>2|ampmean(:,2)>2;

use = goodAll &  inhAll==0 & hasDrift==1& useResp' ;%which cells to include
useInh =goodAll & inhAll==1 & hasDrift==1& useResp';
% %dataAll =goodAll & useResp' & hasDrift==1 ;


%use = goodAll & ~inhAll & hasDrift==1 & useResp' ;%which cells to include
%useInh = goodAll & inhAll==1 & hasDrift==1 & useResp';
layers = {'Layer 2/3','Layer 4','Layer 5','Inh. Units'};
for t=1:2
    figure
    if t==1, set(gcf,'Name','Saline'); else set(gcf,'Name','DOI')
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
            
          %  mnpre = mnpre - repmat(mnpre(1,:),[50 1]);% mnpost = mnpost - repmat(mnpost,[50 1]);

            mnpre= nanmean(prenorm,1); %mnpre = mnpre - repmat(mnpre(:,1),[1 50]);
            mnpre = circshift(mnpre,[0 10]);
            mnpost = nanmean(postnorm,1); %mnpost = mnpost - repmat(mnpost(:,1),[1 50]);
            mnpost = circshift(mnpost,[0 10]); 
            %mnpre = mnpre - repmat(mnpre(1,:),[50 1]);% mnpost = mnpost - repmat(mnpost,[50 1]);
            sempre = nanstd(prenorm)/sqrt(length(data));sempre = circshift(sempre,[0 10]);
            sempost = nanstd(postnorm)/sqrt(length(data));sempost = circshift(sempost,[0 10]);
            ncells = length(data)
        end
        
     %   mnpre = mnpre - repmat(mnpre,[50 1]); mnpost = mnpost - repmat(mnpost,[50 1]);

        x = ((1:length(mnpre)-5)*dt -dt/2); x2 = ((1:length(mnpost)-5)*dt -dt/2);
        subplot(1,4,i)
        preerr=shadedErrorBar(x,mnpre(1:45),sempre(1:45)','b',1); hold on;
        posterr=shadedErrorBar(x2,mnpost(1:45),sempost(1:45)','r',1);%tightfig;
        axis square; set(gca,'fontsize', 18); hold on; %ylim([0 1]);
        xlim([0 2.5]);
        % title(titles{t});
        xlabel('time (s)'); %ylabel('spikes/sec');
        text(1.5 ,.85,max(max((mnpre))+.75), ['n = ' num2str(ncells)],'FontSize',20);
        title(layers{i});
    end
end
%% stats test for temporal components
layers = {'Layer 2/3','Layer 4','Layer 5','Inh. Units'};

for t=1:2
    figure
for i=1:4
    if t==1, set(gcf,'Name','Saline'); else set(gcf,'Name','DOI')
    end
        % figure
        clear mnpre mnpost prenorm postnorm normrange ncells pre post
    clear preTrans postTrans preSus postSus preSpont postSpont 
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
            normrange(c,:)= mean(prerange)% - min(mnPsth(data(c),1,:),[],1); %set range w/pre data to normalize to
            
            prenorm(c,:) =(mnPsth(data(c),1,:)- min(mnPsth(data(c),1,:)))...
                ./normrange(c,1);
            postnorm(c,:) = (mnPsth(data(c),2,:)- min(mnPsth(data(c),2,:)))...
                ./normrange(c,1);
            prenorm(isinf(prenorm))=nan; postnorm(isinf(postnorm)) = nan ;
            pre= circshift(prenorm,[0 10]); post = circshift(postnorm,[0 10]);        
            %              preTrans= pre(c,5:15); postTrans= post(c,5:15);
            %              preSus= pre(c,16:31); postSus = post(c,16:31);
            preSpont(c,:) = pre(c,1:10);postSpont(c,:)=post(c,1:10);
            preTrans(c,:) = pre(c,11:14); postTrans(c,:)= post(c,11:14);
            preSus(c,:)= pre(c,15:39); postSus(c,:) = post(c,15:39);
            [hSpont p] = ttest(preSpont,postSpont);
            [hTrans p] = ttest(preTrans,postTrans);
            [hSus p] = ttest(preSus,postSus);
            %%0 = reject null if no difference
               ncells = length(data);
%            corrTrans = corrcoef(preTrans,postTrans,'rows','pairwise'); corrTrans = corrTrans(2,1);
%            corrSus = corrcoef(preSus,postSus,'rows','pairwise'); corrSus = corrSus(2,1);

        end
        subplot(3,4,i)
        spont = hist(hSpont,0:1:1); bSp=bar(0:1:1,spont/length(hSpont));
        ylim([0 1]); xlim([-.50 1.5]);axis square
        title(layers{i});%ylabel('proportion of cells');
        set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
        
        subplot(3,4,i+4)
        
        trans = hist(hTrans,0:1:1); bT=bar(0:1:1,trans/length(hTrans));
        ylim([0 1]); xlim([-.50 1.5]);axis square
        title(layers{i});%ylabel('proportion of cells');
        set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
        
        subplot(3,4,i+8)
        sus = hist(hSus,0:1:1);
        bS=bar(0:1:1,sus/length(hSus)); ylim([0 1]);xlim([-.5 1.5]);
        axis square;set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
        title(layers{i});%ylabel('proportion of cells');

 end
end

%%
peakresp = squeeze(max(drift_orient,[],2))%-squeeze(drift_spont);%subtract driftspont

ampmean = squeeze(mean(peakresp,2))
useResp = ampmean(:,1)>2&ampmean(:,2)>2;
% use = goodAll &  inhAll==0 & hasDrift==1& useResp' 
% useInh = goodAll &  inhAll==1 & hasDrift==1& useResp' 
layers = {'Layer 2/3','Layer 4','Layer 5','Inh. Units'};


%useResp = peakresp(:,1,1)>0| peakresp(:,1,2)>0;
useOsi_pre = drift_osi(:,1,1)>.2; useOsi_post = drift_osi(:,1,2)>.2; useOsi=(useOsi_pre & useOsi_post)'
use = goodAll &  inhAll==0 & hasDrift==1& useResp' &useOsi==1;
useInh = goodAll &  inhAll==1 & hasDrift==1& useResp' &useOsi==1;

clear used sempre
for t=1:2
    figure
    if t==1, set(gcf,'Name','Saline'); else set(gcf,'Name','DOI'), end
    for i = 1:4
        % figure
        clear mnpre mnpost prenorm postnorm normrange ncells preS postS
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
            mnpre = squeeze(nanmean(pre,1));
            mnpost = squeeze(nanmean(post,1));
            
            sempre = nanstd(pre)/sqrt(length(data));
            sempost = nanstd(post)/sqrt(length(data));
            ncells = length(data);
            preS = circshift(pre,[0 10]); postS= circshift(post,[0 10]);
            
            preSpont(c,:) = preS(c,1:10);postSpont(c,:)=postS(c,1:10);
            preTrans(c,:) = preS(c,11:14); postTrans(c,:)= postS(c,11:14);
            preSus(c,:)= preS(c,15:38); postSus(c,:) = postS(c,15:38);
            preOff(c,:)= preS(c,39:46); postOff(c,:) = postS(c,39:46);

            [hSpont p] = ttest(preSpont,postSpont);
            [hTrans p] = ttest(preTrans,postTrans);
            [hSus p] = ttest(preSus,postSus);
            [hOff p] = ttest(preOff,postOff);

        end
        mnpre = mnpre - repmat(mnpre(:,1),[1 50]);
        mnpre = circshift(mnpre,[0 10]);
        mnpost = mnpost - repmat(mnpost(:,1),[1 50]);
        mnpost = circshift(mnpost,[0 10]);
        sempre = sempre - repmat(sempre(:,1),[1 50]);
        sempre = circshift(sempre,[0 10]);
        sempost = sempost - repmat(sempost(:,1),[1 50]);
        sempost = circshift(sempost,[0 10]);
        
        subplot(5,4,i)
        x = ((1:length(mnpre)-5)*dt -dt/2); x2 = ((1:length(mnpost)-5)*dt -dt/2);
        preerr=shadedErrorBar(x,mnpre(1:45),sempre(1:45)','b',1); hold on;
        posterr=shadedErrorBar(x2,mnpost(1:45),sempost(1:45)','r',1);
        axis square; set(gca,'fontsize', 18); hold on; xlim([0 2.5]);
        title(layers{i}); xlabel('time (s)'); ylabel('spikes/sec');
        text(1.75,max(max((mnpre))+.75), ['n = ' num2str(ncells)],'FontSize',20);
        
        subplot(5,4,i+4)
        spont = hist(hSpont,0:1:1); bSp=bar(0:1:1,spont/length(hSpont));
        ylim([0 1]); xlim([-.50 1.5]);axis square
        title(layers{i});%ylabel('proportion of cells');
        set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
        axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
        
        subplot(5,4,i+8)  
        trans = hist(hTrans,0:1:1); bT=bar(0:1:1,trans/length(hTrans));
        ylim([0 1]); xlim([-.50 1.5]);axis square
        title(layers{i});%ylabel('proportion of cells');
        set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
        axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
        
        subplot(5,4,i+12)
        sus = hist(hSus,0:1:1);
        bS=bar(0:1:1,sus/length(hSus)); ylim([0 1]);xlim([-.5 1.5]);
        axis square;set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
        axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
        
      subplot(5,4,i+16) 
     off = hist(hOff,0:1:1);
     bS=bar(0:1:1,off/length(hOff)); ylim([0 1]);xlim([-.5 1.5]);
     axis square;set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
     axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
     
        
        title(layers{i});%ylabel('proportion of cells');
end
end

%%

layers = {'Layer 2/3','Layer 4','Layer 5','Inh. Units'};
for t=1:2
    figure 
    if t==1, set(gcf,'Name','Saline'); else set(gcf,'Name','DOI')
    end
    for i = 1:4
        clear mnpre mnpost prenorm postnorm normrange ncells pre post
        clear preTrans postTrans preSus postSus preSpont postSpont
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
            pre(c,:) =(mnPsth(data(c),1,:)); post(c,:) =(mnPsth(data(c),2,:));
            
            mnpre = squeeze(nanmean(pre,1));
            mnpost = squeeze(nanmean(post,1));
            
            sempre = nanstd(pre)/sqrt(length(data));
            sempost = nanstd(post)/sqrt(length(data));
            ncells = length(data);
            %  pre = pre - repmat(pre(1,:),[50 1]);%mnpre = circshift(mnpre,[0 10]);  
            
            mnpre = mnpre - repmat(mnpre(:,1),[1 50]);
            mnpre = circshift(mnpre,[0 10]);
            mnpost = mnpost - repmat(mnpost(:,1),[1 50]);
            mnpost = circshift(mnpost,[0 10]);
            sempre = sempre - repmat(sempre(:,1),[1 50]);
            sempre = circshift(sempre,[0 10]);
            sempost = sempost - repmat(sempost(:,1),[1 50]);
            sempost = circshift(sempost,[0 10]);
            %subplot(1,4,i) 
            % pre = pre - repmat(pre(c,:),[50 1]);
            pre = circshift(pre,[0 10]); post= circshift(post,[0 10]);    
            % %              preTrans= pre(c,5:15); postTrans= post(c,5:15);
            % %              preSus= pre(c,16:31); postSus = post(c,16:31);
            preSpont(c,:) = pre(c,2:10);postSpont(c,:)=post(c,2:10);
            preTrans(c,:) = pre(c,11:14); postTrans(c,:)= post(c,11:14);
            preSus(c,:)= pre(c,15:38); postSus(c,:) = post(c,15:38);
            preOff(c,:)= preS(c,39:46); postOff(c,:) = postS(c,39:46);

            [hSpont p] = ttest(preSpont,postSpont);
            [hTrans p] = ttest(preTrans,postTrans);
            [hSus p] = ttest(preSus,postSus);
            [hOff p] = ttest(preOff,postOff);
            %             %%0 = reject null if no difference
            ncells = length(data);
            %            corrTrans = corrcoef(preTrans,postTrans,'rows','pairwise'); corrTrans = corrTrans(2,1);
            %            corrSus = corrcoef(preSus,postSus,'rows','pairwise'); corrSus = corrSus(2,1);
        end
        
        subplot(4,4,i)
        spont = hist(hSpont,0:1:1); bSp=bar(0:1:1,spont/length(hSpont));
        ylim([0 1]); xlim([-.50 1.5]);axis square
        title(layers{i});%ylabel('proportion of cells');
        set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
        axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
        
        subplot(4,4,i+4)
        
        trans = hist(hTrans,0:1:1); bT=bar(0:1:1,trans/length(hTrans));
        ylim([0 1]); xlim([-.50 1.5]);axis square
        title(layers{i});%ylabel('proportion of cells');
        set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
        axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
        
        subplot(4,4,i+8)
        sus = hist(hSus,0:1:1);
        bS=bar(0:1:1,sus/length(hSus)); ylim([0 1]);xlim([-.5 1.5]);
        axis square;set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
        axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
        
        title(layers{i});%ylabel('proportion of cells');
        
    subplot(4,4,i+12) 
     off = hist(hOff,0:1:1);
     bS=bar(0:1:1,off/length(hOff)); ylim([0 1]);xlim([-.5 1.5]);
     axis square;set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
     axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
     
    title(layers{i});%ylabel('proportion of cells');

    end
end


%% orthogonal psth

peakresp = squeeze(max(drift_orient,[],2))%-squeeze(drift_spont);%subtract driftspont
useResp = peakresp(:,1,1)>0 | peakresp(:,1,2)>0;
useOsi_pre = drift_osi(:,1,1)>.3;
useOsi_post = drift_osi(:,1,2)>.3; useOsi=(useOsi_pre & useOsi_post)'

use = goodAll &  inhAll==0 & hasDrift==1& useResp' &useOsi==1;
useInh = goodAll &  inhAll==1 & hasDrift==1& useResp' &useOsi==1;
layers = {'Layer 2/3','Layer 4','Layer 5','Inh. Units'};

clear used sempre
for t=1:2
    figure
    if t==1, set(gcf,'Name','Saline'); else set(gcf,'Name','DOI')
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
            pre(c,:) =(tcourse_orth(data(c),1,1,:));post(c,:) =(tcourse_orth(data(c),1,2,:));
            mnpre = squeeze(nanmean(pre,1));
            mnpost = squeeze(nanmean(post,1));
            
            sempre = nanstd(pre)/sqrt(length(data));
            sempost = nanstd(post)/sqrt(length(data));
            ncells = length(data);
            preS = circshift(pre,[0 10]); postS= circshift(post,[0 10]);
            
            preSpont(c,:) = preS(c,2:10);postSpont(c,:)=postS(c,2:10);
            preTrans(c,:) = preS(c,11:14); postTrans(c,:)= postS(c,11:14);
            preSus(c,:)= preS(c,15:38); postSus(c,:) = postS(c,15:38);
            preOff(c,:)= preS(c,39:46); postOff(c,:) = postS(c,39:46);

            [hSpont p] = ttest(preSpont,postSpont);
            [hTrans p] = ttest(preTrans,postTrans);
            [hSus p] = ttest(preSus,postSus);
            [hOff p] = ttest(preOff,postOff);

        end
        mnpre = mnpre - repmat(mnpre(:,1),[1 50]);
        mnpre = circshift(mnpre,[0 10]);
        mnpost = mnpost - repmat(mnpost(:,1),[1 50]);
        mnpost = circshift(mnpost,[0 10]);
        sempre = sempre - repmat(sempre(:,1),[1 50]);
        sempre = circshift(sempre,[0 10]);
        sempost = sempost - repmat(sempost(:,1),[1 50]);
        sempost = circshift(sempost,[0 10]);
       
        subplot(5,4,i)
        x = ((1:length(mnpre)-5)*dt -dt/2); x2 = ((1:length(mnpost)-5)*dt -dt/2);
        preerr=shadedErrorBar(x,mnpre(1:45),sempre(1:45)','b',1); hold on;
        posterr=shadedErrorBar(x2,mnpost(1:45),sempost(1:45)','r',1);
        axis square; set(gca,'fontsize', 18); hold on; xlim([0 2.5]);
        title(layers{i}); xlabel('time (s)'); ylabel('spikes/sec');
        text(1.75,max(max((mnpre))+.75), ['n = ' num2str(ncells)],'FontSize',20);
        
        subplot(5,4,i+4)
        spont = hist(hSpont,0:1:1); bSp=bar(0:1:1,spont/length(hSpont));
        ylim([0 1]); xlim([-.50 1.5]);axis square
        title(layers{i});%ylabel('proportion of cells');
        set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
        axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);

     subplot(5,4,i+8) 

     trans = hist(hTrans,0:1:1); bT=bar(0:1:1,trans/length(hTrans)); 
     ylim([0 1]); xlim([-.50 1.5]);axis square
     title(layers{i});%ylabel('proportion of cells');
     set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
     axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);

     subplot(5,4,i+12) 
     sus = hist(hSus,0:1:1);
     bS=bar(0:1:1,sus/length(hSus)); ylim([0 1]);xlim([-.5 1.5]);
     axis square;set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
     axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);

     subplot(5,4,i+16) 
     off = hist(hOff,0:1:1);
     bS=bar(0:1:1,off/length(hOff)); ylim([0 1]);xlim([-.5 1.5]);
     axis square;set(gca, 'XTick', [], 'XTickLabel', [],'FontSize',22);
     axis square;set(gca, 'YTick', [], 'YTickLabel', [],'FontSize',22);
     
    title(layers{i});%ylabel('proportion of cells');
end
end

%%

ampmean = squeeze(mean(peakresp,2))
useResp = ampmean(:,1)>2 &ampmean(:,2)>2;
use = goodAll &  inhAll==0 & hasDrift==1& useResp' ;%which cells to include
useN = find(use & (layerAll==4) & treatment==t);
figure
%useN= find(useResp' & goodAll==1 & hasDrift & treatment==doi &layerAll==2)
length(useN)
for c=1:length(useN)
   % set (gcf,'Name','Lyr 2/3 psth/unit') % set up in loop
    subplot(11,11,c)
    plot(squeeze(mnPsth(useN(c),1,:)));hold on;
    plot(squeeze(mnPsth(useN(c),2,:)),'r')
end