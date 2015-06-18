%%% get OSI

%%% get transience

%%% get on/off

%%% get burst/ firing pattern

%%% get crf

%%% compare crf, bursting with running


%%% get grating response with speed

if ~exist('drift_mv_all','var')
    load ('adultData061715.mat')
end
if ~exist('moviedata','var')
    load('C:\wn016alpha1_10hzLg60Hz.mat')
end
%%% main figsspi

%%% supp by contrast

%%% gain across population

%%% bursting across population

%%%for c = 1:length(movement_all);
clear spikebins spikepatterns
close all
moviedata = double(moviedata);
m= (moviedata-mean(moviedata(:)))/mean(moviedata(:));


for c = 1:length(movement_all)
    tic
    c
    
    if ~isempty(movement_all(c).speed)
        spikehist = histc(movement_all(c).spikes,movement_all(c).frameT)*60;
        sp = interp1(movement_all(c).mvlfp_tdata, movement_all(c).speed,movement_all(c).frameT);
        ph= spikehist.*exp(2*pi*sqrt(-1)*movement_all(c).frameNum/600);
        ph = conv(ph,ones(2401,1),'same')./conv(spikehist,ones(2401,1),'same');
        
        %     figure
        %     plot(abs(ph));
        %     hold on
        %     plot(sp/max(sp),'g')
        %
        %     figure
        %     plot(mod(angle(ph),2*pi));
        %     hold on
        %     plot(sp/max(sp),'g')
        
        cycframe = mod(movement_all(c).frameNum,600);
        
        for i = 1:20
            fms = find(cycframe>(i-1)*30+1 & cycframe<i*30 & sp<1);
            length(fms);
            crf(1,i) =   mean(spikehist(fms));
            crferr(1,i) = std(spikehist(fms))/sqrt(length(fms));
            fms = find(cycframe>(i-1)*30+1 & cycframe<i*30 & sp>1);
            crf(2,i) =   mean(spikehist(fms));
            crferr(2,i) = std(spikehist(fms))/sqrt(length(fms));
        end
        figure
        subplot(2,3,1)
        hold on
        errorbar(1:10,crf(1,1:10),crferr(1,1:10),'r','Linewidth',1.5); %plot(crf(1,20:-1:11),'r:','Linewidth',2);
        errorbar(1:10,crf(2,1:10),crferr(2,1:10),'g','Linewidth',1.5);  %plot(crf(2,20:-1:11),'g:','Linewidth',2);
        ylim([0 max(crf(:))])
        xlim([0 10.5])
        
        crfAll(c,:,:) = crf;
        spike_sp = interp1(movement_all(c).mvlfp_tdata, movement_all(c).speed,movement_all(c).spikes);
        
        pre_isi = diff(movement_all(c).spikes(1:end-1))*1000;
        post_isi = diff(movement_all(c).spikes(2:end))*1000;
        
        %     subplot(2,2,2)
        %     figure
        %     use =find(spike_sp<1);
        %     if length(use)>2000
        %         use=use(2:2000);
        %     else use=use(2:end-1);
        %     end
        %     plot(log10(pre_isi(use)),log10(post_isi(use)),'r.','Markersize',12)
        %     hold on
        %       use = find(spike_sp>1);
        %     if length(use)>2000
        %         use=use(2:2000);
        %     else use=use(2:end-1);
        %     end
        %     plot(log10(pre_isi(use)),log10(post_isi(use)),'b.','Markersize',12)
        
        
        
        %%% calculate sta for all spikes, then get transfer function for
        %%% moving vs stationary
        
        histbins{1} = 0:0.25:4;
        histbins{2}=histbins{1};
        
        h=  hist3([log10(pre_isi(spike_sp(2:end-1)<1)) ;log10(post_isi(spike_sp(2:end-1)<1))]','edges',histbins);
        
        spikebins(1,:,:) = h/sum(h(:));
        
        
        
        h = hist3([log10(pre_isi(spike_sp(2:end-1)>1)) ;log10(post_isi(spike_sp(2:end-1)>1))]','edges',histbins);
        
        spikebins(2,:,:) = h/sum(h(:));
        
        subplot(2,3,2);
        imagesc(squeeze(spikebins(1,:,:))); axis square
        subplot(2,3,3)
        imagesc(squeeze(spikebins(2,:,:))); axis square
        
        spikebinsAll(c,:,:,:) = spikebins;
        
        %     spikehist = hist(movement_all(c).spikes,0:0.01:max(movement_all(c).spikes));
        %         figure
        %         subplot(2,1,1)
        %         spect = abs(fft(spikehist));
        %         spect = spect(1:round(end/2));
        %         loglog(spect);
        %
        %         lastT = max(movement_all(c).lfpT);
        %         lfpAll = zeros(2,length(movement_all(c).spikes),501);
        %         lfpAll(:) = NaN;
        %
        %         c
        %         for i = 1:length(movement_all(c).spikes);
        %
        %           spk = movement_all(c).spikes(i);
        %             if spk>1 & spk<lastT-1 & spike_sp(i)<1;
        %              inds = find(movement_all(c).lfpT>spk);
        %
        %               lfpAll(1,i,:) = movement_all(c).lfpV(inds(1)-250:inds(1)+250);
        %             elseif spk>1 & spk<lastT-1 & spike_sp(i)>1;
        %               inds = find(movement_all(c).lfpT>spk);
        %               lfpAll(2,i,:) = movement_all(c).lfpV(inds(1)-250:inds(1)+250);
        %             end
        %         end
        %
        %
        %         subplot(2,1,2)
        %         plot(-250:250,squeeze(nanmean(lfpAll,2))')
        %         xlim([-250 250])
        %
        
        %
        %     figure
        %     plot(log10(pre_isi(1:200)),log10(post_isi(1:200)));
        %     hold on
        %     plot(log10(pre_isi(1:200)),log10(post_isi(1:200)),'.','Markersize',8);
        %
        
        
        % figure
        % hold on
        % plot(log10(pre_isi(1:500)),log10(post_isi(1:500)),'.','Markersize',8);
        %
        % for i = 1:100;
        %
        %    h= plot(log10(pre_isi(i+3)),log10(post_isi(i+3)),'g.','Markersize',24);
        %
        %     plot(log10(pre_isi(i:i+3)),log10(post_isi(i:i+3)));
        %     axis([0 3.5 0 3.5])
        %     mov(i)=getframe(gcf);
        %     delete(h)
        % end
        % figure
        % movie(mov)
        
        
        sta_all= zeros(size(moviedata,1),size(moviedata,2),9);
        %figure
        fnum=movement_all(c).frameNum;
        for lag = 1:9
            for i = (lag+1):length(spikehist);
                sta_all(:,:,lag) =sta_all(:,:,lag)+spikehist(i)*moviedata(:,:,fnum(i-lag));
            end
            
            
            
            %         subplot(3,3,lag);
            %         imagesc(sta,[-0.1 0.1])
            %         axis equal
        end
        
        sta_all = sta_all/sum(spikehist);
        sta_all = (sta_all-128)/128;
        
        obs = reshape(sta_all,size(sta_all,1)*size(sta_all,2),size(sta_all,3));
        
        
        [u s v] = svd(obs');
        
        meanSTA= reshape(v(:,1),size(sta_all,1),size(sta_all,2));
        tcourse = u(:,1);
        if sum(tcourse(1:5))<0;
            tcourse = -tcourse;
            meanSTA=-meanSTA;
        end
        meanSTA = meanSTA*s(1,1)*max(tcourse);
        tcourse = tcourse/max(tcourse);
        
        meanSTAall(c,:,:) = meanSTA;
        
        
        subplot(2,3,4);
        imagesc(meanSTA,[-0.1 0.1]); axis square; axis off
        subplot(2,3,5);
        plot(tcourse); hold on; plot([1 9],[0 0],':'); ylim([-1.05 1.05])
        tcAll(c,:) = tcourse;
        
        
        lag=4;
        %     for i = (lag+1):length(spikehist);
        %         sta =sta+spikehist(i)*moviedata(:,:,movement_all(c).frameNum(i-lag));
        %     end
        %
        %     sta = sta/sum(spikehist);
        %    sta = (sta-128)/128;
        sta = meanSTA;
        %cutoff = max(0.06,max(abs(sta(:)))*0.66);
        cutoff = max(abs(sta(:)))*0.75;
        sta(abs(sta)<cutoff)=0;
        subplot(2,3,4);
        %  imagesc(sta,[-0.1 0.1]); axis square; axis off
        
        staxc=zeros(size(spikehist));
        for i = (lag+1):length(spikehist);
            staxc(i) = sum(sum((m(:,:,movement_all(c).frameNum(i-lag))).*sta));
        end
        staxc = staxc/(sum(abs(sta(:))));
        %         figure
        %         plot(staxc,spikehist,'o');
        
        
        clear xfer err xfer_mv err_mv xfer_stop err_stop
        for i = 1:11;
            matchsp = find(staxc>(i-6.5)*.2 & staxc<(i-5.5)*.2);
            xfer(i)=mean(spikehist(matchsp));
            err(i)=std(spikehist(matchsp))/sqrt(length(matchsp));
            
            matchsp = find(staxc>(i-6.5)*.2 & staxc<(i-5.5)*.2 & sp>1);
            xfer_mv(i)=mean(spikehist(matchsp));
            err_mv(i)=std(spikehist(matchsp))/sqrt(length(matchsp));
            
            matchsp = find(staxc>(i-6.5)*.2 & staxc<(i-5.5)*.2 & sp<1);
            xfer_stop(i)=mean(spikehist(matchsp));
            err_stop(i)=std(spikehist(matchsp))/sqrt(length(matchsp));
        end
        
        subplot(2,3,6)
        errorbar(-10:2:10,xfer,err);
        xlim([-20 20])
        hold on
        plot([0 0],[0 max(xfer)],':')
        errorbar(-10:2:10,xfer_stop,err_stop,'r');
        errorbar(-10:2:10,xfer_mv,err_mv,'g');
        xlim([-10 10])
        if ~isnan(max(xfer_mv))
            ylim([0 max(xfer_mv)+0.1]);
        end
        
        xferAll(c,1,:) = xfer_stop;
        xferAll(c,2,:) = xfer_mv;
        xferSpont(c,1:2) = [xfer_stop(6) xfer_mv(6)];
        xferGain(c,1) = mean([xfer_stop(7)-xfer_stop(6) 0.5*(xfer_stop(8)-xfer_stop(6)) 0.3333*(xfer_stop(9)-xfer_stop(6))]);
        xferGain(c,2) = mean([xfer_mv(7)-xfer_mv(6) 0.5*(xfer_mv(8)-xfer_mv(6)) 0.3333*(xfer_mv(9)-xfer_mv(6))]);
        
        
        R=zeros(8,7,2,3);
        for tf =1:2
            for f=1:3
                R(:,1,tf,f)=drift_all(c).flicker(tf,f);
            end
        end
        R(:,2:7,:,:)=drift_all(c).R;
        R(:,:,:,2:3) = abs(R(:,:,:,2:3));
        osi = squeeze(mean(R(:,2:7,:,:),2)); %osi(:,:,1)=osi(:,:,1)-drift_all(c).spont(1);
        
        sf = squeeze(mean(R,1)); %sf(:,:,1) = sf(:,:,1)-drift_all(c).spont(1);
        dr_spont(c) = drift_all(c).spont(1);
        for tf=1:2
            for f = 1:2
                osi_weight = osi(:,tf,f);
                osi_weight = osi_weight/max(osi_weight);
                osi_weight = 8/sum(osi_weight);
                sf_weight = sf(2:7,tf,f);
                sf_weight = sf_weight/max(sf_weight);
                sf_weight = 6/sum(sf_weight);
                osi(:,tf,f) = osi(:,tf,f)*sf_weight;
                sf(:,tf,f) = sf(:,tf,f)*osi_weight;
            end
        end
        %
        tf_tune = mean(mean(osi(:,:,1:2),3),1);
        [mx tf_pref(c)] = max(tf_tune);
        osi_all(c,:,:) = squeeze(osi(:,tf_pref(c),1:2));
        sf_all(c,:,:) = squeeze(sf(:,tf_pref(c),1:2));
        
        figure  %%% tf,F1,os
        for f = 1:2
            subplot(2,2,f);
            plot(squeeze(osi(:,1:2,f)));
            hold on
            plot([1 7],[drift_all(c).spont(1) drift_all(c).spont(1)],':')
            lower =min(min(osi(:,:,f)));
            upper = max(max(osi(:,:,f))); upper = max(upper, drift_all(c).spont(1));
            ylim([0 1.1*max(upper,5)]);xlim([0.5 8.5])
        end
        
        for f = 1:2
            subplot(2,2,f+2);
            plot(squeeze(sf(:,1:2,f))); hold on
            plot([1 7],[drift_all(c).spont(1) drift_all(c).spont(1)],':')
            lower =min(min(sf(:,:,f)));
            upper = max(max(sf(:,:,f)));upper = max(upper, drift_all(c).spont(1));
            ylim([0 1.1*max(upper,4)]); xlim([0.5 7.5])
        end
        
    end
    
    close all
    toc
end



keyboard

figure
plot(wvform')

staUse = find(max(abs(reshape(meanSTAall,size(meanSTAall,1),size(meanSTAall,2)*size(meanSTAall,3))),[],2)>0.06);
%%%staUse=staUse(staUse<143)


figure
plot(xferSpont(:,1),xferSpont(:,2),'o')
axis([-5 40 -5 40])
axis square
hold on
plot([0 40],[0 40])
title('wn spont')

figure
plot(dr_spont_mv(:,1,1),dr_spont_mv(:,2,1),'o')
axis([-5 40 -5 40])
axis square
hold on
plot([0 40],[0 40])
title('drift spont')

figure
plot(dr_spont_mv(:,1,1),xferSpont(:,1),'o')
axis([-5 40 -5 40])
axis square
hold on
plot([0 40],[0 40])
xlabel('drift stop spont'); ylabel('wnstop spont')

figure
plot(dr_spont_mv(:,2,1),xferSpont(:,2),'o')
axis([-5 40 -5 40])
axis square
hold on
plot([0 40],[0 40])
xlabel('drift mv spont'); ylabel('wn mv spont')




figure
plot(crfAll(:,1,1),crfAll(:,2,1),'o')
axis([0 40 0 40])
axis square
hold on
plot([0 40],[0 40])

figure
plot(xferGain(:,1),xferGain(:,2),'o')
%axis([ 20 0 20])
axis square
hold on
plot([0 20],[0 20])

sbcGain = xferAll(:,:,7)-xferAll(:,:,6);
figure
plot(sbcGain(:,1),sbcGain(:,2),'o')
%axis([ 20 0 20])
axis square
hold on
plot([0 20],[0 20])
plot([-5 20],[0 0]);
plot([0 0],[-5 20])

sbc = find(sbcGain(:,2)<-2 | sbcGain(:,1)<-2);
%sbc = find(xferGain(:,2)<-1 | xferGain(:,1)<-1);
for i = 1:length(sbc)
    figure
    subplot(1,2,1); hold on
    plot(squeeze(xferAll(sbc(i),1,:)),'r')
    plot(squeeze(xferAll(sbc(i),2,:)),'g')
    xlim([1 11])
    
    subplot(1,2,2);hold on
    plot(squeeze(crfAll(sbc(i),1,1:10)),'r')
    plot(squeeze(crfAll(sbc(i),2,1:10)),'g')
    xlim([1 10])
end




for c = 1:size(meanSTAall,1)
    sta= squeeze(meanSTAall(c,:,:)); sta=sta(:);
    [val ind] = max(abs(sta));
    onoff(c)=sign(sta(ind));
end

figure
plot(tcAll(staUse,:)')

figure
hist(min(tcAll(staUse,6:9),[],2)-tcAll(staUse,1))

[coeff score latent] = pca(tcAll(staUse,:));
figure
plot(score(:,1),score(:,2),'o')

id = kmeans(score(:,1:2),2)
figure
plot(score(id==1,1),score(id==1,2),'ro');
hold on
plot(score(id==2,1),score(id==2,2),'go');

tcGroup = zeros(size(score,1),1)
%tcGroup(score(:,1)>0 & score(:,2)>-0.5)=1;
tcGroup(score(:,1)>=score(:,2))=1;
%tcGroup(score(:,1)<0 & score(:,1)>-1.25)=2;
tcGroup(score(:,2)>score(:,1))=2;
tcGroup(score(:,1)<-1.25)=3;
tcGroup(score(:,1)>-0.5 & score(:,2)<-0.5)=4;
tcGroup(tcAll(staUse,8)>0.2)=5;

col = 'rgbcm';
figure
hold on
for i = 1:5
    plot(tcAll(staUse(tcGroup==i),:)',col(i));
end
figure
hist(tcGroup)
figure
hold on
for i = 1:5
    plot(score(tcGroup==i,1),score(tcGroup==i,2),[col(i) 'o'])
end
hold on
plot([-1 1],[-1 1])



F1F0 = mean(osi_all(:,:,2),2)./(mean(osi_all(:,:,1),2) );
figure
h=hist(F1F0,0:0.1:1); bar(0:0.1:1,h)
hold on
h=hist(F1F0(staUse),0:0.1:1); bar(0:0.1:1,h,'r');
xlabel('F1/F0'); xlim([-0.1 1.1]);

ev_osi = squeeze(osi_all(:,:,1)) - repmat(dr_spont',[1 8]);
figure
plot(ev_osi');



evoked = mean(ev_osi,2);
for i = 1:length(ev_osi)
    if max(ev_osi(i,:))>4
        cvosi(i) = abs(calcCVosi(ev_osi(i,:)));
    else
        cvosi(i) = NaN;
    end
end
figure
hist(cvosi,0:0.05:1)
sum(cvosi>0.2)/length(cvosi)
F0OS = find(cvosi>0.2);
figure
plot(ev_osi(F0OS,:)')


figure
h=hist(F1F0,0:0.1:1); bar(0:0.1:1,h)
hold on
h=hist(F1F0(staUse),0:0.1:1); bar(0:0.1:1,h,'r');
h=hist(F1F0(F0OS),0:0.1:1); bar(0:0.1:1,h,'g');
xlabel('F1/F0'); xlim([-0.1 1.1]);

for i = 1:length(F0OS);
    
    figure
    imagesc(squeeze(meanSTAall(F0OS(i),:,:)),[-0.1 0.1])
    axis equal
    
end



for i = 1:length(F0OS)
    figure
    plot(squeeze(osi_all(F0OS(i),:,:)));
    hold on
    plot([1 8], [dr_spont(F0OS(i)) dr_spont(F0OS(i)) ])
    title(sprintf('osi = %0.2f f1f0 = %0.2f sta = %d',cvosi(F0OS(i)),F1F0(F0OS(i)),ismember(F0OS(i),staUse)))
end


onoffOS = find(cvosi>0.2 & F1F0'<0.4);
for i = 1:length(onoffOS);
    figure
    imagesc(squeeze(meanSTAall(onoffOS(i),:,:)),[-0.1 0.1])
    axis equal
end

for i = 1:length(onoffOS)
    figure
    hold on
    plot(squeeze(crfAll(onoffOS(i),1,1:10))','r')
    plot(squeeze(crfAll(onoffOS(i),2,1:10))','g')
end




for i = 1:5
    nanmean(burst_fraction(staUse(tcGroup==i),:),1)
end

titleList = {'sustain','transient','sbc','DS/OS'};
figure
for i = 1:4
    subplot(2,2,i)
    if i ==1
        d= squeeze(nanmean(xferAll(staUse(tcGroup==1),:,:),1));
    elseif i==2
        d= squeeze(nanmean(xferAll(staUse(tcGroup==2),:,:),1));
    elseif i==3
        d= squeeze(nanmean(xferAll(sbc,:,:),1));
    elseif i==4
        d= squeeze(nanmean(xferAll(onoffOS,:,:),1));
    end
    d = d(:,6:10)
    d = (d-min(d(:)))/(max(d(:)-min(d(:))));
    plot(0:0.2:0.8,d(1,:),'r')
    hold on
    plot(0:0.2:0.8,d(2,:),'b')
    xlim([0 0.8]); ylim([0 1]);
    title(titleList{i},'FontSize',16)
    if i==1
        xlabel('contrast','FontSize',12)
        ylabel('norm response','Fontsize',12)
        
        legend('stationary','moving')
    end
end


figure
plot(sbcGain(staUse(tcGroup==1),1),sbcGain(staUse(tcGroup==1),2),'go','LineWidth',2);
hold on
plot(sbcGain(staUse(tcGroup==2),1),sbcGain(staUse(tcGroup==2),2),'bo','LineWidth',2);
plot(sbcGain(sbc,1),sbcGain(sbc,2),'ro','LineWidth',2);
plot(sbcGain(onoffOS,1),sbcGain(onoffOS,2),'ko','LineWidth',2);
axis square
plot([-10 25],[-10 25])
axis([-10 25 -10 25])
xlabel('response gain stationary')
ylabel('response gain moving')


legend('sustained','transient','sbc','OS')


data = [nanmedian(sbcGain(staUse(tcGroup==1),2)) nanmedian(sbcGain(staUse(tcGroup==1),1)); ...
    nanmedian(sbcGain(staUse(tcGroup==2),2)) nanmedian(sbcGain(staUse(tcGroup==2),1)); ...
    nanmedian(sbcGain(sbc,2)) nanmedian(sbcGain(sbc,1)) ; ...
    nanmedian(sbcGain(onoffOS,2)),nanmedian(sbcGain(onoffOS,1))]

%
% data = [nanmean(sbcGain(staUse(tcGroup==1),2)) nanmean(sbcGain(staUse(tcGroup==1),1)); ...
%     nanmean(sbcGain(staUse(tcGroup==2),2)) nanmean(sbcGain(staUse(tcGroup==2),1)); ...
%     nanmean(sbcGain(sbc,2)) nanmean(sbcGain(sbc,1)) ; ...
%     nanmean(sbcGain(onoffOS,2)),nanmean(sbcGain(onoffOS,1))]

err = [nanstd(sbcGain(staUse(tcGroup==1),2))/7 nanstd(sbcGain(staUse(tcGroup==1),1))/7; ...
    nanstd(sbcGain(staUse(tcGroup==2),2))/7 nanstd(sbcGain(staUse(tcGroup==2),1))/7; ...
    nanstd(sbcGain(sbc,2))/4 nanstd(sbcGain(sbc,1))/4 ; ...
    nanstd(sbcGain(onoffOS,2))/4,nanstd(sbcGain(onoffOS,1))/4]
figure
barweb(data(:,[2 1])*10,err(:,[2 1])*10)
set(gca,'xticklabel',{'sustain','transient','sbc','OS'})
ylabel('median gain (sp/sec)')
legend('stationary','moving')
ylim([-40 60])



figure
plot(xferSpont(staUse(tcGroup==1),1),xferSpont(staUse(tcGroup==1),2),'go','LineWidth',2);
hold on
plot(xferSpont(staUse(tcGroup==2),1),xferSpont(staUse(tcGroup==2),2),'bo','LineWidth',2);
plot(xferSpont(sbc,1),xferSpont(sbc,2),'ro','LineWidth',2);
plot(xferSpont(onoffOS,1),xferSpont(onoffOS,2),'ko','LineWidth',2);
axis square
plot([0 30],[0 30])
axis([0 30 0 30])

legend('sustained','transient','sbc','OS')
xlabel('spont stationary (sp/sec)');
ylabel('spont moving (sp/sec)')

data = [nanmedian(xferSpont(staUse(tcGroup==1),2)) nanmedian(xferSpont(staUse(tcGroup==1),1)); ...
    nanmedian(xferSpont(staUse(tcGroup==2),2)) nanmedian(xferSpont(staUse(tcGroup==2),1)); ...
    nanmedian(xferSpont(sbc,2)) nanmedian(xferSpont(sbc,1)) ; ...
    nanmedian(xferSpont(onoffOS,2)),nanmedian(xferSpont(onoffOS,1))]

%
% data = [nanmean(xferSpont(staUse(tcGroup==1),2)) nanmean(xferSpont(staUse(tcGroup==1),1)); ...
%     nanmean(xferSpont(staUse(tcGroup==2),2)) nanmean(xferSpont(staUse(tcGroup==2),1)); ...
%     nanmean(xferSpont(sbc,2)) nanmean(xferSpont(sbc,1)) ; ...
%     nanmean(xferSpont(onoffOS,2)),nanmean(xferSpont(onoffOS,1))]

err = [nanstd(xferSpont(staUse(tcGroup==1),2))/7 nanstd(xferSpont(staUse(tcGroup==1),1))/7; ...
    nanstd(xferSpont(staUse(tcGroup==2),2))/7 nanstd(xferSpont(staUse(tcGroup==2),1))/7; ...
    nanstd(xferSpont(sbc,2))/4 nanstd(xferSpont(sbc,1))/4 ; ...
    nanstd(xferSpont(onoffOS,2))/4,nanstd(xferSpont(onoffOS,1))/4]
figure
barweb(data(:,[2 1]),err(:,[2 1]))
set(gca,'xticklabel',{'sustain','transient','sbc','OS'})
ylabel('median spont rate (sp/sec)')
legend('stationary','moving')
ylim([0 12])



figure
plot(burst_fraction(staUse(tcGroup==1),2),burst_fraction(staUse(tcGroup==1),1),'go','LineWidth',2);
hold on
plot(burst_fraction(staUse(tcGroup==2),2),burst_fraction(staUse(tcGroup==2),1),'bo','LineWidth',2);
plot(burst_fraction(sbc,2),burst_fraction(sbc,1),'ro','LineWidth',2);
plot(burst_fraction(onoffOS,2),burst_fraction(onoffOS,1),'ko','LineWidth',2);
axis square
plot([0 0.3],[0 0.3])
axis([0 0.15 0 0.3])
legend('sustained','transient','sbc','OS')
xlabel('burst fraction moving');
ylabel('burst fraction stationary')

data = [nanmedian(burst_fraction(staUse(tcGroup==1),2)) nanmedian(burst_fraction(staUse(tcGroup==1),1)); ...
    nanmedian(burst_fraction(staUse(tcGroup==2),2)) nanmedian(burst_fraction(staUse(tcGroup==2),1)); ...
    nanmedian(burst_fraction(sbc,2)) nanmedian(burst_fraction(sbc,1)) ; ...
    nanmedian(burst_fraction(onoffOS,2)),nanmedian(burst_fraction(onoffOS,1))]

%
% data = [nanmean(burst_fraction(staUse(tcGroup==1),2)) nanmean(burst_fraction(staUse(tcGroup==1),1)); ...
%     nanmean(burst_fraction(staUse(tcGroup==2),2)) nanmean(burst_fraction(staUse(tcGroup==2),1)); ...
%     nanmean(burst_fraction(sbc,2)) nanmean(burst_fraction(sbc,1)) ; ...
%     nanmean(burst_fraction(onoffOS,2)),nanmean(burst_fraction(onoffOS,1))]

err = [nanstd(burst_fraction(staUse(tcGroup==1),2))/7 nanstd(burst_fraction(staUse(tcGroup==1),1))/7; ...
    nanstd(burst_fraction(staUse(tcGroup==2),2))/7 nanstd(burst_fraction(staUse(tcGroup==2),1))/7; ...
    nanstd(burst_fraction(sbc,2))/6 nanstd(burst_fraction(sbc,1))/6 ; ...
    nanstd(burst_fraction(onoffOS,2))/6,nanstd(burst_fraction(onoffOS,1))/6]


figure
barweb(data,err)
set(gca,'xticklabel',{'sustain','transient','sbc','OS'})
ylabel('median burst fraction')
legend('moving','stationary')
ylim([0 0.08])

figure
hold on
for i = 1:length(score)
    if ~isnan(burst_fraction(staUse(i),1))
        plot(score(i,1),score(i,2),'o','Color',cmapVar(burst_fraction(staUse(i),1),0,0.25))
    end
end

figure
hold on
for i = 1:5
    plot(crfAll(staUse(tcGroup==i),1,1),crfAll(staUse(tcGroup==i),2,1),[col(i) 'o'])
end
axis([0 40 0 40])
axis square
hold on
plot([0 40],[0 40])


transience =tcAll(:,7);



figure
hold on
for i = 1:5
    plot(xferGain(staUse(tcGroup==i),1),xferGain(staUse(tcGroup==i),2),[col(i) 'o'])
end
axis([-10 20 -10 20])
axis square
hold on
plot([0 20],[0 20])
xlabel('xfer gain')

figure
hold on
for i = 1:5
    plot(sbcGain(staUse(tcGroup==i),1),sbcGain(staUse(tcGroup==i),2),[col(i) 'o'])
end
axis([-10 20 -10 20])
axis square
hold on
plot([0 20],[0 20])
xlabel('sbc gain')


figure
hold on
for i = 1:length(score)
    if ~isnan(sustain(staUse(i)))
        plot(score(i,1),score(i,2),'o','Color',cmapVar(sustain(staUse(i)),0,1))
    end
end


figure
hold on
for i = 1:length(score)
    if ~isnan(sustain(staUse(i)))
        plot(sbcGain(staUse(i),1),sbcGain(staUse(i),2),'o','Color',cmapVar(transience(staUse(i)),-1,1))
    end
end
plot([0 20],[0 20])

figure
plot(transience(staUse),sbcGain(staUse,2)-sbcGain(staUse,1),'o')

figure
hold on
for i = 1:length(score)
    if ~isnan(sustain(staUse(i)))
        plot(tcAll(staUse(i),:),'Color',cmapVar(sustain(staUse(i)),0,1))
    end
end

figure
hold on
for i = 1:length(score)
    if ~isnan(onoff(staUse(i)))
        plot(tcAll(staUse(i),:),'Color',cmapVar(onoff(staUse(i)),-1.1,1))
    end
end

figure
hold on
for i = 1:5
    plot(xferGain(staUse(tcGroup==i),1),xferGain(staUse(tcGroup==i),2),[col(i) 'o'])
end
axis([-10 20 -10 20])
axis square
hold on
plot([0 20],[0 20])


obs = reshape(spikebinsAll(:,:,3:end,3:end),[size(spikebinsAll,1)*2 15*15]);
obs(isnan(obs))=0;
figure
for i = 1:2
    subplot(1,2,i);
    imagesc(squeeze(nanmean(spikebinsAll(:,i,:,:),1)));axis equal; axis off
end

[u s v] = svd(obs);
figure

for i = 1:6
    figure
    if mean(u(:,i))<0
        u(:,i)=-u(:,i);
        v(:,i)=-v(:,i)
    end
    plot(u(1:204,i),u(205:end,i),'o')
    axis equal;hold on; plot([-0.1 0.1],[-0.1 0.1])
    title(sprintf('%d',i))
    hold on
    plot(u(staUse,i),u(staUse+204,i),'go')
    plot(u(sbc,i),u(sbc+204,i),'ro')
    plot(u(onoffOS,i),u(onoffOS+204,i),'ko')
end

figure
for i = 1:9
    subplot(3,3,i);
    imagesc(reshape(v(:,i),[15 15]),[-0.5 0.5])
end

figure
plot(diag(s(1:10,1:10)))

%%% find sbc
%%% find OS/DS
%%% find sust/transient

%%% look at change in gain, spont
%%% look at burst pattern
%%% look at lfp?