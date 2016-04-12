close all

load('analysis_032816.mat')
% if ~exist('drift_mv_all','var')
%     %load ('adultData061715.mat')
%    % load ('adultData090815.mat')
%    load('rundata091215.mat')
% end
% if ~exist('movement_all','var')
%     load ('adultData090815.mat')
% end

if ~exist('moviedata','var')
    load('C:\wn016alpha1_10hzLg60Hz.mat')
end


%bad = [ 288 ... % below LGN

%%%for c = 1:length(movement_all);
clear spikebins spikepatterns
close all
moviedata = double(moviedata);
m= (moviedata-mean(moviedata(:)))/mean(moviedata(:));


%%% sta analysis
if ~exist('sta_all','var')
    
    for c = 1:length(movement_all)
        tic
        c
        
        if ~isempty(movement_all(c).speed)
            spikehist = histc(movement_all(c).spikes,movement_all(c).frameT);
            framedur= [diff(movement_all(c).frameT) 1/60];
            spikehist = spikehist./framedur;
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
            
            
            sta_all= zeros(size(moviedata,1),size(moviedata,2),9);
            %figure
            fnum=movement_all(c).frameNum;
            lagframes= find(spikehist>0);
            lagframes=lagframes(lagframes>10);
            for lag = 1:9
                % for i = (lag+1):length(spikehist);
                for i = lagframes
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
end


%%% plot session quality
psfilename = 'test.ps';
dt=1;
close all
if false
    for sess = 1:max(site)
        units = find(site == sess);
        clear sp_hist
        if ~isempty(movement_all(units(1)).spikes)
            for i=1:length(units)
                sp = movement_all(units(i)).spikes;
                sp_hist(i,:) = hist(sp,dt/2:dt:max(movement_all(units(i)).lfpT));
                sp_hist(i,:) = sp_hist(i,:)/std(sp_hist(i,:));
            end
            clear ticklabel
            for j = 1:length(units)
                ticklabel{j} = num2str(units(j));
            end
            
            figure
            imagesc(sp_hist(:,301:450));
            title(sprintf('%d',sess))
            set(gca,'Ytick',1:size(sp_hist,1));
            set(gca,'YTickLabel',ticklabel);
            
            figure
            imagesc(sp_hist);
            title(sprintf('%d',sess))
            set(gca,'Ytick',1:size(sp_hist,1));
            set(gca,'YTickLabel',ticklabel);
            
            
            set(gcf, 'PaperPositionMode', 'auto');
            print('-dpsc',psfilename,'-append');
            
            
            
            figure
            imagesc(imresize(cov(sp_hist'),1),[-1 1]); colormap jet
            set(gca,'Ytick',1:size(sp_hist,1));
            set(gca,'YTickLabel',ticklabel);
            
            set(gcf, 'PaperPositionMode', 'auto');
            print('-dpsc',psfilename,'-append');
            
            
            [coeff score latent] = pca(sp_hist');
            figure
            plot(latent/sum(latent));
            
            figure
            plot(score(:,1),score(:,2));
            
            npca = min(size(score,2),5);
            figure
            for j=1:npca
                subplot(npca+1,1,j)
                plot(score(:,j)); xlim([0 1200])
            end
            subplot(npca+1,1,npca+1);
            plot(movement_all(units(i)).speed); xlim([0 1200]); ylim([0 20])
            
            set(gcf, 'PaperPositionMode', 'auto');
            print('-dpsc',psfilename,'-append');
            
            figure
            range = max(abs(coeff(:)));
            imagesc(coeff,[-range range]); colorbar
            set(gca,'Ytick',1:size(sp_hist,1));
            set(gca,'YTickLabel',ticklabel);
            
            set(gcf, 'PaperPositionMode', 'auto');
            print('-dpsc',psfilename,'-append');
            
        end
    end
    
    [f p] = uiputfile('*.pdf','pdf name');
    ps2pdf('psfile', psfilename, 'pdffile', fullfile(p,f));
    delete(psfilename);
    
end

keyboard

staUse = find(max(abs(reshape(meanSTAall,size(meanSTAall,1),size(meanSTAall,2)*size(meanSTAall,3))),[],2)>0.06);
%%%staUse=staUse(staUse<143)

clear moviedata m v
clear autocorr

clear h
for c = 1:length(movement_all)
    if ~isempty(movement_all(c).mvlfp_tdata)
        sp = interp1(movement_all(c).mvlfp_tdata, movement_all(c).speed,movement_all(c).frameT);
        duration(c) = max(movement_all(c).frameT);
        lastspike(c) = max(movement_all(c).spikes);
        run_fraction(c) = sum(sp>1)/length(sp);
        run_total(c) = sum(sp>1)/60;
        isi = diff(movement_all(c).spikes);
        
        h(c,:) = hist(isi(isi<0.02),0.0005:0.001:0.02)/length(isi);
    end
end

figure
plot(run_total/60);
xlabel('unit')
ylabel('total running time (min)')

refract = h(:,2)./mean(h(:,6:10),2);
figure
plot(run_total,refract,'o')
xlabel('total running time');
ylabel('isi contamination');

figure
plot(refract);
xlabel('unit')
ylabel('isi contaimination')

figure
cdfplot(refract); xlabel('isi contamination')

figure
cdfplot(run_total); xlabel('running time')

goodunit = run_total>300 & refract'<0.25;  %%% most stringent
goodunit = run_total>150 & refract'<0.5;  %%% least stringent
sprintf('%d good out of %d total = %0.2f',sum(goodunit),length(refract),sum(goodunit)/length(refract))



if false
    for c = 1:length(movement_all)
        
        c
        if ~isempty(movement_all(c).mvlfp_tdata)
            % figure
            spike_sp = interp1(movement_all(c).mvlfp_tdata, movement_all(c).speed,movement_all(c).spikes);
            
            for r = 1:2
                if r==1
                    spikehist = hist(movement_all(c).spikes(spike_sp<1),0:0.005:max(movement_all(c).spikes));
                    spikehist = spikehist/sum(spike_sp<1);
                    
                else
                    spikehist = hist(movement_all(c).spikes(spike_sp>1),0:0.005:max(movement_all(c).spikes));
                    spikehist = spikehist/sum(spike_sp>1);
                end
                
                autocorr(c,r,:) = xcorr(spikehist,spikehist,50,'coeff');
                
                xcspect(c,r,:) = abs(fft(squeeze(autocorr(c,r,:))));
                df = 1/(0.005*length(spikehist));
                spect = abs(fft(spikehist));
                spect = spect(1:round(end/2));
                spect = conv(spect,ones(round(0.5/df),1),'valid');
                spect = (spect-0*min(spect))/(max(spect)-0*min(spect));
                
                for i = 1:150
                    f = 0.5*i -0.25; usef = round((f-0.25)/df+2):round((f+0.25)/df);
                    meanspect(c,r,i) = mean(spect(usef));
                end
                meanspect(c,r,:)= meanspect(c,r,:)/max(meanspect(c,r,:),[],3);
                %                         subplot(3,2,r)
                %             plot((1:length(spect))*df+2,spect); xlim([0 75]); ylim([0.5 1])
                %             hold on
                %             plot(0.25:0.5:74.75,squeeze(meanspect(c,r,:)),'g')
            end
            
            spikehist = histc(movement_all(c).spikes,movement_all(c).frameT)*60;
            sp = interp1(movement_all(c).mvlfp_tdata, movement_all(c).speed,movement_all(c).frameT);
            run_fraction(c) = sum(sp>1)/length(sp);
            ph= spikehist.*exp(2*pi*sqrt(-1)*movement_all(c).frameNum/600);
            ph = conv(ph,ones(2401,1),'same')./conv(spikehist,ones(2401,1),'same');
            
            %         subplot(3,2,3:4)
            %         plot(abs(ph),'.');
            %         hold on
            %         plot(sp/max(sp),'g')
            %
            %         subplot(3,2,5:6)
            %         plot(mod(angle(ph),2*pi),'.');
            %         hold on
            %         plot(sp/max(sp)*5,'g')
        end
    end
end

clear lfpTrig lfpTrigFFT lfpTrigFFTcontrol coh

if false %%% lfp coherence analysis
    for c = 1:length(drift_all)
        
        params.Fs = 1/(median(diff(movement_all(c).lfpT)));
        params.tapers = [3 5];
        params.err = 0;
        params.pad =0;
        params.fpass=[0 75];
        
        lfp = movement_all(c).lfpV;
        spk = movement_all(c).spikes;
        if ~isempty(lfp)
            spike_sp = interp1(movement_all(c).mvlfp_tdata, movement_all(c).speed,movement_all(c).spikes);
            figure;
            for r=1:2
                if r==1
                    use_spk = spk(spike_sp<1); col = 'r';
                else
                    use_spk = spk(spike_sp>1);col='g';
                end
                
                [C phi S12 S1,S2,f,zerosp] = coherencysegcpt(lfp',use_spk',1,params,1);
                
                plot(f,C,col); hold on; ylim([0 0.1])
            end
        end
    end
    for c=1:length(drift_all)
        if ~isempty(movement_all(c).mvlfp_tdata)
            lastT = max(movement_all(c).lfpT);
            lfpAll = zeros(2,length(movement_all(c).spikes),1001);
            lfpAll(:) = NaN; lfpAllFFT=lfpAll;
            spike_sp = interp1(movement_all(c).mvlfp_tdata, movement_all(c).speed,movement_all(c).spikes);
            c
            for i = 1:length(movement_all(c).spikes)-1;
                
                spk = movement_all(c).spikes(i);
                if spk>1 & spk<lastT-1 & spike_sp(i)<1;
                    inds = find(movement_all(c).lfpT>spk);
                    w = 500;
                    lfpAll(1,i,:) = (movement_all(c).lfpV(inds(1)-w:inds(1)+w));
                    lfpAllFFT(1,i,:) = fft(movement_all(c).lfpV(inds(1)-w:inds(1)+w));
                elseif spk>1 & spk<lastT-1 & spike_sp(i)>1;
                    inds = find(movement_all(c).lfpT>spk);
                    lfpAll(2,i,:) = (movement_all(c).lfpV(inds(1)-w:inds(1)+w));
                    lfpAllFFT(2,i,:) = fft(movement_all(c).lfpV(inds(1)-w:inds(1)+w));
                end
            end
            
            
            lfpTrig(c,:,:) =squeeze(nanmean(lfpAll,2));
            figure
            subplot(2,2,1)
            plot(squeeze(lfpTrig(c,:,:))')
            
            lfpTrigFFT(c,:,:) = abs(squeeze(nanmean(lfpAllFFT,2)));
            lfpTrigFFTcontrol(c,:,:) = squeeze(nanmean(abs(lfpAllFFT),2));
            
            dt = median(diff(movement_all(c).lfpT));
            df=1/((2*w+1)*dt);
            subplot(2,2,3)
            plot((1:100)*df,squeeze(lfpTrigFFT(c,:,1:100))');
            subplot(2,2,4)
            plot((1:100)*df,squeeze(lfpTrigFFTcontrol(c,:,1:100))')
            
            
            coh(c,:,:) = lfpTrigFFT(c,:,:)./lfpTrigFFTcontrol(c,:,:);
            
            subplot(2,2,2)
            plot((1:100)*df,squeeze(coh(c,1,1:100)),'r')
            hold on; plot((1:100)*df,squeeze(coh(c,2,1:100)),'g')
        end
        
    end
end


for i = 1:length(drift_mv_all);
    for rep = 1:2
        orient_tune_all(:,:,:,rep,i) = drift_mv_all(i,rep).orient_tune;
        sf_tune_all(:,:,:,rep,i) = drift_mv_all(i,rep).sf_tune;
        for f = 1:2
            sf_tune_all(:,f,1,rep,i) =  sf_tune_all(:,f,1,rep,i) - abs(drift_mv_all(i,rep).spont(f));
            dr_spont_all(f,rep,i) = abs(drift_mv_all(i,rep).spont(f));
        end
    end
    
end

%%% orient_tune_all(tf,f,ori,rep,cell)
%%% sf_tune_all(tf,f,sf,rep,cell)

orient_amp = squeeze(mean(mean(orient_tune_all,3),1));

figure
plot(squeeze(orient_amp(1,1,:)),squeeze(orient_amp(1,2,:)),'o')
hold on
plot(squeeze(orient_amp(1,1,goodunit)),squeeze(orient_amp(1,2,goodunit)),'go')
axis equal;hold on;plot([0 20],[0 20],'g')
xlabel('stop F0'); ylabel('move F0');

F0 = squeeze(orient_amp(1,:,:));

driftSBC = find(F0(1,:)<-5 | F0(2,:)<-5);

figure
plot(squeeze(orient_amp(2,1,:)),squeeze(orient_amp(2,2,:)),'o'); hold on
plot(squeeze(orient_amp(2,1,goodunit)),squeeze(orient_amp(2,2,goodunit)),'go')
axis equal;hold on;plot([0 20],[0 20],'g')
xlabel('stop F1'); ylabel('move F1');

mn_orient_tune = squeeze(mean(orient_tune_all,1));
for i = 1:length(mn_orient_tune);
    for rep=1:2
        for f=1:2
            curve =squeeze( mn_orient_tune(f,:,rep,i));
            if max(curve)>3
                cvOSIall(f,rep,i)=abs(calcCVosi(curve));
            else
                cvOSIall(f,rep,i)=NaN;
            end
        end
    end
end

figure
plot(squeeze(cvOSIall(1,1,:)),squeeze(cvOSIall(1,2,:)),'o'); axis equal; hold on
plot(squeeze(cvOSIall(1,1,goodunit)),squeeze(cvOSIall(1,2,goodunit)),'go'); axis equal
hold on; plot([0 0.8],[0 0.8])


figure
hist(squeeze(cvOSIall(1,1,goodunit)),0.05:0.1:1)
xlabel('stop osi')

figure
hist(squeeze(cvOSIall(1,2,goodunit)),0.05:0.1:1)
xlabel('move osi')

cvOSIavg = squeeze(nanmean(cvOSIall(:,:,goodunit),2));
figure
h = hist(squeeze(cvOSIavg(1,:)),0.05:0.1:1);
bar(0.05:0.1:1,h/sum(h))
xlabel('OSI');
ylabel('%')



for i = 1:length(mn_orient_tune);
    for f = 1:2
        
        [alignedOS(f,:,1,i) alignedOS(f,:,2,i)] =coAlignCurve(squeeze(mn_orient_tune(f,:,1,i)),squeeze(mn_orient_tune(f,:,2,i)));
        
    end
end

sbc = find(xferGain(:,2)<-1 | xferGain(:,1)<-1);
OS = find(cvOSIall(1,1,:)>0.25 | cvOSIall(1,2,:)>0.25);
cc= setdiff(staUse,OS);
cc=setdiff(cc,sbc);
sbc=setdiff(sbc,OS)
driftsbc = find(orient_amp(1,1,:)<-5  & orient_amp(1,2,:)<0 | orient_amp(1,2,:)<-5  &  orient_amp(1,1,:)<0 );
driftsbc=setdiff(driftsbc,cc);
OS = setdiff(OS,driftsbc);

driftsbc = intersect(driftsbc,find(goodunit));
OS = intersect(OS,find(goodunit));
cc = intersect(cc,find(goodunit));



figure
plot(squeeze(dr_spont_all(1,1,:)),squeeze(dr_spont_all(1,2,:)),'o');
hold on
plot(squeeze(dr_spont_all(1,1,goodunit)),squeeze(dr_spont_all(1,2,goodunit)),'go');
xlabel('drift spont stop'); ylabel('spont move')
axis equal
plot([0 20],[0 20],'r')

figure
hist(squeeze(dr_spont_all(1,1,goodunit)),1:1:30);
figure
hist(squeeze(dr_spont_all(1,2,goodunit)),1:1:30);

for i = 1:length(drift_mv_all)
    figure
    for tf = 1:2
        for f = 1:2
            subplot(2,3,tf+(f-1)*3);
            plot(squeeze(orient_tune_all(tf,f,:,1,i) + dr_spont_all(f,1,i)),'r');
            hold on
            plot(1:8, dr_spont_all(f,1,i).*ones(1,8),'r:');
            plot(squeeze(orient_tune_all(tf,f,:,2,i) + dr_spont_all(f,2,i)),'g');
            plot(1:8, dr_spont_all(f,2,i).*ones(1,8),'g:');
            yl  = get(gca,'ylim');
            ylim([0 max(20,yl(2))]); xlim([1 8])
        end
    end
    subplot(2,3,1);
    title(sprintf('spont %0.1f %0.1f',dr_spont_all(1,1,i),dr_spont_all(1,2,i)));
    subplot(2,3,2);
    title(sprintf('F0 %0.1f %0.1f',orient_amp(1,1,i),orient_amp(1,2,i)));
    subplot(2,3,4);
    title(sprintf('F1 %0.1f %0.1f',orient_amp(2,1,i),orient_amp(2,2,i)));
    
    subplot(2,3,5);
    if goodunit(i)
        title('good')
    else
        title('bad')
    end
    
    subplot(2,3,3);
    
    plot(squeeze(xferAll(i,1,:)),'r'); xlim([1 11])
    hold on
    plot(squeeze(xferAll(i,2,:)),'g');
    title(sprintf('run=%0.1f isi=%0.2f',run_total(i)/60,refract(i)));
    subplot(2,3,6);
    imagesc(squeeze(meanSTAall(i,:,:)),[-0.1 0.1]); axis square; axis off
    if ismember(i,cc);
        title('cc')
    elseif ismember(i,sbc);
        title('sbc')
    elseif ismember(i,OS)
        title('OS')
    end
    
end

figure
plot(squeeze(orient_amp(1,1,driftsbc)),squeeze(orient_amp(1,2,driftsbc)),'o')
axis equal;hold on
plot([-20 0],[-20 0],'r');

mean(orient_amp(1,2,driftsbc),3)

figure
for r = 1:2
    subplot(1,2,r);
    imagesc(squeeze(nanmean(spikebinsAll(OS,r,:,:),1)));
    axis equal
end



figure
errorbar(1:8,squeeze(mean(alignedOS(1,:,1,OS),4)), squeeze(std(alignedOS(1,:,1,OS),[],4))/sqrt(length(OS)),'r','Linewidth',2  )
hold on
errorbar(1:8,squeeze(mean(alignedOS(1,:,2,OS),4)), squeeze(std(alignedOS(1,:,2,OS),[],4))/sqrt(length(OS)),'b','Linewidth',2  )
ylim([-8 8]); xlim([1 8])
set(gca,'Fontsize',16)
plot([1 8],[0 0],'b:','Linewidth',2)
set(gca,'Xtick',[1 3 5 7]);
set(gca,'Xticklabel',{'-90','0','90','180'})
xlabel('direction');
ylabel('F0 sp/sec')


figure
errorbar(1:8,squeeze(mean(alignedOS(1,:,1,driftsbc),4)), squeeze(std(alignedOS(1,:,1,driftsbc),[],4))/sqrt(length(driftsbc)),'r','Linewidth',2  )
hold on
errorbar(1:8,squeeze(mean(alignedOS(1,:,2,driftsbc),4)), squeeze(std(alignedOS(1,:,2,driftsbc),[],4))/sqrt(length(driftsbc)),'b','Linewidth',2  )
ylim([-8 8]); xlim([1 8])
set(gca,'Fontsize',16)
plot([1 8],[0 0],'b:','Linewidth',2)
set(gca,'Xtick',[1 3 5 7]);
set(gca,'Xticklabel',{'-90','0','90','180'})
xlabel('direction');
ylabel('F0 sp/sec')


figure
errorbar(1:8,squeeze(mean(alignedOS(1,:,1,cc),4)), squeeze(std(alignedOS(1,:,1,cc),[],4))/sqrt(length(cc)),'r','Linewidth',2  )
hold on
errorbar(1:8,squeeze(mean(alignedOS(1,:,2,cc),4)), squeeze(std(alignedOS(1,:,2,cc),[],4))/sqrt(length(cc)),'b','Linewidth',2  )
ylim([-8 8]); xlim([1 8])
set(gca,'Fontsize',16)
plot([1 8],[0 0],'b:','Linewidth',2)
set(gca,'Xtick',[1 3 5 7]);
set(gca,'Xticklabel',{'-90','0','90','180'})
xlabel('direction');
ylabel('F0 sp/sec')


data = [squeeze(mean(alignedOS(1,3,1,cc),4)) squeeze(mean(alignedOS(1,3,2,cc),4)); ...
    squeeze(mean(alignedOS(1,3,1,driftsbc),4)) squeeze(mean(alignedOS(1,3,2,driftsbc),4)) ; ...
    squeeze(mean(alignedOS(1,3,1,OS),4)) squeeze(mean(alignedOS(1,3,2,OS),4))]

err = [squeeze(std(alignedOS(1,3,1,cc),[],4)) squeeze(std(alignedOS(1,3,2,cc),[],4)); ...
    squeeze(std(alignedOS(1,3,1,driftsbc),[],4)) squeeze(std(alignedOS(1,3,2,driftsbc),[],4)) ; ...
    squeeze(std(alignedOS(1,3,1,OS),[],4)) squeeze(std(alignedOS(1,3,2,OS),[],4))]

err(1,:)=err(1,:)/sqrt(length(cc));
err(2,:)=err(2,:)/sqrt(length(driftsbc));
err(3,:)=err(3,:)/sqrt(length(OS));
err

figure
plot(squeeze(alignedOS(1,3,1,cc)),squeeze(alignedOS(1,3,2,cc)),'o'); hold on
plot(squeeze(alignedOS(1,3,1,intersect(cc,find(goodunit)))),squeeze(alignedOS(1,3,2,intersect(cc,find(goodunit)))),'go');
axis equal;
plot([ 0 20],[0 20]);

figure
plot(squeeze(alignedOS(1,3,1,OS)),squeeze(alignedOS(1,3,2,OS)),'o'); hold on
plot(squeeze(alignedOS(1,3,1,intersect(OS,find(goodunit)))),squeeze(alignedOS(1,3,2,intersect(OS,find(goodunit)))),'go');
axis equal;
plot([ 0 20],[0 20]);



figure
barweb(data,err)
ylim([-8 8])
set(gca,'Xticklabel',{'cc','sbc','OS'})
ylabel('peak F0 sp/sec')
legend('stationary','moving')
set(gca,'Fontsize',16)


OS = find(cvOSIall(1,1,:)>0.25 | cvOSIall(1,2,:)>0.25);

figure
plot(squeeze(orient_amp(1,1,:)),squeeze(orient_amp(1,2,:)),'ko')
axis equal;hold on;
xlabel('stop F0'); ylabel('move F0');
plot(squeeze(orient_amp(1,1,staUse)),squeeze(orient_amp(1,2,staUse)),'bo')
plot(squeeze(orient_amp(1,1,OS)),squeeze(orient_amp(1,2,OS)),'go')
plot(squeeze(orient_amp(1,1,sbc)),squeeze(orient_amp(1,2,sbc)),'ro')
legend({'all','sta_use','os','sbc'})
plot([0 20],[0 20],'g')

cc= setdiff(staUse,OS);
cc=setdiff(cc,sbc);

[median(orient_amp(1,1,cc)) median(orient_amp(1,2,cc))]
[median(orient_amp(1,1,OS)) median(orient_amp(1,2,OS))]
[median(orient_amp(1,1,sbc)) median(orient_amp(1,2,sbc))]

[mean(orient_amp(1,1,cc)) mean(orient_amp(1,2,cc))]
[mean(orient_amp(1,1,OS)) mean(orient_amp(1,2,OS))]
[mean(orient_amp(1,1,sbc)) mean(orient_amp(1,2,sbc))]


keyboard

figure
plot(wvform')



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
% for i = 1:length(sbc)
%     figure
%     subplot(1,2,1); hold on
%     plot(squeeze(xferAll(sbc(i),1,:)),'r')
%     plot(squeeze(xferAll(sbc(i),2,:)),'g')
%     xlim([1 11])
%
%     subplot(1,2,2);hold on
%     plot(squeeze(crfAll(sbc(i),1,1:10)),'r')
%     plot(squeeze(crfAll(sbc(i),2,1:10)),'g')
%     xlim([1 10])
% end




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
h=hist(2*F1F0,0.1:0.2:2); bar(0.1:0.2:2,h/sum(h))
hold on
%h=hist(F1F0(staUse),0:0.1:1); bar(0:0.1:1,h,'r');
xlabel('F1/F0'); xlim([0 2]);

ev_osi = squeeze(osi_all(:,:,1)) - repmat(dr_spont',[1 8]);
figure
plot(ev_osi');


figure
data =[length(cc) length(sbc) length(OS) ]/length(ev_osi)
err = [sqrt(length(cc)) sqrt(length(sbc)) sqrt(length(OS))]/length(ev_osi)

bar(data); hold on
errorbar(1:3,data,err,'k.','Linewidth',2)
set(gca,'Xticklabel',{'cc','sbc','OS'})
ylabel('fraction')
set(gca,'Fontsize',16)

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

% for i = 1:length(F0OS);
%
%     figure
%     imagesc(squeeze(meanSTAall(F0OS(i),:,:)),[-0.1 0.1])
%     axis equal
%
% end



% for i = 1:length(F0OS)
%     figure
%     plot(squeeze(osi_all(F0OS(i),:,:)));
%     hold on
%     plot([1 8], [dr_spont(F0OS(i)) dr_spont(F0OS(i)) ])
%     title(sprintf('osi = %0.2f f1f0 = %0.2f sta = %d',cvosi(F0OS(i)),F1F0(F0OS(i)),ismember(F0OS(i),staUse)))
% end


onoffOS = find(cvosi>0.2 & F1F0'<0.4);
onoffOS = OS;
% for i = 1:length(onoffOS);
%     figure
%     imagesc(squeeze(meanSTAall(onoffOS(i),:,:)),[-0.1 0.1])
%     axis equal
% end

% for i = 1:length(onoffOS)
%     figure
%     hold on
%     plot(squeeze(crfAll(onoffOS(i),1,1:10))','r')
%     plot(squeeze(crfAll(onoffOS(i),2,1:10))','g')
% end
%



for i = 1:5
    nanmean(burst_fraction(staUse(tcGroup==i),:),1)
end

titleList = {'sustain','transient','sbc','DS/OS'};
figure
for i = 1:4
    subplot(2,2,i)
    if i ==1
        d= squeeze(nanmean(xferAll(staUse(tcGroup==1 ),:,:),1));
    elseif i==2
        d= squeeze(nanmean(xferAll(staUse(tcGroup==2 ),:,:),1));
    elseif i==3
        d= squeeze(nanmean(xferAll(sbc,:,:),1));
    elseif i==4
        d= squeeze(nanmean(xferAll(OS ,:,:),1));
    end
    d = d(:,6:10)
    d = (d-min(d(:)))/(max(d(:)-min(d(:))));
    plot(0:0.2:0.8,d(1,:),'r','Linewidth',2)
    hold on
    plot(0:0.2:0.8,d(2,:),'b','Linewidth',2)
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
plot(sbcGain(OS,1),sbcGain(OS,2),'ko','LineWidth',2);
axis square
plot([-10 25],[-10 25])
axis([-10 25 -10 25])
xlabel('response gain stationary')
ylabel('response gain moving')


legend('sustained','transient','sbc','OS')


data = [nanmedian(sbcGain(staUse(tcGroup==1),2)) nanmedian(sbcGain(staUse(tcGroup==1),1)); ...
    nanmedian(sbcGain(staUse(tcGroup==2),2)) nanmedian(sbcGain(staUse(tcGroup==2),1)); ...
    nanmedian(sbcGain(sbc,2)) nanmedian(sbcGain(sbc,1)) ; ...
    nanmedian(sbcGain(OS,2)),nanmedian(sbcGain(OS,1))]

%
% data = [nanmean(sbcGain(staUse(tcGroup==1),2)) nanmean(sbcGain(staUse(tcGroup==1),1)); ...
%     nanmean(sbcGain(staUse(tcGroup==2),2)) nanmean(sbcGain(staUse(tcGroup==2),1)); ...
%     nanmean(sbcGain(sbc,2)) nanmean(sbcGain(sbc,1)) ; ...
%     nanmean(sbcGain(onoffOS,2)),nanmean(sbcGain(onoffOS,1))]

err = [nanstd(sbcGain(staUse(tcGroup==1),2))/7 nanstd(sbcGain(staUse(tcGroup==1),1))/7; ...
    nanstd(sbcGain(staUse(tcGroup==2),2))/7 nanstd(sbcGain(staUse(tcGroup==2),1))/7; ...
    nanstd(sbcGain(sbc,2))/4 nanstd(sbcGain(sbc,1))/4 ; ...
    nanstd(sbcGain(OS,2))/4,nanstd(sbcGain(OS,1))/4]
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
plot(burst_fraction(staUse(tcGroup==1),1),burst_fraction(staUse(tcGroup==1),2),'go','LineWidth',2);
hold on
plot(burst_fraction(staUse(tcGroup==2),1),burst_fraction(staUse(tcGroup==2),2),'bo','LineWidth',2);
plot(burst_fraction(sbc,1),burst_fraction(sbc,2),'ro','LineWidth',2);
plot(burst_fraction(onoffOS,1),burst_fraction(onoffOS,2),'ko','LineWidth',2);
axis square
plot([0 0.3],[0 0.3])
axis([0 0.25 0 0.25])
legend('sustained','transient','sbc','OS')
xlabel('burst fraction stationary');
ylabel('burst fraction moving')

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
barweb(data(:,2:-1:1),err(:,2:-1:1))
set(gca,'xticklabel',{'sustain','transient','sbc','OS'})
ylabel('median burst fraction')
legend('stationary','moving')
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