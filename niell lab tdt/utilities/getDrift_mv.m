function drift_mv = getWn_mv(clustfile,afile, block,redo)
%%% read in spike data for drifting gratings and sort by stimuli
%%% works for gratings stim 6sf x 12 orient + spont + flicker

load(clustfile,'Block_Name','Tank_Name','stimEpocs');
blocknum = find(strcmp(block,Block_Name));
Block_Name = Block_Name{blocknum}

load(afile,'drift_mv');

if ~exist('drift_mv','var') | length(drift_mv)<blocknum  | isempty(drift_mv(blocknum).tuning)| ~isfield(drift_mv,'cv_osi') | redo
    %try
    display('getting spikes')
    spikes=getSpikes(clustfile,afile,block,redo);
    display('getting speed')
    spd = getSpeed(clustfile,afile,block,redo);
    cond = stimEpocs{blocknum}(1,:);
    condT = stimEpocs{blocknum}(2,:);
    
    frameSpd = interp1(spd.t,spd.v,condT+0.5,'nearest');
    dt= 0.05;
    display('doing hist')
    for c =1 :length(spikes.sp);
        c
        for i= 1:length(cond)-1;
            s = find(spikes.sp{c}>=condT(i) & spikes.sp{c}<condT(i)+2.5);
            trialPsth(c,i,:) = hist(spikes.sp{c}(s)-condT(i),dt/2:dt:2.5)/dt;
        end
    end
    
    for i = 1:74;
        ori = ceil(i/6); sf = mod(i-1,6)+1;
        for m = 1:2;
            
            if m==1
                trials = find(cond==i & frameSpd<1);
            else
                trials = find(cond==i & frameSpd>1);
            end
            trials= trials(trials<=size(trialPsth,2));
            mn= mean(mean(trialPsth(:,trials,2:31),3),2); err = std(mean(trialPsth(:,trials,2:31),3),[],2)/sqrt(length(trials));
            if i<=72
                drift_mv(blocknum).tuning(:,ori,sf,m)=mn; drift_mv(blocknum).tuning_err(:,ori,sf,m)=err;
                drift_mv(blocknum).trialOrient(cond==i)=ori;
                drift_mv(blocknum).trialSF(cond==i)=sf;
                
            elseif i==73;
                drift_mv(blocknum).spont(:,m)=mn; drift_mv(blocknum).spont_err(:,m)=err;
                drift_mv(blocknum).trialOrient(cond==i)=NaN;
                drift_mv(blocknum).trialSF(cond==i)=NaN;
            elseif i==74;
                drift_mv(blocknum).flicker(:,m)=mn; drift_mv(blocknum).flicker_err(:,m)=err;
                drift_mv(blocknum).trialOrient(cond==i)=NaN;
                drift_mv(blocknum).trialSF(cond==i)=NaN;
            end
            
        end
    end
    
    
    drift_mv(blocknum).interSpont(:,1) = mean(mean(trialPsth(:,frameSpd(1:end-1)<1,end-10:end),3),2);
    drift_mv(blocknum).interSpont_err(:,1) =std(mean(trialPsth(:,frameSpd(1:end-1)<1,end-10:end),3),[],2)/sqrt(sum(frameSpd(1:end-1)<1));
    drift_mv(blocknum).interSpont(:,2) = mean(mean(trialPsth(:,frameSpd(1:end-1)>1,end-10:end),3),2);
    drift_mv(blocknum).interSpont_err(:,2) =std(mean(trialPsth(:,frameSpd(1:end-1)>1,end-10:end),3),[],2)/sqrt(sum(frameSpd(1:end-1)>1));
    
    drift_mv(blocknum).orient_tune(:,:,:) =squeeze(nanmean(drift_mv(blocknum).tuning,3));
    drift_mv(blocknum).sf_tune(:,2:7,:) = squeeze(nanmean(drift_mv(blocknum).tuning,2));
    drift_mv(blocknum).sf_tune(:,1,:) = drift_mv(blocknum).flicker;
    for i = 1:size(drift_mv(blocknum).orient_tune,1)
        for mv = 1:2
            drift_mv(blocknum).cv_osi(i,mv) = calcCVosi(squeeze(drift_mv(blocknum).orient_tune(i,:,mv)) - drift_mv(blocknum).interSpont(i,mv));
        end
    end
    
    %     for c = 1:length(spikes.sp);
    %         figure
    %         subplot(2,3,1)
    %         imagesc(squeeze(drift_mv(blocknum).tuning(c,:,:,1)));
    %         subplot(2,3,2)
    %         imagesc(squeeze(drift_mv(blocknum).tuning(c,:,:,2)));
    %         subplot(2,3,3);
    %         d1 = drift_mv(blocknum).orient_tune(c,:,1); d2 = drift_mv(blocknum).orient_tune(c,:,2);
    %         plot(d1(:)-drift_mv(blocknum).interSpont(c,1),d2(:)-drift_mv(blocknum).interSpont(c,1),'o'); hold on; plot([0 5], [0 5]);
    %
    %         [r m b] = regression(d1(:)'- drift_mv(blocknum).interSpont(c,1),d2(:)'-drift_mv(blocknum).interSpont(c,1));
    %         plot([0 10],[0 10]*m+b);
    %         subplot(2,3,4);
    %         plot(drift_mv(blocknum).orient_tune(c,:,1),'r'); hold on; plot([1 12],[drift_mv(blocknum).interSpont(c,1) drift_mv(blocknum).interSpont(c,1)],'r:');
    %         plot(drift_mv(blocknum).orient_tune(c,:,2),'b'); hold on; plot([1 12],[drift_mv(blocknum).interSpont(c,2) drift_mv(blocknum).interSpont(c,2)],'b:');
    %         xlim([0.5 12.5])
    %
    %         subplot(2,3,5);
    %         plot(drift_mv(blocknum).sf_tune(c,:,1),'r'); hold on; plot([1 7],[drift_mv(blocknum).interSpont(c,1) drift_mv(blocknum).interSpont(c,1)],'r:');
    %         plot(drift_mv(blocknum).sf_tune(c,:,2),'b'); hold on; plot([1 7],[drift_mv(blocknum).interSpont(c,2) drift_mv(blocknum).interSpont(c,2)],'b:');
    %         xlim([0.5 7.5])
    %     end
    
    figure
    plot(drift_mv(blocknum).interSpont(:,1),drift_mv(blocknum).interSpont(:,2),'o'); hold on
    plot([0 10],[0 10]);
    xlabel('stationary inter spont'); ylabel('mv inter spont');
    
    figure
    plot(drift_mv(blocknum).interSpont(:,2),drift_mv(blocknum).spont(:,2),'o'); hold on
    plot([0 10],[0 10]);
    xlabel('stationary inter spont'); ylabel('stationary drift spont');
    
    
    drift_mv(blocknum).frameSpd = frameSpd(1:end-1);
    drift_mv(blocknum).trialPsth = trialPsth;
    
%     catch
%         drift_mv(blocknum).tuning = [];
%         drift_mv(blocknum).trialOrient = [];
%         drift_mv(blocknum).frameSpd = [];
%         drift_mv(blocknum).trialSF =[];
%         drift_mv(blocknum).R = [];
%         drift_mv(blocknum).cv_osi = [];
% end
    save(afile,'drift_mv','-append')
end

drift_mv = drift_mv(blocknum);
