
GreenCapBatch

for f = 1:length(files)
    ID{f}=files(f).subj;
    lighting(f) = files(f).lighting;
    Type(f)=files(f).type;
    contrast(f) = files(f).contrast;
    group(f)=files(f).group;
    sex(f)=files(f).sex;
    %latency(f) = files(f).CapTime; %files(f).FrameEnd-files(f).FrameS;
    
     %when auto Cricket tracking is working use the following function...
     fn = [pathname files(f).Tfiles];
     load (fn)
%      fn = [pathname files(f).trackpts];
%     if ismac
%         fn(fn=='\')='/';
%     end
    
    if ~isempty(files(f).Tfiles)
   % [ location(:,:,f), targ{f}, body{f}, r{f}, thetaR{f}, thetaSac{f}, vhead{f}, vbody{f}, abody{f}, vtarg{f}, atarg{f}, head{f}, sub{f}] = analyzeCapture(fn,files(f).fps,files(f).scale,files(f).subj);
     %dist(f) = mean(r{f});   
     %when auto Cricket tracking is working use the following function...
     [targ{f}, head{f}, r{f}, thetaR{f}, vhead{f}, vtarg{f} ] = Chase(Left, Right, CentroidB, CentroidCNT,files(f).scale,files(f).fps,files(f).subj);

    end
    
    %title(['animal ',num2str(files(f).subj)])
end

%%% location = histogram of location relative to target
%%% targ = target position (x,y)
%%% body = body position (x,y)
%%% r = distance from head to target
%%% thetaR = angle between head and target
%%% thetaSac = angle of head relative to body
%%% vhead = angular velocity of head
%%% vbody = angular velocity of body
%%% abody = angular acceleration of body
%%% vtarg = angular velocity of target relative to body
%%% atarg = angular velocity of target relative to body
%time to capture


clear rhist speedhist thetahist tdhist tdhist_mv midtheta midtheta_mv spdist basin dRhist approachR
dbins = 1.5:2.5:35; spdbins = 1.25:2.5:50; thetabins = linspace(-150, 150, 13);
dbinsMid=5:2.5:35;
grouplabels = {'Ntsr1_hM4+CNO','GRP_hM4+CNO','PV_hM4+CNO','cre_mCherry+CNO','cre_hM4_saline'};
lightlabels = {'light','dark'};
close all
makeMov=0;

% trials = find(group ==1);
% trials = find(group ==2 & lighting==0);
% trials = find(group ==5 & lighting==1);
% trials = find(group ==5 & lighting==0);

for tr = 1:length(head);   
    %%% loop over all trials
    %%% put extraction of values here
    
    if ~isempty(files(tr).Tfiles)
        
    if makeMov
        close all
    end
    
    group(tr)=files(tr).group;
    
    clear vid;
    if makeMov;
        fn = [pathname files(tr).Moviefile];
        
        %%% load raw avi, downsampled in space and time by factor 'bin'
        bin=5;
        if ~isempty(files(tr).Moviefile)
            vid = loadMouseMovie(fn,bin,files(tr).FrameS,files(tr).FrameEnd);% files(tr).FrameS,files(tr).FrameEnd
            %frameRange = 1:min(files(tr).FrameEnd, max(find(~isnan(thetaR{tr}))))/bin;%sometimes the size of the framerange is off??
            frameRange=1:(files(tr).FrameEnd-files(tr).FrameS)/bin;
            fn = [pathname files(tr).Moviefile];
            %            figure; imshow(vid(:,:,i))
        else
            vid = [];
            frameRange = 1: max(find(~isnan(thetaR{tr})))/bin;
            fn = [pathname files(tr).trackpts];
            
            
        end
        
        figName = sprintf('trial %d %s %s',tr,grouplabels{group(tr)+1},lightlabels{lighting(tr)+1});
        mov = trackMovie(vid,head{tr},targ{tr},thetaR{tr},frameRange,bin,figName);
        
        %%% save movie
        vidObj = VideoWriter([fn(1:end-4) '_tracked.avi']);
        vidObj.FrameRate = 30/bin;
        open(vidObj);
        writeVideo(vidObj,mov);
        close(vidObj);
        
    end
    
    
    if group(tr)~=0 %%% to exclude bad data
       dur = length(r{tr});
       TrackEndT(tr) = dur/files(tr).fps;%time in seconds of track
       %captureT(tr)=files(tr).CapTime; %time in seconds 
      
        % range over whole session
        range = r{tr}(1:dur);
        rangeS=conv(range,ones(1,3),'same')/3;
        rangeSmooth{tr}=rangeS;
        
        rhist(:,tr) = hist(rangeS,dbins)/dur;
        
        %distribution of ranges greater than 5cm
        rhist_mid(:,tr)=hist(rangeS(rangeS>4),dbinsMid)/length(rangeS(rangeS>4));
        %         figure
        %         hist(rangeS(rangeS>4),dbinsMid)
        
        %mouse speeds over whole session
        filt = ones(15,1); filt = filt/length(filt);
        vx = filter(filt,1,diff(head{tr}(:,1))); vy = filter(filt,1,diff(head{tr}(:,2))); speed = sqrt(vx.^2 + vy.^2) * files(tr).fps / files(tr).scale;
        speed = speed(1:dur-1); speed(1:end-6) = speed(7:end); %sp=NaN(size(speed,1)+1,1); sp(1:end-1)=speed(1:end);sp(end)=speed(end);speed=sp;%%% account for phase shift
        
        mSpeedS{tr}=speed;
        speedhist(:,tr) = hist(speed,spdbins)/dur;
        
       
        %target speeds over whole session
        
        filt = ones(5,1); filt = filt/length(filt);
        vx = filter(filt,1,diff(targ{tr}(:,1))); vy = filter(filt,1,diff(targ{tr}(:,2))); speedtarg = sqrt(vx.^2 + vy.^2) * files(tr).fps / files(tr).scale;
        speedtarg = speedtarg(1:dur-1); speedtarg(1:end-2) = speedtarg(3:end);
        
        tSpeedS{tr}=speedtarg;
        speedhistT(:,tr) = hist(speedtarg,spdbins)/dur;

        % thetaS over whole session, mouse theta relative to cricket
        rangeMid=rangeS(1:dur-1);
        
        th = mod(thetaR{tr}(1:dur)+180,360)-180;
        thetaS=conv(th,ones(1,3),'same')/3;
        
        d=diff(thetaS); dd= diff(rangeS);
        thetaSC = thetaS;
        thetaSC(abs(d)>180)=NaN;
        
        thetaSmoothC{tr}=thetaSC;
        
        
        %         figure
        %         hist(thetaS(rangeMid>5& speed>5),thetabins)
        %
        thetahist(:,tr) = hist(thetaSC(1:dur),thetabins)/dur;
        
        thetahist_mid(:,tr)=hist(thetaSC(rangeS>2),thetabins)/length(thetaSC(rangeS>2));
        thetahist_midMove(:,tr) = hist(thetaSC(rangeMid>2 & speed>2),thetabins)/length(thetaSC( rangeMid>2 & speed>2)); 
        tdhist_mv(:,:,tr) = myHist2(r{tr}(rangeMid>2 & speed>2),th(rangeMid>2& speed>2),dbinsMid,thetabins)/sum(rangeMid>2& speed>2);
        %tdhist_mv_mid(:,:,tr) = myHist2(r{tr}(rangeMid),th(rangeMid),dbinsMid,thetabins)/length(rangeMid);
        %[tdhist_mv(:,:,tr),avg_T(tr,:)] = myHist2Avg(r{tr}(rangeMid>5 & speed>5),thetaS(rangeMid>5 & speed>5),dbins,thetabins);
        spdist(:,:,tr) = myHist2(rangeS(1:dur-1),speed,dbins,spdbins)/(dur-1);
        
        figure
        hold on
        plot(speed/25,'g');hold on
        plot(r{tr}(1:dur)/10,'b');
        %plot(abs(th)/180,'r');
        legend('speed','dist','theta'); xlim([0 100])
        title(sprintf('trial %d group %d light %d',tr,group(tr),lighting(tr)));
        
    end
    
    
    contact = range<3.5;
    contact = medfilt2(contact,[10 1]);
    contact(1)=0; contact(end)=0;
    
    contactStart = find(diff(contact)>0);
    contactEnd = find(diff(contact)<0);
    
    %     for j = 1:length(contactEnd)-1;
    %         if (contactStart(j+1)-contactEnd(j))<30
    %             contact(contactEnd(j):contactStart(j+1))=1;
    %         end
    %     end
    
    for j = 1:length(contactStart);
        if j==1
            peak(j)=max(rangeS(1:contactStart(j)));
            if peak(j)<6
                contact(1:contactStart(j))=1;
            end
            
        else
            peak(j)=max(rangeS(contactEnd(j-1):contactStart(j)));
            if peak(j)<6
                contact(contactEnd(j-1):contactStart(j))=1;
            end
        end
    end
    
    contactStart = find(diff(contact)>0)/files(f).fps;
    contactEnd = find(diff(contact)<0)/files(f).fps;
    
    %sometime cicket is already on contact with mouse before trial starts,
    %this gets rid of that "contact"
    
    if length(contactEnd) > length(contactStart);
        contactEnd=contactEnd(2:end);
    end
    
    
    if length(contactStart)>0
        firstContact(tr) = contactStart(1);
    else
        firstContact(tr)=NaN;
    end
    
    nContacts(tr) = length(contactStart);
   % captureProb(tr)=1/length(contactStart);
   captureProb(tr)=1/length(contactStart)  
   
   if captureProb(tr)==0;%cap time took longer than video and tracking
       numIntercepts=(length(contactStart)/TrackEndT(tr))
       if numIntercepts==0
           numIntercepts=1/TrackEndT(tr);
       end
       captureProb(tr)=1/numIntercepts;
       
   end
   
   if length(contactStart)>1
       interContact(tr) = mean(contactStart(2:end) - contactEnd(1:end-1));
       
   else
       interContact(tr)=NaN;    
   end
    
    %%
    if length(contactStart)>0
        contactDuration(tr) = mean(contactEnd-contactStart);
        
        %%%% states
        caught =1; chase =2; sit = 3; roam=4;
        range=range(1:dur-1);
        notCaught = sum(range>2 | speed>4);
        state(caught,tr)= sum(range<2 & speed<4)/length(range);
        state(chase,tr) = sum(range<10 & ~(range<2 & speed<4))/notCaught;
        state(sit,tr) = sum(range>10 & speed<4)/notCaught;
        state(roam,tr)= sum(range>10 & speed>4)/notCaught;
    end
    
    close(gcf)
    
    deltaR = diff(rangeS(1:dur));
    filt = ones(15,1);
    deltaR = filter(filt,1,deltaR); deltaR(1:end-6) = deltaR(7:end);
    dRhist(:,:,tr) = myHist2(speed,deltaR, 2.5:5:50,-19:2:20);
    %figure;imagesc(dRhist(:,:,tr))
    
    %% define approaches based on mouse behavior
    approach = (deltaR<-1 & speed>1); % | contact(1:end-1); %approach = medfilt2(approach,[30 1]);
    approach(1)=0; approach(end)=0;
    starts = find(diff(approach)>0); ends = find(diff(approach)<0);
    
    %     for i = 2:length(ends);
    %         if (starts(i) - ends(i-1))< 35 && (rangeS(starts(i))-rangeS(ends(i-1)))<4 %%%eliminate false starts
    %             approach(starts(i):ends(i-1))=0;
    %         end
    %     end
    % starts = find(diff(approach)>0); ends = find(diff(approach)<0);
    
    
    for i = 2:length(ends);
        if (starts(i) - ends(i-1))< 35 && (rangeS(starts(i))-rangeS(ends(i-1)))<5  %%% stitch together ones that are close in time and continual approach
            approach(ends(i-1):starts(i))=1;
        else if (starts(i) - ends(i-1))<35 && (rangeS(starts(i))-rangeS(ends(i-1)))>5  %%% can't immediately start a new approach
                approach(starts(i):ends(i))=0;
            end
        end
    end
    
    
    % index starts and preApproach frames (~500 ms before start)
    starts = find(diff(approach)>0); ends = find(diff(approach)<0);
   
    for i = 1:length(starts)
        if (rangeS(starts(i)) - rangeS(ends(i)))<10 %| range(starts(i))<15
            approach(starts(i):ends(i)+1)=0;
        end
    end   
   
    
    starts = find(diff(approach)>0);ends = find(diff(approach)<0);
    %% calc number of approaches & frequency
    nApproach(tr) = sum(diff(approach)>0);
    if nApproach(tr)>0
        freqApproach(tr) = (nApproach(tr)/(sum(~contact)/files(tr).fps))*60;%number of events per minute
    else
        freqApproach(tr) = NaN;
    end
    
    %% interapproach Interval
     IAI(tr)=nanmean(diff(starts))/files(tr).fps;
     
     if isnan(IAI(tr));
         IAI(tr)=starts(1)/files(tr).fps;
     end
    
    %% calc approach succes
    clear endDist
    
    if nApproach(tr)>0
        for i = 1:length(ends)
            finalPt = min(ends(i)+10,length(rangeS));
            endDist(i) = min(rangeS(starts(i):finalPt));
        end
        
        outcome = endDist<3;
        approachSuccess(tr) = sum(outcome==1)/length(outcome);
    else
        outcome=NaN;
        approachSuccess(tr)=NaN;
        
    end
    
    %% calc cricket speed prior to approach for successful approaches
    clear preApp
   if  nApproach(tr)>0
       
    if starts(1)<31
        preApp=starts-30;
        preApp(1)=1;
    elseif starts(1)>=31;
    preApp=starts-30;    
    end
    
    preApproach=zeros(size(approach));
    
    for i = 1:length(starts)       
             preApproach(preApp(i):starts(i))=1;
    end
   
    
    approachRange = rangeS; approachRange(~approach)=NaN;
    thetaRange=thetaSC;thetaRange(~approach)=NaN;
    preAppRange=rangeS;preAppRange(~preApproach)=NaN;
    preAppTheta=thetaSC;preAppTheta(~preApproach)=NaN;
    PreAppTargSp=speedtarg;PreAppTargSp(~preApproach)=NaN;
    speedtargApp=speedtarg;speedtargApp(~approach)=NaN;
    
    clear avCsp
    for i=1:length(starts)
    avCsp(i)=nanmean(PreAppTargSp(preApp(i):starts(i)));
    end
    
    avgCspeedSucc{tr}=avCsp(outcome==1);
    avgCspeedFail{tr}=avCsp(outcome==0);
    
    mSpeedRange=speed; mSpeedRange(~approach)=NaN;
    
   end
    %% calc tracks of all approaches
    clear goodAppR goodAppTheta
    
    if nApproach(tr)>0
        for i = 1:length(ends)
            goodAppR{i}=approachRange(starts(i)+1:ends(i));
            goodAppTheta{i} = thetaRange(starts(i)+1:ends(i)); %make thetas symetric
        end
        
    else
        goodAppR{i} = [];
        goodAppTheta{i} = [];
        
    end
    
    AppR{tr}=goodAppR;
    AppT{tr}=goodAppTheta;
   
    %%  plot peri approach data for a given session
    
  %if lighting(tr)==0 & (group(tr)==2 | group(tr)==5);
      figure
        hold on; plot((1:length(deltaR))/30,deltaR,'g'); hold on; plot((1:length(rangeS))/30,rangeS,'k')
    
        hold on
        plot(starts(outcome==0)/30,ones(1,sum(outcome==0))*40,'r*')
        hold on
        plot(starts(outcome==1)/30,ones(1,sum(outcome==1))*40,'g*')
        plot((1:length(approachRange))/30,approachRange,'r');
        xrange = max(30,length(approach)/30);
        xlim([0 xrange]); ylim([-10 60])
        title(sprintf('trial %d group %d light %d n %d freq %0.2f success %0.2f',tr,group(tr),lighting(tr),nApproach(tr),freqApproach(tr),approachSuccess(tr)));
  %end 
        
    end
end

%%
close all
Gr={'ntsr1_CNO','GRP_CNO','PV_CNO','Control_mcherry_CNO','hM4Di_saline'};
colorlabel = 'bgrkcm';
for cond = 1:4
    clear tr
    if cond ==1
        trials = find(group==1);  %%%   or should be group==2 & lighting==1
    elseif cond ==2
        trials = find(group ==2);
    elseif cond ==3
        trials = find(group ==3);
    elseif cond ==4
        trials = find(group ==4 );
    
    end
    
    
%% bar graph of approach frequency and success
    clear s

    R{cond}=vertcat(rangeSmooth{trials});
    T{cond}=vertcat(thetaSmoothC{trials});
    Msp{cond}=vertcat(mSpeedS{trials});
    
%%
CapTime{cond}=TrackEndT(trials);
medianCapTime(cond)=nanmedian(TrackEndT(trials)');
errCT(cond)=semedian(TrackEndT(trials)');

IntAppI{cond}=IAI(trials);
medianIAI(cond)=nanmedian(IAI(trials)');
errIAI(cond)=semedian(IAI(trials)');

data_appFreq{cond}=freqApproach(trials);
dataMu(cond,1) = nanmedian(freqApproach(trials)');
err(cond,1) = semedian(freqApproach(trials)');

data_Int{cond}=approachSuccess(trials);
dataMu(cond,2) = nanmedian(approachSuccess(trials));
err(cond,2)=semedian(approachSuccess(trials));


data_ContDur{cond}=contactDuration(trials);
ContactDataMu(cond,1)=nanmedian(contactDuration(trials)');
Contacterr(cond,1) = semedian(contactDuration(trials)');
%
data_CapProb{cond}=captureProb(trials);
ContactDataMu(cond,2)=nanmedian(captureProb(trials)');
Contacterr(cond,2) = semedian(captureProb(trials)');
    
    
    %% change to account for eccentricity over middle distance
    
    PreyE_mid{cond}=thetahist_midMove(:,trials);
    figure
    plot(thetabins,nanmean(PreyE_mid{cond},2),colorlabel(cond),'LineWidth',3);
    ylim([0 0.25])
    title(sprintf('Prey eccentricity dist >5 speed >4 %s',Gr{cond}));
    
    %range over whole session
    Rsession{cond}=rhist(:,trials);
    figure
    plot(dbins,nanmean(Rsession{cond},2),colorlabel(cond),'LineWidth',3);
    title(sprintf('range whole session %s',Gr{cond}));
    ylim([0 0.3])
    
    %speed over whoel session
    
    Spsession{cond}=speedhist(:,trials);
    figure
    plot(spdbins,nanmean(Spsession{cond},2),colorlabel(cond),'LineWidth',3);
    title(sprintf('range whole session %s',Gr{cond}));
    ylim([0 0.3])
       
    
    clear avg_theta thetaRangeHist
    
    %% to hist distribution to determine natural distance cutoff
    
    dbin=1:2:40;
    ebins=0:15:130;
    
    for tr=1:length(trials);
        clear dist Te
        dist=rangeSmooth{trials(tr)};%extract the trial's ranges
        Te=abs(thetaSmoothC{trials(tr)});%extract the trial's theta error
        [thetaRangeHist(:,:,tr),avg_theta(:,tr)]=myHist2Avg(dist,Te,dbin,ebins);
    end
    clear n
    %extract out N based on whether trace has a theta at a given distance
    for i=1:length(trials)
        e=(avg_theta(:,i));
        n(:,:,i)=~isnan(e);
    end
    
    nAll=sum(n,3);
    stdE=nanstd(avg_theta,[],2)./sqrt(nAll);
    avg_theta_full{cond} = avg_theta;
    
    figure; hold on
    plot(dbin,nanmean(avg_theta,2),colorlabel(cond),'LineWidth',2);hold on
    plot(dbin,(nanmean(avg_theta,2)-stdE),'k','LineWidth',1.5);
    plot(dbin,(nanmean(avg_theta,2)+ stdE),'k','LineWidth',1.5);
    set(gca,'xdir','reverse','XLim',[0 40],'YLim',[0 180])
    title(sprintf('avg theta whole session %s',Gr{cond}));
    
    %% approach paths 
    clear goodApp AppEcc AppErr
    goodApp=0;
    
    for tr = 1:length(trials);
        
        for i = 1:length(AppR{trials(tr)});
            
            clear x y maxX preTouch dist mindist
            x=AppR{trials(tr)}{:,i};
            y=AppT{trials(tr)}{:,i};
            preTouch=find(x<5);%find places where mouse is within 5cm
            
            if ~isempty (preTouch);
                dist=preTouch(1);
            else
                mindist=min(x);
                dist=find(x==mindist);
            end
            
            lookAway =90-abs(y)<0;%find first point where prey azimuth is <90  %%% replaced function nearest, unknown?
            startList=find(diff(lookAway)<0);
            
            if ~isempty(startList) && startList(1)<dist;
            start=startList(1);
            else
            start=1;
            end
            
            
            if ~isempty(x)
                goodApp=goodApp+1;
                
                AppEcc{goodApp}=[x(start:dist),abs(y(start:dist))];
                  
                %AppErr(goodApp)=abs(y(dist));
                
                if ~isempty(x(dist));
                    AppErr(goodApp)=abs(y(dist));
                else
                    AppErr(goodApp)=NaN;
                end
                
            end
        end
    end
    
    A{cond}=vertcat(AppEcc{:});
    
    preContT{cond}=AppErr;
    
    preContT_mu(cond,1)=nanmean(AppErr);
    preContErr(cond,1)=nanstd(AppErr)/sqrt(length(AppErr));
    
    binD=4:3:28;
    binE=-60:5:60;
    
    clear thetaRangeHist avg_theta
    
    for i=1:length(AppEcc);
        if ~isempty(AppEcc{i})
            [thetaRangeHist(:,:,i),avg_theta(i,:)]=myHist2Avg(AppEcc{i}(:,1),AppEcc{i}(:,2),binD,binE);
        end
    end
    
    
    figure; hold on
    col=[0,0,1,.5]
    for i=1:length(AppEcc)
        if ~isempty(AppEcc{i})
            plot(AppEcc{i}(:,1),AppEcc{i}(:,2),'color',col,'LineWidth',0.5)
        end
    end
    set(gca,'xdir','reverse','XLim',[0 25],'YLim',[0 120]);
    n=sum(~isnan(avg_theta));
    stdE=nanstd(avg_theta,[],1)./sqrt(n);
    plot(binD,nanmean(avg_theta),colorlabel(cond),'LineWidth',2);hold on
    plot(binD,(nanmean(avg_theta)-stdE),'k','LineWidth',1.5);
    plot(binD,(nanmean(avg_theta)+ stdE),'k','LineWidth',1.5);
    title(sprintf('average theta approaches %s',Gr{cond}));
    set(gca,'xdir','reverse','XLim',[0 25],'YLim',[0 100])
    avg_theta_approach{cond} = avg_theta;
    stdE_app{cond}=stdE;
    
end
    
figure; hold on
       shadedErrorBar(binD,nanmean(avg_theta_approach{1}),stdE_app{1},'b'); hold on
       shadedErrorBar(binD,nanmean(avg_theta_approach{2}),stdE_app{2},'k'); hold on
       shadedErrorBar(binD,nanmean(avg_theta_approach{3}),stdE_app{3},'g'); hold on
       shadedErrorBar(binD,nanmean(avg_theta_approach{4}),stdE_app{4},'r'); hold on
       set(gca,'xdir','reverse','XLim',[0 25],'YLim',[0 60])
       
       


%% outcome bar graphs figure generation

%% capture times

figure
bar(medianCapTime(1:4)); hold on; errorbar(1:4,medianCapTime(1:4),errCT(1:4),'o')
ylabel('CaptureTime')
set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'NTSR1_C','GRP_C','PV_C','control'});


%% Interapproach Interval
figure
bar(medianIAI(1:4)); hold on; errorbar(1:4,medianIAI(1:4),errIAI(1:4),'o')
ylabel('Interapproach interval')
set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'NTSR1_C','GRP_C','PV_C','control'});

hold on
for ii=1:4
    tmp=IntAppI{ii}; %temporarily store data in variable "tmp"
    x = repmat(ii,1,length(tmp)); %the x axis location
    x = (x+(rand(size(x))-0.6)*0.1-.1); %add a little random "jitter" to aid visibility

    plot(x,tmp,'ok')
end

hold off

%%  approach Frequency plots
figure
bar(dataMu(1:4,1)); hold on; errorbar(1:4,dataMu(1:4,1),err(1:4,1),'o')
ylabel('approach Frequency')
set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'NTSR1_C','GRP_C','PV_C','control'});

hold on
for ii=1:4
    tmp=data_appFreq{ii}; %temporarily store data in variable "tmp"
    x = repmat(ii,1,length(tmp)); %the x axis location
    x = (x+(rand(size(x))-0.6)*0.1-.1); %add a little random "jitter" to aid visibility

    plot(x,tmp,'ok')
end
hold off

% AF_light=horzcat(data_appFreq{1}',repmat(1,length(data_appFreq{1}),1));
% AF_dark=horzcat(data_appFreq{2}',repmat(2,length(data_appFreq{2}),1));
% AF_EPL=horzcat(data_appFreq{3}',repmat(3,length(data_appFreq{3}),1));
% AF_EPD=horzcat(data_appFreq{4}',repmat(4,length(data_appFreq{4}),1));
% AF_ESL=horzcat(data_appFreq{5}',repmat(5,length(data_appFreq{5}),1));
% %AF_ESD=horzcat(data_appFreq{6}',repmat(6,length(data_appFreq{6}),1));
% 
% AppFreq=vertcat(AF_light,AF_dark,AF_EPL,AF_EPD,AF_ESL,AF_ESD);
% 
% [p,tbl,stats] = kruskalwallis(AppFreq(:,1),AppFreq(:,2))
% c=multcompare(stats);



%% prey interception rates
figure
bar(dataMu(1:4,2)); hold on; errorbar(1:4,dataMu(1:4,2),err(1:4,2),'o')
ylabel('prey intercept Success')
set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'NTSR1','GRP','PV','ctl'});
title('figure 3D bottom')

 hold on
for ii=1:4
    tmp=(data_Int{ii}+(rand(size(data_Int{ii}))+0.2)*0.05);
    x = repmat(ii,1,length(tmp)); %the x axis location
    x = (x+(rand(size(x))-0.5)*0.1-.1); %add a little random "jitter" to aid visibility

    plot(x,tmp,'ok')
end
hold off

% PI_light=horzcat(data_Int{1}',repmat(1,length(data_Int{1}),1));
% PI_dark=horzcat(data_Int{2}',repmat(2,length(data_Int{2}),1));
% PI_EPL=horzcat(data_Int{3}',repmat(3,length(data_Int{3}),1));
% PI_EPD=horzcat(data_Int{4}',repmat(4,length(data_Int{4}),1));
% preyInt=vertcat(PI_light,PI_dark,PI_EPL,PI_EPD);
% 
% [p,tbl,stats] = kruskalwallis(preyInt(:,1),preyInt(:,2))
% c=multcompare(stats);

%% capture Probability   
figure
bar(ContactDataMu(1:4,2)); hold on; errorbar(1:4,ContactDataMu(1:4,2),Contacterr(1:4,2),'o')
ylabel('Cap Probability')
set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'NTSR1','GRP','PV','ctl'});
hold on
for ii=1:4
    
    tmp=(data_CapProb{ii}+(rand(size(data_CapProb{ii}))+0.1)*0.05)
    x = repmat(ii,1,length(tmp)); %the x axis location
    x = (x+(rand(size(x))-0.5)*0.1-.1); %add a little random "jitter" to aid visibility

    plot(x,tmp,'ok')
end
hold off
% CP_light=horzcat(data_CapProb{1}',repmat(1,length(data_CapProb{1}),1));
% CP_dark=horzcat(data_CapProb{2}',repmat(2,length(data_CapProb{2}),1));
% CP_EPL=horzcat(data_CapProb{3}',repmat(3,length(data_CapProb{3}),1));
% CP_EPD=horzcat(data_CapProb{4}',repmat(4,length(data_CapProb{4}),1));
% CapP=vertcat(CP_light,CP_dark,CP_EPL,CP_EPD);
% 
% [p,tbl,stats] = kruskalwallis(CapP(:,1),CapP(:,2))
% c=multcompare(stats);

%% contact duration
  
% figure
% bar(ContactDataMu(1:4,1)); hold on; errorbar(1:4,ContactDataMu(1:4,1),Contacterr(1:4,1),'o')
% ylabel('Contact Duration')
% set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'light','dark','deaf','dark + deaf'});
% hold on
% for ii=1:4
%     
%     tmp=(data_ContDur{ii}+(rand(size(data_ContDur{ii}))+0.1)*0.05)
%     x = repmat(ii,1,length(tmp)); %the x axis location
%     x = (x+(rand(size(x))-0.5)*0.1-.1); %add a little random "jitter" to aid visibility
% 
%     plot(x,tmp,'ok')
% end
% hold off
% CP_light=horzcat(data_ContDur{1}',repmat(1,length(data_ContDur{1}),1));
% CP_dark=horzcat(data_ContDur{2}',repmat(2,length(data_ContDur{2}),1));
% CP_EPL=horzcat(data_ContDur{3}',repmat(3,length(data_ContDur{3}),1));
% CP_EPD=horzcat(data_ContDur{4}',repmat(4,length(data_ContDur{4}),1));
% ContD=vertcat(CP_light,CP_dark,CP_EPL,CP_EPD);
% 
% [p,tbl,stats] = kruskalwallis(ContD(:,1),ContD(:,2))
% c=multcompare(stats);

%% Pre Contact approach error
figure
bar(preContT_mu(1:4,1)); hold on; errorbar(1:4,preContT_mu(1:4,1),preContErr(1:4,1),'o')
ylabel('Pre Contact prey azimuth')
set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'NTSR1','GRP','PV','ctl'});

hold on
for ii=1:4
    tmp=(preContT{ii}+(rand(size(preContT{ii}))+ 0.1)*0.05); %temporarily store data in variable "tmp"
    x = repmat(ii,1,length(tmp)); %the x axis location
    x = (x+(rand(size(x))-0.5)*0.1-.1); %add a little random "jitter" to aid visibility

    plot(x,tmp,'ok')
end
hold off

% PA_light=horzcat(preContT{1}',repmat(1,length(preContT{1}),1));
% PA_dark=horzcat(preContT{2}',repmat(2,length(preContT{2}),1));
% PA_EPL=horzcat(preContT{3}',repmat(3,length(preContT{3}),1));
% PA_EPD=horzcat(preContT{4}',repmat(4,length(preContT{4}),1));
% PA_ESL=horzcat(preContT{5}',repmat(5,length(preContT{5}),1));
% PA_ESD=horzcat(preContT{6}',repmat(6,length(preContT{6}),1));
% 
% PAmain=vertcat(PA_light,PA_dark,PA_EPL,PA_EPD, PA_ESL,PA_ESD);
% 
% [P,T,STATS,TERMS]=anovan(PAmain(:,1),PAmain(:,2)) ;%test whether there are differences by factor/groups
%  results=multcompare(STATS)



%% hist theta, range and speed

figure
plot(thetabins,nanmedian(PreyE_mid{1},2),'b','LineWidth',3); hold on
plot(thetabins,nanmedian(PreyE_mid{4},2),'k','LineWidth',3);
figure
plot(thetabins,nanmedian(PreyE_mid{2},2),'g','LineWidth',3);hold on
plot(thetabins,nanmedian(PreyE_mid{4},2),'k','LineWidth',3);

figure
plot(thetabins,nanmedian(PreyE_mid{3},2),'r','LineWidth',3);hold on
plot(thetabins,nanmedian(PreyE_mid{4},2),'k','LineWidth',3);

title 'Prey eccentricity all conditions >5cm mouse moving'

figure
plot(thetabins,nanmedian(PreyE_mid{1},2),'b','LineWidth',3); hold on
plot(thetabins,nanmedian(PreyE_mid{2},2),'g','LineWidth',3);
plot(thetabins,nanmedian(PreyE_mid{3},2),'r','LineWidth',3);
plot(thetabins,nanmedian(PreyE_mid{4},2),'k','LineWidth',3); 


% PLOT ALL CONDITIONS range over trial
figure
plot(nanmedian(Rsession{1},2),'b','LineWidth',3);hold on
plot(nanmedian(Rsession{2},2),'g','LineWidth',3);hold on
plot(nanmedian(Rsession{3},2),'r','LineWidth',3);hold on
plot(nanmedian(Rsession{4},2),'k','LineWidth',3);
ylabel('RangeOverall')

% speed plots
figure
plot(nanmedian(Spsession{1},2),'b','LineWidth',3);hold on
plot(nanmedian(Spsession{2},2),'g','LineWidth',3);hold on
plot(nanmedian(Spsession{3},2),'r','LineWidth',3);hold on
plot(nanmedian(Spsession{4},2),'k','LineWidth',3);

ylabel('speed')



