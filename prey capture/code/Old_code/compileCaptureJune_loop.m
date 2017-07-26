%close all
%clear all

captureBatchNew

if ismac
    pathname = '/Users/crisniell/Dropbox/Prey capture/Prey_tracks/'
end


for f = 1:length(files)
    ID{f}=files(f).subj;
    lighting(f) = files(f).lighting;
    contrast(f) = files(f).contrast;
    group(f)=files(f).group;
    sex(f)=files(f).sex;
    latency(f) = files(f).CapTime; %files(f).FrameEnd-files(f).FrameS;
    fn = [pathname files(f).trackpts];
    if ismac
        fn(fn=='\')='/';
    end
    if ~isempty(files(f).trackpts)
    [ location(:,:,f) targ{f} body{f} r{f} thetaR{f} thetaSac{f} vhead{f} vbody{f} abody{f} vtarg{f} atarg{f} head{f} sub{f}] = analyzeCapture(fn,files(f).fps,files(f).scale,files(f).subj);
     dist(f) = mean(r{f});   
    end
    
    title(['animal ',num2str(files(f).subj)])
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


d = [mean(latency(lighting==1 & group==1)),...
    mean(latency(lighting==0 & group==2)),...
    mean(latency(lighting==1 & group==3)),...
    mean(latency(lighting==0 & group==3)),...
    mean(latency(lighting==1 & group==4)),...
    mean(latency(lighting==0 & group==4)),...
    mean(latency(lighting==1 & group==5)),...
    mean(latency(lighting==0 & group==5))]/60;

dMain = [mean(latency(lighting==1 & group==1)),...
    mean(latency(lighting==0 & group==2)),...
    mean(latency(lighting==1 & group==5)),...
    mean(latency(lighting==0 & group==5))]/60;

err = [std(latency(lighting==1 & group==1))/sqrt(sum(lighting==1 & group==1)),...%light premanipulation group
    std(latency(lighting==0 & group==2))/sqrt(sum(lighting==0 & group==2)),...% dark exposure for the first time group, alternated with light trials
    std(latency(lighting==1 & group==3))/sqrt(sum(lighting==1 & group==3)),...%light eye suture group
    std(latency(lighting==0 & group==3))/sqrt(sum(lighting==0 & group==3)),...% dark eye suture group
    std(latency(lighting==1 & group==4))/sqrt(sum(lighting==1 & group==4)),...%light, Blk background group
    std(latency(lighting==0 & group==4))/sqrt(sum(lighting==0 & group==4)),...%dark Blk Background group
    std(latency(lighting==1 & group==5))/sqrt(sum(lighting==1 & group==5)),...%light, EP'd group
    std(latency(lighting==0 & group==5))/sqrt(sum(lighting==0 & group==5))]/60;%dark EP group

errMain = [std(latency(lighting==1 & group==1))/sqrt(sum(lighting==1 & group==1)),...%light premanipulation group
    std(latency(lighting==0 & group==2))/sqrt(sum(lighting==0 & group==2)),...% dark exposure for the first time group, alternated with light trials
    std(latency(lighting==1 & group==5))/sqrt(sum(lighting==1 & group==5)),...%light, EP'd group
    std(latency(lighting==0 & group==5))/sqrt(sum(lighting==0 & group==5))]/60;%dark EP group

dEye = [mean(latency(lighting==1 & group==1)),...
    mean(latency(lighting==0 & group==2)),...
    mean(latency(lighting==1 & group==3)),...
    mean(latency(lighting==0 & group==3))]/60;

errEye = [std(latency(lighting==1 & group==1))/sqrt(sum(lighting==1 & group==1)),...
    std(latency(lighting==0 & group==2))/sqrt(sum(lighting==0 & group==2)),...% dark exposure for the first time group, alternated with light trials
    std(latency(lighting==1 & group==3))/sqrt(sum(lighting==1 & group==3)),...%light eye suture group
    std(latency(lighting==0 & group==3))/sqrt(sum(lighting==0 & group==3))]/60;% dark eye suture group



% dmed = [median(latency(lighting==1 & group==1)),...
%     median(latency(lighting==0 & group==2)),...
%     median(latency(lighting==1 & group==3)),...
%     median(latency(lighting==0 & group==3)),...
%     median(latency(lighting==1 & group==4)),...
%     median(latency(lighting==0 & group==4)),...
%     median(latency(lighting==1 & group==5)),...
%     median(latency(lighting==0 & group==5))]/60;

figure
barweb(d,err)
ylabel('time to caputure (s)')
legend('Light','dark', 'ES light','ES Dark','BlkBgd light','BlkBgd dark','EP light', 'EP dark')

figure
barweb(dMain,errMain,'b')
ylabel('time to capture (s)')


figure
% col=[1,1,1,0.2]
barweb(dEye,errEye,'b')
ylabel('time to capture (s)')
legend('dark','ES light', 'ES dark');

clear rhist speedhist thetahist tdhist tdhist_mv midtheta midtheta_mv spdist basin dRhist approachR
dbins = 1.5:2.5:35; spdbins = 1.25:2.5:50; thetabins = linspace(-150, 150, 13);
dbinsMid=5:2.5:35;
grouplabels = {'null','dark','suture','black bgd','earplug','','','','',''};
lightlabels = {'light','dark'};
close all
makeMov=0;

% trials = find(group ==1);
% trials = find(group ==2 & lighting==0);
% trials = find(group ==5 & lighting==1);
% trials = find(group ==5 & lighting==0);

for tr = 1:length(head); % trials(example)%44 is the start of the newly added Ear plug videos
    
    %%% loop over all trials
    %%% put extraction of values here
    
    if ~isempty(files(tr).trackpts)
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
        dur = max(find(~isnan(r{tr})));
       TrackEndT(tr) = dur/files(tr).fps;%time in seconds of track
       captureT(tr)=files(tr).CapTime; %time in seconds
        
        %% plot head orientation on top of tracks
        %                 if (group(tr)==1 & lighting(tr)==1) | (group(tr)==4 & lighting(tr)==1)
        %                     figure
        %                     hold on
        %
        %
        %                     plot(targ{tr}(:,1),targ{tr}(:,2),'k');
        %                     for t = 1:3:dur
        %
        %                         %plot(targ{tr}(t,1),targ{tr}(t,2),'go');
        %                         %plot(head{tr}(t,1),head{tr}(t,2),'ko');
        %                         rangeangle = atan2( targ{tr}(t,2) - head{tr}(t,2), targ{tr}(t,1) - head{tr}(t,1))*180/pi;
        %                         headangle = rangeangle -thetaR{tr}(t);
        %                         if r{tr}(t)<10
        %                             plot([head{tr}(t,1) head{tr}(t,1) + cosd(headangle)*50 ], [head{tr}(t,2) head{tr}(t,2) + sind(headangle)*50],'r');
        %                             %  plot(head{tr}(t,1),head{tr}(t,2),'r.','Markersize',8)
        %                         else
        %                             plot([head{tr}(t,1) head{tr}(t,1) + cosd(headangle)*50 ], [head{tr}(t,2) head{tr}(t,2) + sind(headangle)*50],'b');
        %                             % plot(head{tr}(t,1),head{tr}(t,2),'b.','Markersize',8)
        %                         end
        %                         plot([head{tr}(t,1) head{tr}(t,1) + cosd(rangeangle)*50 ], [head{tr}(t,2) head{tr}(t,2) + sind(rangeangle)*50],'g');
        %                         %             plot(head{tr}(1:t,1),head{tr}(1:t,2))
        %                         %             plot(targ{tr}(1:t,1),targ{tr}(1:t,2),'r');
        %                         %             axis([0 2000 0 1200])
        %                         %             mov(t) = getframe(gcf);
        %                     end
        %                     axis([0 2000 0 1200])
        %                     title(sprintf('trial %d %s %s',tr,grouplabels{group(tr)+1},lightlabels{lighting(tr)+1}));
        %                 end
        
        
        %create a list of sub IDs to align with all instances of a measure
        s=sub{tr};
        sID={s};
        ID{tr}=repmat(sID,[dur,1]);
        
        sx=sex(tr);
        Sx{tr}=repmat(sx,[dur,1]);
        
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
        
        thetahist_mid(:,tr)=hist(thetaSC(rangeS>4),thetabins)/length(thetaSC(rangeS>4));
        thetahist_midMove(:,tr) = hist(thetaSC(rangeMid>4 & speed>4),thetabins)/length(thetaSC( rangeMid>4 & speed>4));
        
        tdhist_mv(:,:,tr) = myHist2(r{tr}(rangeMid>4 & speed>4),th(rangeMid>4& speed>4),dbinsMid,thetabins)/sum(rangeMid>4& speed>4);
        %tdhist_mv_mid(:,:,tr) = myHist2(r{tr}(rangeMid),th(rangeMid),dbinsMid,thetabins)/length(rangeMid);
        
        %[tdhist_mv(:,:,tr),avg_T(tr,:)] = myHist2Avg(r{tr}(rangeMid>5 & speed>5),thetaS(rangeMid>5 & speed>5),dbins,thetabins);
        
        spdist(:,:,tr) = myHist2(rangeS(1:dur-1),speed,dbins,spdbins)/(dur-1);
        
        figure
        hold on
        plot(speed/25,'g');
        plot(r{tr}(1:dur)/10,'b');
        %plot(abs(th)/180,'r');
        legend('speed','dist','theta'); xlim([0 3000])
        
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
            if peak(j)<8
                contact(1:contactStart(j))=1;
            end
            
        else
            peak(j)=max(rangeS(contactEnd(j-1):contactStart(j)));
            if peak(j)<8
                contact(contactEnd(j-1):contactStart(j))=1;
            end
        end
    end
    
    contactStart = find(diff(contact)>0)/files(f).fps;
    contactEnd = find(diff(contact)<0)/files(f).fps;
    
    %sometime cicket is already on contact with mouse before trial starts,
    %this gets rid of that "contact"
    
    if length(contactEnd) > length(contactStart)
        contactEnd=contactEnd(2:end)
    end
    
    
    if length(contactStart)>0
        firstContact(tr) = contactStart(1);
    else
        firstContact(tr)=NaN;
    end
    
    nContacts(tr) = length(contactStart);
   % captureProb(tr)=1/length(contactStart);
    if files(tr).CapTime>70;%cap time took longer than video and tracking
    numIntercepts=(length(contactStart)/TrackEndT(tr))*(files(tr).CapTime)
    captureProb(tr)=1/numIntercepts;
    else
    
    captureProb(tr)=1/length(contactStart)

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
    
    %% define appraoches based on mouse behavior
    approach = (deltaR<-1 & speed>5); % | contact(1:end-1); %approach = medfilt2(approach,[30 1]);
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
    
    %     for i = 2:length(ends);
    %         if (starts(i) - ends(i-1))< 35 & (rangeS(starts(i))-rangeS(ends(i-1)))<5  %%% stitch together ones that are close in time and continual approach
    %             approach(ends(i-1):starts(i))=1;
    %         else if (starts(i) - ends(i-1))<60 & (rangeS(starts(i))-rangeS(ends(i-1)))>5  %%% can't immediately start a new approach
    %                 approach(starts(i):ends(i))=0;
    %             end
    %         end
    %     end
    
    
    starts = find(diff(approach)>0); ends = find(diff(approach)<0);
    
    for i = 1:length(starts)
        if (rangeS(starts(i)) - rangeS(ends(i)))<10 %| range(starts(i))<15
            approach(starts(i):ends(i)+1)=0;
        end
    end
    
    
    
    approachRange = rangeS; approachRange(~approach)=NaN;
    thetaRange=thetaSC;thetaRange(~approach)=NaN;
    
    
    mSpeedRange=speed; mSpeedRange(~approach)=NaN;
    speedtargRange=speedtarg;speedtargRange(~approach)=NaN;
    
    nApproach(tr) = sum(diff(approach)>0);
    if nApproach(tr)>0
        freqApproach(tr) = nApproach(tr)/(sum(~contact)/files(tr).fps);
    else
        freqApproach(tr) = NaN;
    end
    
    starts = find(diff(approach)>0);
    ends = find(diff(approach)<0);
    
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
    
    
    clear goodAppR goodAppTheta
    
    if nApproach(tr)>0
        for i = 1:length(ends)
            %appR{tr}(i,1:length(starts(i)+1:ends(i))) = NaN; appTheta{tr}(i,1:length(starts(i)+1:ends(i))) = NaN;  %%% fill with Nans for short segements
            
            goodAppR{i}=approachRange(starts(i)+1:ends(i));
            goodAppTheta{i} = thetaRange(starts(i)+1:ends(i)); %make thetas symetric
            
        end
        
    else
        goodAppR{i} = [];
        goodAppTheta{i} = [];
        
    end
    %
    AppR{tr}=goodAppR;
    AppT{tr}=goodAppTheta;
    %%
    % plot peri approach data for a given session
    
  if lighting(tr)==0 & (group(tr)==2 | group(tr)==5)
      figure
        hold on; plot((1:length(deltaR))/60,deltaR,'g'); hold on; plot((1:length(rangeS))/60,rangeS,'k')
    
        hold on
        plot(starts(outcome==0)/60,ones(1,sum(outcome==0))*40,'r*')
        hold on
        plot(starts(outcome==1)/60,ones(1,sum(outcome==1))*40,'g*')
        plot((1:length(approachRange))/60,approachRange,'r');
        xrange = max(60,length(approach)/60);
        xlim([0 xrange]); ylim([-10 50])
        title(sprintf('trial %d group %d light %d n %d freq %0.2f success %0.2f',tr,group(tr),lighting(tr),nApproach(tr),freqApproach(tr),approachSuccess(tr)));
  end 
    %     behavFig = figure;
    %
    % % plot range between mouse and cricket
    % figure(behavFig);
    % subplot(1,3,1);plot((1:length(rangeSmooth{tr}(2:end)))/60,rangeSmooth{tr}(2:end),'k','LineWidth',3); axis([0 length(rangeSmooth{tr}(2:end))/60 0 40]);hold on
    %     plot(starts(outcome==0)/60,ones(1,sum(outcome==0))*40,'r*');hold on
    %     plot(starts(outcome==1)/60,ones(1,sum(outcome==1))*40,'g*'); hold on
    %         plot((1:length(approachRange))/60,approachRange,'r','LineWidth',2);
    %
    %  % plot bearing between mouse and cricket(prey eccentricity)
    % subplot(1,3,2);plot((1:length(rangeSmooth{tr}(2:end)))/60,thetaSmooth{tr}(2:end),'k','LineWidth',3); axis([0 length(rangeSmooth{tr}(2:end))/60 -180 180]); hold on
    % plot(starts(outcome==0)/60,ones(1,sum(outcome==0))*40,'r*');hold on
    %     plot(starts(outcome==1)/60,ones(1,sum(outcome==1))*40,'g*');hold on
    %     plot((1:length(thetaRange))/60,thetaRange,'r','LineWidth',2);
    %
    %  % plot mouse speed
    % subplot(1,3,3);plot((1:length(rangeSmooth{tr}(2:end)))/60, mSpeedS{tr},'k','LineWidth',3); axis([0 length(rangeSmooth{tr}(2:end))/60 0 55]); hold on
    % plot(starts(outcome==0)/60,ones(1,sum(outcome==0))*40,'r*');hold on
    %     plot(starts(outcome==1)/60,ones(1,sum(outcome==1))*40,'g*');hold on
    %     plot((1:length(mSpeedRange))/60,mSpeedRange,'r','LineWidth',2);
    
    %subplot(4,1,4);plot(tSpeedS{i},'k','LineWidth',3); axis([0 length(rangeSmooth{i}) 0 100]);
    
    %     figure
    %             plot(targ{tr}(:,1),targ{tr}(:,2),'k');hold on
    %             plot(head{tr}(:,1),head{tr}(:,2),'b');
    %             axis([0 2000 0 1200]); axis off;
    %             %col=colormap('jet')
    %             %col = 'rgbcmk';
    % %
    %             for i=1:length(starts)
    %             k=rand(1,3)
    %             hold on;plot(head{tr}(starts(i):ends(i),1),head{tr}(starts(i):ends(i),2),'color',k,'Linewidth',2); hold on;
    %             plot(head{tr}(starts(i),1),head{tr}(starts(i),2),'go','Linewidth',2);hold on;
    %             plot(targ{tr}(starts(i):ends(i),1),targ{tr}(starts(i):ends(i),2),'color',k,'Linewidth',2);hold on;
    %             plot(targ{tr}(starts(i),1),targ{tr}(starts(i),2),'r*','Linewidth',2);hold on;
    %
    % %             rangeangle = atan2( targ{trials(tr)}(starts(i):ends(i),2) - head{trials(tr)}(starts(i):ends(i),2), targ{trials(tr)}(starts(i):ends(i),1) - head{trials(tr)}(starts(i):ends(i),1))*180/pi;
    % %             headangle = rangeangle -thetaR{trials(tr)}(starts(i):ends(i));
    % %             plot([head{trials(tr)}(starts(i):ends(i),1) head{trials(tr)}(starts(i):ends(i),1) + cosd(headangle)*50] , [head{trials(tr)}(starts(i):ends(i),2) head{trials(tr)}(starts(i):ends(i),2) + sind(headangle)*50],'b','Linewidth',2);
    %             end
    %col(mod(i,6)+1)
end
end

%%
close all
Gr={'Light','Dark','EPL','EPD','ESL','ESD'};
colorlabel = 'bkgrcm';
for cond = 1:6
    clear tr
    if cond ==1
        trials = find(group==1);  %%%   or should be group==2 & lighting==1
    elseif cond ==2
        trials = find(group ==2 & lighting==0);
    elseif cond ==3
        trials = find(group ==5 & lighting==1);
    elseif cond ==4
        trials = find(group ==5 & lighting==0);
    elseif cond ==5
        trials = find(group ==3 & lighting==1);
    elseif cond ==6
        trials = find(group ==3 & lighting==0);
    end
    
    
    %% bar graph of approach frequency and success
    clear s
%     for j =1:length(trials) 
%     s(j)=Sx{trials(j)}(1);
%     end

%     FM{cond}=s; %whether subject in each trial is female or male
    
    R{cond}=vertcat(rangeSmooth{trials});
    T{cond}=vertcat(thetaSmoothC{trials});
    Msp{cond}=vertcat(mSpeedS{trials});
    
    timeToCap{cond} = latency(trials);  
    timeToCapMu(cond)=nanmedian(latency(trials)');
    errTime(cond)=semedian(latency(trials)')
%     
    data_appFreq{cond}=freqApproach(trials);
    dataMu(cond,1) = nanmedian(freqApproach(trials)');
    err(cond,1) = semedian(freqApproach(trials)');
    
    data_Int{cond}=approachSuccess(trials);
    dataMu(cond,2) = nanmedian(approachSuccess(trials));
    err(cond,2)=semedian(approachSuccess(trials));
    
%     data_Cont{cond}=firstContact(trials);
%     ContactDataMu(cond,1)=nanmean(firstContact(trials)');
%     Contacterr(cond,1) = nanstd(firstContact(trials)')/sqrt(length(trials));
%     
%     ContactDataMu(cond,2)=nanmean(interContact(trials)');
%     Contacterr(cond,2) = nanstd(interContact(trials)')/sqrt(length(trials));
%     
%     ContactDataMu(cond,3)=nanmean(nContacts(trials)');
%     Contacterr(cond,3) = nanstd(nContacts(trials)')/sqrt(length(trials));
    
     data_ContDur{cond}=contactDuration(trials);
     ContactDataMu(cond,4)=nanmedian(contactDuration(trials)');
     Contacterr(cond,4) = semedian(contactDuration(trials)');
%     ContactDataMu(cond,4)=nanmean(contactDuration(trials)');
%     Contacterr(cond,4) = nanstd(contactDuration(trials)')/sqrt(length(trials));
%   
    data_CapProb{cond}=captureProb(trials);
    ContactDataMu(cond,5)=nanmedian(captureProb(trials)');
    Contacterr(cond,5) = semedian(captureProb(trials)');

%     Contacterr(cond,5) = nanstd(captureProb(trials)')/sqrt(length(trials));
%     ContactDataMu(cond,5)=nanmean(captureProb(trials)');
%     Contacterr(cond,5) = nanstd(captureProb(trials)')/sqrt(length(trials));
%    
    
    
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
            preTouch=find(x<4);%find places where mouse is within 5cm
            
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
%     preContT_mu(cond,1)=nanmedian(AppErr);
%     preContErr(cond,1)=semedian(AppErr);
    
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
bar(timeToCapMu(1,1:4)); hold on; errorbar(1:4,timeToCapMu(1,1:4),errTime(1,1:4),'o')
ylabel('time to capture')
set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'light','dark','deaf','dark + deaf'});

hold on
for ii=1:4
    tmp=timeToCap{ii}; %temporarily store data in variable "tmp"
    x = repmat(ii,1,length(tmp)); %the x axis location
    x = (x+(rand(size(x))-0.6)*0.1-.1); %add a little random "jitter" to aid visibility

    plot(x,tmp,'ok')
end
hold off

TC_light=horzcat(timeToCap{1}',repmat(1,length(timeToCap{1}),1));
TC_dark=horzcat(timeToCap{2}',repmat(2,length(timeToCap{2}),1));
TC_EPL=horzcat(timeToCap{3}',repmat(3,length(timeToCap{3}),1));
TC_EPD=horzcat(timeToCap{4}',repmat(4,length(timeToCap{4}),1));
TC_ESL=horzcat(timeToCap{5}',repmat(5,length(timeToCap{5}),1));
TC_ESD=horzcat(timeToCap{6}',repmat(6,length(timeToCap{6}),1));

TimeToCapDat=vertcat(TC_light,TC_dark,TC_EPL,TC_EPD,TC_ESL,TC_ESD);

[p,tbl,stats] = kruskalwallis(TimeToCapDat(:,1),TimeToCapDat(:,2))
c=multcompare(stats);

% [P,T,STATS,TERMS]=anovan(TimeToCapDat(:,1),TimeToCapDat(:,2)) ;%test whether there are differences by factor/groups
% results=multcompare(STATS,'CType','bonferroni')


%% Capture times blind controls

figure
bar(timeToCapMu(1,[1 2 5 6])); hold on; errorbar(1:4,timeToCapMu(1,[1 2 5 6]),errTime(1,[1 2 5 6]),'o')
ylabel('time_toCap_blindcontrols')
set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'light','dark','ES light','ES dark'});

hold on
trial=[1 2 5 6];
for ii=1:4
    tmp=timeToCap{trial(ii)}; %temporarily store data in variable "tmp"
    x = repmat(ii,1,length(tmp)); %the x axis location
    x = (x+(rand(size(x))-0.5)*0.1-.1); %add a little random "jitter" to aid visibility

    plot(x,tmp,'ok')
end
hold off


%% 
figure
bar(dataMu(1:4,1)); hold on; errorbar(1:4,dataMu(1:4,1),err(1:4,1),'o')
ylabel('approach frequency')
set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'light','dark','deaf','dark + deaf'});

hold on
for ii=1:4
    tmp=data_appFreq{ii}; %temporarily store data in variable "tmp"
    x = repmat(ii,1,length(tmp)); %the x axis location
    x = (x+(rand(size(x))-0.6)*0.1-.1); %add a little random "jitter" to aid visibility

    plot(x,tmp,'ok')
end
hold off

AF_light=horzcat(data_appFreq{1}',repmat(1,length(data_appFreq{1}),1));
AF_dark=horzcat(data_appFreq{2}',repmat(2,length(data_appFreq{2}),1));
AF_EPL=horzcat(data_appFreq{3}',repmat(3,length(data_appFreq{3}),1));
AF_EPD=horzcat(data_appFreq{4}',repmat(4,length(data_appFreq{4}),1));
AF_ESL=horzcat(data_appFreq{5}',repmat(5,length(data_appFreq{5}),1));
AF_ESD=horzcat(data_appFreq{6}',repmat(6,length(data_appFreq{6}),1));

AppFreq=vertcat(AF_light,AF_dark,AF_EPL,AF_EPD,AF_ESL,AF_ESD);

[p,tbl,stats] = kruskalwallis(AppFreq(:,1),AppFreq(:,2))
c=multcompare(stats);



%% prey intercept
figure
bar(dataMu(1:4,2)); hold on; errorbar(1:4,dataMu(1:4,2),err(1:4,2),'o')
ylabel('prey intercept')
set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'light','dark','deaf','dark + deaf'});
title('figure 3D bottom')

 hold on
for ii=1:4
    tmp=(data_Int{ii}+(rand(size(data_Int{ii}))+0.2)*0.05);
    x = repmat(ii,1,length(tmp)); %the x axis location
    x = (x+(rand(size(x))-0.5)*0.1-.1); %add a little random "jitter" to aid visibility

    plot(x,tmp,'ok')
end
hold off


PI_light=horzcat(data_Int{1}',repmat(1,length(data_Int{1}),1));
PI_dark=horzcat(data_Int{2}',repmat(2,length(data_Int{2}),1));
PI_EPL=horzcat(data_Int{3}',repmat(3,length(data_Int{3}),1));
PI_EPD=horzcat(data_Int{4}',repmat(4,length(data_Int{4}),1));
preyInt=vertcat(PI_light,PI_dark,PI_EPL,PI_EPD);

[p,tbl,stats] = kruskalwallis(preyInt(:,1),preyInt(:,2))
c=multcompare(stats);

% figure
% bar(ContactDataMu(1:4,1)); hold on; errorbar(1:4,ContactDataMu(1:4,1),Contacterr(1:4,1),'o')
% ylabel('FirstContactTime')
% set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'light','dark','deaf','dark + deaf'});
% title('figure 3E top')
% 
% hold on
% for ii=1:4
%     tmp=data_Cont{ii}; %temporarily store data in variable "tmp"
%     x = repmat(ii,1,length(tmp)); %the x axis location
%     x = (x+(rand(size(x))-0.5)*0.1-.1); %add a little random "jitter" to aid visibility
% 
%     plot(x,tmp,'ok')
% end
% hold off

% figure
% bar(ContactDataMu(1:4,4)); hold on; errorbar(1:4,ContactDataMu(1:4,4),Contacterr(1:4,4),'o')
% ylabel('Contact duration')
% set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'light','dark','deaf','dark + deaf'});
% title('fig3')
% 
% hold on
% for ii=1:4
%     tmp=data_ContDur{ii}; %temporarily store data in variable "tmp"
%     x = repmat(ii,1,length(tmp)); %the x axis location
%     x = (x+(rand(size(x))-0.5)*0.1-.1); %add a little random "jitter" to aid visibility
% 
%     plot(x,tmp,'ok')
% end
% hold off
% 
% CD_light=horzcat(data_ContDur{1}',repmat(1,length(data_ContDur{1}),1));
% CD_dark=horzcat(data_ContDur{2}',repmat(2,length(data_ContDur{2}),1));
% CD_EPL=horzcat(data_ContDur{3}',repmat(3,length(data_ContDur{3}),1));
% CD_EPD=horzcat(data_ContDur{4}',repmat(4,length(data_ContDur{4}),1));
% ContDur=vertcat(CD_light,CD_dark,CD_EPL,CD_EPD);
% 
% [p,tbl,stats] = kruskalwallis(ContDur(:,1),ContDur(:,2))
% c=multcompare(stats);


%% contact duration    
figure
bar(ContactDataMu(1:4,5)); hold on; errorbar(1:4,ContactDataMu(1:4,5),Contacterr(1:4,5),'o')
ylabel('Cap Probability 1/numConacts')
set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'light','dark','deaf','dark + deaf'});
title('fig3')
hold on
for ii=1:4
    
    tmp=(data_CapProb{ii}+(rand(size(data_CapProb{ii}))+0.1)*0.05)
    x = repmat(ii,1,length(tmp)); %the x axis location
    x = (x+(rand(size(x))-0.5)*0.1-.1); %add a little random "jitter" to aid visibility

    plot(x,tmp,'ok')
end
hold off
CP_light=horzcat(data_CapProb{1}',repmat(1,length(data_CapProb{1}),1));
CP_dark=horzcat(data_CapProb{2}',repmat(2,length(data_CapProb{2}),1));
CP_EPL=horzcat(data_CapProb{3}',repmat(3,length(data_CapProb{3}),1));
CP_EPD=horzcat(data_CapProb{4}',repmat(4,length(data_CapProb{4}),1));
CapP=vertcat(CP_light,CP_dark,CP_EPL,CP_EPD);

[p,tbl,stats] = kruskalwallis(CapP(:,1),CapP(:,2))
c=multcompare(stats);


%% PA
figure
bar(preContT_mu(1:4,1)); hold on; errorbar(1:4,preContT_mu(1:4,1),preContErr(1:4,1),'o')
ylabel('Pre Contact prey azimuth')
set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'light','dark','deaf','dark + deaf'});

hold on
for ii=1:4
    tmp=(preContT{ii}+(rand(size(preContT{ii}))+ 0.1)*0.05); %temporarily store data in variable "tmp"
    x = repmat(ii,1,length(tmp)); %the x axis location
    x = (x+(rand(size(x))-0.5)*0.1-.1); %add a little random "jitter" to aid visibility

    plot(x,tmp,'ok')
end
hold off

PA_light=horzcat(preContT{1}',repmat(1,length(preContT{1}),1));
PA_dark=horzcat(preContT{2}',repmat(2,length(preContT{2}),1));
PA_EPL=horzcat(preContT{3}',repmat(3,length(preContT{3}),1));
PA_EPD=horzcat(preContT{4}',repmat(4,length(preContT{4}),1));
PA_ESL=horzcat(preContT{5}',repmat(5,length(preContT{5}),1));
PA_ESD=horzcat(preContT{6}',repmat(6,length(preContT{6}),1));

PAmain=vertcat(PA_light,PA_dark,PA_EPL,PA_EPD, PA_ESL,PA_ESD);

[P,T,STATS,TERMS]=anovan(PAmain(:,1),PAmain(:,2)) ;%test whether there are differences by factor/groups
 results=multcompare(STATS)



    %extract out non NAN values for ranges and theata during approaches

    H1=find(~isnan(A{1}(:,2))& (A{1}(:,1))<25 & (A{1}(:,2))<100);
    
    g=A{1}(H1,1);
    h=A{1}(H1,2);
    [r p]=corrcoef(g,h)

    P = polyfit(g,h,1);
    f = polyval(P,g); 
    [r2 rmse] = rsquare(h,f); 
    figure; plot(g,h,'b*'); 
    hold on; plot(g,f,'y-','LineWidth',2); 
    title(strcat(['R2 = ' num2str(r2) '; RMSE = ' num2str(rmse)])) 
    
   H2=find(~isnan(A{2}(:,2))& (A{2}(:,1))<25& (A{2}(:,2))<100 );
    
    g2=A{2}(H2,1);
    h2=A{2}(H2,2);
    [r p]=corrcoef(g2,h2)

    P = polyfit(g2,h2,1);
    f = polyval(P,g2); 
    [r2 rmse] = rsquare(h2,f); 
    figure; plot(g2,h2,'k*'); 
    hold on; plot(g2,f,'y-','LineWidth',2); 
    title(strcat(['R2 = ' num2str(r2) '; RMSE = ' num2str(rmse)])) 
    
   H3=find(~isnan(A{3}(:,2))& (A{3}(:,1))<25 & (A{3}(:,2))<60);
    
    g3=A{3}(H3,1);
    h3=A{3}(H3,2);
    [r p]=corrcoef(g3,h3)

    P = polyfit(g3,h3,1);
    f = polyval(P,g3); 
    [r2 rmse] = rsquare(h3,f); 
    figure; plot(g3,h3,'g*'); 
    hold on; plot(g3,f,'y-','LineWidth',2); 
    title(strcat(['R2 = ' num2str(r2) '; RMSE = ' num2str(rmse)])) 
    
    H4=find(~isnan(A{4}(:,2))& (A{4}(:,1))<25 & (A{4}(:,2))<100);
    g4=A{4}(H4,1);
    h4=A{4}(H4,2);
    [r p]=corrcoef(g4,h4)

    P = polyfit(g4,h4,1);
    f = polyval(P,g4); 
    [r2 rmse] = rsquare(h4,f); 
    figure; plot(g4,h4,'r*'); 
    hold on; plot(g4,f,'y-','LineWidth',2); 
    title(strcat(['R2 = ' num2str(r2) '; RMSE = ' num2str(rmse)]))
    
%corrcoef(A,B)
%fitlm(patients,modelspec)

% [p,tbl,stats] = anova(PAmain(:,1),PAmain(:,2))
% c=multcompare(stats);

%% PA blind
figure
bar(preContT_mu([1 2 5 6],1)); hold on; errorbar(1:4,preContT_mu([1 2 5 6],1),preContErr([1 2 5 6],1),'o')
ylabel('Pre Contact prey azimuth eye suture controls')
set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'light','dark','sutured','dark + sutured'});

hold on
trial=[1 2 5 6];
for ii=1:4
    tmp=(preContT{trial(ii)}+(rand(size(preContT{trial(ii)}))+0.2)*0.05); %temporarily store data in variable "tmp"
    x = repmat(ii,1,length(tmp)); %the x axis location
    x = (x+(rand(size(x))-0.5)*0.1-.1); %add a little random "jitter" to aid visibility

    plot(x,tmp,'ok')
end
hold off


% figure
% bar(ContactData(:,3)); hold on; errorbar(1:4,ContactData(:,3),Contacterr(:,3),'o')
% ylabel('Num of contacts')
% set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'light','dark','deaf','dark + deaf'});
% title('figure 3F top')
%plot all condition theta

figure
plot(thetabins,nanmedian(PreyE_mid{2},2),'k','LineWidth',3); hold on
plot(thetabins,nanmedian(PreyE_mid{1},2),'b','LineWidth',3);
plot(thetabins,nanmedian(PreyE_mid{3},2),'g','LineWidth',3);
plot(thetabins,nanmedian(PreyE_mid{4},2),'r','LineWidth',3);

title 'Prey eccentricity all conditions >5cm mouse moving'

figure
plot(thetabins,nanmedian(PreyE_mid{1},2),'b','LineWidth',3); hold on%eye suture light
plot(thetabins,nanmedian(PreyE_mid{2},2),'k','LineWidth',3);%eye suture light
plot(thetabins,nanmedian(PreyE_mid{5},2),'c','LineWidth',3); %eye suture light
plot(thetabins,nanmedian(PreyE_mid{6},2),'m','LineWidth',3); %eye suture


% PLOT ALL CONDITIONS range over trial
figure
plot(nanmedian(Rsession{2},2),'k','LineWidth',3);hold on
plot(nanmedian(Rsession{1},2),'b','LineWidth',3);hold on
plot(nanmedian(Rsession{3},2),'g','LineWidth',3);hold on
plot(nanmedian(Rsession{4},2),'r','LineWidth',3);

figure
plot(nanmedian(Spsession{2},2),'k','LineWidth',3);hold on
plot(nanmedian(Spsession{1},2),'b','LineWidth',3);hold on
plot(nanmedian(Spsession{3},2),'g','LineWidth',3);hold on
plot(nanmedian(Spsession{4},2),'r','LineWidth',3);


for i=1:6
moving=Msp{i}>2;
ag(i)=nanmean(Msp{i}(moving));
end

for i=1:6
moving=Msp{i}>2;
spMed(i)=nanstd(Msp{i}(moving))/sqrt(20);
end

figure
bar(ag(1:4)); hold on; errorbar(1:4,ag(1:4),spMed(1:4),'o')
ylabel('Avg Speed per  session')



keyboard
close all
keyboard

%% plot first contact time, number of contacts and inter-contactinterval

WT=1; Dark=2;EyeSuture=3; BlkBack = 4; EarPlug = 5; Whisk=6; Plexi=7; Reopen=8; Unplug=9; rd=10; rdhet=11;
thetaFig = figure; set(gcf,'Name','theta');
distFig = figure; set(gcf,'Name','distance');
speedFig = figure; set(gcf,'Name','speed');
dthistFig = figure; set(gcf,'Name','dist and theta');
dthist_mvFig = figure; set(gcf,'Name','dist and theta moving');
spdistFig = figure; set(gcf,'Name','speed vs dist')
stateFig = figure; set(gcf,'Name','state')
dt_curveFig = figure; set(gcf,'Name','width vs distance')
ncond=4;

clear midthetaAll midtheta_mvAll speedAll distAll tdhistAll tdhist_mvAll spdistAll indAll
for m=1:ncond;
    for light=1:2;
        trials =find(group == m & lighting == light-1);
        figure; set(gcf,'Name',sprintf('%s %s',grouplabels{m+1},lightlabels{light}));
        for tr =1:length(trials)
            
            subplot(8,8,tr);
            hold on
            %             for i = 1:8:length(targ{trials(tr)})-8
            %                 plot([ targ{trials(tr)}(i,1) targ{trials(tr)}(i+8,1)] ,[ targ{trials(tr)}(i,2) targ{trials(tr)}(i+8,2) ],'LineWidth',2,'Color',cmapVar(i,1,length(targ{trials(tr)}),jet ));
            %                 plot(head{trials(tr)}(i,1),head{trials(tr)}(i,2),'o','Color',cmapVar(i,1,length(targ{trials(tr)}),jet ) );
            %             end
            
            plot(targ{trials(tr)}(:,1),targ{trials(tr)}(:,2),'g');
            plot(head{trials(tr)}(:,1),head{trials(tr)}(:,2),'b');
            axis([0 2000 0 1200]); axis off
            
        end
        
%         timeRange = 120;
%         hit_rg=[]; hit_th=[]; miss_rg=[]; miss_th=[];
% %         figure;  set(gcf,'Name',sprintf('capture %s %s',grouplabels{m+1},lightlabels{light}));
%         for tr = 1:length(trials)
%             rg = basin{trials(tr)}(:,1); thet = basin{trials(tr)}(:,2); nxt = basin{trials(tr)}(:,3);
%             
%             hold on
%             plot(rg(nxt>timeRange),abs(thet(nxt>timeRange)),'r.')
%             plot(rg(nxt<timeRange),abs(thet(nxt<timeRange)),'g.')
%             hit_rg  = [hit_rg ; rg(nxt<timeRange)]; hit_th = [hit_th ; abs(thet(nxt<timeRange))];
%             miss_rg  = [miss_rg ; rg(nxt>timeRange)]; miss_th = [miss_th ; abs(thet(nxt>timeRange))];
%         end
%         
%         hits = myHist2(hit_rg,hit_th,2.5:5:50,15:30:180);
%         miss =  myHist2(miss_rg,miss_th,2.5:5:50,15:30:180);
%         %         figure
%         %         imagesc(hits)
%         %         figure
%         %         imagesc(miss)
%         %
%         figure
%         imagesc(hits./(hits+miss),[0 1]);
%         axis xy
%         
%         
%         approachhist=0; inds =0;
%         figure
%         hold on
%         col = 'rgbcmk'; set(gcf,'Name',sprintf('capture %s %s',grouplabels{m+1},lightlabels{light}));
%         ntr=0;
%         for tr= 1:length(trials)
%             for j = 1:size(approachR{trials(tr)},1)
%                 
%                 for i = 1:length(approachTheta{trials(tr)}(j,:));
%                     straight(i) = all(abs(approachTheta{trials(tr)}(j,i:end))<60);
%                 end
%                 
%                 [y ind] = max(approachR{trials(tr)}(j,:));
%                 if y>5
%                     ntr=ntr+1;
%                     plot(approachR{trials(tr)}(j,ind), approachTheta{trials(tr)}(j,ind),'bo');
%                     for k = ind:59
%                         plot([approachR{trials(tr)}(j,k) approachR{trials(tr)}(j,k+1)], [approachTheta{trials(tr)}(j,k) approachTheta{trials(tr)}(j,k+1)] ,'Color',cmapVar(k,1,59,jet));
%                         %plot(approachSpeed{trials(tr)}(i,:));
%                     end
%                     approachD = find( approachR{trials(tr)}(j,ind:end)>5 & approachR{trials(tr)}(j,ind:end)<10)+ind-1;
%                     bins = 10:20:180;
%                     approachhist(ntr,1:length(bins)) = hist(abs(approachTheta{trials(tr)}(j,approachD)),bins)/length(approachD);
%                 end
%                 
%             end
%         end
%         axis([ 2.5 30 -180 180])
        
        
        bins = 10:20:180;
        
        
        %%% get distribution of starting points
%         inds =0;
%         ntr=0;
%         for tr= 1:length(trials)
%             for j = 1:size(approachR{trials(tr)},1)
%                 for i = 1:length(approachTheta{trials(tr)}(j,:));
%                     straight(i) = all(abs(approachTheta{trials(tr)}(j,i:end))<60);
%                 end
%                 [y ind] = max(approachR{trials(tr)}(j,:).*straight); %%% find most distant point where always aimed at target (abs(theta)<90)
%                 if y>5
%                     ntr=ntr+1;
%                     inds(ntr)=y;
%                 end
%             end
%         end
%         
%         
%         
%         indAll(m,light,1:length(dbins)) = hist(inds,dbins)/length(inds);
%         
%         figure
%         figure; set(gcf,'Name',sprintf('approach %s %s',grouplabels{m+1},lightlabels{light}));
%         plot(dbins,squeeze(indAll(m,light,:))); xlim([ 0 30]); xlabel('distance of approach initiation (cm)');
%         
%         
%         
%         approachhist = nanmean(approachhist,1);
%         
%         figure; set(gcf,'Name',sprintf('approach %s %s',grouplabels{m+1},lightlabels{light}));
%         plot(bins,approachhist); ylim([0 1]); xlim([0 90])
%         approachHistAll(m,light,:) = approachhist; xlabel('theta 5- 10 cmn')
%         %
%         %         escapehist=0;
        %         figure
        %         hold on
        %         col = 'rgbcmk'; set(gcf,'Name',sprintf('escape %s %s',grouplabels{m+1},lightlabels{light}));
        %         ntr=0;
        %         for tr= 1:length(trials)
        %             for j = 1:size(escapeR{trials(tr)},1)
        %                 ntr=ntr+1;
        %                 for k = 1:59
        %                     plot([escapeR{trials(tr)}(j,k) escapeR{trials(tr)}(j,k+1)], [escapeTheta{trials(tr)}(j,k) escapeTheta{trials(tr)}(j,k+1)] ,'Color',cmapVar(k,1,59,jet));
        %
        %                 end
        %                 escapeD = escapeR{trials(tr)}(j,:);
        %                 bins = 15:30:180;
        %                 escapehist(ntr,1:length(bins))= hist(abs(escapeTheta{trials(tr)}(j,escapeD>5 & escapeD<10)),bins)/sum(escapeD>5 & escapeD<10);
        %
        %             end
        %         end
        %         axis([ 2.5 30 -180 180])
        %         escapehist = nanmean(escapehist,1);
        %
        %         figure; set(gcf,'Name',sprintf('escape %s %s',grouplabels{m+1},lightlabels{light}));
        %         plot(escapehist); ylim([0 0.6])
        %
        %         figure
        %         hold on
        %         col = 'rgbcmk'; set(gcf,'Name',sprintf(' capture %s %s',grouplabels{m+1},lightlabels{light}));
        %         ntr=0;
        %         for tr= 1:length(trials)
        %             for j = 1:size(approachSpeedTarg{trials(tr)},1)
        %                 plot(approachSpeed{trials(tr)}(j,:),'g'); plot(approachSpeedTarg{trials(tr)}(j,:),'r');
        %             end
        %         end
        
        %         figure
        %         hold on
        %         col = 'rgbcmk'; set(gcf,'Name',sprintf(' capture %s %s',grouplabels{m+1},lightlabels{light}));
        %         ntr=0;
        %         for tr= 1:length(trials)
        %             for j = 1:size(approachSpeedTarg{trials(tr)},1)
        %                 [y ind] = max( approachSpeedTarg{trials(tr)}(j,:));
        %                 plot(ind,y,'o')
        %             end
        %         end
        %         xlabel('time'); ylabel('max target speed')
        %
        %         figure
        %         hold on
        %         col = 'rgbcmk'; set(gcf,'Name',sprintf(' escape %s %s',grouplabels{m+1},lightlabels{light}));
        %         ntr=0;
        %         for tr= 1:length(trials)
        %             for j = 1:size(escapeTarg{trials(tr)},1)
        %                 plot(escapeSpeed{trials(tr)}(j,:),'g'); plot(escapeTarg{trials(tr)}(j,:),'r');
        %             end
        %         end
        %
        %         figure
        %         hold on
        %         col = 'rgbcmk'; set(gcf,'Name',sprintf(' escape %s %s',grouplabels{m+1},lightlabels{light}));
        %         ntr=0;
        %         for tr= 1:length(trials)
        %             for j = 1:size(escapeTarg{trials(tr)},1)
        %                 plot(escapeSpeed{trials(tr)}(j,10), escapeTarg{trials(tr)}(j,10),'o');
        %             end
        %         end
        %         plot([0 50],[0 50]); axis square; axis([ -5 75 -5 75])
        %
        
        
        %% % if data has been extracted for each trial, compile or average
        %%% here. List of trials in this group is 'trials'
        
        %%% e.g. nanmean(data(trials))  or  for tr=1:length(trials)
        
        
        
        contactDurAll(m,light) = nanmean(contactDuration(trials));
        interContactAll(m,light) = nanmean(interContact(trials));
        nContactAll(m,light) = nanmean(nContacts(trials));
        firstContactAll(m,light) = nanmean(firstContact(trials));
        
        
        ntr(m,light)=length(trials);
        %tAll(m,light) = mean(captureT(trials))
        %midthetaAll(m,light,:) = nanmean(midtheta(:,trials),2);
        %midtheta_mvAll(m,light,:) = nanmean(midtheta_mv(:,trials),2);
        speedAll(m,light,:) = nanmean(speedhist(:,trials),2);
        distAll(m,light,:) = nanmean(rhist(:,trials),2);
        %tdhistAll(m,light,:,:) = nanmean(tdhist(:,:,trials),3);
        %tdhist_mvAll(m,light,:,:) = nanmean(tdhist_mv(:,:,trials),3);
        spdistAll(m,light,:,:) = nanmean(spdist(:,:,trials),3);
        stateAll(m,light,:) = nanmean(state(:,trials),2);
        
%         td = squeeze(tdhist_mvAll(m,light,:,:));
%         for d = 1:size(td,2);
%             tdcurve = td(:,d); tdcurve(tdcurve<max(tdcurve)*0.1)=0;
%             meantheta(m,light,d) = sum(tdcurve.*thetabins')/sum(tdcurve);
%             stdtheta(m,light,d) = sqrt(sum(tdcurve.*(thetabins.^2)')/sum(tdcurve));
%         end
        
        
        figure(distFig); subplot(ncond,2,2*(m-1)+(3-light)); plot(dbins,squeeze(distAll(m,light,:))); ylim([0 0.5])
       % figure(thetaFig); subplot(ncond,2,2*(m-1)+(3-light)); plot(thetabins,squeeze(midthetaAll(m,light,:))); hold on; plot(thetabins,squeeze(midtheta_mvAll(m,light,:)),'g'); ylim([0 0.5]); xlim([-180 180])
        figure(speedFig); subplot(ncond,2,2*(m-1)+(3-light)); plot(spdbins,squeeze(speedAll(m,light,:))); ylim([0 0.5])
        %figure(dthistFig); subplot(ncond,2,2*(m-1)+(3-light)); imagesc(squeeze(tdhistAll(m,light,:,:)),[0 0.05]);
        %figure(dt_curveFig); subplot(ncond,2,2*(m-1)+(3-light));  plot(squeeze(stdtheta(m,light,:))); ylim([0 120]); xlim([0 20])
        %figure(dthist_mvFig); subplot(ncond,2,2*(m-1)+(3-light));   imagesc(squeeze(tdhist_mvAll(m,light,:,:)),[0 0.05]);
        figure(spdistFig); subplot(ncond,2,2*(m-1)+(3-light)); imagesc(squeeze(spdistAll(m,light,:,:)),[0 0.1]);
        figure(stateFig); subplot(ncond,2,2*(m-1)+(3-light));  bar(squeeze(stateAll(m,light,:))); ylim([0 1])
        
        
    end
end

figure
hold on
plot(bins,squeeze(approachHistAll(1,2,:)),'b'); sqrt(sum((bins.^2)'.*squeeze(approachHistAll(1,2,:))))

plot(bins,squeeze(mean(approachHistAll(2,:,:),2)),'r');  sum(bins'.*squeeze(mean(approachHistAll(2,:,:),2)))
plot(bins,squeeze(approachHistAll(1,1,:)),'r'); sqrt(sum((bins.^2)'.*squeeze(approachHistAll(1,1,:))))
plot(bins,squeeze(approachHistAll(4,2,:)),'g'); sqrt(sum((bins.^2)'.*squeeze(approachHistAll(4,2,:))))
legend('light','suture','earplug'); axis([0 90 0 0.75])





keyboard



