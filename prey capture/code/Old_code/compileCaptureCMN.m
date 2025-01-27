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
    latency(f) = files(f).FrameEnd-files(f).FrameS;
    fn = [pathname files(f).trackpts];
    if ismac
        fn(fn=='\')='/';
    end
    [ location(:,:,f) targ{f} body{f} r{f} thetaR{f} thetaSac{f} vhead{f} vbody{f} abody{f} vtarg{f} atarg{f} head{f} sub{f}] = analyzeCapture(fn,files(f).fps,files(f).scale,files(f).subj);
    dist(f) = mean(r{f});
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

dEye = [mean(latency(lighting==0 & group==2)),...
    mean(latency(lighting==1 & group==3)),...
    mean(latency(lighting==0 & group==3))]/60;  

errEye = [std(latency(lighting==0 & group==2))/sqrt(sum(lighting==0 & group==2)),...% dark exposure for the first time group, alternated with light trials
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
legend('Light','dark', 'EP light','EP Dark');

figure
% col=[1,1,1,0.2]
barweb(dEye,errEye,'b')
ylabel('time to capture (s)')

clear rhist speedhist thetahist tdhist tdhist_mv midtheta midtheta_mv spdist basin dRhist approachR
dbins = 1.25:2.5:45; spdbins = 1.25:2.5:50; thetabins = linspace(-150, 150, 13); theta2bins = -180:30:150;
dbinsMid=5:2.5:45;
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
    
    if makeMov
        close all
    end
    
    group(tr)=files(tr).group
    
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
        captureT(tr) = dur/files(tr).fps;
        
        %%% plot head orientation on top of tracks
        %         if (group(tr)==2 & lighting(tr)==1) | (group(tr)==4 & lighting(tr)==1)
        %             figure
        %             hold on
        %
        %
        %             plot(targ{tr}(:,1),targ{tr}(:,2),'k');
        %             for t = 1:3:dur
        %
        %                 %plot(targ{tr}(t,1),targ{tr}(t,2),'go');
        %                 %plot(head{tr}(t,1),head{tr}(t,2),'ko');
        %                 rangeangle = atan2( targ{tr}(t,2) - head{tr}(t,2), targ{tr}(t,1) - head{tr}(t,1))*180/pi;
        %                 headangle = rangeangle -thetaR{tr}(t);
        %                 if r{tr}(t)<10
        %                     plot([head{tr}(t,1) head{tr}(t,1) + cosd(headangle)*50 ], [head{tr}(t,2) head{tr}(t,2) + sind(headangle)*50],'r');
        %                     %  plot(head{tr}(t,1),head{tr}(t,2),'r.','Markersize',8)
        %                 else
        %                     plot([head{tr}(t,1) head{tr}(t,1) + cosd(headangle)*50 ], [head{tr}(t,2) head{tr}(t,2) + sind(headangle)*50],'b');
        %                     % plot(head{tr}(t,1),head{tr}(t,2),'b.','Markersize',8)
        %                 end
        %                 plot([head{tr}(t,1) head{tr}(t,1) + cosd(rangeangle)*50 ], [head{tr}(t,2) head{tr}(t,2) + sind(rangeangle)*50],'g');
        %                 %             plot(head{tr}(1:t,1),head{tr}(1:t,2))
        %                 %             plot(targ{tr}(1:t,1),targ{tr}(1:t,2),'r');
        %                 %             axis([0 2000 0 1200])
        %                 %             mov(t) = getframe(gcf);
        %             end
        %             axis([0 2000 0 1200])
        %             title(sprintf('trial %d %s %s',tr,grouplabels{group(tr)+1},lightlabels{lighting(tr)+1}));
        %         end
        %
        
        %create a list of sub IDs to align with all instances of a measure
        s=sub{tr}
        sID={s}
        ID{tr}=repmat(sID,[dur,1]);
        
        sx=sex(tr)
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
        
        %target speeds over whole session
        mSpeedS{tr}=speed;    
        speedhist(:,tr) = hist(speed,spdbins)/dur;
        
        filt = ones(5,1); filt = filt/length(filt);
        vx = filter(filt,1,diff(targ{tr}(:,1))); vy = filter(filt,1,diff(targ{tr}(:,2))); speedtarg = sqrt(vx.^2 + vy.^2) * files(tr).fps / files(tr).scale;
        speedtarg = speedtarg(1:dur-1); speedtarg(1:end-2) = speedtarg(3:end);
        
        tSpeedS{tr}=speedtarg; 
        
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
    
    
    contact = range<3.75;
    
    
    contact = medfilt2(contact,[10 1]);
    contact(1)=0; contact(end)=0;
    
    
    %plot(contact*0.2+0.1,'k');
    
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
            if peak(j)<10
            contact(1:contactStart(j))=1;
            end
        
        else 
        peak(j)=max(rangeS(contactEnd(j-1):contactStart(j)));
            if peak(j)<10
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
    if length(contactStart)>1
        interContact(tr) = mean(contactStart(2:end) - contactEnd(1:end-1));
        
    else
        interContact(tr)=NaN;
        
    end
    
    %% Cris's contact definition
%     if length(contactStart)>0
%         for i = 1:length(contactStart)
%             grabFrames = 60; %1 sec prior to contact
%             contactFrame = round(contactStart(i)*files(f).fps);
%             nFrames = min(grabFrames,contactFrame);
%             startFrame = max(1,contactFrame-grabFrames+1);
%             %define approach range, theta, predator speed and prey speed in
%             %1 sec prior to contact
%             approachR{tr}(i,1:grabFrames) = NaN; approachTheta{tr}(i,1:grabFrames) = NaN; approachSpeed{tr}(i,1:grabFrames) = NaN; approachSpeedTarg{tr}(i,1:grabFrames) = NaN;  %%% fill with Nans for short segements
%             approachR{tr}(i,(grabFrames-nFrames+1) : grabFrames) = rangeS(startFrame:contactFrame);
%             approachTheta{tr}(i,(grabFrames-nFrames+1) : grabFrames) = thetaSC(startFrame:contactFrame);
%             approachSpeed{tr}(i,(grabFrames-nFrames+1) : grabFrames) = speed(startFrame:contactFrame);
%             approachSpeedTarg{tr}(i,(grabFrames-nFrames+1) : grabFrames) = speedtarg(startFrame:contactFrame);
%             
%             
%             startFrame = round(contactEnd(i)*files(f).fps);
%             if startFrame~=length(range)-1
%                 nFrames = min(grabFrames,length(speed)-startFrame+1);
%                 
%                 escapeR{tr}(i,1:grabFrames) = NaN; escapeTheta{tr}(i,1:grabFrames) = NaN; escapeSpeed{tr}(i,1:grabFrames) = NaN;  escapeTarg{tr}(i,1:grabFrames) = NaN; %%% fill with Nans for short segements
%                 escapeR{tr}(i,(grabFrames-nFrames+1) : grabFrames) = rangeS(startFrame:(startFrame+nFrames-1));
%                 escapeTheta{tr}(i,(grabFrames-nFrames+1) : grabFrames) = thetaSC(startFrame:(startFrame+nFrames-1));
%                 escapeSpeed{tr}(i,(grabFrames-nFrames+1) : grabFrames) = speed(startFrame:(startFrame+nFrames-1));
%                 escapeTarg{tr}(i,(grabFrames-nFrames+1) : grabFrames) = speedtarg(startFrame:(startFrame+nFrames-1));
%             end
%             
%             
%             if makeMov
%                 frameRange = ceil((contactStart(i)*files(tr).fps -120 : bin : contactEnd(i)*files(tr).fps + 120) /bin);
%                 frameRange = frameRange(frameRange>=1);
%                 if ~isempty(vid)
%                     frameRange = frameRange(frameRange<= size(vid,3));
%                 end
%                 frameRange = frameRange(frameRange*bin < length(thetaR{tr}));
%                 mov = trackMovie(vid,head{tr},targ{tr},thetaR{tr},frameRange ,bin,figName);
%                 
%                 %%% save movie
%                 vidObj = VideoWriter(sprintf('%s_contact%d.avi', fn(1:end-4),i));
%                 vidObj.FrameRate = 60/bin;
%                 vidObj.Quality = 100;
%                 open(vidObj);
%                 writeVideo(vidObj,mov);
%                 close(vidObj);
%             end
%             
%         end
%     else
%         approachR{tr} = [];
%         approachTheta{tr} = [];
%         approachSpeed{tr} = [];
%         
%         escapeR{tr} = [];
%         escapeTheta{tr} = [];
%         escapeSpeed{tr} = [];
%     end
% 

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
    
%     near = range<3;
%     nearStarts = find(diff(near)>0)+1;
%     nextContact=zeros(length(range),1);
%     for i = 1:length(range);
%         if ~near(i)
%             nexts = nearStarts-i; nexts=nexts(nexts>0);
%             m = min(nexts);
%             if ~isempty(m)
%                 nextContact(i) = m;
%             else
%                 nextContact(i)=NaN;
%             end
%         end
%     end
%     %
%     %     figure
%     %     plot(range(nextContact<60),th(nextContact<60),'g.')
%     %     hold on
%     %         plot(range(nextContact>60),th(nextContact>60),'r.')
%     
%     basin{tr}(:,1) = range; basin{tr}(:,2) = th(1:length(range)); basin{tr}(:,3)=nextContact;
    
    close(gcf)
    
    deltaR = diff(rangeS(1:dur));
    filt = ones(15,1);
    deltaR = filter(filt,1,deltaR); deltaR(1:end-6) = deltaR(7:end);
    
  
    
    dRhist(:,:,tr) = myHist2(speed,deltaR, 2.5:5:50,-19:2:20);
    %figure;imagesc(dRhist(:,:,tr))
   
    %% define appraoches based on mouse behavior 
    approach = (deltaR<-1 & speed>4); % | contact(1:end-1); %approach = medfilt2(approach,[30 1]);
    approach(1)=0; approach(end)=0;
    starts = find(diff(approach)>0); ends = find(diff(approach)<0);
    
    for i = 2:length(ends);
        if (starts(i) - ends(i-1))<60 & (rangeS(starts(i))-rangeS(ends(i-1)))<5  %%% stitch together ones that are close in time and continual approach
            approach(ends(i-1):starts(i))=1;
        else if (starts(i) - ends(i-1))<60 & (rangeS(starts(i))-rangeS(ends(i-1)))>5  %%% can't immediately start a new approach
                approach(starts(i):ends(i))=0;
            end
        end
    end
    
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

clear tr
%% light trials plot figures
trials = find(group ==1);

figure
for tr = 1:length(trials)
           %figure
            subplot(4,6,tr);
            hold on
            %             for i = 1:8:length(targ{trials(tr)})-8
            %                 plot([ targ{trials(tr)}(i,1) targ{trials(tr)}(i+8,1)] ,[ targ{trials(tr)}(i,2) targ{trials(tr)}(i+8,2) ],'LineWidth',2,'Color',cmapVar(i,1,length(targ{trials(tr)}),jet ));
            %                 plot(head{trials(tr)}(i,1),head{trials(tr)}(i,2),'o','Color',cmapVar(i,1,length(targ{trials(tr)}),jet ) );
            %             end
            
            plot(targ{trials(tr)}(:,1),targ{trials(tr)}(:,2),'k');
            plot(head{trials(tr)}(:,1),head{trials(tr)}(:,2),'c');
            axis([0 2000 0 1200]); axis off  
end
 title 'light trials'

 %%
% %save out data for stats into excel files
% 
% %theta and ranges over whole trial by session type
% Light(:,1)=vertcat(thetaS{trials});
% Light(:,2)=vertcat(rangeS{trials});
% Light(:,3) = repmat(1,length(Light(:,2)),1);
% 
% %speeds over whole trial by session type
% sLight(:,1)=vertcat(tSpeedS{trials});
% sLight(:,2)=vertcat(mSpeedS{trials});
% 
% xlswrite('allThetaRangeSpeedsLight', Light,1,'A2');
% xlswrite('allThetaRangeSpeedsLight',vertcat(ID{trials}),1,'D2');
% xlswrite('allThetaRangeSpeedsLight',vertcat(Sx{trials}) ,1,'E2');
% 
% name = {'allTheta','allRange','Group','ID','Sex'};
% xlswrite('allThetaRangeSpeedsLight', sLight,2,'A2');
% name2={'targSpeed','mouseSpeed'};
% xlswrite('allThetaRangeSpeedsLight', name,1,'A1' );
% xlswrite('allThetaRangeSpeedsLight',name2,2,'A1');
% 
clear err
%% bar graph of approach frequency and success

data(1,1) = nanmean(freqApproach(trials)');
err(1,1) = nanstd(freqApproach(trials)')/sqrt(length(trials));
data(1,2) = nanmean(approachSuccess(trials));
err(1,2) = nanstd(approachSuccess(trials))/sqrt(length(trials));

ContactData(1,1)=nanmean(firstContact(trials)');
Contacterr(1,1) = nanstd(firstContact(trials)')/sqrt(length(trials));
ContactData(1,2)=nanmean(interContact(trials)');
Contacterr(1,2) = nanstd(interContact(trials)')/sqrt(length(trials));
ContactData(1,3)=nanmean(nContacts(trials)');
Contacterr(1,3) = nanstd(nContacts(trials)')/sqrt(length(trials));

%% change to account for eccentricity over middle distance
% prey eccentricity over whole session
PreyEsession=thetahist(:,trials);
figure
plot(nanmean(PreyEsession,2),'b','LineWidth',3);
title 'Prey eccentricity whole session with eharing and vision-Light'

% PreyE_mid_L=thetahist_mid(:,trials);
% figure
% plot(thetabins,nanmean(PreyE_mid_L,2),'b','LineWidth',3);
% title 'Prey eccentricity whole session with eharing and vision-Light'

PreyE_mid_L=thetahist_midMove(:,trials);
figure
plot(thetabins,nanmean(PreyE_mid_L,2),'b','LineWidth',3);
title 'Prey eccentricity whole session with eharing and vision-Light'

%range over whole session
Rsession_L=rhist(:,trials);
figure
plot(dbins,nanmean(Rsession_L,2),'b','LineWidth',3);
title 'range whole session with hearing and vision-Light'

% theta x range between mice, moving 
TR=tdhist_mv(:,:,trials);
figure
imagesc(squeeze(nanmean(TR,3)));
title ' bearing(theta) by range when mouse is moving with hearing and vision-Light'

%approach paths for light
clear goodApp AppEcc
goodApp=0;

for tr = 1:length(trials);
   
            for i = 1:length(AppR{trials(tr)});
                
                       clear x y maxX
                       x=AppR{trials(tr)}{:,i};
                       y=AppT{trials(tr)}{:,i};
                       preTouch=find(x<3);
                       if ~isempty (preTouch);
                       dist=preTouch(1);
                       else 
                           mindist=min(x);
                           dist=find(x==mindist);
                       end
                       %minX=min(x);dist=find(x==minX);
                       
                       goodApp=goodApp+1;
                       
                       AppEcc{goodApp}=[x(1:dist),abs(y(1:dist))];               
                  %[approachR{trials(tr)}(i,dist:end)',approachTheta{trials(tr)}(i,dist:end)']
            end 
end


figure; hold on
col=[0,0,1,0.2]
for i=1:length(AppEcc)
    if ~isempty(AppEcc{i})
    plot(AppEcc{i}(:,1),AppEcc{i}(:,2),'color',col,'LineWidth',0.5)
    end
end
set(gca,'xdir','reverse','XLim',[1 35],'YLim',[0 180]);


binD=3:3:24;
binE=-30:5:30; 

clear thetaRangeHist avg_theta

for i=1:length(AppEcc);
     if ~isempty(AppEcc{i})
 [thetaRangeHist(:,:,i),avg_theta(i,:)]=myHist2Avg(AppEcc{i}(:,1),AppEcc{i}(:,2),binD,binE);
     end
end
 
%plot stad deviation plot
 figure; hold on
 
      stdD=nanstd(avg_theta,[],1)%/sqrt(length(AppEcc));
      shadedErrorBar(binD,nanmean(avg_theta),stdD,'b'); hold on
      stdD=nanstd(avg_theta,[],1)/sqrt(length(AppEcc));
      shadedErrorBar(binD,nanmean(avg_theta),stdD,'y'); hold on
      title 'approachError with vision'
      set(gca,'xdir','reverse','XLim',[0 25],'YLim',[-15 110])
%overlay plot of all the tracks and make them transparent

% col=[0,0,1,0.2] %4th entry makes line transparent
% 
% for i=1:length(AppEcc)
%     
%     plot(AppEcc{i}(:,1),AppEcc{i}(:,2),'color',col);
% end
% 
% hold on;
plot(binD,nanmean(avg_theta),'b','LineWidth',2)
    set(gca,'xdir','reverse','XLim',[0 20],'YLim',[-100 100])

    

%% dark trials
trials = find(group ==2 & lighting==0);
% 
figure

for tr =1:length(trials)
           
            subplot(4,4,tr);
            hold on
            %             for i = 1:8:length(targ{trials(tr)})-8
            %                 plot([ targ{trials(tr)}(i,1) targ{trials(tr)}(i+8,1)] ,[ targ{trials(tr)}(i,2) targ{trials(tr)}(i+8,2) ],'LineWidth',2,'Color',cmapVar(i,1,length(targ{trials(tr)}),jet ));
            %                 plot(head{trials(tr)}(i,1),head{trials(tr)}(i,2),'o','Color',cmapVar(i,1,length(targ{trials(tr)}),jet ) );
            %             end
          
            plot(targ{trials(tr)}(:,1),targ{trials(tr)}(:,2),'k');hold on
            plot(head{trials(tr)}(:,1),head{trials(tr)}(:,2),'m');
            axis([0 2000 0 1200]); axis off  
end
 title 'dark trials'
% 
% %theta and ranges over whole trial by session type
% Dark(:,1)=vertcat(thetaS{trials});
% Dark(:,2)=vertcat(rangeS{trials});
% Dark(:,3) = repmat(2,length(Dark(:,2)),1); % 2= group2
% 
% %speeds over whole trial by session type
% sDark(:,1)=vertcat(tSpeedS{trials});
% sDark(:,2)=vertcat(mSpeedS{trials});
% 
% xlswrite('allThetaRangeSpeedsDark', Dark,1,'A2');
% xlswrite('allThetaRangeSpeedsDark',vertcat(ID{trials}),1,'D2');
% xlswrite('allThetaRangeSpeedsDark',vertcat(Sx{trials}) ,1,'E2');
% 
% name = {'allTheta','allRange','Group','ID','Sex'};
% xlswrite('allThetaRangeSpeedsDark', sLight,2,'A2');
% name2={'targSpeed','mouseSpeed'};
% xlswrite('allThetaRangeSpeedsDark', name,1,'A1' );
% xlswrite('allThetaRangeSpeedsDark',name2,2,'A1');

data(2,1) = nanmean(freqApproach(trials));
err(2,1) = nanstd(freqApproach(trials))/sqrt(length(trials));
data(2,2) = nanmean(approachSuccess(trials));
err(2,2) = nanstd(approachSuccess(trials))/sqrt(length(trials));

ContactData(2,1)=nanmean(firstContact(trials)');
Contacterr(2,1) = nanstd(firstContact(trials)')/sqrt(length(trials));
ContactData(2,2)=nanmean(interContact(trials)');
Contacterr(2,2) = nanstd(interContact(trials)')/sqrt(length(trials));
ContactData(2,3)=nanmean(nContacts(trials)');
Contacterr(2,3) = nanstd(nContacts(trials)')/sqrt(length(trials));

%Prey eccentricity over whole trial
clear PreyEsession PreyE_mid Rsession TR TR_mid h ht stdD men med

% PreyEsession=thetahist(:,trials);
% figure
% plot(nanmean(PreyEsession,2),'k','LineWidth',3);
% title 'Prey eccentricity whole session with hearing only'
% 
% PreyE_mid_D=thetahist_mid(:,trials);
% figure
% plot(thetabins,nanmean(PreyE_mid_D,2),'k','LineWidth',3);
% title 'Prey eccentricity session with hearing only'

PreyE_mid_D=thetahist_midMove(:,trials);
figure
plot(thetabins,nanmean(PreyE_mid_D,2),'k','LineWidth',3);
title 'Prey eccentricity miving middle dist hearing only'

%range over whole session
Rsession_D=rhist(:,trials);
figure
plot(dbins,nanmean(Rsession_D,2),'k','LineWidth',3);
title 'range whole session with hearing only-Dark'

% theta x range between mice moving 
TR=tdhist_mv(:,:,trials);
figure
imagesc(squeeze(nanmean(TR,3)));
title ' bearing(theta) by range when mouse is moving hearing only-Dark '

%approach paths for Dark

clear goodApp AppEcc
goodApp=0;

for tr = 1:length(trials);
   
             for i = 1:length(AppR{trials(tr)});
                
                       clear x y maxX
                       x=AppR{trials(tr)}{:,i};
                       y=AppT{trials(tr)}{:,i};
                       preTouch=find(x<3);
                       if ~isempty (preTouch);
                       dist=preTouch(1);
                       else 
                           mindist=min(x);
                           dist=find(x==mindist);
                       end
                       %minX=min(x);dist=find(x==minX);
                       
                       goodApp=goodApp+1;
                       
                       AppEcc{goodApp}=[x(1:dist),abs(y(1:dist))];               
                  %[approachR{trials(tr)}(i,dist:end)',approachTheta{trials(tr)}(i,dist:end)']
            end 
end


figure; hold on
for i=1:length(AppEcc)
    if ~isempty(AppEcc{i})
    plot(AppEcc{i}(:,1),AppEcc{i}(:,2))
    end
end

clear thetaRangeHist avg_theta

for i=1:length(AppEcc);
     if ~isempty(AppEcc{i})
 [thetaRangeHist(:,:,i),avg_theta(i,:)]=myHist2Avg(AppEcc{i}(:,1),AppEcc{i}(:,2),binD,binE);
     end
end
 
%plot stad deviation plot
 figure; hold on
 
      stdD=nanstd(avg_theta,[],1)%/sqrt(length(AppEcc));
      shadedErrorBar(binD,nanmean(avg_theta),stdD,'k'); hold on
      stdD=nanstd(avg_theta,[],1)/sqrt(length(AppEcc));
      shadedErrorBar(binD,nanmean(avg_theta),stdD,'y'); hold on
      title 'approachError Dark'
      set(gca,'xdir','reverse','XLim',[0 25],'YLim',[-15 110]);

%overlay plot of all the tracks and make them transparent

% col=[1,0,1,0.2] %4th entry makes line transparent
% 
% for i=1:length(AppEcc)
%     
%     plot(AppEcc{i}(:,1),AppEcc{i}(:,2),'color',col);
% end
% 
% hold on;
% plot(binD,nanmean(avg_theta),'k','LineWidth',2)
% title 'bearing deviation with vision and hearing'
% set(gca,'xdir','reverse','XLim',[0 20],'YLim',[-100 100])  

    

%% deaf and light
trials = find(group ==5 & lighting==1);

figure

for tr =1:length(trials)
           
            subplot(4,4,tr);
            hold on
            %             for i = 1:8:length(targ{trials(tr)})-8
            %                 plot([ targ{trials(tr)}(i,1) targ{trials(tr)}(i+8,1)] ,[ targ{trials(tr)}(i,2) targ{trials(tr)}(i+8,2) ],'LineWidth',2,'Color',cmapVar(i,1,length(targ{trials(tr)}),jet ));
            %                 plot(head{trials(tr)}(i,1),head{trials(tr)}(i,2),'o','Color',cmapVar(i,1,length(targ{trials(tr)}),jet ) );
            %             end
          
            plot(targ{trials(tr)}(:,1),targ{trials(tr)}(:,2),'k');hold on
            plot(head{trials(tr)}(:,1),head{trials(tr)}(:,2),'g');
            axis([0 2000 0 1200]); axis off  
end
 title 'EP_light_trials'
% 
% %theta and ranges over whole trial by session type
% EP_L(:,1)=vertcat(thetaS{trials});
% EP_L(:,2)=vertcat(rangeS{trials});
% EP_L(:,3) = repmat(3,length(EP_L(:,2)),1); % 3= EP light group
% 
% %speeds over whole trial by session type
% sEP_L(:,1)=vertcat(tSpeedS{trials});
% sEP_L(:,2)=vertcat(mSpeedS{trials});
% 
% xlswrite('allThetaRangeSpeedsEPlight', EP_L,1,'A2');
% xlswrite('allThetaRangeSpeedsEPlight',vertcat(ID{trials}),1,'D2');
% xlswrite('allThetaRangeSpeedsEPlight',vertcat(Sx{trials}) ,1,'E2');
% 
% name = {'allTheta','allRange','Group','ID','Sex'};
% xlswrite('allThetaRangeSpeedsEPlight', sLight,2,'A2');
% name2={'targSpeed','mouseSpeed'};
% xlswrite('allThetaRangeSpeedsEPlight', name,1,'A1' );
% xlswrite('allThetaRangeSpeedsEPlight',name2,2,'A1');

data(3,1) = nanmean(freqApproach(trials));
err(3,1) = nanstd(freqApproach(trials))/sqrt(length(trials));
data(3,2) = nanmean(approachSuccess(trials));
err(3,2) = nanstd(approachSuccess(trials))/sqrt(length(trials));

ContactData(3,1)=nanmean(firstContact(trials)');
Contacterr(3,1) = nanstd(firstContact(trials)')/sqrt(length(trials));
ContactData(3,2)=nanmean(interContact(trials)');
Contacterr(3,2) = nanstd(interContact(trials)')/sqrt(length(trials));
ContactData(3,3)=nanmean(nContacts(trials)');
Contacterr(3,3) = nanstd(nContacts(trials)')/sqrt(length(trials));

clear PreyEsession PreyE_mid Rsession TR TR_mid h ht stdD men med

PreyEsession=thetahist(:,trials);
figure
%plot(nanmean(PreyEsession,3)); hold on
plot(nanmean(PreyEsession,2),'g','LineWidth',3);
title 'Prey eccentricity whole session with NO hearing-ear plug + light'

PreyE_mid_EPL=thetahist_mid(:,trials);
figure
plot(thetabins,nanmean(PreyE_mid_EPL,2),'g','LineWidth',3);
title 'Prey eccentricity session with hearing only'

PreyE_mid_EPL=thetahist_midMove(:,trials);
figure
plot(thetabins,nanmean(PreyE_mid_EPL,2),'g','LineWidth',3);
title 'Prey eccentricity EP light'

%range over whole session
Rsession_EPL=rhist(:,trials);
figure
plot(dbins,nanmean(Rsession_EPL,2),'g','LineWidth',3);hold on
%plot(nanmean(Rsession,3));hold on
title 'range whole session with with NO hearing-ear plug + light'

% theta x range between mice moving 
TR=tdhist_mv(:,:,trials);
figure
imagesc(squeeze(nanmean(TR,3)));
title ' bearing(theta) by range when mouse is moving, with NO hearing-ear plug + light '

%approach paths for ear plug Light

clear goodApp AppEcc
goodApp=0;

for tr = 1:length(trials);
   
            for i = 1:length(AppR{trials(tr)});
                
                       clear x y maxX
                       x=AppR{trials(tr)}{:,i};
                       y=AppT{trials(tr)}{:,i};
                       preTouch=find(x<3);
                       if ~isempty (preTouch);
                       dist=preTouch(1);
                       else 
                           mindist=min(x);
                           dist=find(x==mindist);
                       end
                       %minX=min(x);dist=find(x==minX);
                       
                       goodApp=goodApp+1;
                       
                       AppEcc{goodApp}=[x(1:dist),abs(y(1:dist))];               
                  %[approachR{trials(tr)}(i,dist:end)',approachTheta{trials(tr)}(i,dist:end)']
            end 
end


figure; hold on
for i=1:length(AppEcc)
    if ~isempty(AppEcc{i})
    plot(AppEcc{i}(:,1),AppEcc{i}(:,2))
    end
end


clear thetaRangeHist avg_theta

for i=1:length(AppEcc);
     if ~isempty(AppEcc{i})
 [thetaRangeHist(:,:,i),avg_theta(i,:)]=myHist2Avg(AppEcc{i}(:,1),AppEcc{i}(:,2),binD,binE);
     end
end
 
%plot stad deviation plot
 figure; hold on
 
      stdD=nanstd(avg_theta,[],1)%/sqrt(length(AppEcc));
      shadedErrorBar(binD,nanmean(avg_theta),stdD,'g'); hold on
      stdD=nanstd(avg_theta,[],1)/sqrt(length(AppEcc));
      shadedErrorBar(binD,nanmean(avg_theta),stdD,'y'); hold on
      title 'bearing deviation EP light'
      set(gca,'xdir','reverse','XLim',[0 25],'YLim',[-15 110])

%overlay plot of all the tracks and make them transparent

% col=[0,1,0,0.3] %4th entry makes line transparent
% 
% for i=1:length(AppEcc)
%     
%     plot(AppEcc{i}(:,1),AppEcc{i}(:,2),'color',col);
% end
% 
% hold on;
% plot(binD,nanmean(avg_theta),'g','LineWidth',2)
% title 'bearing deviation EP light'
% set(gca,'xdir','reverse','XLim',[0 20],'YLim',[-100 100])


%% deaf and dark
trials = find(group ==5 & lighting==0);

figure

for tr =1:length(trials)
           
            subplot(4,4,tr);
            hold on
            %             for i = 1:8:length(targ{trials(tr)})-8
            %                 plot([ targ{trials(tr)}(i,1) targ{trials(tr)}(i+8,1)] ,[ targ{trials(tr)}(i,2) targ{trials(tr)}(i+8,2) ],'LineWidth',2,'Color',cmapVar(i,1,length(targ{trials(tr)}),jet ));
            %                 plot(head{trials(tr)}(i,1),head{trials(tr)}(i,2),'o','Color',cmapVar(i,1,length(targ{trials(tr)}),jet ) );
            %             end
          
            plot(targ{trials(tr)}(:,1),targ{trials(tr)}(:,2),'k');hold on
            plot(head{trials(tr)}(:,1),head{trials(tr)}(:,2),'r');
            axis([0 2000 0 1200]); axis off  
end
 title 'EP_dark_trials'
% 
% %theta and ranges over whole trial by session type
% EP_D(:,1)=vertcat(thetaS{trials});
% EP_D(:,2)=vertcat(rangeS{trials});
% EP_D(:,3) = repmat(4,length(EP_D(:,2)),1); % 4= EP dark group
% 
% %speeds over whole trial by session type
% sEP_D(:,1)=vertcat(tSpeedS{trials});
% sEP_D(:,2)=vertcat(mSpeedS{trials});
% 
% xlswrite('allThetaRangeSpeedsEPdark', EP_D,1,'A2');
% xlswrite('allThetaRangeSpeedsEPdark',vertcat(ID{trials}),1,'D2');
% xlswrite('allThetaRangeSpeedsEPdark',vertcat(Sx{trials}) ,1,'E2');
% 
% name = {'allTheta','allRange','Group','ID','Sex'};
% xlswrite('allThetaRangeSpeedsEPdark', sLight,2,'A2');
% name2={'targSpeed','mouseSpeed'};
% xlswrite('allThetaRangeSpeedsEPdark', name,1,'A1' );
% xlswrite('allThetaRangeSpeedsEPdark',name2,2,'A1');

data(4,1) = nanmean(freqApproach(trials));
err(4,1) = nanstd(freqApproach(trials))/sqrt(length(trials));
data(4,2) = nanmean(approachSuccess(trials));
err(4,2) = nanstd(approachSuccess(trials))/sqrt(length(trials));

ContactData(4,1)=nanmean(firstContact(trials)');
Contacterr(4,1) = nanstd(firstContact(trials)')/sqrt(length(trials));
ContactData(4,2)=nanmean(interContact(trials)');
Contacterr(4,2) = nanstd(interContact(trials)')/sqrt(length(trials));
ContactData(4,3)=nanmean(nContacts(trials)');
Contacterr(4,3) = nanstd(nContacts(trials)')/sqrt(length(trials));

clear PreyEsession PreyE_mid Rsession TR TR_mid h ht stdD men med

PreyEsession=thetahist(:,trials);
figure
plot(nanmean(PreyEsession,2),'r','LineWidth',3);
title 'Prey eccentricity whole session with NO hearing NO vision- ear plug & dark'

PreyE_mid_EPD=thetahist_mid(:,trials);
figure
plot(thetabins,nanmean(PreyE_mid_EPD,2),'r','LineWidth',3);
title 'Prey eccentricity ear plugs no light'

PreyE_mid_EPD=thetahist_midMove(:,trials);
figure
plot(thetabins,nanmean(PreyE_mid_EPD,2),'r','LineWidth',3);
title 'Prey eccentricity EP dark'

%range over whole session
Rsession=rhist(:,trials);
figure
plot(dbin,nanmean(Rsession,2),'r','LineWidth',3);
title 'range whole session with NO hearing NO vision- ear plug & dark'

% theta x range between mice moving 
TR=tdhist_mv(:,:,trials);
figure
imagesc(squeeze(nanmean(TR,3)));
title ' bearing(theta) by range when mouse is moving, with NO hearing NO vision- ear plug & dark'

%approach paths for ear plug Light

clear goodApp AppEcc
goodApp=0;

for tr = 1:length(trials);
   
            for i = 1:length(AppR{trials(tr)});
                
                       clear x y maxX
                       x=AppR{trials(tr)}{:,i};
                       y=AppT{trials(tr)}{:,i};
                       preTouch=find(x<3);
                       if ~isempty (preTouch);
                       dist=preTouch(1);
                       else 
                           mindist=min(x);
                           dist=find(x==mindist);
                       end
                       %minX=min(x);dist=find(x==minX);
                       
                       goodApp=goodApp+1;
                       
                       AppEcc{goodApp}=[x(1:dist),abs(y(1:dist))];               
                  %[approachR{trials(tr)}(i,dist:end)',approachTheta{trials(tr)}(i,dist:end)']
            end 
end


figure; hold on
for i=1:length(AppEcc)
    if ~isempty(AppEcc{i})
    plot(AppEcc{i}(:,1),AppEcc{i}(:,2))
    end
end
 
%plot stad deviation plot
 clear thetaRangeHist avg_theta

for i=1:length(AppEcc);
     if ~isempty(AppEcc{i})
 [thetaRangeHist(:,:,i),avg_theta(i,:)]=myHist2Avg(AppEcc{i}(:,1),AppEcc{i}(:,2),binD,binE);
     end
end
 
%plot stad deviation plot
 figure; hold on
 
      stdD=nanstd(avg_theta,[],1)%/sqrt(length(AppEcc));
      shadedErrorBar(binD,nanmean(avg_theta),stdD,'r'); hold on
      stdD=nanstd(avg_theta,[],1)/sqrt(length(AppEcc));
      shadedErrorBar(binD,nanmean(avg_theta),stdD,'y'); hold on
      title 'bearing deviation EP light'
      set(gca,'xdir','reverse','XLim',[0 25],'YLim',[-15 110])

%% outcome bar graphs figure generation
figure
bar(data(:,1)); hold on; errorbar(1:4,data(:,1),err(:,1),'o')
ylabel('approach frequency')
set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'light','dark','deaf','dark + deaf'});


figure
bar(data(:,2)); hold on; errorbar(1:4,data(:,2),err(:,2),'o')
ylabel('approach success')
set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'light','dark','deaf','dark + deaf'});
title('figure 3D bottom')

figure
bar(ContactData(:,1)); hold on; errorbar(1:4,ContactData(:,1),Contacterr(:,1),'o')
ylabel('FirstContactTime')
set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'light','dark','deaf','dark + deaf'});
title('figure 3E top')

figure
bar(ContactData(:,2)); hold on; errorbar(1:4,ContactData(:,2),Contacterr(:,2),'o')
ylabel('InterContactInterval')
set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'light','dark','deaf','dark + deaf'});
title('figure 3E bottom')

% figure
% bar(ContactData(:,3)); hold on; errorbar(1:4,ContactData(:,3),Contacterr(:,3),'o')
% ylabel('Num of contacts')
% set(gca,'Xtick',1:4); set(gca,'XtickLabel',{'light','dark','deaf','dark + deaf'});
% title('figure 3F top')
%plot all condition theta

figure
plot(thetabins,nanmean(PreyE_mid_D,2),'k','LineWidth',3); hold on
plot(thetabins,nanmean(PreyE_mid_L,2),'b','LineWidth',3);
plot(thetabins,nanmean(PreyE_mid_EPL,2),'g','LineWidth',3);
plot(thetabins,nanmean(PreyE_mid_EPD,2),'r','LineWidth',3);

title 'Prey eccentricity all conditions >5cm mouse moving'

% PLOT ALL CONDITIONS range over trial
figure
plot(nanmean(Rsession_D,2),'k','LineWidth',3);hold on
plot(nanmean(Rsession_L,2),'b','LineWidth',3);hold on
plot(nanmean(Rsession_EPL,2),'g','LineWidth',3);hold on
plot(nanmean(Rsession,2),'r','LineWidth',3);

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
            
            subplot(4,4,tr);
            hold on
            %             for i = 1:8:length(targ{trials(tr)})-8
            %                 plot([ targ{trials(tr)}(i,1) targ{trials(tr)}(i+8,1)] ,[ targ{trials(tr)}(i,2) targ{trials(tr)}(i+8,2) ],'LineWidth',2,'Color',cmapVar(i,1,length(targ{trials(tr)}),jet ));
            %                 plot(head{trials(tr)}(i,1),head{trials(tr)}(i,2),'o','Color',cmapVar(i,1,length(targ{trials(tr)}),jet ) );
            %             end
            
            plot(targ{trials(tr)}(:,1),targ{trials(tr)}(:,2),'g');
            plot(head{trials(tr)}(:,1),head{trials(tr)}(:,2),'b');
            axis([0 2000 0 1200]); axis off
            
        end
        
        timeRange = 120;
        hit_rg=[]; hit_th=[]; miss_rg=[]; miss_th=[];
        figure;  set(gcf,'Name',sprintf('capture %s %s',grouplabels{m+1},lightlabels{light}));
        for tr = 1:length(trials)
            rg = basin{trials(tr)}(:,1); thet = basin{trials(tr)}(:,2); nxt = basin{trials(tr)}(:,3);
            
            hold on
            plot(rg(nxt>timeRange),abs(thet(nxt>timeRange)),'r.')
            plot(rg(nxt<timeRange),abs(thet(nxt<timeRange)),'g.')
            hit_rg  = [hit_rg ; rg(nxt<timeRange)]; hit_th = [hit_th ; abs(thet(nxt<timeRange))];
            miss_rg  = [miss_rg ; rg(nxt>timeRange)]; miss_th = [miss_th ; abs(thet(nxt>timeRange))];
        end
        
        hits = myHist2(hit_rg,hit_th,2.5:5:50,15:30:180);
        miss =  myHist2(miss_rg,miss_th,2.5:5:50,15:30:180);
        %         figure
        %         imagesc(hits)
        %         figure
        %         imagesc(miss)
        %
        figure
        imagesc(hits./(hits+miss),[0 1]);
        axis xy
        
        
        approachhist=0; inds =0;
        figure
        hold on
        col = 'rgbcmk'; set(gcf,'Name',sprintf('capture %s %s',grouplabels{m+1},lightlabels{light}));
        ntr=0;
        for tr= 1:length(trials)
            for j = 1:size(approachR{trials(tr)},1)
                
                for i = 1:length(approachTheta{trials(tr)}(j,:));
                    straight(i) = all(abs(approachTheta{trials(tr)}(j,i:end))<60);
                end
                
                [y ind] = max(approachR{trials(tr)}(j,:));
                if y>5
                    ntr=ntr+1;
                    plot(approachR{trials(tr)}(j,ind), approachTheta{trials(tr)}(j,ind),'bo');
                    for k = ind:59
                        plot([approachR{trials(tr)}(j,k) approachR{trials(tr)}(j,k+1)], [approachTheta{trials(tr)}(j,k) approachTheta{trials(tr)}(j,k+1)] ,'Color',cmapVar(k,1,59,jet));
                        %plot(approachSpeed{trials(tr)}(i,:));
                    end
                    approachD = find( approachR{trials(tr)}(j,ind:end)>5 & approachR{trials(tr)}(j,ind:end)<10)+ind-1;
                    bins = 10:20:180;
                    approachhist(ntr,1:length(bins)) = hist(abs(approachTheta{trials(tr)}(j,approachD)),bins)/length(approachD);
                end
                
            end
        end
        axis([ 2.5 30 -180 180])
        
        
        bins = 10:20:180;

        
        %%% get distribution of starting points
        inds =0;
        ntr=0;
        for tr= 1:length(trials)
            for j = 1:size(approachR{trials(tr)},1)
                for i = 1:length(approachTheta{trials(tr)}(j,:));
                    straight(i) = all(abs(approachTheta{trials(tr)}(j,i:end))<60);
                end
                [y ind] = max(approachR{trials(tr)}(j,:).*straight); %%% find most distant point where always aimed at target (abs(theta)<90)
                if y>5
                    ntr=ntr+1;
                    inds(ntr)=y;
                end
            end
        end
        
        
        
        indAll(m,light,1:length(dbins)) = hist(inds,dbins)/length(inds);
        
        figure
        figure; set(gcf,'Name',sprintf('approach %s %s',grouplabels{m+1},lightlabels{light}));
        plot(dbins,squeeze(indAll(m,light,:))); xlim([ 0 30]); xlabel('distance of approach initiation (cm)');
        
        
        
        approachhist = nanmean(approachhist,1);
        
        figure; set(gcf,'Name',sprintf('approach %s %s',grouplabels{m+1},lightlabels{light}));
        plot(bins,approachhist); ylim([0 1]); xlim([0 90])
        approachHistAll(m,light,:) = approachhist; xlabel('theta 5- 10 cmn')
        %
        %         escapehist=0;
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
        tAll(m,light) = mean(captureT(trials))
        midthetaAll(m,light,:) = nanmean(midtheta(:,trials),2);
        midtheta_mvAll(m,light,:) = nanmean(midtheta_mv(:,trials),2);
        speedAll(m,light,:) = nanmean(speedhist(:,trials),2);
        distAll(m,light,:) = nanmean(rhist(:,trials),2);
        tdhistAll(m,light,:,:) = nanmean(tdhist(:,:,trials),3);
        tdhist_mvAll(m,light,:,:) = nanmean(tdhist_mv(:,:,trials),3);
        spdistAll(m,light,:,:) = nanmean(spdist(:,:,trials),3);
        stateAll(m,light,:) = nanmean(state(:,trials),2);
        
        td = squeeze(tdhist_mvAll(m,light,:,:));
        for d = 1:size(td,2);
            tdcurve = td(:,d); tdcurve(tdcurve<max(tdcurve)*0.1)=0;
            meantheta(m,light,d) = sum(tdcurve.*thetabins')/sum(tdcurve);
            stdtheta(m,light,d) = sqrt(sum(tdcurve.*(thetabins.^2)')/sum(tdcurve));
        end
        
        
        figure(distFig); subplot(ncond,2,2*(m-1)+(3-light)); plot(dbins,squeeze(distAll(m,light,:))); ylim([0 0.5])
        figure(thetaFig); subplot(ncond,2,2*(m-1)+(3-light)); plot(thetabins,squeeze(midthetaAll(m,light,:))); hold on; plot(thetabins,squeeze(midtheta_mvAll(m,light,:)),'g'); ylim([0 0.5]); xlim([-180 180])
        figure(speedFig); subplot(ncond,2,2*(m-1)+(3-light)); plot(spdbins,squeeze(speedAll(m,light,:))); ylim([0 0.5])
        figure(dthistFig); subplot(ncond,2,2*(m-1)+(3-light)); imagesc(squeeze(tdhistAll(m,light,:,:)),[0 0.05]);
        figure(dt_curveFig); subplot(ncond,2,2*(m-1)+(3-light));  plot(squeeze(stdtheta(m,light,:))); ylim([0 120]); xlim([0 20])
        figure(dthist_mvFig); subplot(ncond,2,2*(m-1)+(3-light));   imagesc(squeeze(tdhist_mvAll(m,light,:,:)),[0 0.05]);
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



