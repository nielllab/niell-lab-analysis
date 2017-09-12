function [latency, tracks, Time1stApp IAppI mouseTouch, cricketEdge,mouseTouchEd,MouseTouchL, Rstarts, hits, range, tdhist, cricketTouch] = analyzePlexi(fname,fps,scale, virtual)

% [f p] = uigetfile('*.*','behav data');
% fname = fullfile(p,f);
% scale = 35;
% fps = 60

in=importdata(fname);
data = in.data; %%% numeric portion
data = data/scale;

box = data(1,:);
wallL=data(2,1:2);
wallR=data(2,3:4);

figure
plot(box([1 3 7 5 1]),box([2 4 8 6 2])); hold on
plot(wallL(1),wallL(2),'ko'); hold on
plot(wallR(1),wallR(2),'ko');
hold on

contact = ~isnan(data(4:end,4));
contact(1) = 0; contact(end)=0;
starts = find(diff(contact)>0)+1;
Time1stApp=starts(1)/fps;
ends = find(diff(contact)<0);
n = length(starts);

%calculate inter-approach intervals

if n>1
    for i=1:(n-1)
        gap(i)=(starts(i+1)-ends(i))/fps;
    end
    IAppI=nanmean(gap);
    
else
    IAppI=NaN;
end

%generate approach paths in terms of lateral error
clear ear1 ear2 body targ head
nApproach=0;  %%% approaches that satisfy constraints. start from > 5 cm away, no rearing

% if Tnum==1
for i = 1:n
    ear1{i}(:,1:2) = data(starts(i):ends(i),1:2);
    ear2{i}(:,1:2) = data(starts(i):ends(i),3:4);
    body{i}(:,1:2) = data(starts(i):ends(i),5:6);
    targ{i}(:,1:2) = data(starts(i):ends(i),7:8);
    
    head{i} = 0.5*(ear1{i} + ear2{i});
         
    %eliminate NaN from tracks, data cleaning
    if isnan(head{i}(end,1));
      num=find(~isnan(head{i}(:,1))); 
      last=num(end);
      
      head{i}=head{i}(1:last,:);
      targ{i}=targ{i}(1:last,:);
      ear1{i}=ear1{i}(1:last,:);
      ear2{i}=ear2{i}(1:last,:);
      
    end
     
%     pt1 = ear2{i};%left ear coordinates
%     pt2 = ear1{i};%right ear
%     headvec = pt2-pt1;
%     headtheta = getSmoothAngle(headvec)-pi/2
%    
%     targvec = targ{i}-head{i};
%     targtheta = getSmoothAngle(targvec)
%     
%     Rtheta=headtheta-targtheta;
%     thetaR{i}=mod(Rtheta*180/pi,360);
%     
%     th{i} = mod(thetaR{i}+180,360)-180
    
    %r{i} = sqrt((head{i}(:,1) - targ{i}(:,1)).^2 + (head{i}(:,2) - targ{i}(:,2)).^2);

    plot(head{i}(:,1),head{i}(:,2),'g')
    axis equal
    plot(targ{i}(:,1),targ{i}(:,2),'ro')
   
    % where is the target?
        if mean(targ{i}(:,1))>20
        targleft=0;
        else targleft=1;
        end
    % where does the mouse touch
    if virtual==0;
        if head{i}(end,1) < wallL(1) +10;
            MouseTouchLeft=1;
        elseif head{i}(end,1) > wallR(1)-35;% touch point of track is within 10cm of right wall
            MouseTouchLeft=0;
        end
    end
    if virtual==1
    MouseTouchLeft=0;
    targleft=0;
    end
        
    % does the mouse go to the correct side?
        if targleft && MouseTouchLeft; %mouse touched left wall
            hit=1;
        elseif ~targleft && ~MouseTouchLeft;
            hit=1;
        else
            hit=0;
        end
%         if targleft
%         wall(i) = min(head{i}(:,1));
%         else
%         wall(i) = max(head{i}(:,1));
%         end
    hits(i)=hit;

  if hit==1
      
       %track relative to cricket x coordinates
        if targleft && MouseTouchLeft
            d{i}= head{i}(:,1)-wallL(1);
            bottom = box(4);
            top=box(2);
            track{i}(:,1) = head{i}(:,1) - targ{i}(:,1);
            track{i}(:,2) = -(head{i}(:,2)-targ{i}(:,2));
            
        elseif targleft==0 && MouseTouchLeft==0;
            d{i}= wallR(1)-head{i}(:,1);
            track{i}(:,1) = targ{i}(:,1) - head{i}(:,1) ;
            bottom = box(8);
            top=box(6);
            track{i}(:,2) = head{i}(:,2)-targ{i}(:,2);
      
        end
    
    
    % track y coordinates
    
    [dist startApproach] = max(d{i});
    
    pathA=d{i}(startApproach:end)
    closest=min(pathA)
    firstContact = (find(pathA==closest))+ (startApproach-1); %length(d{i})
    
  else
        % save out error trial data, distance to run to wall etc. 
        track{i}=NaN; % can't make a distance track if they touch the wrong wall
        
        if targleft==0 && MouseTouchLeft
        d{i}=head{i}(:,1)-wallL(1)
        [dist startApproach] = max(d{i});
        bottom = box(4);
        top=box(2);
        
        pathA=d{i}(startApproach:end);
        closest=min(pathA);
        firstContact = (find(pathA==closest))+ (startApproach-1);%length(d{i})
        
        elseif  targleft && MouseTouchLeft==0
        d{i}=wallR(1)-head{i}(:,1);
        [dist startApproach] = max(d{i});
        bottom = box(8);
        top=box(6);
        
        pathA=d{i}(startApproach:end)
        closest=min(pathA);
        firstContact = (find(pathA==closest))+ (startApproach-1); %length(d{i})
        end
        
 end
    plot(head{i}(firstContact,1),head{i}(firstContact,2),'go')
    plot(head{i}(startApproach:firstContact,1),head{i}(startApproach:firstContact,2),'g', 'Linewidth',2);
    plot(targ{i}(firstContact,1),targ{i}(firstContact,2),'ro')
    plot(targ{i}(startApproach:firstContact,1),targ{i}(startApproach:firstContact,2),'r', 'Linewidth',2);
    
%     figure
     bind=3:3:30;
     bint=-30:2:30;
%      hrt=hist2(d{i},track{i}(:,2),bind,bint);
%      hrt= hrt/sum(hrt(:));
%     imagesc(flipud(hrt));
%     title 'dist to targ by diff-Y to Target'
    
  % exclude wall touch where cricket was in corner
    cricketCorner= targ{i}(firstContact,2) < bottom+3 | targ{i}(firstContact,2)> top-3;

    if dist>3 & hit==1 & cricketCorner==0
   
    nApproach = nApproach+1;
    %CK_speed(nApproach)= calc targ speed 6 frames (100ms) before
    %StartApproach and 6 frames (100 ms after startApproach)
    latency(nApproach) = (firstContact +starts(i))/fps;
    cricketTouch(nApproach) = targ{i}(firstContact,2);
    mouseTouch(nApproach) = head{i}(firstContact,2);
    cricketEdge(nApproach)=cricketTouch(nApproach) < bottom+2.5 | cricketTouch(nApproach)> top-2.5;
    mouseTouchEd(nApproach)=mouseTouch(nApproach) < bottom+2.5 | mouseTouch(nApproach)> top-2.5;
    MouseTouchL(nApproach)=MouseTouchLeft;
    cricketTouch(nApproach) = targ{i}(firstContact,2)- bottom;
    mouseTouch(nApproach) = head{i}(firstContact,2)-bottom;
    

    Rstarts(nApproach)=dist;
    
   
    tracks{nApproach} = track{i}(startApproach:firstContact,:);
    range{nApproach} = d{i}(startApproach:firstContact);
    
    tdhist(:,:,nApproach) = myHist2(d{i}(d{i}>.5),track{i}(:,2),bind,bint)/length(d{i});
        
end

title(fname);

% end
    if dist>3 & hit==1 & cricketCorner==0

figure
hold on
for i = 1:nApproach
    plot(range{i});
end

figure
hold on
for i = 1:nApproach
    plot(tracks{i}(:,1),tracks{i}(:,2),'.');
    axis([0 50 -30 30])
end
title(fname);
    end 
close all
end




