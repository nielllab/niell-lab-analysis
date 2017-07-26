function [latency, tracks, mouseTouch, mouseTouchEd, Rstarts, range, tdhist,cricketTouch] = analyzePlexi(fname,fps,scale,Tnum)


%import data from excel file 

in=importdata(fname);
data = in.data; %%% numeric portion
data = data/scale;

box = data(1,:);
wallL=data(2,1:2);
wallR=data(2,3:4)

figure
plot(box([1 3 7 5 1]),box([2 4 8 6 2])); hold on
plot(wallL(1),wallL(2),'ko'); hold on
plot(wallR(1),wallR(2),'ko');
hold on

% define contact times
contact = ~isnan(data(:,8));
contact(1) = 0; contact(end)=0;
starts = find(diff(contact)>0)+1; %index of start of approaches
ends = find(diff(contact)<0); % index end of approach (contact!!)
n = length(starts);

clear ear1 ear2 body targ head
nApproach=0;  %%% approaches that satisfy constraints. start from > 5 cm away, no rearing

% analysis of each approach for each file
for i = 1:n
    ear1{i}(:,1:2) = data(starts(i):ends(i),1:2);
    ear2{i}(:,1:2) = data(starts(i):ends(i),3:4);
    body{i}(:,1:2) = data(starts(i):ends(i),5:6);
    targ{i}(:,1:2) = data(starts(i):ends(i),7:8);
    
    if Tnum==2
    targ2{i}(:,1:2) = data(starts(i):ends(i),9:10);
    end
    
    head{i} = 0.5*(ear1{i} + ear2{i});
    plot(head{i}(:,1),head{i}(:,2),'g')
    axis equal
    plot(targ{i}(:,1),targ{i}(:,2),'r')
    
    if Tnum==2
    plot(targ2{i}(:,1),targ{i}(:,2),'m')
    end
    
    if isnan(head{i}(end,1));
      num=find(~isnan(head{i}(:,1))); 
      last=num(end);
      head{i}=head{i}(1:last,:);
      targ{i}=targ{i}(1:last,:);
      if Tnum==2
         targ2{i}=targ2{i}(1:last,:); 
      end
    end
    
    if head{i}(end,1) < (wallL(1)+5);
    Lapp=1
    else
    Lapp=0
    end

       
    if Tnum==1;
        
        if mean(targ{i}(:,1))>20
        targleft=0;
        else targleft=1;
        
        end
           
    
    %track App coordinates relative to target
    
        if Lapp & targleft %target is on left and mouse went left
            d{i}= head{i}(:,1)-wallL(1);
            bottom = box(4);
            top=box(2)
            track{i}(:,1) = head{i}(:,1) - targ{i}(:,1);
            track{i}(:,2) = head{i}(:,2)- targ{i}(:,2);
            
        elseif Lapp & targleft==0 %target is on the right and mouse went left
            d{i}= wallR(1)-head{i}(:,1);
            bottom = box(8);
            top=box(6);
            %track{i}(:,1) = targ{i}(:,1) - head{i}(:,1);
            
            
        elseif Lapp==0 & targleft %target is on left and mouse went right
            d{i}= wallR(1)-head{i}(:,1);
            bottom = box(8);
            top=box(6);        
            
        elseif Lapp==0 & targleft==0 %target is on right and mouse went right
            d{i}= wallR(1)-head{i}(:,1);
            bottom = box(8);
            top=box(6);
            track{i}(:,1) = targ{i}(:,1) - head{i}(:,1);
            track{i}(:,2) = targ{i}(:,2)- head{i}(:,2);
           
            
        end
            
%             x=find(diff(d{i}) < -0.1)           
%             dist  = d{i}(x(1+1)); 
%             startApproach= x(1+1);%based on slope of change in range towards correct wall
%             firstContact = min(find(d{i}<1)); %index first contact
            
                [dist startApproach] = max(d{i});
                firstContact = min(find(d{i}<1.5))      
    end
    
    if Tnum==2
        targleft=2;
        if Lapp==1 %going towards left hand target
             d{i}= head{i}(:,1)- wallL(1); %distance from target L
             track{i}(:,1) = head{i}(:,1) - targ{i}(:,1);
             track{i}(:,2) = targ{i}(:,2)- head{i}(:,2);
             bottom = box(4);
             top=box(2);
             
        elseif Lapp==0;% going towards the right target
            d{i}= wallR(1)-head{i}(:,1);
            track{i}(:,1) = targ2{i}(:,1) - head{i}(:,1);
            track{i}(:,2) = targ2{i}(:,2)- head{i}(:,2);
            bottom = box(8);
            top=box(6);
        end
         [dist startApproach] = max(d{i});
         firstContact = min(find(d{i}<1.5))
        
    end
    
    plot(head{i}(firstContact,1),head{i}(firstContact,2),'go')
    plot(head{i}(startApproach:firstContact,1),head{i}(startApproach:firstContact,2),'g', 'Linewidth',2)
    plot(targ{i}(firstContact,1),targ{i}(firstContact,2),'ro')
    plot(targ{i}(startApproach:firstContact,1),targ{i}(startApproach:firstContact,2),'r', 'Linewidth',2)
      
     bind=3:3:30;
     bint=-30:2:30;
    
    if dist>4 %& targ{i}(firstContact,2)> bottom+2 & targ{i}(firstContact,2)< top+2 % starting distance is greater than 5cm away and cricket is not in corner
    
    nApproach = nApproach+1;
    %CK_speed(nApproach)= calc targ speed 6 frames (100ms) before
    %StartApproach and 6 frames (100 ms after startApproach)
    latency(nApproach) = (firstContact +starts(i))/fps;
    cricketTouch(nApproach) = targ{i}(firstContact,2)-bottom;
    mouseTouch(nApproach) = head{i}(firstContact,2)-bottom;
    mouseTouchEd(nApproach)=mouseTouch(nApproach) < bottom+2 | mouseTouch(1)> top-2;
    sideApp(nApproach)=Lapp;
    sideTarg(nApproach)=targleft;
    
    tracks{nApproach} = track{i}(startApproach:firstContact,:);
    range{nApproach} = d{i}(startApproach:firstContact);
    Rstarts{nApproach}=dist;
    
    tdhist(:,:,nApproach) = myHist2(d{i}(d{i}>1),track{i}(:,2),bind,bint)/length(d{i});
  
    end
    
end
title(fname);

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




