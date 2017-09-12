function [range_smooth contact move_away hits approachStart pre_approach_start escape time_recovery contactDwell_time numEscape_trial]=rangeAnalysis(range,TH,range_end,fps)

%% define "contact" events and determine peri-approach moments
 %range is the list of distances at each frame of the movie concatenated by
 %group of interest
 
 %TH is the distance threshold that defines close contact between mouse and
 %target
 
 %range_end is the frame number the defines the end of a given trial within
 %a group, when the cricket is caught for good
 
 %fps is the frame rate of the movie so we can convert frames to a time
 %metric
sz=5; % smoothing factor
range_smooth=conv(range,ones(1,sz),'same')/sz;

isTouch = range_smooth<=TH; 
contact   = find(diff(isTouch)==1) + 1;%generate indices where tranistion from NOT touching (0) to touching (1) occurs
move_away  = find(diff(isTouch)==-1) + 1; %generate indices where tranistion from touching (1) to NOT touching (0) occurs
 
  
%   figure 
%   plot(range_smooth); hold on;
%   plot(contact,0,'g*'); hold on;
%   plot(move_away,0,'rs'); hold on;
%   plot(range_end, 45,'ko'); hold on;
%   legend ('range','contacts','escapes','trial end')
%   title 'range with approaches and escapes'

%determine the frames that define the pre-approach period
clear n
n=1;hit=[];pre_approach_start=[];approachStart=[]

for i=1:length(contact);
    clear cont peri_approach pre_dist change
    cont=contact(i);
    pre_contact_ranges=range_smooth(n:cont);
    TH_dist=pre_contact_ranges>7.5; %10 cm is threshold for having been at a "distance" prior to contact
    hit(i,:)=sum(TH_dist);% # of frames mouse was more than TH distance away from cricket
    
    if hit(i,:)>0;   
    approaching=diff(pre_contact_ranges)<0;
    pre=find(diff(approaching)==1)-12; %include 200 ms prior to when mouse is consistently decreasing it's range relative to cricket
    pre_approach_start(i,:)=pre(end)+n;
    approachStart(i,:)=pre(end)+12+n;

    else
    pre_approach_start(i,:)=NaN;
    approachStart(i,:)=NaN
   
    end
    
    edge_check=range_end(range_end>pre_approach_start(i)& range_end<cont)
    if  edge_check>0
        pre_approach_start(i,:)=edge_check + 1;
        approachStart(i,:)=edge_check+1;
    end    
     n=cont+1;
end

 %% histogram of times between events:(1) contact to contact=intercontact_dist (2)escape to next contact=recovery (3) contact dwell time =likelihood of escape after contact 
%find frames/indices between contact(hit) and escape events
clear hits n i trial_hit secondary_contact recovery trial_escapes intercontact_dist escapeEnd contactDwell
hits=contact(hit>0);n=0;
trial_hit={};secondary_hit={};recovery={};trial_escapes={};contactEnd={};contactDwell={};intercontact_dist={}

for i=1:length(range_end)% size of range_end equals number of trials
   
    trial_hit{i}=hits(hits>=n & hits<range_end(i));
    
    if range_smooth(n+1)<5
        trial_hit{i}(2:end+1,1)=trial_hit{i};
        trial_hit{i}(1)=n+1;
    end
%     intercontact_dist{i}=diff(trial_hit{i});%dist in terms of frame number between bonifide contacts 
    clear idx_esc away_trials away_escapes j 
    
 if  length(trial_hit{i})>1;%indicates that there was an escape and recovery before final capture in a trial
    secondary_hit{i}=trial_hit{i}(2:end);%contacts that follow the first contact which is then followed by an escape
    away_trials=move_away(move_away>trial_hit{i}(1) & move_away<trial_hit{i}(end));
  
        for j=1:length(trial_hit{i})-1
        idx_esc=away_trials(away_trials>trial_hit{i}(j)& away_trials<trial_hit{i}(j+1));
        away_escapes(j,1)=idx_esc(end);
        end 
     
        trial_escapes{i}=away_escapes;
    
        contactEnd{i}=trial_escapes{i};
        contactEnd{i}(end+1,1)=range_end(i);%make sure it is a row vector
    
        contactDwell{i}=contactEnd{i}-trial_hit{i};
    
 else
        trial_escapes{i}=NaN;% there were no escapes in the trial
        secondary_hit{i}=NaN;
    
        if ~isempty(trial_hit{i});
        contactDwell{i}=range_end(i)-trial_hit{i}(1);% get dwell time mouse was in close contact with the target from first contact
        else
        contactDwell{i}=range_end(i)-(n+1);%spent the whole trial in contact with target
        end
 end
    
    if ~isnan (trial_escapes{i});
        recovery{i}=secondary_hit{i}-trial_escapes{i};
    else
        recovery{i}=NaN;
    end
   
    n=range_end(i)+1;
end
Fx = @(recovery) any(isnan(recovery))
idx_nan=cellfun(Fx,recovery);

escape=vertcat(trial_escapes{idx_nan==0});
time_recovery=vertcat(recovery{idx_nan==0})/fps;
contactDwell_time=vertcat(contactDwell{idx_nan==0})/fps;

%intercontactTime_L=vertcat(intercontact_dist{idx_nan_L==0})/fps;
  figure
  plot(range_smooth); hold on;
  plot(hits,0,'g*'); hold on;
  plot(pre_approach_start(~isnan(pre_approach_start)),0,'k*');
  plot(range_end, 40,'ko'); hold on;
%   plot(away_L,2.5,'rs');
  plot(escape,2.5,'rs');
  
Fl=@(trial_escapes) length(trial_escapes);
numEscape_trial=cellfun(Fl,trial_escapes);
  
  