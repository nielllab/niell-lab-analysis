function appending_drift_depth(afile, L4_tip1, L4_tip2,angle,L5_end)


if ~exist ('afile','var');
[afname, apname] = uigetfile('*.mat','analysis data'); %load data set that you want to append layer data to
afile = fullfile(apname,afname);
end
load(afile)

%histology variables
if ~exist ('tip_loc_1','var')
L4_tip1 = input('Channel/site nearest the lower border of layer 4 Anterior shank: '); 
L4_tip2 = input('Channel/site nearest the lower border of layer 4 Posterior shank: ');
angle = input('angle on electrode penetration : ');
end

if ~exist ('L5_end','var')
L5_end = input('Does electrode end in layer 5 or deeper? 1=yes, 0=no ');
end

L4_end  = [L4_tip1 L4_tip2];

peakChan = peakchan'; 


% 375 lower depth of layer4
D = zeros(size(peakChan));
D(peakChan<33) = 500-((L4_end(1)-peakChan(peakChan<33))*25*sin(angle*pi/180));
D(peakChan>32) = 500-((L4_end(2)-peakChan(peakChan>32))*25*sin(angle*pi/180));
%350/375 is around the top of layer4

if L5_end==0;
    
L4_offset(1) = input('Num sites superficial to L4 bottom anterior shank '); %apply standard depth in uM for location of tip in top, middle or bottom third of any layer
L4_offset(2) = input('Num sites superficial to L4 bottom posterior shank  ');
  
D(peakChan<33) = 500-(((L4_end(1)+ L4_offset(1))-peakChan(peakChan<33))*25*sin(angle*pi/180));
D(peakChan>32) = 500-(((L4_end(2)+ L4_offset(2))-peakChan(peakChan>32))*25*sin(angle*pi/180));

end

depth=D;    
     
save(afile, 'depth','-append')
end
