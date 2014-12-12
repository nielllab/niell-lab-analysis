%load analysis file using input function to manually load relevant files

% tip_loc = input('depth in uM of electrode tip : ');
% angle = input('angle on electrode penetration : ');

clear all
close all

[afname, apname] = uigetfile('*.mat','analysis data'); %load data set that you want to append layer data to
noisepname = apname;
afile = fullfile(apname,afname);
load(afile)

%histology variables
tip_loc_1 = input('depth in uM of electrode tip first shank : '); %apply standard depth in uM for location of tip in top, middle or bottom third of any layer
tip_loc_2 = input('depth in uM of electrode tip second shank : ');
angle = input('angle on electrode penetration : ');

 %extracts first site of all tetrodes
 
 tet = cells(:,1);
    
 %converts tetrode location into phsyical distance from tip of electrode
    
    Distance = zeros(length(cells),1);
   
    Distance(tet==1)= 775;
    Distance(tet==5)= 675;
    Distance(tet==9)= 575;
    Distance(tet==13)= 475;
    Distance(tet==17)= 375;
    Distance(tet==21)= 275;
    Distance(tet==25)= 175;
    Distance(tet==29)= 75;
    
    Distance(tet==33)= 775;
    Distance(tet==37)= 675;
    Distance(tet==41)= 575;
    Distance(tet==45)= 475;
    Distance(tet==49)= 375;
    Distance(tet==53)= 275;
    Distance(tet==57)= 175;
    Distance(tet==61)= 75;
    
vert_dist = sind(angle)*Distance(:); % distance of top tetrode site relative to the t"top" of cortex
layer = zeros(length(cells),1);

if tet<=29
    depth_site = (tip_loc_1 - vert_dist);
    if depth_site <= 0
        error('layer less than 0')
    else
        
        layer(depth_site<=50)=1;
        layer((depth_site>50)&(depth_site<=125))=2;
        layer((depth_site>125)&(depth_site<=300))=3;
        layer((depth_site>300)&(depth_site<=425))=4;
        layer((depth_site>425)&(depth_site<=625))=5;
        layer(depth_site>625)=6;
    end
    
else
    depth_site = (tip_loc_2 - vert_dist);
    if depth_site <= 0
        error('layer less than 0')
    else
        
        layer(depth_site<=50)=1;
        layer((depth_site>50)&(depth_site<=125))=2;
        layer((depth_site>125)&(depth_site<=300))=3;
        layer((depth_site>300)&(depth_site<=425))=4;
        layer((depth_site>425)&(depth_site<=625))=5;
        layer(depth_site>625)=6;
        
    end
end
    

% n=zeros(length(cells),1);
% [drift(n).layer]='layer'; %add new field 'layer' to structure "drift"


for i = 1:length(layer); % for loop to add layer info of each unit to each cell of the new field 'layer' in the structure 'drift'
    [drift(i,1).layer]=layer(i,1);
    [drift(i,2).layer]=layer(i,1);
end

save(afile, 'drift','-append')


%if you need to remove a corrupted field, 'layer' for example, evoke the following code:
 field = 'layer';
 drift = rmfield(drift,field);
% % save(afile, 'drift','-append')
