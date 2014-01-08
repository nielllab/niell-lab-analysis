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
tip_loc = input('depth in uM of electrode tip : '); %apply standard depth in uM for location of tip in top, middle or bottom third of any layer
angle = input('angle on electrode penetration : ');

 %extracts first site of all tetrodes
 
 tet = cells(:,1);
    
 %converts tetrode location into phsyical distance from tip of electrode
    
    Distance = zeros(length(tet),1);
   
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
    
vert_dist = sind(angle)*Distance(:);
    
depth_site = (tip_loc - vert_dist);
if depth_site <= 0
    error('layer less than 0')
else
    
layer = zeros(length(tet),1);
layer(depth_site<=50)=1;
layer((depth_site>50)&(depth_site<=125))=2;
layer((depth_site>125)&(depth_site<=250))=3;
layer((depth_site>250)&(depth_site<=375))=4;
layer((depth_site>375)&(depth_site<=650))=5;
layer(depth_site>650)=6;


%layer = calculate_layer(tip_loc,angle);% run function to calculate layer info if you have run the old drift

[drift(length(cells)).layer]='layer'; %add new field 'layer' to structure "drift"


for i = 1:length(layer); % for loop to add layer info of each unit to each cell of the new field 'layer' in the structure 'drift'
    [drift(i,1).layer]=layer(i,1);
    [drift(i,2).layer]=layer(i,1);
end
save(afile, 'drift','-append')

end

%if you need to remove a corrupted field, 'layer' for example, evoke the following code:
%  field = 'layer';
%  drift = rmfield(drift,field);
%  save(afile, 'drift','-append')
