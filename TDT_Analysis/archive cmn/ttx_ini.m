% Matlab codes for reading from TTank
% and plot histgrams and rasters with each orientation
% Jianhua Cang 06-27-03; 10-27/03  
% connect to the tank and read the block
% temporary code, use GUI in the future

%cd C:\MATLAB6p5\TDT

close all
%Tank_Name='080805'

Select_Duration = [1330 1630];
%Select_Duration = [620 1160];
%channel_no = 4;
Orientation_Angle = 0; % for contstim
max_events =10000;

plot_duration=8; %in second
hist_range=[0:0.1:9];
axis_range=[0 plot_duration 0 4];

Event_Name_Snip='Snip'
Event_Name_Wave='PDec'
Sample_Interval=0.04096 % 24414.0625Hz
Sample_Number_Snip=64
Dec_Factor=32; % 
T_Snip=[0:Sample_Interval:(Sample_Number_Snip-3)*Sample_Interval]; 
T_Wave=[Select_Duration(1): Sample_Interval/1000*Dec_Factor :Select_Duration(2)];
Save_File_Snip=[Tank_Name Block_Name '_s.mat']
Save_File_Wave=[Tank_Name Block_Name '_w.mat']

TTX = actxcontrol('ttank.x', [1 1 1 1]);
if (invoke(TTX, 'ConnectServer', 'Local', 'jc') ~= 1)
  err = 'error connecting to server'
end
if (invoke(TTX, 'OpenTank', Tank_Name, 'r') ~= 1)
  err = 'error opening tank'
end
if (invoke(TTX, 'SelectBlock', Block_Name) ~= 1)
  err = 'error selecting block'
end
