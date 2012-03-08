% Matlab codes for gettting stim parameters
% in order to plot histgrams and rasters
% Jianhua Cang 11-05-03
function Bar_Time = episostim_params (Orientation_Angle)
%only works for sweeping bar now

Screen_Distance = 25; % cms
Screen_Width =40;   % cms
Screen_Height = 30; % cms
Stripe_Size_Deg = 5; %in deg

postshutter = 100; % ms; from PCV after trigger
sweeplength = 75; %deg
sweepspeed = 25; % deg/sec

deg_per_cm = atan(1/25)*180/pi; % or 80deg/37cms (Val) 
% if (Orientation_Angle == 0 | Orientation_Angle == 180)
%     Screen_in_Degree = Screen_Width * deg_per_cm;
% elseif (Orientation_Angle == 90 | Orientation_Angle == 270)
%     Screen_in_Degree = Screen_Height * deg_per_cm;
% end
Screen_in_Degree = Screen_Width * deg_per_cm;

start_time = postshutter/1000;
duration = sweeplength/sweepspeed;
if (sweeplength > Screen_in_Degree)
    start_time = start_time + -1*(Screen_in_Degree - sweeplength)/2/sweepspeed; % Assume the bar centers in the screen
    duration = Screen_in_Degree/sweepspeed;
end

Screen_in_Degree;
sweepspeed;
duration;
Bar_Time = [start_time start_time + duration];
       