function [tstamp] = read_camTTL(TTX);
W=0;
max_events=10^6;
N = invoke(TTX, 'ReadEventsV', max_events, 'CAM_',1, 0, 0,0,'All');
%            N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_EEG,tet, 0, ...
%             0,0,'All');
if (N==max_events)
    warning('max number of events aquired');
end
laser = invoke(TTX, 'ParseEvV', 0, N);
size(laser)
laser=laser(:);

Wave_TS = invoke(TTX, 'ParseEvInfoV', 0, N, 6);  %   6  = Time Stamp
size(Wave_TS)
Wave_TS = Wave_TS(:);

    bl_length = 256;


dt=mean(diff(Wave_TS))/bl_length;
tstamp = zeros(length(W),1);
for i = 1:length(Wave_TS);
    tstamp((i-1)*bl_length+1:i*bl_length) = Wave_TS(i)+(0:(bl_length-1))*dt;
end
tstamp = tstamp(1:length(laser));

rising = find(diff(laser)>2.5);
tstamp = tstamp(rising);

