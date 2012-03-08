function [tstamp laser] = read_laser(TTX);
W=0;
max_events=10^6;
N = invoke(TTX, 'ReadEventsV', max_events, 'TTL_',1, 0, 0,0,'All');
%            N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_EEG,tet, 0, ...
%             0,0,'All');
if (N==max_events)
    warning('max number of events aquired');
end
laser = invoke(TTX, 'ParseEvV', 0, N);
laser=laser(:);
Wave_TS = invoke(TTX, 'ParseEvInfoV', 0, N, 6);  %   6  = Time Stamp
Wave_TS = Wave_TS(:);

dt=mean(diff(Wave_TS))/64;
tstamp = zeros(length(W),1);
for i = 1:length(Wave_TS);
    tstamp((i-1)*64+1:i*64) = Wave_TS(i)+(0:63)*dt;
end
tstamp = tstamp(1:length(laser));
