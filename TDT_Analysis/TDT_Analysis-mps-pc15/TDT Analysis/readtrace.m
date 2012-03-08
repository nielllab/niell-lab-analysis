function [trace tstamp] = readtrace(TTX,trace_name,ch);

max_events=10^6;
invoke(TTX,'CreateEpocIndexing');

MyEpocs = invoke(TTX, 'GetEpocsV', 'xTrg', 0,0, 5000);
nEpocs = size(MyEpocs,2)

N = invoke(TTX, 'ReadEventsV', max_events, trace_name,ch, 0, ...
    0,0,'All');
N
if (N==max_events)
    warning('max number of events aquired');
end
trace = invoke(TTX, 'ParseEvV', 0, N);
trace_TS = invoke(TTX, 'ParseEvInfoV', 0, N, 6);  %   6  = Time Stamp
size(trace)
size(trace_TS)

TS_all = zeros(size(trace));
size(TS_all)

TS_all(1,:)=trace_TS;

block_size = size(trace,1)

trace = trace(:);
trace_TS = trace_TS(:);

dt=mean(diff(trace_TS))/32;

tstamp = zeros(length(trace),1);


for i = 2:block_size;
    i
    TS_all(i,:) = TS_all(1,:)+(i-1)*dt;
end

tstamp=TS_all(:);
tstamp = tstamp(1:length(trace));