function [ W tstamp ] = readWaveAll(TTX,event_code, max_events, max_t )
%reads in continuous data from TDT tank
%and assigns timestamps cmn 2011
ch=32;
% max_events=10^12;
% max_t=10^9;
tic
% N = TTX.ReadEventsV(   max_events ,event_code,ch, 0, ...
%     0,max_t,'All')
N = TTX.ReadEventsV(   max_events ,event_code,ch, 0, ...
    0,max_t,'All')
toc
if (N==max_events)
    warning('max number of events aquired');
end
W = invoke(TTX, 'ParseEvV', 0, N);
block_length = size(W,1);
Wave_TS = invoke(TTX, 'ParseEvInfoV', 0, N, 6);  %   6  = Time Stamp
W=W(:);
W_EV = W';
Wave_TS = Wave_TS(:);

dt=mean(diff(Wave_TS))/block_length;

tstamp = zeros(length(W),1);
for i = 1:length(Wave_TS);
    tstamp((i-1)*block_length+1:i*block_length) = Wave_TS(i)+(0:(block_length-1))*dt;
end
tstamp = tstamp(1:length(W));

ch
TTX.SetGlobalV('WavesMemLimit',2.1e9);
TTX.SetGlobalV('MaxReturn',1e7);
TTX.SetGlobalV('Channel', 0);
TTX.SetGlobalV('T2', 600);
TTX.SetGlobalV('T1', 0);

% event_code
tic
W= TTX.ReadWavesV(event_code)';
toc
size(W)

%size(W)
TTX.SetGlobals('FillItem=Time')
TTX.SetGlobalV('Channel', 1)
tic
t = TTX.ReadWavesV(event_code);
toc

inds = [1 ;find(t>0)];
tstamp = interp1(inds,double(t(inds)),1:max(inds));
W = W(:,1:length(tstamp));


% 
% figure
% plot(W_EV(1:100000));
% hold on
% plot(W(32,1:100000),'g')

figure
plot(W_EV(1:1000000),W(32,1:1000000),'o')
drawnow