function [smdata MaxIdx MinIdx Amp Vel mmIdx]=StateTransition(data,smooth,dt,baseline,IPI,A)

% figure
% subplot(1,2,1)
% plot(data);
% axis square

%% set values for parameters for findpeaks function
y=data;
time=length(data);
x=[1:1:time];


 sz=smooth; % smoothing factor of 3 seconds match to smootwisth parameter above to look at your smoothened series
 smdata=conv(data,ones(1,sz),'same')/sz;
 
 
 
% [...] = FINDPEAKS(X,'MINPEAKHEIGHT',MPH) finds only those peaks that
%   are greater than MINPEAKHEIGHT MPH. Specifying a minimum peak height
%   may help in reducing the processing time. MPH is a real valued scalar.
%   The default value of MPH is -Inf.
%
%   [...] = FINDPEAKS(X,'MINPEAKDISTANCE',MPD) finds peaks that are at
%   least separated by MINPEAKDISTANCE MPD. MPD is a positive integer
%   valued scalar. This parameter may be specified to ignore smaller peaks
%   that may occur in close proximity to a large local peak. For example,
%   if a large local peak occurs at index N, then all smaller peaks in the
%   range (N-MPD, N+MPD) are ignored. If not specified, MPD is assigned a
%   value of one. 
%
%   [...] = FINDPEAKS(X,'THRESHOLD',TH)finds peaks that are at least
%   greater than their neighbors by the THRESHOLD TH. TH is real valued
%   scalar greater than or equal to zero. The default value of TH is zero.



[Maxima,MaxIdx] = findpeaks(smdata,'minpeakdistance',IPI,'minpeakheight',A); 
DataInv = 1.01*max(smdata) - smdata;

[Minima,MinIdx] = findpeaks(DataInv,'minpeakdistance',IPI,'minpeakheight',A);

%  figure
%  plot(DataInv-43)

%The true minima will then be:
Minima = smdata(MinIdx);

 
%  figure
%  plot(smdata);hold on
%  plot(MaxIdx,180,'g*');hold on 
%  plot(MinIdx,180,'ro')
 
 mmIdx=cat(1,MaxIdx,MinIdx);mmIdx=sort(mmIdx);
 
 figure
 plot(smdata);hold on
 plot(mmIdx,baseline,'k*');hold on
 plot(MaxIdx,4,'c*')
 title 'peak' %and valley indices'
 
 Maxmins=smdata(mmIdx);
 Amp=diff(Maxmins);
 numframe=diff(mmIdx);
 timeSac=numframe*dt;
 Vel=Amp./timeSac;
 
 mx=smdata(MaxIdx);
 
%  figure
%  nhist(Amp);
%  title 'distribution of Amp between local max and min'
%  
%  figure
%  hist(Vel);
%  title 'distribution of Vel between local max and min'
%  
 figure
 plot(abs(Amp),abs(Vel),'ko')
 title '(x) local peak to valley or valley to peak Amplitude vs. (Y) angular Velocity VtoP or P to V'

%%
 


