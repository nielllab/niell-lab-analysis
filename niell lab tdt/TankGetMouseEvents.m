function [dTS, dData] = TankGetMouseEvents(strTank, strBlock)

strUDP_Tag = 'UDP_';
nChannels = 5;

nEventsMax = 2000000; % max number of events to read
iAllChnls   = 0;    % All Channels
iAllSortCodes  = 0;    % All Sorted Spikes 
dTimeStartAll = 0.0;
dTimeEndAll = 0.0;

dTS = [];
dData = [];


[hTank hFig] = TankOpen( strTank, strBlock);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Count Mouse UDP events (read channel 1)
nEvents = hTank.ReadEventsV( nEventsMax, strUDP_Tag, 1, iAllSortCodes, ... 
                             dTimeStartAll, dTimeEndAll, 'NODATA');
if nEvents >0 
  fprintf('Got %d Mouse UDP events\n', nEvents);
else 
  fprintf('Could not find any Mouse UDP events\n');
  hTank.ReleaseServer;
  delete(hFig);
  return;
end

% Get timestamps
dTS = hTank.ParseEvInfoV( 0, nEvents, 6); % 0 = startIdx, nEvents, 6 = get timestamp

% Get data
tic
% for iEvent = 1:nEvents
%   nEv = hTank.ReadEventsV( nChannels, strUDP_Tag, iAllChnls, iAllSortCodes, ...  % nChannels = read N events
%                                dTS(iEvent), dTimeEndAll, 'ALL');   % dTS() - start time
%   dData(:,iEvent) = hTank.ParseEvV( 0, nChannels); % 0 = startIdx, nEvents
% end
% 
  nEv = hTank.ReadEventsV( nEventsMax, strUDP_Tag, iAllChnls, iAllSortCodes, ...  % nChannels = read N events
                               dTimeStartAll, dTimeEndAll, 'ALL')   % dTS() - start time
  dData = hTank.ParseEvV( 0, nEv); % 0 = startIdx, nEvents
  dData = reshape(dData,5,nEv/5);

toc
hTank.ReleaseServer;
delete(hFig);

