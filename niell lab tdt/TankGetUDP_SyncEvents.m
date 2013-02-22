function [dTS, dData] = TankGetUDP_SyncEvents(strTank, iBlock)

strUDP_Tag = 'UDP_';
nChannels = 2;

nEventsMax = 1000000; % max number of events to read
iAllChnls   = 0;    % All Channels
iAllSortCodes  = 0;    % All Sorted Spikes 
dTimeStartAll = 0.0;
dTimeEndAll = 0.0;

dTS = [];
dData = [];


[hTank hFig] = TankOpen( strTank, iBlock);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Count UDP events (read channel 1)
nEvents = hTank.ReadEventsV( nEventsMax, strUDP_Tag, 1, iAllSortCodes, ... 
                             dTimeStartAll, dTimeEndAll, 'NODATA');
if nEvents >0 
  fprintf('Got %d UDP Sync events\n', nEvents);
else 
  fprintf('Could not find any UDP Sync events\n');
  hTank.ReleaseServer;
  delete(hFig);
  return;
end

% Get timestamps
dTS = hTank.ParseEvInfoV( 0, nEvents, 6); % 0 = startIdx, nEvents, 6 = get timestamp

% Get data
for iEvent = 1:nEvents
  nEv = hTank.ReadEventsV( nChannels, strUDP_Tag, iAllChnls, iAllSortCodes, ...  % nChannels = read N events
                               dTS(iEvent), dTimeEndAll, 'ALL');   % dTS() - start time
  dData(:,iEvent) = hTank.ParseEvV( 0, nChannels); % 0 = startIdx, nEvents
end
hTank.ReleaseServer;
delete(hFig);

