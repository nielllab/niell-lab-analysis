function T =getSpikeTimes(Tank_Name,Block_Name,ch);

TTX =openTTX(Tank_Name,Block_Name);
invoke(TTX,'CreateEpocIndexing');

N = invoke(TTX, 'ReadEventsV', 10^6, 'Snip', ch, 0,0,0,'ALL');
T= invoke(TTX, 'ParseEvInfoV', 0, N, 6);

invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');