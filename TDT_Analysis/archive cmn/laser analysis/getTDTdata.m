function tdtData= getTDTdata(Tank_Name, Block_Name, nChan, flags);
%%% function to read in all relevant data from a block
%%% and convert to matlab data structure
%%% cmn 06/2011


tdtData=struct('streamV',[],'streamT',[],'MUspikeT',[], 'spikeT',[], 'lfpT' ,[],...
    'lfpData',[], 'spectT',[], 'spectF',[], 'spectData' ,[], ...
    'frameEpocs' ,[],'stimEpocs' ,[], ...
    'mouseT',[], 'mouseV' ,[],'laserT',[], 'laserTTL',[])

TTX = openTTX(Tank_Name,Block_Name);
invoke(TTX,'CreateEpocIndexing');
ep = invoke(TTX, 'GetEpocsV', 'xTrg', 0,0, 1000)


max_events = 10^6;
max_time = 10^9;
if flags.MUspike
    Event_Name_Snip='Snip'
    
    for ch = 1:nChan
        N= invoke(TTX, 'ReadEventsV', max_events, Event_Name_Snip, ch, 0, 0,max_time,'ALL')
        tdtData.MUspikeT{ch} = invoke(TTX, 'ParseEvInfoV', 0, N, 6)   %   3  = event codes
    end
end

if flags.stream
    event_code = 'pAll';
    max_events = 10^6;
    max_t = 10^9;
for ch = 1:nChan
        [tdtData.streamV{ch} tdtData.streamT{ch}] = readWave(TTX,ch, event_code,max_events,max_t);
    end
end

if flags.oldCluster
    tdtData.spikeT= getBlockSpikes;
elseif flags.newCluster
    [afname, apname] = uigetfile('*.mat','analysis data');
    if afname~=0
        afile = fullfile(apname,afname);
        load(afile,spikeT);
        tdtData.spikeT=spikeT
    end
    
end

if flags.lfpTseries | flags.lfpSpectra
    [tdtData.lfpT tdtData.lfpData tdtData.spectT tdtData.spectF tdtData.spectData] =...
        analyzeLFP_chronux(TTX,nChan,flags.lfpTseries,flags.lfpSpectra);
end

if flags.visStim
    tdtData.frameEpocs = invoke(TTX, 'GetEpocsV', 'fTrg',0,0,max_events);
    tdtData.stimEpocs = invoke(TTX, 'GetEpocsV', 'xTrg', 0,0,max_events);
end

if flags.laser
    [tdtData.laserT tdtData.laserTTL] = read_laser(TTX);
end

if flags.mouse
    [tdtData.mouseT tdtData.mouseV groomstate] =getBlockVelocity_groom(Tank_Name,Block_Name);
end

invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');