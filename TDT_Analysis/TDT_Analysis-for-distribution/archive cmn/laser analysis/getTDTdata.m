function tdtData= getTDTdata(Tank_Name, Block_Name, nChan, MUspikeOn, oldClusterOn, newClusterOn, lfpTseriesOn, lfpSpectrumOn,visStimOn, mouseOn, laserOn);
%%% function to read in all relevant data from a block
%%% and convert to matlab data structure
%%% cmn 06/2011


tdtData=struct('MUspikeT',[], 'spikeT',[], 'lfpT' ,[],...
    'lfpData',[], 'spectT',[], 'spectF',[], 'spectData' ,[], ...
    'frameEpocs' ,[],'stimEpocs' ,[], ...
    'mouseT',[], 'mouseV' ,[],'laserT',[], 'laserTTL',[])

TTX = openTTX(Tank_Name,Block_Name);
invoke(TTX,'CreateEpocIndexing');
ep = invoke(TTX, 'GetEpocsV', 'xTrg', 0,0, 1000)

MUspikeOn=1
nChan=16
max_events = 10^6;
max_time = 10^9;
if MUspikeOn
    Event_Name_Snip='Snip'

    for ch = 1:nChan
 N= invoke(TTX, 'ReadEventsV', max_events, Event_Name_Snip, ch, 0, 0,max_time,'ALL')
        tdtData.MUspikeT{ch} = invoke(TTX, 'ParseEvInfoV', 0, N, 6)   %   3  = event codes
    end
end

if oldClusterOn
    tdtData.spikeT= getBlockSpikes;
elseif newClusterOn
    [afname, apname] = uigetfile('*.mat','analysis data');
    if afname~=0
        afile = fullfile(apname,afname);
        load(afile,spikeT);
        tdtData.spikeT=spikeT
    end

end

if lfpTseriesOn | lfpSpectrumOn
  [tdtData.lfpT tdtData.lfpData tdtData.spectT tdtData.spectF tdtData.spectData] = analyzeLFP_chronux(TTX,nChan,lfpTseriesOn,lfpSpectrumOn);   
end

if visStimOn
    tdtData.frameEpocs = invoke(TTX, 'GetEpocsV', 'fTrg',0,0,max_events);
    tdtData.stimEpocs = invoke(TTX, 'GetEpocsV', 'xTrg', 0,0,max_events);
end

if laserOn
    [tdtData.laserT tdtData.laserTTL] = read_laser(TTX);
end

if mouseOn
    [tdtData.mouseT tdtData.mouseV groomstate] =getBlockVelocity_groom(Tank_Name,Block_Name);
end
