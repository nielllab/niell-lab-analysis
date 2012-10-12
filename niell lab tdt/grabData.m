pname = uigetdir('C:\data\tdt tanks','block data')
delims = strfind(pname,'\');
selected_path = pname(1 :delims(length(delims))-1)
Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
Block_Name = pname(delims(length(delims))+1 :length(pname))

nChan=input('# channels : ');
tic
flags =struct('lfpTseries',1,'lfpSpectra',1,'visStim',1,'mouseOn',1,'newCluster',1);
tdtData= getTDTdata(Tank_Name, Block_Name, 1:nChan, flags);
toc
