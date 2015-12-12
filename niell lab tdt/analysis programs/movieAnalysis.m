%%% load cluster file
[fname, pname] = uigetfile('*.mat','cluster data');
clustfile=fullfile(pname,fname);
load(clustfile);

%%% select block and load analysis file
blocknum = listdlg('ListString',Block_Name,'SelectionMode','single');
[afname, apname] = uigetfile('*.mat','analysis data');
afile = fullfile(apname,afname);
load(afile);
[pname fname] = fileparts(afile);
Block_Name = Block_Name{blocknum}

%%% get mouse data and stim timing
flags = struct('mouseOn',1,'visStim',1);
tdtData= getTDTdata(Tank_Name, Block_Name, 1, flags);
tsampDark = tdtData.mouseT;  %%% time points for speed measurements
vsmoothDark = tdtData.mouseV;  %%% speed
framenums = tdtData.frameEpocs;  %%% stim timing (frame # and time)

%%% get spikes for this cell and block
for c = 1:length(spikeT)   
    sp = spikeT{c};
    sp = sp-(blocknum-1)*10^5;
    sp = sp(sp>0 & sp<10^5);
    spikes{c} = sp;   
end

%%% get eye data as well?

%%% now you have spike times, mouse speed, and stim times