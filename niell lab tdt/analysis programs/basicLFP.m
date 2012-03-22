function basicLFP

pname = uigetdir('C:\data\tdt tanks','block data')
delims = strfind(pname,'\');
selected_path = pname(1 :delims(length(delims))-1)
Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
Block_Name = pname(delims(length(delims))+1 :length(pname))

nChan=input('# channels : ');
tic
flags =struct('lfpTseries',1,'lfpSpectra',1);
tdtData= getTDTdata(Tank_Name, Block_Name, 1:nChan, flags);
toc

% figure
% lfp =tdtData.lfpData{1};
% lfpT = tdtData.lfpT{1};
% plot(lfpT(1:1000),lfp(1:1000));
% lfp =tdtData.lfpData{16};
% hold on
% plot(lfpT(1:1000),lfp(1:1000),'g');


for ch=1:nChan;
    lfp = tdtData.spectData{ch};
    normalizer = 1:size(lfp,2);
    normalizer = repmat(normalizer,size(lfp,1),1);
    lfpnorm = lfp.*normalizer;
   
    figure   
    imagesc(lfpnorm',[0 prctile(lfpnorm(:),95)]);
    axis xy
    df = median(diff(tdtData.spectF{ch}));
    dt = median(diff(tdtData.spectT{ch}));
    set(gca,'YTick',(10:10:80)/df);
    set(gca,'YTickLabel',{'10','20','30','40','50','60','70','80'})
    title(sprintf('channel = %d',ch));
end