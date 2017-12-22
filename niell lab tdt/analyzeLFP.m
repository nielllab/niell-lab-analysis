pname = uigetdir('C:\data\tdt tanks','block data')
delims = strfind(pname,'\');
selected_path = pname(1 :delims(length(delims))-1)
Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
Block_Name = pname(delims(length(delims))+1 :length(pname))

% Tank_Name='06222011_mlr_stim_awake'
% Block_Name='wn3_72hz'
[fname pname] =uiputfile('*.ps'); psfilename=fullfile(pname,fname);  %%% get ps filename
%psfilename = 'c:/test.ps';   %%% default location
if exist(psfilename,'file')==2;delete(psfilename);end %%% check for previous file


nChan=64;
tic
flags =struct('lfpTseries',1,'lfpSpectra',1);

tdtData= getTDTdata(Tank_Name, Block_Name, 1:nChan, flags);
toc
for ch=1:nChan;
    lfp = tdtData.spectData{ch};
    normalizer = 1:size(lfp,2);
    normalizer = repmat(normalizer,size(lfp,1),1);
    lfpnorm = lfp.*normalizer;
   
    figure   
subplot(2,1,1)
    imagesc(lfpnorm',[0 prctile(lfpnorm(:),95)]);
    axis xy
    df = median(diff(tdtData.spectF{ch}));
    dt = median(diff(tdtData.spectT{ch}));
    set(gca,'YTick',(10:10:80)/df);
    set(gca,'YTickLabel',{'10','20','30','40','50','60','70','80'})
    title(sprintf('channel = %d',ch));
    
    subplot(2,1,2)
    plot(lfpnorm)
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
end

close all


ps2pdf('psfile', psfilename, 'pdffile', [psfilename(1:(end-2)) 'pdf']);
delete(psfilename);
