clear all
close all

pname = uigetdir('C:\data\TDT tanks','block data')
delims = strfind(pname,'\');
selected_path = pname(1 :delims(length(delims))-1)
Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
Block_Name = pname(delims(length(delims))+1 :length(pname))


nchan = input('# chans : ');
uselaser = 1
chans = 1:4:nchan;
SU = input('multiunit (0) or single-unit (1) : ');

if SU
    [fname, pname] = uigetfile('*.mat','cluster data');
    load(fullfile(pname,fname));
    for i =1:length(Block_Name);
        sprintf('%d : %s ',i,Block_Name{i})
    end
    block = input('which block to analyze ? ');
    Block_Name = Block_Name{block}
    [afname, apname] = uigetfile('*.mat','analysis data');
    noisepname = apname;
    afile = fullfile(apname,afname);
    load(afile);
    use_afile=1;
    cells
end

if SU
    flags = struct('laserOn',1)
else
    flags = struct('laserOn',1,'MUspike',1)
end

tdtData= getTDTdata(Tank_Name, Block_Name, chans, flags);

laserT = tdtData.laserT;
dV = diff(tdtData.laserTTL);
edges = laserT(find(dV>2.5)+1);


if SU
    cell_range = 1:size(cells,1)
else
    cell_range=1:4:nchan
end


dt=0.001;
histbins = -0.1:dt:0.1;
psth=zeros(length(cell_range),length(hist(0,histbins)));
nfig=0;
for c = 1:length(cell_range)
    cell_n=cell_range(c);
    ch= cell_n;
    if SU
        channel_no = cells(cell_n,1)
        clust_no = cells(cell_n,2)
        channel_times =spikeT{cell_n} - (block-1)*10^5;
        times = channel_times(channel_times>0 & channel_times<10^5);
    else
        clust_no = [];
        channel_no = cell_n;
        times=tdtData.MUspikeT{cell_n};
    end
    
    histfig(c)=figure;
   subplot(2,1,1);
   plot(laserT,tdtData.laserTTL,'g');
   hold on
   plot(0:1:max(times),hist(times,0:1:max(times)));
   
   title(sprintf('ch %d ',cell_n))
  
    for t = 1:length(edges);
        tdiff = times-edges(t);
        tdiff = tdiff(abs(tdiff)<max(histbins));
        psth(c,:) = psth(c,:) + hist(tdiff,histbins);
    end
    psth(c,:)=psth(c,:)/(dt*length(edges));
    subplot(2,1,2);
    plot(histbins*1000,psth(c,:));
    hold on
    plot([0 0],[0 max(psth(c,:))],'g');
    xlabel('msec')
    ylabel('sp/sec')
end
mainfig=figure;
plot(histbins*1000,psth);
plot([0 0],[0 max(psth(:))],'g');
xlabel('msec')
    ylabel('sp/sec')

if SU
    save(fullfile(apname,afname),'psth','histbins','-append');
end


[fname pname] =uiputfile('*.ps'); psfname=fullfile(pname,fname);
if exist(psfname,'file')==2;delete(psfname);end

figure(mainfig)
set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',fullfile(pname,fname),'-append');


for i = 1:length(histfig)
    figure(histfig(i))
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',fullfile(pname,fname),'-append');
end
if SU
    ps2pdf('psfile', psfname, 'pdffile', [psfname(1:(end-3)) 'SU.pdf']);
else
   ps2pdf('psfile', psfname, 'pdffile', [psfname(1:(end-3)) 'MU.pdf']); 
end
