function OptoTag(clustfile, afile, pdfFile,Block_Name,blocknum,nchan, sessionNum);

if ~exist('Block_Name','var');
    SU = menu('recording type','multi-unit','single unit')-1;
    useArgin=0;
else
    SU=1;
    useArgin=1;
end

if SU
    if ~useArgin
        [fname, pname] = uigetfile('*.mat','cluster data');
        clustfile=fullfile(pname,fname);
    end
    
    load(clustfile);
    
    if ~useArgin
        blocknum = listdlg('ListString',Block_Name,'SelectionMode','single');
        [afname, apname] = uigetfile('*.mat','analysis data');
        tagpname = apname;
        afile = fullfile(apname,afname);
    end
    load(afile);
    afile
    if ~exist ('sessionNum','var')
        sessionNum=input('Stim Session number, 1-5?')
    end
    
    [tagpname tagfname] = fileparts(afile);
    Block_Name = Block_Name{blocknum}
    use_afile=1;
    cells
else
    pname = uigetdir('C:\data\','block data')
    delims = strfind(pname,'\');
    selected_path = pname(1 :delims(length(delims))-1)
    Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
    Block_Name = pname(delims(length(delims))+1 :length(pname))
    nchan = input('number of channels : ');
    flags = struct('visStim',1,'MUspike',1);
    data = getTDTdata(Tank_Name,Block_Name,1:4:nchan,flags);
    
end

Block_Name

if useArgin
    psfilename = [pdfFile(1:end-4) Block_Name '.ps'];
else
    [fname pname] =uiputfile('*.ps'); psfilename=fullfile(pname,fname);
end

if exist(psfilename,'file')==2;delete(psfilename);end
% pname = uigetdir('C:\data\TDT tanks','block data')
% delims = strfind(pname,'\');
% selected_path = pname(1 :delims(length(delims))-1)
% Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
% Block_Name = pname(delims(length(delims))+1 :length(pname))

if ~exist('nchan', 'var')
nchan = input('# chans : ');
end

uselaser = 1;


% if SU
%     [fname, pname] = uigetfile('*.mat','cluster data');
%     load(fullfile(pname,fname));
%     for i =1:length(Block_Name);
%         sprintf('%d : %s ',i,Block_Name{i})
%     end
%     block = input('which block to analyze ? ');
%     Block_Name = Block_Name{block}
%     [afname, apname] = uigetfile('*.mat','analysis data');
%     noisepname = apname;
%     afile = fullfile(apname,afname);
%     load(afile);
%     use_afile=1;
%     cells
% end

if SU
 chans=1:nchan;
else
 chans = 1:4:max(nchan,1);
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
edgeDiff=zeros(size(edges));
edgeDiff(2:end)=diff(edges);
edgeDiff(1)=100;

trainEdges=edges(edgeDiff>1);


if SU
    cell_range = 1:size(cells,1)
else
    cell_range=1:4:nchan
end


dt=0.001;
histbins = -0.05:dt:0.05;
longdt=2;
longBins=-17:longdt:34;
psth=zeros(length(cell_range),length(hist(0,histbins)));
longPsth = zeros(length(cell_range),length(hist(0,longBins)));
tdiff_total=zeros(length(cell_range),length(edges));

nfig=0;
for cell_n = cell_range;

  cell_n

    if SU
        channel_no = cells(cell_n,1);
        clust_no = cells(cell_n,2);
        channel_times =spikeT{cell_n} - (blocknum-1)*10^5;
        times = channel_times(channel_times>0 & channel_times<10^5);
    else
        clust_no = [];
        channel_no = cell_n;
        times=tdtData.MUspikeT{cell_n};
    end
    
   OTagfig(cell_n)=figure;
   set(gcf,'position',[200 200 600 700]);
   subplot(3,2,1);
   plot(laserT,tdtData.laserTTL,'g');
   hold on
   plot(0:1:max(times),hist(times,0:1:max(times)));
   title(sprintf('unit %d %d',channel_no,clust_no))
  
  
 %figure  
 subplot(3,2,2)
 hold on
    for t = 1:length(edges);
        tdiff = times-edges(t);
        tdiff = tdiff(abs(tdiff)<max(histbins));
        plot(tdiff*1000,t*ones(size(tdiff)),'ks','MarkerSize',3)
        psth(cell_n,:) = psth(cell_n,:) + hist(tdiff,histbins);
    end
    
    hold on
    plot([0 0],[0 1000],'g','linewidth',2);
    subplot(3,2,4);
    psth(cell_n,:)=psth(cell_n,:)/(dt*length(edges));
   
    
    plot(histbins*1000,psth(cell_n,:),'k','linewidth',2);
    hold on
    plot([0 0],[0 max(psth(cell_n,:))],'g','linewidth',2);
    
    xlabel('msec')
    ylabel('sp/sec')
  
   subplot(3,2,3)
   hold on
   for t = 1:length(trainEdges);
       tdiff = times-trainEdges(t);
       tdiff = tdiff(tdiff<max(longBins) & tdiff>min(longBins));
       plot([tdiff;tdiff],[t*ones(size(tdiff));t*ones(size(tdiff))-1],'k');
       longPsth(cell_n,:) = longPsth(cell_n,:) + hist(tdiff,longBins);
   end
    longPsth(cell_n,:)=longPsth(cell_n,:)/(longdt*length(trainEdges));
    xlim([min(longBins) max(longBins)]);
    subplot(3,2,5);
    plot(longBins,longPsth(cell_n,:));
    hold on
    plot([0 0],[0 max(longPsth(cell_n,:))],'g');
    xlabel('sec');
    ylabel('sp/sec');
    xlim([min(longBins) max(longBins)]);
  
end

%save(afile, 'base','Sdev','PinpFR','PINPed','PINP_sup','-append')
histbins=histbins*1000;

mainfig=figure;
plot(histbins,psth(:,:));hold on
plot([0 0],[0 (max(psth(:)))],'g');
xlabel('msec')
ylabel('sp/sec')

if SU
    load(afile,'psthblocks');
    psthblocks{sessionNum} = psth;
    save(afile,'psth','psthblocks','histbins','-append')
end


% if SU && sessionNum==1;
%      psthStim1=psth;
%    save(afile,'psth','histbins','-append');
% elseif SU && sessionNum==2;
%     psthStim2=psth;
%     save(afile,'psthStim2','-append');
% elseif SU && sessionNum==3;
%     psthStim3=psth;
%     save(afile,'psthStim3','-append');
% elseif SU && sessionNum==4;
%     psthStim4=psth;
%     save(afile,'psthStim4','-append');
% elseif SU && sessionNum==5;
%     psthStim5=psth;
%     save(afile,'psthStim5','-append');
% end


figure(mainfig)
set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');


for i = 1:length(OTagfig)
    figure(OTagfig(i))
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
end
if SU
   ps2pdf('psfile', psfilename, 'pdffile', [psfilename(1:(end-2)) 'pdf']);
    delete(psfilename);
else
   ps2pdf('psfile',psfilename , 'pdffile', [psfilename(1:(end-3)) 'pdf']);
   
end
close all
end
