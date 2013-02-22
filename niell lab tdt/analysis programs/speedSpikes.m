clear all

pname = uigetdir('C:\data\TDT tanks','block data')
delims = strfind(pname,'\');
selected_path = pname(1 :delims(length(delims))-1)
Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
Block_Name = pname(delims(length(delims))+1 :length(pname))

% Tank_Name='06222011_mlr_stim_awake'
% Block_Name='wn3_72hz'

nchan = input('# chans : ');

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


    flags = struct('mouseOn',1,'MUspike',1)


tdtData= getTDTdata(Tank_Name, Block_Name, chans, flags);
%tdtData= getTDTdata(Tank_Name, Block_Name, nChan, 0, 0, 0, 1,1, 0, 0, 0);

histsize = 2;
histbins = histsize/2:histsize:max(tdtData.mouseT);
%histbins = histsize/2:histsize:1000;
clear mousehist
    for i = 1:length(histbins);
        
        mousehist(i) = mean(tdtData.mouseV(find(tdtData.mouseT>histbins(i)-histsize/2 & tdtData.mouseT<histbins(i)+histsize/2)));
    end

if SU
    cell_range = 1:size(cells,1)
else
    cell_range=1:4:nchan
end
maxV = max(tdtData.mouseV);
close all

   d=0;
     p=0;
     
for cell_n = cell_range                     %runs through each unit
    ch= cell_n;
    if SU
        channel_no = cells(cell_n,1);
        clust_no = cells(cell_n,2);
        channel_times =spikeT{cell_n} - (block-1)*10^5;
        times = channel_times(channel_times>0 & channel_times<10^5);
  
    else
        clust_no = [];
        channel_no = cell_n;
       
        times=tdtData.MUspikeT{cell_n};
    end
   
    R = hist(times,histbins);
    speedfig(cell_n) = figure;
    subplot(2,1,1)
    plot(histbins,R);                       %firing rate vs time
    hold on
    plot(histbins, mousehist,'g');          %speed vs time
    legend('rate','cm/sec')
    title(sprintf('ch = %d',ch));
    %mousehist = interp1(tdtData.mouseT,tdtData.mouseV,histbin)
%     figure
%     plot(mousehist,mousehist2,'o')
    %subplot(2,1,2)
   subplot(2,1,2)
    plot(mousehist,R,'o');                  %firing rate vs speed
   hold on
    xlabel('cm/sec')
    ylabel('hz')
    R = R(~isnan(mousehist));
    mousehist_nonan = mousehist(~isnan(mousehist));
    
    [c,pvals] = corrcoef(mousehist_nonan,R);
     title(sprintf('ch %d - corr = %f',ch,c(1,2)))
     co(cell_n) = c(1,2);
     allPvals(cell_n)=pvals(1,2);
     
     
     [shuffle_c,shuffle_p] = corrcoef(mousehist_nonan(randperm(length(mousehist_nonan))),R);
     shuffle_co(cell_n) = shuffle_c(1,2);
     allShuffle_p(cell_n)=shuffle_p(1,2);
     
     interval = [0 0.5 1 2 4 8 16 32 64];
  
     for i = 1:length(interval)-1;
         d(cell_n,i) = mean(R(find(mousehist_nonan>interval(i) & mousehist_nonan<interval(i+1))));
         p(cell_n,i) = mean(interval(i:i+1));
     end
 
   plot(p(cell_n,:),d(cell_n,:),'g');    
     
end


histfig=figure
h = hist([co ; shuffle_co]', -0.4:0.2:0.8)
bar( -0.4:0.2:0.8, h)
legend({'raw','shuffled'})

wvfig = figure
hold on
for i = 1:length(co);
    
if co(i)>0.5;
    plot(wv(:,i),'g');
elseif co(i)<-0.5
    plot(wv(:,i),'r');
    else
         plot(wv(:,i),'b');
end
end


if SU
    save(fullfile(apname,afname),'co','shuffle_co' ,'p','d','-append');
end


[fname pname] =uiputfile('*.ps'); psfname=fullfile(pname,fname);
if exist(psfname,'file')==2;delete(psfname);end

figure(histfig)
set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',fullfile(pname,fname),'-append');
figure(wvfig)
set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',fullfile(pname,fname),'-append');

for i = 1:length(speedfig)
    figure(speedfig(i))
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',fullfile(pname,fname),'-append');
end
if SU
    ps2pdf('psfile', psfname, 'pdffile', [psfname(1:(end-3)) 'SU.pdf']);
else
   ps2pdf('psfile', psfname, 'pdffile', [psfname(1:(end-3)) 'MU.pdf']); 
end
