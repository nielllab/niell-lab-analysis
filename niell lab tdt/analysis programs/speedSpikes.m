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
histbins = histsize/2:histsize:1000;
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
for cell_n = cell_range
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
    figure
    %subplot(2,1,1)
    plot(histbins,R);
    %plot(histbins,R/max(R));
    hold on
    %plot(histbins, mousehist/max(mousehist),'g');
    plot(histbins, mousehist,'g');
    legend('rate','cm/sec')
    title(sprintf('ch = %d',ch));
    %mousehist = interp1(tdtData.mouseT,tdtData.mouseV,histbins);

    
%     figure
%     plot(mousehist,mousehist2,'o')
    %subplot(2,1,2)
    figure
    plot(mousehist,R,'o');
   hold on
    xlabel('cm/sec')
    ylabel('hz')
    R = R(~isnan(mousehist));
    mousehist = mousehist(~isnan(mousehist));
    
    c = corrcoef(mousehist,R);
     title(sprintf('ch %d - corr = %f',ch,c(1,2)))
     co(cell_n) = c(1,2);
     
     shuffle_c = corrcoef(mousehist(end:-1:1),R);
     shuffle_co(cell_n) = shuffle_c(1,2);
     
     interval = [0 0.5 1 2 4 8 16 32 ];
     d=0;
     p=0;
     for i = 1:length(interval)-1;
         d(i) = mean(R(find(mousehist>interval(i) & mousehist<interval(i+1))));
         p(i) = mean(interval(i:i+1));
     end
 
   plot(p,d,'g');    
     
end

figure
h = histc(co,[-0.5 -0.3 0.3 1]);
bar(h(1:3)/sum(h))


figure
h = hist([co_all ; shuffle_co_all]', -0.4:0.2:0.8)
bar( -0.4:0.2:0.8, h/31)

figure
hold on
for i = 1:length(co);
    
if co(i)>0.3;
    plot(wv(:,i),'g');
else
    plot(wv(:,i),'b');
end
end




