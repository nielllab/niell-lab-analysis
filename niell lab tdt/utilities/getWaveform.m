function [inh mid] = getWaveform(clustfile,afile,redo)
%%% calculate cell type based on waveform
%%% this is pretty simple but nice to do it just in one line!

load(afile,'inh','mid','trough2peak','trough_depth','peak_height');

if ~exist('inh','var')  | isempty(inh) | redo
    try
        alldata( :,1) = trough2peak;
        alldata( :,2) = -trough_depth./peak_height;
        inh = alldata(:,2)<1.75 & alldata(:,1)<8.5 ;
        mid = alldata(:,2)>1.75 & alldata(:,1)<10;
        
        figure
        plot(alldata(find(inh),1),alldata(find(inh),2),'ro');
        hold on
        plot(alldata(find(~inh),1),alldata(find(~inh),2),'ko');
        plot(alldata(find(mid),1),alldata(find(mid),2),'go');
        
        
    catch
        inh = [];
        mid = [];
    end
    save(afile,'inh','mid','-append')
end

