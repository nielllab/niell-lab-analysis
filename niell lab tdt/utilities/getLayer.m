function layerInfo = getLayer(clustfile,afile,tip1,tip2, angle, redo )
%%% calculate layer info based on geometry;
%%% tip = depth of tip; angle = angle relative to horizontal
%%% if only one shank, set tip2=[];
%%% checks to see if data already exists
%%% otherwise reads it in and saves to analysis file, with entry for each block
%%% returns data for specified block
%%% layerInfo.sites = all recording sites, layerInfo.units = selected untis

load(afile,'layerInfo','peakchan');

if ~exist('layerInfo','var')  | isempty(layerInfo.units) | redo
  try
        %%% geometry of sites
        dist(1:32) = 775:-25:0;
        if ~isempty(tip2)
            dist(33:64) = dist(1:32);
        end
        
        vert_dist = sind(angle)*dist;
        vert_dist(1:32) = tip1 - vert_dist(1:32);
        if ~isempty(tip2)
            vert_dist(33:64) = tip2 - vert_dist(33:64);
        end
        
        layer = zeros(size(vert_dist));       
        layer(vert_dist<0)=0;
        layer(vert_dist>0 & vert_dist<=50)=1;
        layer(vert_dist>50 & vert_dist<=125) = 2;
        layer(vert_dist>125 & vert_dist<=300)=3;
        layer(vert_dist>300 & vert_dist<=425)=4;
        layer(vert_dist>425 & vert_dist<=625) = 5;
        layer(vert_dist>625)=6;
        
        for i= 1:length(peakchan);
            unitlayer(i) = layer(peakchan(i));
        end
        layerInfo.sites = layer;
        layerInfo.units = unitlayer;
        
  catch
        layerInfo.sites=[];
        layerInfo.units=[];
  end
    save(afile,'layerInfo','-append')
end
