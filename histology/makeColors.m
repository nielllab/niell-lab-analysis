function [labels cmap clim] = makeColors(data,n,type,color)
min_x=min(data);
max_x=max(data);
clim = [min_x max_x];
labels = zeros(length(data),3);
if exist('n','var')
    labels(:)=n;
end
inds = find(~isnan(data));

if exist('type','var');
    if strcmp(type,'qual')
        cmap = cbrewer(type,color,length(min_x:max_x));
    else
        cmap= cbrewer(type,color,64);
    end
else
    cmap=colormap(jet);
end
labels(inds,:) = scalar2color(data(inds),min_x,max_x,cmap);