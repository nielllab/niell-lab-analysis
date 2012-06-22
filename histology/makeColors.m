function labels = makeColors(data,n)
min_x=min(data);
max_x=max(data);
labels = zeros(length(data),3);
inds = find(~isnan(data));
labels(inds,:) = scalar2color(data(inds),min_x,max_x);