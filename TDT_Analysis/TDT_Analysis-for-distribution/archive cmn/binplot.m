function [xbin ybin] = binplot(x,y,bins);
for i = 1:length(bins)-1;
    pts = find(x>bins(i) & x<bins(i+1));
    xbin(i) = (bins(i) + bins(i+1))/2;
    ybin(i) = mean(y(pts));
end
