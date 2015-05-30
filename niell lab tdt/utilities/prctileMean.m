function out = prctileMean(data,prct);
n = size(data,1);
ndrop = round(n*prct/100);
if ndrop<1
    ndrop=1;
end
if ndrop>floor(n/2)
    ndrop=floor(n/2);
end


for i = 1:size(data,2);
    sortData = sort(data(:,i));
    out(1,i) = mean(sortData(ndrop+1:end-ndrop));
end
