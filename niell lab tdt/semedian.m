function s= semedian(data)
s = nanstd(bootstrp(1000,@(x) nanmedianMW(x,1),data));
