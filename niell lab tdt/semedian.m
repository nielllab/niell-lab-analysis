function s= semedian(data)
s = nanstd(bootstrp(1000,@(x) nanmedian(x,1),data));
