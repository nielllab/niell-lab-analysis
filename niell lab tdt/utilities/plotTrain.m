figure
hold on
col = 'rgbmk'
n=0;
for i = 1:length(spikeT);
    sp = spikeT{i}-2*10^5;
    sp = sp(sp >0 &sp<120);
    length(sp)
    if length(sp)<500
        n= n+1;
    for s = 1:length(sp);
        plot([sp(s) sp(s)],[n-0.2 n+0.2],col(mod(n,5)+1));
    end
    end
end
    