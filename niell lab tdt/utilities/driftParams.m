tflist = [2 8];
sflist = [0.01 0.02 0.04 0.08 0.16 0.32];
orientlist = 0:45:315;
for rep = 1:2
    for s = 1:6
        for o=1:8;
            cond = (o-1)*12+(s-1)*2+rep;
            tf(cond)= tflist(rep);
            sf(cond)=sflist(s);
            orient(cond)=orientlist(o);
        end
    end
end