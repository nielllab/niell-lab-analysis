clist = [-1 1];
for c = 1:2
    for y = 1:6;
        for x=1:9
            cond = (x-1)*12 +(y-1)*2 + c;
            contrast(cond) = clist(c);
            xpos(cond)=x;
            ypos(cond)=y;
        end
    end
end