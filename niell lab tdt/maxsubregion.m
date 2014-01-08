function img = maxsubregion(sta,width,c)


if length(sta)<width
    img = sta;
    return;
end

width = width-1;

if strcmp(c,'neg');
    test = -sta;
elseif strcmp(c,'abs')
    test = abs(sta);
else
    test=sta;
end
[linemax y] = max(test,[],1);
linemax

[m x] = max(linemax);

x
y=y(x)

x = max(x,width/2+1);
x = min(x,size(sta,2)-width/2);


y = max(y,width/2+1);
y = min(y,size(sta,1)-width/2);

x
y



img = sta(y-width/2 :y+width/2, x-width/2: x+ width/2); 

