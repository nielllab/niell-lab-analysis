function c=scalar2color(x0, min_x, max_x,map);

x0(x0<min_x)=min_x;
x0(x0>max_x)=max_x;

xnew = round((length(map)-1)*(x0-min_x)/(max_x-min_x))+1
c = map(xnew,:);
