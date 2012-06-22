function c=scalar2color(x0, min_x, max_x,cmap);

if ~exist('cmap','var')
    %cmap = 'cool';
    cmap=zeros(64,3);
    cmap(33:64,1)=(0:31)/31;
    cmap(1:32,2) = (0:31)/31;
    cmap(33:64,2) = (31:-1:0)/31;
    cmap(1:32,3) = (31:-1:0)/31;
end
map = colormap(cmap);

x0(x0<min_x)=min_x;
x0(x0>max_x)=max_x;

xnew = round((length(map)-1)*(x0-min_x)/(max_x-min_x))+1
c = map(xnew,:);
