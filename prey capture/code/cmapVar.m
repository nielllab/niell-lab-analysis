function cm=cmapVar(val,mn,mx,cmap)
%%% map values into a colormap
%%%cm=cmapVar(val,mn,mx,cmap)
%%% val = values; mn = minimum of range; mx = max of range; cmap = colormap to use (should have 64 entries, matlab default)
if ~exist('cmap','var')
cmap = jet;
end
v=(val-mn)/(mx-mn);
if v<0
    v=0;
elseif v>1
    v=1;
end
v = ceil(v*63+1);
cm = cmap(v,:);