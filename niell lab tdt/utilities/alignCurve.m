function out = alignCurve(in)
mid = floor(length(in)/2);
[y m ] = max(in);
if size(in,2)>size(in,1);
    out = circshift(in,[0,mid-m]);
else
     out = circshift(in,[mid-m,0]);
end