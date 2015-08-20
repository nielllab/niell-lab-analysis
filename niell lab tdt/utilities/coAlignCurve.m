function [out1 out2] = alignCurve(in1,in2)
in = max(in1,in2)

mid = ceil(length(in)/4)+1
cv = sum(in.*exp((1:length(in))*sqrt(-1)*2*pi/(length(in)/2)));
ph = mod(angle(cv),2*pi);
m = round(0.5*length(in)*ph/(2*pi));

if m==0
    m=4;
end
if in(m+4)>in(m)
    m=m+4;
end




if size(in,2)>size(in,1);
    out1 = circshift(in1,[0,mid-m]);
    out2 = circshift(in2,[0,mid-m]);
else
     out1 = circshift(in1,[mid-m,0]);
     out2 = circshift(in2,[mid-m,0]);
end