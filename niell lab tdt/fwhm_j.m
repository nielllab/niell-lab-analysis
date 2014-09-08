function[Fullw  peak]=fwhm_j(x, y)
% This function determines  full width at half 
% maximum of a peak if the input data  has two columns:
% Column 1 = x
% Column 2 = y
%Coded by Ebo Ewusi-Annan
%University of Florida 
%August 2012
% x = data(:,1);
% y= data(:,2);
maxy = max(y); 



f = find(y==maxy); 
range=max(x);
if f>1 & f<range
peak = x(f);% ignore Matlabs suggestion to fix!!!
elseif f==1 
    peak=x(f)+2
elseif f==max(x)
    peak=x(f)-2
end

%y1=y-(min(y))
%[new_min idx_new_min]=min(y(x<peak))
base=min(y)
y1=y-base % new y values with local minimum subtracted off
y1(y1<0)=nan
y2= y1./max(y1);

ydatawr(:,1) = y2;
ydatawr(:,2) = x;
newFit1=find(x>= peak);
newFit2=find(x < peak);



ydatawr2 = ydatawr(min(newFit1):max(newFit1),:);
ydatawr3 = ydatawr(min(newFit2):max(newFit2),:);

sp1 = spline(ydatawr2(:,1),ydatawr2(:,2),0.4);
sp2 = spline(ydatawr3(:,1),ydatawr3(:,2),0.4);
Fullw = sp1-sp2;

 end