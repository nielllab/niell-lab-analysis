function [OSI width peak] = calculate_tuning(A1,A2,B,w);


if A1<0
    A1=-A1;
    A2=-A2;
    B=-B;
    inverted=1
else
    inverted=0;
end

if w<0.1
    w=0.1;
end
if B<0
    B=0;
end
peak = B + A1 + A2*exp(-0.5*(pi/w)^2);
null = B + A1*exp(-0.5*(0.5*pi/w)^2) + A2*exp(-0.5*(0.5*pi/w)^2);

OSI = (peak-null)/(peak+null)
if OSI>0.3333
    theta = 0:.01:(pi/2);
    curve = A1*exp(-0.5*(theta/w).^2) + A2*exp(-0.5*((theta-pi)/w).^2)
    curve = curve/max(curve);
    width = interp1(curve,theta,0.5);
else
    width = 0;
end
if inverted
    peak=-peak;
end
[OSI width peak]
