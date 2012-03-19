function [osi preftheta]= calcOSI(R,dsi);
theta = 0:(2*pi/length(R)):(2*pi);
if dsi
     theta=theta(1:end-1);
else
    theta=2*theta(1:end-1);  %%% double for orientation
end
R(R<0)=0;

osi = sum(R.*exp(sqrt(-1)*theta'))/sum(R);
preftheta=0.5*angle(osi);
osi = abs(osi);
