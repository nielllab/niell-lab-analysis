function [osi preftheta]= calcCVOSI_new(R,dsi);
theta = 0:(2*pi/length(R)):(2*pi);
if dsi
     theta=theta(1:end-1);
else
    theta=2*theta(1:end-1);  %%% double for orientation
end
R(R<0)=0;
%osi = sum(R.*exp(sqrt(-1)*theta'))/sum(R);
osi = sum(R'.*exp((1:length(R))*sqrt(-1)*2*pi/(0.5*length(R)))); %cv


    th = angle(osi);
    th(th<0) = th(th<0)+2*pi;
    if dsi
        preftheta=th;
    else
        preftheta=0.5*th;
    end

osi = abs(osi);
tuning = R;
tuning = tuning-min(tuning);
end
