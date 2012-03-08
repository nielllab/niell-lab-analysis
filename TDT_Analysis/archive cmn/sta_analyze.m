function [wpref ,theta_pref, A1, w, null,b] = sta_analyze(subfield);

sta_fft =fftshift(abs(fft2(subfield)));

wmax = size(subfield,2);

dx = floor(wmax/2)+1;
dy =dx;


[wx wy ] = meshgrid((1:wmax)-dx,(1:wmax)-dy);
i=sqrt(-1);
sta_fft(wx<-dx/3) =0;
sta_fft(wx>dx/3) =0;
sta_fft(wy<-dx/3) =0;
sta_fft(wy>dx/3) =0;

%sta_fft(wx==0 & wy==0)=0;

mu = sum(sum( sta_fft.*exp(2*i*atan2(wy,wx))))./sum(sum(sta_fft));

theta_pref = angle(mu)/2

windex = 0:0.1:dx/3;

w_tuning = interp2(wx,wy,sta_fft,windex*cos(theta_pref),windex*sin(theta_pref),'spline');

[peak_amp index] = max(w_tuning);
peak_amp
wpref = windex(index)

wxpref = cos(theta_pref)*wpref
wypref = sin(theta_pref)*wpref

figure
imagesc(sta_fft);
hold on;
plot(wxpref+dx,wypref+dy,'g*');


dtheta = 2*pi/16;
theta_index = 0:dtheta:2*pi-dtheta;

theta_tuning = interp2(wx,wy,sta_fft,cos(theta_index)*wpref,sin(theta_index)*wpref);


[theta_pref OSI A1 A2 w b null yfit] = fit_tuningcurve(theta_tuning, theta_index);





