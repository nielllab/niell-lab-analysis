function [ph amp t] = phaseFromTimeseries(lfp,tsig,f0,sigma)
% close all
% lfp = rand(5000,1)-0.5;
% tsig = (1:length(lfp))*.002;
% f0 = 100;
% sigma = 20;

dt = .001;
tsamp = min(tsig):dt:max(tsig);
if round(length(tsamp)/2)==length(tsamp)/2 %%% odd # of samples
    tsamp = tsamp(1:length(tsamp)-1);
end


sampfreq = 1/dt;
lfp_samp = interp1(tsig,lfp,tsamp);
lfp_fft = fftshift(fft(lfp_samp));

% figure
% plot(abs(lfp_fft))
% 
% figure
% plot(tsamp,lfp_samp,'*');
% hold on
% plot(tsamp,lfp_samp);
% plot(tsig,lfp,'go')

dF = 1/max(tsamp);
halflength = floor(0.5*length(lfp_fft));

fsamp = -halflength:1:halflength;
fsamp = fsamp*dF;
filterenvelope = exp(-0.5*(abs(fsamp)-f0).^2 /sigma^2);
figure
hold on
plot(abs(lfp_fft));
plot(filterenvelope*0.5 / max(abs(lfp_fft)),'g');

lfp_filtered = ifft(ifftshift(lfp_fft.*filterenvelope));

% figure
% plot(tsig,lfp);
% hold on
% plot(tsamp,lfp_filtered,'g');

i = 2:length(tsamp)-1;
peaks = find(lfp_filtered(i)>lfp_filtered(i-1) & lfp_filtered(i)>lfp_filtered(i+1)) + 1;
amp = lfp_filtered(peaks);
t = tsamp(peaks);
ph = (1:length(peaks))*2*pi;

figure
plot(tsig,lfp);
hold on
plot(tsamp,lfp_filtered,'g');
plot(t,zeros(size(t)),'*');

axis([0 0.5 -1*10^-4 1*10^-4]);




