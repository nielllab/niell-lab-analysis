close all
clear all
dt = .01;
tmax = 20;
t = 0:dt:tmax;
% R=(t>3 &t<3.5)*10;
% R(t>5 & t<5.5)=5;
% R(t>10 & t<10.5) = 20;
% R(t>15 & t<15.5) =4;
R=(t>10 & t<10.5)*2;
R(t>9 & t<9.5) =2;
R(t>11& t<11.5) =2;


%R(t>2& t<2.5)=10;
rep =1;

spike = poissrnd(R*dt);
N(rep) = sum(spike(t>10 & t<10.5));
figure
plot(spike)

tau = 2;
impulse = exp(-1*(0:dt:tmax)/tau);

baselineF=5000;   %%% photons per microsecond at rest 1200 = 2500@ 5hz =.02 rms noise
dF = .02;
F = baselineF + baselineF*dF*conv(impulse,spike) ;
% figure
% plot(F)

sample_rate = 5;
sampT = 1/sample_rate;
dwelltime = 10 / sample_rate;  %%% dwell time in usec

image_t = 0:sampT:tmax;
Fsamp = interp1(t,F(1:length(t)),image_t);
figure
plot(Fsamp);


Fout =  poissrnd(Fsamp*dwelltime);

baseFout = baselineF*dwelltime;
deltaF = (Fout - baseFout)/baseFout;
figure
plot(deltaF);
% 
% drift = (1:length(deltaF))/length(deltaF);
% drift = 0.1*drift;
% deltaF=deltaF+drift;
% figure
% plot(deltaF);


 
    df_interp = interp1(image_t,deltaF,t);
%       figure
%       plot(interp1(image_t,deltaF,t));
%     hold on
%     plot(spike.*max(deltaF),'g');
    est(1,rep) = sum(df_interp(1000:1200));
    title('raw data');
 
     impulse_samp = interp1(t,impulse,image_t);
    impulse_samp = impulse_samp(length(impulse_samp):-1:1);
    
    matched_filt = conv(deltaF,impulse_samp);
    matched_filt = matched_filt(length(deltaF):length(matched_filt));

   
%     figure
%       plot(interp1(0:sampT:sampT*(length(matched_filt)-1),matched_filt,t));
%      hold on
%      plot(spike.*max(deltaF),'g');
   matched_filt_interp = interp1(0:sampT:sampT*(length(matched_filt)-1),matched_filt,t);
       est(5,rep) = sum(matched_filt_interp(1000:1050));
       
% 
%      figure
%    tic
%    plot(conv(interp1(0:sampT:(length(deltaF)-1)*sampT,deltaF,t),spike(length(spike):-1:1)))
%    toc 
   
for j =1:3;
    SNR = 10^(j+1);
    

    impulse_samp = interp1(t,impulse,image_t);
    
    
    impulse_snip = impulse_samp(1:min(find(image_t>5*tau)));
    
    df_est = dF;
    tau_est=tau;
    impulse  = df_est * exp(-image_t/tau_est);
    
    pulsewidth = 0.001;
    t_samp = image_t-image_t/2;
    s = exp(-t_samp.^2/pulsewidth.^2);
    s = s./sum(s);
    S = (fft(s));
    figure
    plot(abs(S))
    figure
    plot(ifft(1./S))
    
    plot(conv(ifft(1./S),s));
    
    S2= S.^2;

    H = fft(impulse);

    SNR = 10;
    N = mean(S2.*H.*conj(H))/SNR;
    
    G = conj(H).*S2./((abs(H).^2).*S2 + N);
    figure
    plot(abs(G))

    figure
    plot((ifft((G))));

 Fdecon = ifft(G.*fft(deltaF));

figure
plot(image_t,deltaF/df_est);
hold on
plot(image_t,Fdecon,'g');
plot(image_t,Fdecon>0.5,'r')

    impulse_snip = impulse_samp(1:min(find(image_t>5*tau)));

    tau0=tau;
   
   f = deltaF;
    s = round(Fdecon);
    s(s<0)=0;
    p = impulse;
    
    for i = 1:10
        f_est = conv(s,p);
        f_est = f_est(1:length(s));
        figure
        plot(f);
        hold on
        plot(f_est,'g');
        err = f - f_est;
        p_norm = sum(p.^2);
        delta = conv(p(length(p):-1:1)/p_norm,err);
        delta = delta(length(p):length(delta));
        s = s+delta;
     %   for t = shuffle(1:length(s));
            
        figure
        plot(s)
        hold on
        s = round(s);
        s(s<0)=0;
        plot(s,'g');
    end
    
        

% 

    figure
    plot(interp1(0:sampT:(length(Fdecon)-1)*sampT,Fdecon,t));
    hold on
        plot(interp1(0:sampT:(length(Fdecon)-1)*sampT,deltaF./max(deltaF),t),'r');
    plot(spike/4,'g');
    
%     F = (fft(Fdecon));
%     Fmeas(rep,j+1) =abs(F(25))/abs(F(1));
%     angmeas(rep,j+1) =angle(F(25));
%     figure
%     plot(abs(fft(Fdecon)));
decon_interp = interp1(0:sampT:(length(Fdecon)-1)*sampT,Fdecon,t);
est(j+1,rep) = sum(decon_interp(1000:1050));

%    figure
%    tic
%    plot(conv(interp1(0:sampT:(length(Fdecon)-1)*sampT,Fdecon,t),spike(length(spike):-1:1)))
%    toc
title(sprintf('SNR = %f',SNR));
end



figure
plot(N,est(1,:),'o');
title('averaging')
figure
plot(N,est(2,:),'o');
title('snr 100');
figure
plot(N,est(3,:),'o');
title('snr 1000');
figure
plot(N,est(4,:),'o');
title('snr 10000');
figure
plot(N,est(5,:),'o');
title('matched filter');

for i = 1:5
    corr(N',est(i,:)')
end

