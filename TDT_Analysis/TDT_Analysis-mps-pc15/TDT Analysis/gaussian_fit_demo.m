%%% set up datga

x = -5:5;
y = 5+3.5*exp(-0.5*((x-2)/2).^2) +0.5*rand(1,11);


%%% estimate parameters for initial guess

baseline_est = min(y);
[peak_est x0_est]=max(y);
x0_est = x(x0_est);
peak_est=peak_est-baseline_est;
sigma_est = 2;  %%% sigma is the hardest to estimate - I just use a reasonable a priori value

fit_coeff = nlinfit(x,y,@gauss_fit,[ baseline_est peak_est x0_est sigma_est])

%%% parse out results
baseline = fit_coeff(1)
peak = fit_coeff(2)
x0=fit_coeff(3)
sigma_est=fit_coeff(4)

%%% plot raw data and fit
figure
plot(x,y)
hold on
plot(x,gauss_fit(fit_coeff,x),'g')