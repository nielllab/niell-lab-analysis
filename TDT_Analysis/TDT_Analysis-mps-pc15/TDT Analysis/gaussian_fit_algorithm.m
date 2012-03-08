%%% fit firing rate from sweeping bars to a gaussian
%%% written by cmn, extracted from single unit analysis code 02/2008

spont = R(cell_n,size(bar_orients,2)+1);   % spontaneous rate

%%% range to fit over ... only use spikes in here
fit_range = 0:0.1:duration +0.5;
fit_int = fit_range(2)-fit_range(1);
Spike_Timing = Spike_Timing(find((Spike_Timing>min(fit_range))&(Spike_Timing<max(fit_range))));


%%% create histogram of spike timing
obs = hist(Spike_Timing, fit_range)/(fit_int*numtrials);



%%% linearly interpolated to avoid overfitting
fit_range_interp = 0:.02:Bar_Time(2)+0.5;
obs_interp = interp1(fit_range,obs,fit_range_interp);
obs =obs_interp - spont;  %% subtract off spontaneous rate
baseline(cell_n,orientation+1) = spont;

%%% coefficients are amplitude, location, width
peak_guess = median(fit_range_interp(find(obs> 0.75*max(obs))))
fit_coeff = nlinfit(fit_range_interp,obs,@rf_fit_nobaseline,[ max(obs) peak_guess duration/10])
if isnan(fit_coeff(1))
    fit_coeff(1)=0;
end

amp(cell_n,orientation+1) = fit_coeff(1);
if isnan(fit_coeff(2))
    amp(cell_n,orientation+1)=0;
end
if orientation<4
    x0(cell_n,orientation+1) = fit_coeff(3) + bar_width_time/2;
else
    x0(cell_n,orientation+1) = stim_duration-fit_coeff(3) - bar_width_time/2;
end
width(cell_n,orientation+1) = abs(fit_coeff(3));

%%% look for aberrant results, and set all values to zero
if abs(fit_coeff(3)>(Bar_Time(2)-Bar_Time(1))) | (fit_coeff(1)<0) | fit_coeff(3)<.045 | sum(isnan(fit_coeff))>0
    fit_coeff(2)=0;
    amp(cell_n,orientation+1)=0;
    x0(cell_n,orientation+1) = 0;
    width(cell_n,orientation+1) = 0;
    baseline(cell_n,orientation+1) = mean(obs);
    fit_coeff(1)=mean(obs);
end


hold on
plot(fit_range_interp, rf_fit_nobaseline(fit_coeff,fit_range_interp)+spont,'g','LineWidth',1.5);
