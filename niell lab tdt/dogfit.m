function [w_pref_dog w_bw_dog peak_value] = dogfit(drift_sftuning);

data = drift_sftuning([1 5 9 13 17 21 25]);

data(8)=0;
% if min(data)<0
%     data(8) = min(data);
% end
if abs(min(data))>max(data)  %%% more negative than positive
    data = -data;
    inverted=1
else
    inverted=0
end

data(data<-1)=0;


sf = [0 .01 .02 .04 .08 .16 .32 .64];
fit_sfs = [0 .005 .01 .015 .02 .03 .04 .06 .08 .12 .16 .24 .32 .48 .64];
fit_sfs = sf;
interpdata = interp1(sf,data,fit_sfs,'spline');
%     figure
%     plot(sf,data);
%     hold on
%     plot(fit_sfs,interpdata,'g');

[A1 w1 A2 w2] = fit_dog(interpdata,fit_sfs)

w_pref_dog = dog_peak([A1 w1 A2 w2])

sf_plot = [0 .01 .02 .04 .08 .16 .32];

% figure
% intp=1:.1:7;
% yfit = diff_of_gauss([A1 w1 A2 w2],.01*2.^(intp-2));
% plot(data(1:7),'o');
% hold on
% plot(data(1:7),'b-');
% %plot(intp,interp1(1:7,data(1:7),intp,'spline'),'b-');
% %plot(intp,yfit,'g--');
% hold on
% yfit = diff_of_gauss([A1 w1 A2 w2],sf_plot);
% plot(intp,interp1(1:7,yfit(1:7),intp,'spline'),'g--');


yfit = diff_of_gauss([A1 w1 A2 w2],fit_sfs)
[y m] = max(yfit)
if ~isreal(w_pref_dog) | m==1   %%peak at 0
    w_pref_dog=0;
end
peak_value = diff_of_gauss([A1 w1 A2 w2],w_pref_dog)
if (w1==0 | w2 ==0 | peak_value<=0.1)
    w_pref_dog=-1;
    peak_value=0;
    w_bw_dog=0;
else
    f0 = diff_of_gauss([A1 w1 A2 w2],0)
    if f0>(peak_value/2)
        w_bw_dog = 0;
    else
        yfit = diff_of_gauss([A1 w1 A2 w2],fit_sfs);

        yfit(yfit>peak_value*0.6)=0;
        %         figure
        %         plot(yfit);
        [y lo] = max(yfit(find(fit_sfs<w_pref_dog)));
        upper= find(fit_sfs>w_pref_dog);
        [y hi] = max(yfit(upper));
        fit_sfs(lo)
        fit_sfs(upper(hi))

        w_low = fzero(@(w_low) solve_half_dog(w_low,[A1 w1 A2 w2 peak_value]),w_pref_dog/4)
        if w_low<.005
            w_low=.005;
        end
        w_hi = fzero(@(w_hi) solve_half_dog(w_hi,[A1 w1 A2 w2 peak_value]),w_pref_dog*4)
        w_bw_dog = log2(w_hi/w_low)/2
    end
end

if inverted
    peak_value=-peak_value;
end

%#ok<NOPRT>