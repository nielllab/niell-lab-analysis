function se = stderr(x);
se = nanstd(x)/sqrt(length(~isnan(x)));