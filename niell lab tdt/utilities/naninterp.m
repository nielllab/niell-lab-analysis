function out = naninterp(data);

use = find(~isnan(data));
out = interp1(use,data(use),1:length(data));
