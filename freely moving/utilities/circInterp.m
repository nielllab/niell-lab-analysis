function y = circInterp(x0,y0,x);
%%% perform interpolation on circular variables (from -pi to pi)
%%% uses cumulative sum to unwrap values

y0 = mod(y0+pi,2*pi)-pi; %%% get centered between -pi to pi. Some are larger (e.g. azimuth)
y0raw = y0;
nanList = isnan(y0);
if isnan(y0(1))
    y0(1) =0;
end

for i =2:length(y0);
    if isnan(y0(i))
        y0(i)= y0(i-1);
    end
end

dy0 = diff(y0);
dy0(dy0>pi)=dy0(dy0>pi)-2*pi;
dy0(dy0<-pi)=dy0(dy0<-pi)+2*pi;
y0unwrap = cumsum([y0(1) dy0]);

y0unwrap(nanList) = NaN;
% 
% figure; 
% plot(y0raw); hold on; plot(y0unwrap+0.1);
% use = (~isnan(y0(1:length(x))));
yUnwrap = interp1(x0,y0unwrap,x);

% yUnwrap = interp1(x0(use),y0unwrap(use),x(use));
y = mod(yUnwrap + pi,2*pi)-pi;

% figure
% plot(x0,y0raw);
% hold on;
% plot(x0,y0unwrap+0.1);
% plot(x,y+0.2);
% 
% figure
% plot(diff(y)); hold on; plot(diff(y0raw))


