function [p g] = fitLGNrf(z);
%%% fit gaussian rf cmn 2012
%%% peak has to be reasonably high above background for this to work
%%%% fits to gauss2d
%%% p = [A B x0 y0 sigx sigy], where ...
%%% A = amplitude of baseline
%%% B = baseline
%%% x0,y0 = center coordinates
%%% sigx, sigy = width of gaussian
%%%%
%%% g = STA from fit



%%% estimate initial values
baseline = 0;
[amp ind] = max(abs(z(:)));
amp=z(ind);
[max_x max_y] = ind2sub(size(z),ind);;

peakarea = sum(abs(z(:))>0.4*abs(amp));
sigx = sqrt(peakarea/pi);
sigy= sigx;

[y x ] = meshgrid(1:size(z,2),1:size(z,1));

%%% vector of xy coordinates
xy = [x(:) y(:)]; 

%%% put initial estimates into parameter vector
p0(1) = amp;
p0(2) = baseline;
p0(3) = max_x;
p0(4) = max_y;
p0(5)= sigx;
p0(6)=sigy;

%%% perform fit
p = nlinfit(xy,z(:),@gauss2d,p0)

%%% calculate estimated RF from fit parameters
fit= gauss2d(p,xy);
g= reshape(fit,size(z,1),size(z,2));
