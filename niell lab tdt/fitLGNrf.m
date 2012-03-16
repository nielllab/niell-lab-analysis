function [p g] = fitLGNrf(z);
%%% fit gaussian rf cmn 2012
%%% peak has to be reasonably high above background
%%% p = [A B x0 y0 sigx sigy]
%%%% fits to gauss2d

baseline = 0
[amp ind] = max(abs(z(:)))
amp=z(ind);
[max_x max_y] = ind2sub(size(z),ind);;

peakarea = sum(abs(z(:))>0.6*abs(amp));
sigx = sqrt(peakarea/pi);
sigy= sigx;

[y x ] = meshgrid(1:size(z,2),1:size(z,1));

xy = [x(:) y(:)];
p0(1) = amp;
p0(2) = baseline;
p0(3) = max_x;
p0(4) = max_y;
p0(5)= sigx;
p0(6)=sigy;

% figure 
% subplot(2,1,1)
% imagesc(z)
% axis equal

% guess = gauss2d(p0,xy);
% figure
% imagesc(reshape(guess,size(z,1),size(z,2)))
% axis equal

p = nlinfit(xy,z(:),@gauss2d,p0)

fit= gauss2d(p,xy);
g= reshape(fit,size(z,1),size(z,2));

% subplot(2,1,2);
% imagesc(g)
% axis equal
