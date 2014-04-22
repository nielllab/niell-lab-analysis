%%%Load STA data for all cells

STA = wn.sta;

%%%Dtermine time point with maximial response
[m ind] = max(abs(STA(:)-127));
[x y t_lag] = ind2sub(size(STA),ind)

STA1 = STA(:,:,t_lag)-128;
figure
colormap(redblue)
imagesc(STA1)

%%%crop image
width = 41;
img = maxsubregion(STA1,width,'neg');
stafig1 = figure;
colormap(redblue);
imagesc(img)
axis off;
saveas(stafig1,'Group2.png','png'); %come back and give unique, updating name

opts.errorbars = 'none';
opts.tilted = true;
opts.iso = false;
[x y] = meshgrid(1:width,1:width);
results = autoGaussianSurf(x,y,img-128,opts);

figure
subplot(1,2,1);imagesc(x(:),y(:),img);axis equal
subplot(1,2,2);imagesc(x(:),y(:),results.G);axis equal


%%%gabor fitting, applies gabor filtering to an image where I can determine
%%%input parameters

lambda  = 8;
theta   = 0;
psi     = [0 pi/2];
gamma   = 0.5;
bw      = .1; %%%likely SF estimates??
N       = 1; %%possible orientations??

img_in = im2double(imread('Group2.png'));
img_in(:,:,2:3) = [];   % discard redundant channels, make gray
img_out = zeros(size(img_in,1), size(img_in,2), N);
for n=1:N
    gb = gabor_fn(bw,gamma,psi(1),lambda,theta)...
        + 1i * gabor_fn(bw,gamma,psi(2),lambda,theta);
    % gb is the n-th gabor filter
    img_out(:,:,n) = imfilter(img_in, gb, 'symmetric');
    % filter output to the n-th channel
    theta = theta + 2*pi/N;
    % next orientation
end
figure(1);
imshow(img_in);
title('input image');
figure(2);
img_out_disp = sum(abs(img_out).^2, 3).^0.5;
% default superposition method, L2-norm
img_out_disp = img_out_disp./max(img_out_disp(:));
% normalize
imshow(img_out_disp);
title('gabor output, L-2 super-imposed, normalized');