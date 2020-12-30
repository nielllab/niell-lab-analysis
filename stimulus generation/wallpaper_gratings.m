ncyc = 32;
sz = 512;
period = sz/ncyc;
mag = 4

[x y] = meshgrid(1:sz,1:sz);
im = 0.5 + 0.5*cos(2*pi*x/period);
figure
imshow(im)

[f p] = uiputfile('*.tif')
imwrite(imresize(im,mag),fullfile(p,f));


