%%% simple script to make grating icons to use in paper figures

figure
[x y] = meshgrid(0:5:5*360,0:5:5*360);
thetas = [-90 -45 0 45 90];
for i = 1:5
    subplot(2,3,i);
    im = cosd(cosd(thetas(i))*x + sind(thetas(i))*y);
    imagesc(im);
    colormap gray;
    axis equal
    xlim([0 360])
    axis off
end

    