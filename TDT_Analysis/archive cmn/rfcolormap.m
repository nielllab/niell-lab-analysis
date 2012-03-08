function map=rfcolormap(gamma)
%%%%  intersting colormaps for STA rfs
caxis([-60 60])
baseline =1;
gamma =2;
map =zeros(64,3)+baseline;

gammamap= ((1:32)/32).^gamma;
gammamap2 = ((1:32)/32).^(2*gamma);

map(33:64,1) = baseline + (1-baseline)*gammamap;
map(33:64,2) = baseline - baseline*gammamap;
map(33:64,3) = baseline - baseline*gammamap;

map(32:-1:1,3) = baseline + (1-baseline)*gammamap;
map(32:-1:1,2) =baseline - baseline*gammamap;
map(32:-1:1,1) = baseline - baseline*gammamap;
colormap(map)
axis square

baseline =0;
gamma = 2;
gammamap= ((1:32)/32).^gamma;
gammamap2 = ((1:32)/32).^(2*gamma);
map =zeros(64,3)+baseline;

map(33:64,1) = baseline + (1-baseline)*gammamap;
map(33:64,2) = baseline + 0.5*(1-baseline)*gammamap2;
map(33:64,3) = baseline + 0.5*(1-baseline)*gammamap2;

map(32:-1:1,3) = baseline + (1-baseline)*gammamap;
map(32:-1:1,2) =baseline + 0.5*(1-baseline)*gammamap2;
map(32:-1:1,1) = baseline + 0.5*(1-baseline)*gammamap2;
colormap(map)