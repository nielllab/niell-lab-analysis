
nii = load_nii('C:\data\allen\volume_LGd.nii');
x = nii.img>100;
p = patch(isosurface(x));
set(p,'FaceAlpha',.3,'FaceColor','b','EdgeColor','none')
material dull
lighting phong
camlight


