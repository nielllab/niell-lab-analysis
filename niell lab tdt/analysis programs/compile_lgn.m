clear all
afile = {'C:\data\ephys matlab data\030912_LGN\recording3\analysis3.mat' ...
    'C:\data\ephys matlab data\031512_LGN\analysis2.mat'}
n=0;
for i = 1:length(afile);
    load(afile{i});
    cell_range = n+(1:length(drift));
    n=n+length(drift);
    if size(wn,2)==2
        wn_all(cell_range,:)=wn;
    else
        wn_all(cell_range,1)=wn;
    end
    drift_all(cell_range)=drift;
    mv_all(cell_range)=mv;
    fl_all(cell_range)=fl;
end

clear p

for i = 15:20
    sta = double(squeeze(wn(i).svd_xy(1,:,:)));    
    p{i} = fitLGNrf(sta);
end