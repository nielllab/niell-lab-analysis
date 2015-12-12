







afiles = {'NR2A\KO\11_2_15\11_2_5_analysis_2.mat',...
    'NR2A\KO\11_3_15\11_3_15_analysis_2.mat',...
    'NR1\4_28_15\4_28_15_analysis_2.mat',...
    'NR1\5_1_15\5_1_15_analysis_2.mat',...
    'NR1\11_9_15\11_9_15_analysis_2.mat',...
    'NR1\11_10_15\11_10_15_analysis_2.mat',...
    'NR1\11_12_15\rec1\11_12_15_rec1_analysis_2.mat',...
    'NR1\11_12_15\11_12_15_rec2_analysis_2.mat'}

% afiles = {'NR2A\KO\9_14_15\analysis_9_14_15.mat',...
%     'NR2A\KO\9_16_15\analysis_2.mat'}


apath = 'D:\Jen_analysis\';
       
        for i = 1:length(afiles)
            afiles{i}
            fit2dgabor([apath afiles{i}]);
        end