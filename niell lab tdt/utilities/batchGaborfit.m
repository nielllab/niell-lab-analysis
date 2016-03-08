







afiles = {'01_08_16\analysis_010816.mat',...
    '01_29_16\analysis_012916'}

% afiles = {'NR2A\KO\9_14_15\analysis_9_14_15.mat',...
%     'NR2A\KO\9_16_15\analysis_2.mat'}


apath = 'D:\Angie_analysis\DOI_experiments\';
       
        for i = 1:length(afiles)
            afiles{i}
            fit2dgabor([apath afiles{i}]);
        end