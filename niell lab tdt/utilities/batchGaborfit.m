






afiles = {'NR2A\KO\9_14_15\analysis_9_14_15.mat',...
    'NR2A\KO\9_16_15\analysis_2.mat'}

apath = 'E:\jennifer\Jennifer_analysis\';
       % for i = 1:length(afiles)
        for i = 1:length(afiles)
            afiles{i}
            fit2dgabor([apath afiles{i}]);
        end