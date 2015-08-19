






afiles = {'NR2A\KO\7_29_15\analysis_2.mat',...
    'NR5A_Pinping\6_25_15\analysis_2.mat',...
    'NR2A\KO\7_17_15\7_17_15\analysis.mat'}

apath = 'D:\Jen_analysis\';
       % for i = 1:length(afiles)
        for i = 1:length(afiles)
            afiles{i}
            fit2dgabor([apath afiles{i}]);
        end