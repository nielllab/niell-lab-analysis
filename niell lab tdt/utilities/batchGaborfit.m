







afiles = {'1_18_16_analysis_2A.mat',... 
    }

% afiles = {'NR2A\KO\9_14_15\analysis_9_14_15.mat',...
%     'NR2A\KO\9_16_15\analysis_2.mat'}


apath = 'D:\Jen_analysis\NR5A_Pinping\1_18_16\';
       
        for i = 1:length(afiles)
            afiles{i}
            fit2dgabor([apath afiles{i}]);
        end