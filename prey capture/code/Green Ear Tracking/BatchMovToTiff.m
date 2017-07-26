dbstop if error

% m_file = {   'E:\Dropbox\mucsimol_V1\Trial videos for Tristans tracking\Iphone5\live_2_LT_1hr_postCNO.mov',...
%             'E:\Dropbox\mucsimol_V1\Trial videos for Tristans tracking\Iphone5\live_3_LT_1hr_postCNO.mov'};                       

%m_file = {   'F:\Dropbox\behavior videos\iPhone\5_25_17\ALGS0101.mov'};                       
 m_file = {   'F:\GreenEarTracking\iPhone\5_25_17\ALGS0101.mov '}      
%cd F:\GreenEarTracking\iPhone\5_25_17

for i = 1:length(m_file);
   i
   m_file(i)
     MovToTiff(m_file{i}) ;
end%
    