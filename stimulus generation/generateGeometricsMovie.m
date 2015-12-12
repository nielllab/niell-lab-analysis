%PRLP 10/29/2015 Niell Lab
%This code is used to generate simple unidirectional stimuli based on
%images made in Adobe Illustrator
%run after running generateGeometrics

cnt = 1;
for n = 1:size(ordstimlib,4)
    moviedata(:,:,cnt:cnt+29) = uint8(ordstimlib(:,:,:,n));
    moviedata(:,:,cnt+30:cnt+89) = uint8(127);
    cnt = cnt + 90;
end

dir = 'C:\Users\nlab\Documents\MATLAB';
nam = 'GeomStim';
save(fullfile(dir,nam),'moviedata','ordleg');


% 
%                     
% figure
% for j = 1:128
%     imshow(ordstimlib(:,:,j,6));
% end