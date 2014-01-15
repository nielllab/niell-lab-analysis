close all
for i = 1:length(ds)
   figure
    sta = wn_all(ds(i),1).sta_final
   if ~isnan(sta)
      
   imagesc(sta,1.0*[-max(abs(sta(:))) max(abs(sta(:)))])
    % imagesc(sta)
   colormap(gray)
   end
   
   title(sprintf('cell %d',ds(i)))
   axis equal
end


figds = [124 157 182];

for i = 1:length(figds); figure; imagesc(wn_all(figds(i),1).sta_final,[-0.1 0.1]); colormap(gray);end
