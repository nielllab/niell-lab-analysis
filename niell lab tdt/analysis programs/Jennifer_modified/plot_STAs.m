    colormap redblue
    
    ss=find(exp>0.55)
    
    for j=1:length(all_img)
    if ismember(j,ss)
    figure
    if ~isempty(all_img(1,j))
    colormap redblue
    imagesc(all_img{1,j}',[-40 40]); axis equal
   
    end
    end
    end
    
    
    for j=1:length(all_fit);
        if ismember(j,ss);
    figure
    if ~isempty(all_fit(j));
    colormap redblue
    imagesc(all_fit{j}',[-40 40]); axis equal
    end
        end
    end
    
    for w = 1:length(wn)
    STA = wn(w).sta;
    
    %%%Dtermine time point with maximial response
    [m ind] = max(abs(STA(:)-127));
    [x y t_lag] = ind2sub(size(STA),ind);
    
    STA1{w} = STA(:,:,t_lag)-128;

 figure
 imagesc(STA1{1,w}',[-64 64]); axis equal
    end