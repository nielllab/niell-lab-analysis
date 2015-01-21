function [onset_hist resp err_hist baseline baseresp onset_bins] = doOnOffFlash(moviedata,sz_mov,eps,epSpikes,ny);

for yi = 1:ny;
    
    y = yi*2;
    tseries = single(squeeze(moviedata(y,:)));
    size_series = single(squeeze(sz_mov(y,:)));
    
    %%%%% histrograms relative to onset/offset
    
    onset = (tseries(1:end-1)==127 & diff(tseries)==0 );
    
    [baseline(1,yi,:) baseresp(1,yi)]= getFlashHist(epSpikes,eps,onset);
  
   
    sizes = [1 2 4 8 16 255];
    for rep = 1:2
        for sz = 1:length(sizes);
            n=0;
            ontime=0;
            
            if rep ==1
                   % onset = (tseries(1:end-2)==127 & tseries(2:end-1)==0 & tseries(3:end)==127 &size_series(2:end-1)==sizes(sz));
                   onset = (tseries(1:end-2)==127 & tseries(2:end-1)==0 &size_series(2:end-1)==sizes(sz));

            elseif rep ==2
                    %onset = (tseries(1:end-2)==127 & tseries(2:end-1)==255 & tseries(3:end)==127 &size_series(2:end-1)==sizes(sz));
                    onset = (tseries(1:end-2)==127 & tseries(2:end-1)==255  &size_series(2:end-1)==sizes(sz));
      end           
            [onset_hist(1,yi,rep,sz,:) resp(1,yi,rep,sz) err_hist(1,yi,rep,sz,:) onset_bins] = getFlashHist(epSpikes,eps,onset);                        
        end
    end   
   
end