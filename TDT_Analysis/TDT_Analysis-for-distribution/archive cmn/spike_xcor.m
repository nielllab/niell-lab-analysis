function spike_xcor(event_times_all,idx_all,cells)
%%%  calculate x correlation
%%% this gives times of cell 2 relative to cell 1;
%%% i.e. if cell 2 fires 5msec after cell 1, see peak at +5

% ch1 = 5 ; cl1=3;
% ch2=9; cl2=3;
for c =1:size(cells,1);
    for c2 = c:size(cells,1);
        ch1 = cells(c,1);
        ch2= cells(c2,1);
        cl1 =cells(c,2);
        cl2 = cells(c2,2);

        t1all= squeeze(event_times_all(ch1,find(idx_all(ch1,:) == cl1)));
        t2all= squeeze(event_times_all(ch2,find(idx_all(ch2,:) == cl2)));
        b1=floor(t1all/10^5);
        b2 = floor(t2all/10^5);
        t1all = t1all-b1*10^5;
        t2all = t2all-b2*10^5;

        xc_all=0;

        for block=0:max(b1);
            block
            t1 = t1all(b1==block);
            t2 = t2all(b2==block);
            size(t1)
            size(t2)

            if length(t1)==0
                t1=zeros(size(t2));
            elseif length(t2)==0
                t2 = zeros(size(t1));
            end
            maxt = max(max(t1),max(t2));


            tic
            hist1 = hist(t1,0:0.002:maxt);
            hist2 = hist(t2,0:0.002:maxt);
            toc

            clear xc
            mint = -25;
            maxt = 25;
            tic
            for t= mint:maxt;
                if t<0
                    % xc(t-mint+1) = sum(circshift(hist1,[0 t]).*hist2);
                    xc(t-mint+1) = sum(hist1(abs(t)+1 : length(hist1)).*hist2(1:length(hist2)-abs(t)));
                elseif t==0
                    xc(t-mint+1) = sum(hist1.*hist2);
                elseif t>0
                    xc(t-mint+1) = sum(hist2(abs(t)+1 : length(hist1)).*hist1(1:length(hist1)-abs(t)));
                end
            end
            toc
            xc_all = xc_all+xc;
        end

        t0 = 1-mint;

        figure
        plot(mint:maxt,xc_all,'o');
        hold on
        plot(mint:maxt,xc_all);
        axis([mint maxt 0 max(xc_all)])
        title(sprintf('%d %d - %d %d',ch1,cl1,ch2,cl2));

    end
end



%     %dt = zeros(size(t1,2),size(t2,2));
%
%      xc = zeros(size(mint:maxt));
%
%      tic
%      for i = 1:size(t1,2);
%          dt = t2-t1(i);
%          xc = xc+hist(dt*1000,mint:maxt);
%      end
%      toc
%
%      tic
%      xc=hist(dt(:)*1000,mint:maxt);
%      toc
%      clear dt
%      pack