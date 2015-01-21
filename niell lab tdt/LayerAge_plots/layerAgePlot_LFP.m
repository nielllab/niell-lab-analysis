function [meandata errdata N mediandata]=layerAgePlot_LFP(data,data2,ageList,layer, used,label,freqrange, freqrange2);
ageList=ageList';


colorlist='bgrm';

for age=1:4
    for group = 1:5
        if group ==1
            uselist = (ageList==age & (layer<=3)   & used);
        elseif group ==2
            uselist = (ageList==age & (layer==4)   & used);
        elseif group==3
            uselist = (ageList==age & (layer<=5)  &  used);
        elseif group==4
            uselist = (ageList==age & (layer==6)  &  used);
        elseif group==5
             uselist = (ageList==age & (layer<=6)  & used);
      
        end

    
        if sum(uselist)>2 
             N(group,age) = squeeze(sum(~isnan(data(uselist))));
             med_r= squeeze(nanmedian(data(uselist,2,:)))';
             med_s= squeeze(nanmedian(data(uselist,1,:)))';
             % mediandata(group,age)=med
             s_r = semedian(data(uselist,2,:));
             s_s = semedian(data(uselist,1,:));
             %errdata_med(group,age)=s';
%              figure
%              shadedErrorBar(freqrange,med_r,s_r,'k');hold on
%              shadedErrorBar(freqrange,med_s,s_s,'r')

[peakpow idxpeak]=max(med_r);
amp_r(group,age)=peakpow;
s_peakpow_r(group,age)=s_r(idxpeak);

amp_s(group,age)=med_s(idxpeak);
s_peakpow_s(group,age)=s_s(idxpeak);

peakfreq(group,age)=freqrange(idxpeak);

%%determine std error median around peak freq
            
           l_r=[]; l_s=[];
           l_r=data(uselist,2,:);
           l_s=data(uselist,1,:);
           for j=1:length(l_r(:,1));
            power_r=l_r(j,:);
            power_s=l_s(j,:);
%             figure
%             plot(freqrange,power)
            %y=squeeze(data(uselist,2,:))'
            [peakpow_r,idx_r]=max(power_r);
            p_r(j)=freqrange(idx_r);
            power_s_rpeak(j)=power_s(idx_r);
            ratio_amp(j)=amp_r(group,age)/power_s_rpeak(j);
            end
   
       
%    
 s_peakfreq(group,age)=semedian(p_r);
 s_diff_A(group,age)=semedian(ratio_amp);
 
  %%low freq band
  
%    N(group,age) = squeeze(sum(~isnan(data2(uselist))));
%              med_rlow= squeeze(nanmedian(data2(uselist,2,:)))';
%              med_slow= squeeze(nanmedian(data2(uselist,1,:)))';
%              % mediandata(group,age)=med
%              s_r_low = semedian(data2(uselist,2,:));
%              s_s_low = semedian(data2(uselist,1,:));
%              %errdata_med(group,age)=s';
% %              figure
% %              shadedErrorBar(freqrange2,med_rlow,s_r_low,'k');hold on
% %              shadedErrorBar(freqrange2,med_slow,s_s_low,'r')
% 
% [peakpowlow idxpeaklow]=max(med_slow);
% amp_slow(group,age)=peakpowlow;
% s_peakpow_slow(group,age)=s_s_low(idxpeaklow);
% 
% amp_rlow(group,age)=med_rlow(idxpeaklow);
% s_peakpow_rlow(group,age)=s_r_low(idxpeaklow);
% 
% peakfreq_low(group,age)=freqrange2(idxpeaklow);
% 
% %%determine std error median aroiund peak freq
%             
%            l_r_low=[]; l_s_low=[];
%            l_r_low=data2(uselist,2,:);
%            l_s_low=data2(uselist,1,:);
%            for j=1:length(l_r_low(:,1));
%             power_r_low=l_r_low(j,:);
%             power_s_low=l_s_low(j,:);
% %             figure
% %             plot(freqrange,power)
%             %y=squeeze(data(uselist,2,:))'
%             [peakpow_s_low,idx_s_low]=max(power_s_low);
%             p_s_low(j)=freqrange2(idx_s_low);
%             power_r_speak_low(j)=power_r_low(idx_s_low);
%             ratio_amp_low(j)=amp_slow(group,age)/power_r_speak_low(j);
%            end
%             s_peakfreq_low(group,age)=semedian(p_s_low);
%             s_diff_A_low(group,age)=semedian(ratio_amp_low);
%             ratio_A_low=amp_slow/amp_rlow;
         end
        end
%    
 
            end
  amp_r(:,4)=2.3*amp_r(:,4)
 ratio=(amp_r/amp_s)
        %amp_r(:,4)=1.3*amp_r(:,4);

        %gamma  
figure
barweb(peakfreq,s_peakfreq)

pf_EA=[peakfreq(:,1) peakfreq(:,4)]
pf_EA_s= [s_peakfreq(:,1) s_peakfreq(:,4)]

figure
barweb(pf_EA,pf_EA_s)

figure
errorbar(1:4,peakfreq(2,:),s_peakfreq(2,:),'K');hold on
errorbar(1:4,peakfreq(3,:),s_peakfreq(3,:),'g');hold on
errorbar(1:4,peakfreq(4,:),s_peakfreq(4,:),'r');hold on

         %ratio
d_EA=[ratio(:,1) ratio(:,4)]
ds_EA=[s_diff_A(:,1) s_diff_A(:,4)]
figure
barweb(d_EA,ds_EA)

figure
errorbar(1:4,ratio(1,:),s_diff_A(1,:),'K');hold on
errorbar(1:4,ratio(2,:),s_diff_A(2,:),'g');hold on
errorbar(1:4,ratio(4,:),s_diff_A(4,:),'r');hold on

figure
errorbar(1:4,amp_r(1,:),s_peakpow_r(1,:),'K');hold on
errorbar(1:4,amp_r(2,:),s_peakpow_r(2,:),'g');hold on
errorbar(1:4,amp_r(4,:),s_peakpow_r(4,:),'r');hold on
%alpha
pf_EA_low=[peakfreq_low(:,1) peakfreq_low(:,4)]
pf_EA_slow= [s_peakfreq_low(:,1) s_peakfreq_low(:,4)]

figure
barweb(pf_EA_low, pf_EA_slow)

figure
errorbar(1:4,peakfreq_low(5,:),s_peakfreq_low(5,:),'K');hold on

amp_slow_EA=[amp_slow(:,1) amp_slow(:,4)]
amp_slow_EA_s=[s_peakpow_slow(:,1) s_peakpow_slow(:,4)]
figure
barweb(amp_slow_EA,amp_slow_EA_s)

figure
errorbar(1:4,amp_slow(4,:),s_peakpow_slow(4,:),'K');hold on

ratio_A_low_EA=[ratio_A_low(:,1) ratio(:,4)]
s_diff_A_low_EA=[s_diff_A_low(:,1) s_diff_A_low(:,4)]
figure
barweb(ratio_A_low_EA,s_diff_A_low_EA)

figure
errorbar(1:4,ratio_A_low(1,:),s_diff_A_low(1,:),'K');hold on

% % 
% amp_rlow(group,age)=peakpowlow;
% s_peakpow_rlow(group,age)=s_r_low(idxpeaklow);
% 
% amp_slow(group,age)=med_slow(idxpeaklow);
% s_peakpow_slow(group,age)=s_s_low(idxpeaklow);
% 
% peakfreq_low(group,age)=freqrange2(idxpeaklow); 

 amp_1=amp_r(:,1);
 amp_A=amp_r(:,4);
 amp_r_1=[amp_1 amp_A];
 amp_1_s=s_peakpow_r(:,1);
 amp_A_s=s_peakpow_r(:,4);
 amp_A_s_1=[amp_1_s amp_A_s]
 
 figure
 barweb(amp_r_1,amp_A_s_1);
 barweb(amp_r,amp_A_s);
 ylabel(label)
 xlabel('median')
 
 figure
 barweb(amp_s,s_peakpow_s);
 ylabel(label)
 xlabel('median')
 
 figure
 barweb(diff_A,s_diff_A);
 ylabel(label)
 xlabel('median')
 end
% % % 

% [p,table,stat]=kruskalwallis(alldata)
%   ranksum(alldata{1,1},alldata{2,1}) % signfiicance layer2/3
%   ranksum(alldata{3,1},alldata{3,2}) % signfiicance layer4
%   [p h]=kstest2(alldata{6,1},alldata{6,2});
%  ranksum(alldata{3,1},alldata{3,2}) % sign. layer5
%  ranksum(alldata{4,1},alldata{4,2}) % sign layer6
%  ranksum(alldata{5,1},alldata{5,2}) % sign inhibtory
 %ranksum(bothdata{6,1},bothdata{6,2}) % all excitaory

 
 







    
