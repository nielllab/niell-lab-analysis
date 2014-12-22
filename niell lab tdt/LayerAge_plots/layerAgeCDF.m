function [meandata errdata N]=layerAgeCDF(data,ageList,layer,inh,used,label);
ageList=ageList';




    for group = 1:6
        
        if group ==1
                figure
                hold on
                colorlist='bmrg';
            for age=1:2
                uselist = (ageList==age & (layer==2  | layer==3) & ~inh & used);
            if sum(uselist)>2
                
            CDF = data(uselist);
            [f,x,Flo,Fup]= ecdf(CDF);
            stairs(x,f,'lineWidth',4,'color',colorlist(age));
            hold on
            plot(x,Flo,'color',colorlist(age));plot(x,Fup,'color',colorlist(age)); 
            axis xy
            xlabel(label)
            title 'layer 2/3'
            bothdata{group,age}= data(uselist);           
            end
        
            end
         
        elseif group ==2
                figure
                hold on
                colorlist='bmrg';
            for age=1:4
                uselist = (ageList==age & (layer==4) & ~inh & used);
                
                if sum(uselist)>2
                    
            CDF = data(uselist);
            [f,x,Flo,Fup]= ecdf(CDF);
            stairs(x,f,'lineWidth',4,'color',colorlist(age));
            hold on
            plot(x,Flo,'color',colorlist(age));plot(x,Fup,'color',colorlist(age)); 
            axis xy
            xlabel(label)
            title 'layer 4'
            bothdata{group,age}= data(uselist);
                end
            end
%             
            elseif group ==3
                figure
                hold on
                colorlist='bmrg';
            for age=1:4
                uselist = (ageList==age & (layer==5) & ~inh & used);
                if sum(uselist)>2
            CDF = data(uselist);
            [f,x,Flo,Fup]= ecdf(CDF);
            stairs(x,f,'lineWidth',4,'color',colorlist(age));
            hold on
            plot(x,Flo,'color',colorlist(age));plot(x,Fup,'color',colorlist(age));   
            axis xy
            xlabel(label)
            title 'layer 5'
            bothdata{group,age}= data(uselist);
                end
            end
%             
            elseif group ==4
                figure
                hold on
                colorlist='bmrg';
            for age=1:4
                uselist = (ageList==age & (layer==6) & ~inh & used);
                if sum(uselist)>2
            CDF = data(uselist);
            [f,x,Flo,Fup]= ecdf(CDF);
            stairs(x,f,'lineWidth',4,'color',colorlist(age));
            hold on
            plot(x,Flo,'color',colorlist(age));plot(x,Fup,'color',colorlist(age)); 
            axis xy
            xlabel(label)
            title 'layer 6'
            bothdata{group,age}= data(uselist);
                end
            end
%             
%             elseif group ==5
%                 figure
%                 hold on
%                 colorlist='bg';
%             for age=1:2
%                 uselist = (ageList==age & inh & used);  
%                 if sum(uselist)>0
%                 
%             CDF = data(uselist);
%             [f,x,Flo,Fup]= ecdf(CDF);
%             stairs(x,f,'lineWidth',4,'color',colorlist(age));
%             hold on
%             plot(x,Flo,'color',colorlist(age));plot(x,Fup,'color',colorlist(age)); 
%             axis xy
%             xlabel(label)
%             title 'inhibitory';
%            
%                 end
%             end
            
            elseif group ==5
                figure
                hold on
                colorlist='bmrg';
            for age=1:4
                uselist = (ageList==age & inh & used);
            
            if sum(uselist)>2
            CDF = data(uselist);
            [f,x,Flo,Fup]= ecdf(CDF);
            stairs(x,f,'lineWidth',4,'color',colorlist(age));
            hold on
            plot(x,Flo,'color',colorlist(age));plot(x,Fup,'color',colorlist(age)); 
            axis xy
            %xlim([0 2]);
            xlabel(label)
            title 'total population'
            
            bothdata{group,age}= data(uselist);
            end
            end
            
            elseif group ==6
                figure
                hold on
                colorlist='bmrg';
            for age=1:4
                uselist = (ageList==age & (layer<=6) & ~inh & used);
            
            if sum(uselist)>2
            CDF = data(uselist);
            [f,x,Flo,Fup]= ecdf(CDF);
            stairs(x,f,'lineWidth',4,'color',colorlist(age));
            hold on
            plot(x,Flo,'color',colorlist(age));plot(x,Fup,'color',colorlist(age)); 
            axis xy
            %xlim([0 2]);
            xlabel(label)
            title 'total population'
            
            bothdata{group,age}= data(uselist);
            end
            end
            
              
            
               [h,p]= kstest2(bothdata{1,1},bothdata{1,2})
               [h1,p1]= kstest2(bothdata{2,1},bothdata{2,2})
               [h2,p2]= kstest2(bothdata{3,1},bothdata{4,2})
               [h3,p3]= kstest2(bothdata{4,1},bothdata{4,2})
               [h4,p4]= kstest2(bothdata{5,1},bothdata{5,2})
               
        end
    end         
        
%         elseif group==3
%             uselist = (ageList==age & (layer==5) & ~inh & used);
%         elseif group==4
%             uselist = (ageList==age & (layer==6) & ~inh & used);
%         elseif group==5
%             uselist = (ageList==age & inh & used);
%         elseif group == 6
%             uselist = (ageList==age & (layer<=6)& ~inh & used);
       
%         if sum(uselist)>0
%                 figure
%                 hold on
%                 colorlist='bg';
%                 
%             CDF = data(uselist);
%             [f,x,Flo,Fup]= ecdf(CDF);
%             stairs(x,f,'lineWidth',4,'color',colorlist(age));
%             hold on
%             plot(x,Flo,'color',colorlist(age));plot(x,Fup,'color',colorlist(age)); 
%         end
            
            
    

% figure
% barweb(meandata,errdata);
% ylabel(label);
% set(gca,'Xtick',1:6);
% set(gca,'Xticklabel',{'2/3','4','5','6','inh','mid'});
% legend('EO1','adult');






% hold on
% OSItuning_adult= ?;
% [f,x,Flo,Fup]= ecdf(OSItuning_adult);
% stairs(x,f,'lineWidth',4,'color','b');
% hold on
% plot(x,Flo,'b');plot(x,Fup,'b');
% hold on
% title 'OSI tuning width'



