function [meandata errdata N]=layerAgeNxNy(data1,data2,ageList,layer,inh,used,label,titlestr);
ageList=ageList';

    for group = 1:5
        
        if group ==1
            figure
            
            for age=1:2
              
            colors = 'bmrg'
            uselist = (ageList==age & (layer<=6) & ~inh & used');
            
                    % jitterAmount = 0.05;
                     jitterValuesX = 2*(rand(size(data1(uselist)))-0.5)*0.02;   % +/-jitterAmount max
                     jitterValuesY = 2*(rand(size(data2(uselist)))-0.5)*0.02;   % +/-jitterAmount max
                     
                     scatter(data1(uselist)+jitterValuesX, data2(uselist)+ jitterValuesY,[colors(age) 'o']);hold on
%                    % plot(nanmedian(data1(uselist)),nanmedian(data2(uselist)), '--rs','MarkerSize',10,'MarkerEdgeColor',colors(age),'LineWidth',2);
                    
                     if sum(uselist)>2
                         
                     mx=nanmedian(data1(uselist));
                     sx= semedian(data1(uselist));
                     
                     my=nanmedian(data2(uselist));
                     sy=semedian(data2(uselist));
%                      
                     mediandata_x(group,age)=mx;
                     errdata_med_x(group,age)=sx;
                     
                     mediandata_y(group,age)=my;
                     errdata_med_y(group,age)=sy;
                     
                     plot(mx,my, '--rs','MarkerSize',10,'MarkerEdgeColor',colors(age),'LineWidth',2);
                     plot ([(mx-sx) (mx-sx)],my,'-s','MarkerSize',2,'MarkerEdgeColor',colors(age),'LineWidth',2);
                     plot ([(mx+sx) (mx+sx)],my,'-s','MarkerSize',2,'MarkerEdgeColor',colors(age),'LineWidth',2);
                     plot (mx,(my-sy),'-s','MarkerSize',2,'MarkerEdgeColor',colors(age),'LineWidth',2);
                     plot (mx,(my+sy),'-s','MarkerSize',2,'MarkerEdgeColor',colors(age),'LineWidth',2);
                     
                     
                     end
            end
                   
                    legend('EO1','Adult');
                    xlabel(label{1}); ylabel(label{2});
                    axis square
                    axis equal
                    plot([0 1],[0 1])
                    
%                     
%                     
%                     
                    
                    
                    
                    
%         elseif group ==2
%            figure
%            
%            for age=1:4
%            
%            uselist = (ageList==age & (layer==4) & ~inh & used');
%            
%                     % jitterAmount = 0.05;
%                      jitterValuesX = 2*(rand(size(data1(uselist)))-0.5)*.02;   % +/-jitterAmount max
%                      %jitterValuesY = 2*(rand(size(data2(uselist)))-0.5)*jitterAmount;   % +/-jitterAmount max
%                      scatter(data1(uselist)+jitterValuesX, data2(uselist),[colors(age) 'o']);hold on
%                     %  plot(data1(uselist),data2(uselist),[colors(age) 'o']);hold on
%                      plot(nanmedian(data1(uselist)),nanmedian(data2(uselist)), '--rs','MarkerSize',10,'MarkerEdgeColor',colors(age),'LineWidth',2);
%                     
%                      if sum(uselist)>2
%                      mx=nanmedian(data1(uselist));
%                      sx= semedian(data1(uselist));
%                      
%                      my=nanmedian(data2(uselist));
%                      sy=semedian(data2(uselist));
% %                      
%                      plot ([(mx-sx) (mx-sx)],my,'-s','MarkerSize',2,'MarkerEdgeColor',colors(age),'LineWidth',2);
%                      plot ([(mx+sx) (mx+sx)],my,'-s','MarkerSize',2,'MarkerEdgeColor',colors(age),'LineWidth',2);
%                      plot (mx,[(my-sy) (my-sy)],'-s','MarkerSize',2,'MarkerEdgeColor',colors(age),'LineWidth',2);
%                      plot (mx,[(my+sy) (my+sy)],'-s','MarkerSize',2,'MarkerEdgeColor',colors(age),'LineWidth',2);
%                      end
%             end
%                     title 'layer4'
%                     legend('EO1','Adult');
%                     xlabel(label{1}); ylabel(label{2});
%                     axis square
%                     axis equal
%                     plot([0 1],[0 1])
%         elseif group==3
%             figure
%             
%             for age=1:4
%             uselist = (ageList==age & (layer==5) & ~inh & used');
%           
%                     % jitterAmount = 0.05;
%                      jitterValuesX = 2*(rand(size(data1(uselist)))-0.5)*.02;   % +/-jitterAmount max
%                      %jitterValuesY = 2*(rand(size(data2(uselist)))-0.5)*jitterAmount;   % +/-jitterAmount max
%                      scatter(data1(uselist)+jitterValuesX, data2(uselist),[colors(age) 'o']);hold on
%                     %  plot(data1(uselist),data2(uselist),[colors(age) 'o']);hold on
%                      plot(nanmedian(data1(uselist)),nanmedian(data2(uselist)), '--rs','MarkerSize',10,'MarkerEdgeColor',colors(age),'LineWidth',2);
%                     
%                      if sum(uselist)>2
%                      mx=nanmedian(data1(uselist));
%                      sx= semedian(data1(uselist));
%                      
%                      my=nanmedian(data2(uselist));
%                      sy=semedian(data2(uselist));
% %                      
%                      plot ([(mx-sx) (mx-sx)],my,'-s','MarkerSize',2,'MarkerEdgeColor',colors(age),'LineWidth',2);
%                      plot ([(mx+sx) (mx+sx)],my,'-s','MarkerSize',2,'MarkerEdgeColor',colors(age),'LineWidth',2);
%                      plot (mx,[(my-sy) (my-sy)],'-s','MarkerSize',2,'MarkerEdgeColor',colors(age),'LineWidth',2);
%                      plot (mx,[(my+sy) (my+sy)],'-s','MarkerSize',2,'MarkerEdgeColor',colors(age),'LineWidth',2);
%                      end 
%             end
%                     title 'layer5'
%                     legend('EO1','Adult');
%                     xlabel(label{1}); ylabel(label{2});
%                     axis square
%                     axis equal
%                     plot([0 1],[0 1])
% %         elseif group==4
% % %             figure
% %             
% %             for age=1:4
% %             uselist = (ageList==age & (layer==6) & ~inh &  used');
% %               
% %                     % jitterAmount = 0.05;
% %                      jitterValuesX = 2*(rand(size(data1(uselist)))-0.5)*.02;   % +/-jitterAmount max
% %                      %jitterValuesY = 2*(rand(size(data2(uselist)))-0.5)*jitterAmount;   % +/-jitterAmount max
% %                      scatter(data1(uselist)+jitterValuesX, data2(uselist),[colors(age) 'o']);hold on
%                     %  plot(data1(uselist),data2(uselist),[colors(age) 'o']);hold on
%                     
%                      plot(nanmedian(data1(uselist)),nanmedian(data2(uselist)), '--rs','MarkerSize',10,'MarkerEdgeColor',colors(age),'LineWidth',2);
%                     
%                      if sum(uselist)>2
%                      mx=nanmedian(data1(uselist));
%                      sx= semedian(data1(uselist));
%                      
%                      my=nanmedian(data2(uselist));
%                      sy=semedian(data2(uselist));
% %                      
%                      plot ([(mx-sx) (mx-sx)],my,'-s','MarkerSize',2,'MarkerEdgeColor',colors(age),'LineWidth',2);
%                      plot ([(mx+sx) (mx+sx)],my,'-s','MarkerSize',2,'MarkerEdgeColor',colors(age),'LineWidth',2);
%                      plot (mx,[(my-sy) (my-sy)],'-s','MarkerSize',2,'MarkerEdgeColor',colors(age),'LineWidth',2);
%                      plot (mx,[(my+sy) (my+sy)],'-s','MarkerSize',2,'MarkerEdgeColor',colors(age),'LineWidth',2);
%                      end
%                      
%             end
%             
%                     title 'layer6'
%                     legend('EO1','Adult');
%                     xlabel(label{1}); ylabel(label{2});
%                     axis square
%                     axis equal
%                     plot([0 1],[0 1])
% %         elseif group==5
% %             uselist = (ageList==age & inh & used');
         elseif group==2
            figure
            for age=1:4
            
            uselist = (ageList==age & (layer<=6)& ~inh & used');
                  
                    % jitterAmount = 0.05;
                     jitterValuesX = 2*(rand(size(data1(uselist)))-0.5)*0.02;   % +/-jitterAmount max
                     jitterValuesY = 2*(rand(size(data2(uselist)))-0.7)*0.1;   % +/-jitterAmount max
                     scatter(data1(uselist), data2(uselist),[colors(age) 'o']);hold on
                    %  plot(data1(uselist),data2(uselist),[colors(age) 'o']);hold on
    
                     plot(nanmedian(data1(uselist)),nanmedian(data2(uselist)), '--rs','MarkerSize',10,'MarkerEdgeColor',colors(age),'LineWidth',2);
                    
                     if sum(uselist)>2
                     mx=nanmedian(data1(uselist));
                     sx= semedian(data1(uselist));
                     
                     my=nanmedian(data2(uselist));
                     sy=semedian(data2(uselist));
%                      
                     plot ([(mx-sx) (mx-sx)],my,'-s','MarkerSize',2,'MarkerEdgeColor',colors(age),'LineWidth',2);
                     plot ([(mx+sx) (mx+sx)],my,'-s','MarkerSize',2,'MarkerEdgeColor',colors(age),'LineWidth',2);
                     plot (mx,[(my-sy) (my-sy)],'-s','MarkerSize',2,'MarkerEdgeColor',colors(age),'LineWidth',2);
                     plot (mx,[(my+sy) (my+sy)],'-s','MarkerSize',2,'MarkerEdgeColor',colors(age),'LineWidth',2);
                     end
            end
                    title 'All layers'
                    legend('EO1','Adult');
                    xlabel(label{1}); ylabel(label{2});
                    axis square
                    axis equal
                    plot([0 0.75],[0 0.75])
    end
    
  
            end
        end

             
             
      
    


      

   




    
