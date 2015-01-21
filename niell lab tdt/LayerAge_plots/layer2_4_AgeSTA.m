function [meandata errdata N]=layerAgeActivity(data1,data2,ageList,layer,inh,used,label,titlestr);
ageList=ageList';

    for group = 1:5
        
        if group ==1
            figure
            
            for age=1:2
            colors = 'bg'
            uselist = (ageList==age & (layer==2 | layer==3) & ~inh & used');
            
                    
                    % jitterAmount = 0.05;
                     jitterValuesX = 2*(rand(size(data1(uselist)))-0.5)*.02;   % +/-jitterAmount max
                     %jitterValuesY = 2*(rand(size(data2(uselist)))-0.5)*jitterAmount;   % +/-jitterAmount max
                     scatter(data1(uselist)+jitterValuesX, data2(uselist),[colors(age) 'o']);hold on
                    %  plot(data1(uselist),data2(uselist),[colors(age) 'o']);hold on
            end
                    title 'layer2/3'
                    legend('EO1','Adult');
                    xlabel(label{1}); ylabel(label{2});
                    axis square
                    axis equal
                    plot([0 1],[0 1])
                    
        elseif group ==2
           figure
           
           for age=1:2
           
           uselist = (ageList==age & (layer==4) & ~inh & used');
           
                    % jitterAmount = 0.05;
                     jitterValuesX = 2*(rand(size(data1(uselist)))-0.5)*.02;   % +/-jitterAmount max
                     %jitterValuesY = 2*(rand(size(data2(uselist)))-0.5)*jitterAmount;   % +/-jitterAmount max
                     scatter(data1(uselist)+jitterValuesX, data2(uselist),[colors(age) 'o']);hold on
                    %  plot(data1(uselist),data2(uselist),[colors(age) 'o']);hold on
            end
                    title 'layer4'
                    legend('EO1','Adult');
                    xlabel(label{1}); ylabel(label{2});
                    axis square
                    axis equal
                    plot([0 1],[0 1])
        elseif group==3
            figure
            
            for age=1:2
            uselist = (ageList==age & (layer==5) & ~inh & used');
          
                    % jitterAmount = 0.05;
                     jitterValuesX = 2*(rand(size(data1(uselist)))-0.5)*.02;   % +/-jitterAmount max
                     %jitterValuesY = 2*(rand(size(data2(uselist)))-0.5)*jitterAmount;   % +/-jitterAmount max
                     scatter(data1(uselist)+jitterValuesX, data2(uselist),[colors(age) 'o']);hold on
                    %  plot(data1(uselist),data2(uselist),[colors(age) 'o']);hold on
            end
                    title 'layer5'
                    legend('EO1','Adult');
                    xlabel(label{1}); ylabel(label{2});
                    axis square
                    axis equal
                    plot([0 1],[0 1])
        elseif group==4
            figure
            
            for age=1:2
            uselist = (ageList==age & (layer==6) & ~inh &  used');
              
                    % jitterAmount = 0.05;
                     jitterValuesX = 2*(rand(size(data1(uselist)))-0.5)*.02;   % +/-jitterAmount max
                     %jitterValuesY = 2*(rand(size(data2(uselist)))-0.5)*jitterAmount;   % +/-jitterAmount max
                     scatter(data1(uselist)+jitterValuesX, data2(uselist),[colors(age) 'o']);hold on
                    %  plot(data1(uselist),data2(uselist),[colors(age) 'o']);hold on
            end
            
                    title 'layer6'
                    legend('EO1','Adult');
                    xlabel(label{1}); ylabel(label{2});
                    axis square
                    axis equal
                    plot([0 1],[0 1])
%         elseif group==5
%             uselist = (ageList==age & inh & used');
        elseif group==5
            figure
            for age=1:2
            
            uselist = (ageList==age & (layer<=6)& ~inh & used');
                  
                    % jitterAmount = 0.05;
                     jitterValuesX = 2*(rand(size(data1(uselist)))-0.5)*.02;   % +/-jitterAmount max
                     %jitterValuesY = 2*(rand(size(data2(uselist)))-0.5)*jitterAmount;   % +/-jitterAmount max
                     scatter(data1(uselist)+jitterValuesX, data2(uselist),[colors(age) 'o']);hold on
                    %  plot(data1(uselist),data2(uselist),[colors(age) 'o']);hold on
    
            end
                    title 'All layers'
                    legend('EO1','Adult');
                    xlabel(label{1}); ylabel(label{2});
                    axis square
                    axis equal
                    plot([0 1],[0 1])
    end
    
  
            end
        end

             
             
      
    


      

   




    
