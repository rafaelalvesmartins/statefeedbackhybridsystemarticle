function [] = plotWithArea(time,mean,std,axisVector,imageName)
    close all;
    hold on;
    plot(time,mean(:,1:end)+std(:,1:end),'LineWidth', 2,'Color',[0.502 0.502 0.502]);
    plot(time,mean(:,1:end)-std(:,1:end),'LineWidth', 2,'Color',[0.502 0.502 0.502]);
    fill([time;time(end:-1:1)],[mean(:,1:2)+std(:,1:2);mean(end:-1:1,1:2)-std(end:-1:1,1:2)],[0.860 0.860 0.860]);
    plot(time, mean(:,1:end), 'LineWidth', 2, 'Color', [0.5 0.6 0.9]);
 
    
    axis(axisVector);
    box on
   
    
    %Graphic Configuration
    ax = gca;
    ax.XLabel.Interpreter = 'latex';
    ax.XLabel.String = '$t(s)$'; %Name of X axes in Latex format

    ax.YLabel.Interpreter = 'latex';
    ax.YLabel.String = '$x_1(t),x_2(t)$'; %Name of Y axes in Latex format

    ax.FontName = 'Helvetica';
    ax.FontSize = 14;

    %Handle the size of the graph
    ax.PlotBoxAspectRatio = [1 0.7696 0.7696];
    %Handle the size of the window
    af = gcf;
    x0 = 403;
    y0 = 246;
    width = 560;
    height = 420;
    af.Position = [x0 y0 width height];

    %More useful configurations
    %ax.TickLabelInterpreter = 'latex';
    %ax.XLim = [0 10];
    %ax.YLim = [-0.5 2.0];
    %ax.XTick = [0 2.5 5.0 7.5 10];
    %ax.YTick = [-0.5 0 0.5 1 1.5 2];
    ax.YGrid = 'on';
    ax.XGrid = 'on';

%     legend('step','vertex 1','vertex 2','vertex 3','vertex 4','vertex 5','vertex 6','vertex 7','vertex 8')

    %Saving as matlab figure and as eps
    savefig([imageName '.fig']);
    saveas(af,imageName,'epsc');

end