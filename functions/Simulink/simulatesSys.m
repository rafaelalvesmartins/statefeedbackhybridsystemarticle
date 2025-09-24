function [outPut] = simulatesSys(combPoly,K,param)
        
    eigVec = zeros(length(combPoly.A.alphaVecs),length(combPoly.A.polytopicMatrices{1}));
    if(nargin == 3 && isfield(param,'disc'))
         for j=1:length(combPoly.A.alphaVecs)
            ACl = combPoly.A.polytopicMatrices{j}+ combPoly.B2.polytopicMatrices{j}*K;
            BCl = combPoly.B1.polytopicMatrices{j};
            CCl  =combPoly.C.polytopicMatrices{j}+combPoly.D2.polytopicMatrices{j}*K;
            DCl = combPoly.D1.polytopicMatrices{j};
            sys = ss(ACl,BCl,CCl,DCl,param.h)*(1/sqrt(param.h));
            [hinfJ(j),FPEAK(j)] = hinfnorm(sys);
            eigVec(j,:) = (eig(ACl))';
         end
         close all;
         scatter(real(eigVec(:)),imag(eigVec(:)),100,[0.45 0.45 0.45],'*','LineWidth',1);
          hold on;
          
    %     Plot unity circle
            theta = [0:0.001:2*pi];
            f = exp(i*theta);
            plot(f, 'k');
            axis square;

            line(xlim, [0 0],'color','k');
            line([0 0],  ylim,'color','k');
            incr = 0.1;
            if(min(real(eigVec(:)))>-1)
                minX = -incr-1;
            else
                minX = -0.1+min(real(eigVec(:)));
            end
            
            if(max(real(eigVec(:)))<1)
                maxX = incr+1;
            else
                maxX = 0.1+max(real(eigVec(:)));
            end
            
            if(min(imag(eigVec(:)))>-1)
                minY = -incr-1;
            else
                minY = -0.1+min(real(eigVec(:)));
            end
            
            if(max(imag(eigVec(:)))<1)
                maxY = incr+1;
            else
                maxY = 0.1+max(imag(eigVec(:)));
            end
            
            axisVector = [minX maxX minY maxY];
              axis(axisVector);
            box on;


            %Graphic Configuration
            ax = gca;
            ax.XLabel.Interpreter = 'latex';
            ax.XLabel.String = 'Real'; %Name of X axes in Latex format

            ax.YLabel.Interpreter = 'latex';
            ax.YLabel.String = 'Imaginary'; %Name of Y axes in Latex format

            ax.FontName = 'Helvetica';
            ax.FontSize = 14;

            %Handle the size of the graph
%             ax.PlotBoxAspectRatio = [1 0.7696 0.7696];
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
            set(findall(gca, 'Type', 'Line'),'LineWidth',1);
            imageNameFile = param.imageNameFile;
            %Saving as matlab figure and as eps
            savefig([imageNameFile '.fig']);
            saveas(af,imageNameFile,'epsc');

         
    else
        for j=1:length(combPoly.A.alphaVecs)
            ACl = combPoly.A.polytopicMatrices{j}+ combPoly.B2.polytopicMatrices{j}*K;
            BCl = combPoly.B1.polytopicMatrices{j};
            CCl  =combPoly.C.polytopicMatrices{j}+combPoly.D2.polytopicMatrices{j}*K;
            DCl = combPoly.D1.polytopicMatrices{j};
            sys = ss(ACl,BCl,CCl,DCl);
            [hinfJ(j),FPEAK(j)] = hinfnorm(sys);
         end
    end
    meanNorm = mean(hinfJ);
    stdNorm = std(hinfJ);
    [maxNorm, indexMaxJ] = max(hinfJ);
    
    
    % Save parameters to outPut
    outPut.meanNorm = meanNorm;
    outPut.stdNorm = stdNorm;
    outPut.maxNorm = maxNorm;
    outPut.indexMaxNorm = indexMaxJ;
    outPut.fpeak = FPEAK(indexMaxJ);
 
end