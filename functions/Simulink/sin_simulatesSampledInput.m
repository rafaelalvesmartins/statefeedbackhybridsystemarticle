function [outPut] = sin_simulatesSampledInput(combPoly,h,K,imageName,axisVector,w)
    % This function simulates a system with sampled input and calculates quadratic
    % cust
    
    model = 'sin_simulinkSampledInput';    
    options = simset('SrcWorkspace','current');
    load_system(model);
    simulationTime = 30;
    set_param(model, 'StopTime', int2str(simulationTime));
    for i=1:length(combPoly.A.alphaVecs)

        A = combPoly.A.polytopicMatrices{i};
        B2 = combPoly.B2.polytopicMatrices{i};
        B1 = combPoly.B2.polytopicMatrices{i};
        C = combPoly.C.polytopicMatrices{i};
        D2  =combPoly.D2.polytopicMatrices{i};
        D1  =combPoly.D1.polytopicMatrices{i};
        
        B = [B1 B2];
       
       sim(model,[],options);
       custVector(i) = simout.Data(end);
       
% Time
       if(i==1)
          time = simout.Time;
          lenTime = length(time);
       end
       
      states{i} =  normalizeSizeVec(simoutStates.Data,lenTime);          
      input = normalizeSizeVec(simoutInput.Data,lenTime);

%        disp([i custVectorCenter(i) ]);
    end
    
    meanCost = mean(custVector);
    stdCost = std(custVector);
    [maxCost, indexMaxCost] = max(custVector);

    
    meanMat = [];
    stdCenterMat = [];
    xLabel = [];
    for i=1:size(states{1},2)
        
        state = cell2mat(cellfun(@(X) X(:,i),states,'UniformOutput',false));
       
        meanCustVecState = mean(state,2);
        
        stdCustVecState = std(state,0,2);
        
        meanMat = [meanMat meanCustVecState];
        stdCenterMat = [stdCenterMat stdCustVecState];
        if(i~=size(states{1},2))
            xLabel = [xLabel '$x_' int2str(i) '(t)$,'];
        else
            xLabel = [xLabel '$x_' int2str(i) '(t)$'];
        end
    
    end
    
   
    plotWithArea(time,meanMat,stdCenterMat,axisVector,input,xLabel,imageName);
    
    % Save parameters to outPut
    outPut.custVector = custVector;
    outPut.meanCost = meanCost;
    outPut.stdCost = stdCost;
    outPut.maxCost = maxCost;
    outPut.indexMaxCost = indexMaxCost; 
    outPut.time = time;
%     outPut.states = states;
end