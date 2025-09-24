function [outPut] = simulatesSampledInput(combPoly,h,K,imageName,axisVector,delta,tol,isInt)
                                      
    % This function simulates a system with sampled input and calculates quadratic
    % cust
    
    model = 'simulinkSampledInput';    
    options = simset('SrcWorkspace','current');
    load_system(model);
    simulationTime = 30;
    set_param(model, 'StopTime', int2str(simulationTime));
    nx = size(combPoly.A.polytopicMatrices{1},1);
    nu = size(combPoly.B.polytopicMatrices{1},2);
    nw = size(combPoly.E.polytopicMatrices{1},2);
    
    for i=1:length(combPoly.A.alphaVecs)

        A = combPoly.A.polytopicMatrices{i};
        B = combPoly.B.polytopicMatrices{i};
        E = combPoly.E.polytopicMatrices{i};
        C = combPoly.C.polytopicMatrices{i};
        D  =combPoly.D.polytopicMatrices{i};
        D1  =combPoly.D1.polytopicMatrices{i};
        
        
        ACal = [A B ; zeros(nu,nx) zeros(nu,nu)];
        ECal = [E ; zeros(nu,nw)];
        CCal = [C D];
        KCal = [eye(nx) zeros(nx,nu); K zeros(nu)];
        
        if isInt == 1
             saidaEstHInfLMILab = valEstHInfLMILabSemInt(ACal,ECal,CCal,KCal,h,delta,tol);
             gammaLMI(i) = saidaEstHInfLMILab.gamma;
             B = [E B];

               sim(model,[],options);
               custVector(i) = simout.Data(end);

        % Time
               if(i==1)
                  time = simout.Time;
                  lenTime = length(time);
               end

              states{i} =  normalizeSizeVec(simoutStates.Data,lenTime);          
              input = normalizeSizeVec(simoutInput.Data,lenTime);
        else
           poly.APoly{1} = ACal;
           poly.EPoly{1} = ECal;
           poly.CPoly{1} = CCal;
           
           saidaEstHInfLMILab = estHInfAnaPolyLMILabOriginal(poly,KCal,h,delta,tol);
           gammaLMI(i) = saidaEstHInfLMILab.gamma;
        end
        

%        disp([i custVectorCenter(i) ]);
    end
    if isInt == 0
         meanGamma = mean(gammaLMI);
        stdGamma = std(gammaLMI);
        [maxGamma, indexMaxGamma] = max(gammaLMI);
        gamma.meanGamma = meanGamma;
        gamma.stdGamma = stdGamma;
        gamma.maxGamma = maxGamma;
        gamma.indexMaxGamma = indexMaxGamma;
    
    
        outPut.gamma = gamma;
        return
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
  
    meanGamma = mean(gammaLMI);
    stdGamma = std(gammaLMI);
    [maxGamma, indexMaxGamma] = max(gammaLMI);
    gamma.meanGamma = meanGamma;
    gamma.stdGamma = stdGamma;
    gamma.maxGamma = maxGamma;
    gamma.indexMaxGamma = indexMaxGamma;
    
   
   
    plotWithArea(time,meanMat,stdCenterMat,axisVector,input,xLabel,imageName);
    
    % Save parameters to outPut
    outPut.custVector = custVector;
    outPut.meanCost = meanCost;
    outPut.stdCost = stdCost;
    outPut.maxCost = maxCost;
    outPut.indexMaxCost = indexMaxCost; 
    outPut.time = time;
    outPut.gamma = gamma;
    
%     outPut.states = states;
end