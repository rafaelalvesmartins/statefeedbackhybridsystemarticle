function outPut = masterFunction(paraSys, folderExample,opt)
    
   %% ------------------Load Parameters from var-------------------------
    % sys
    A = paraSys.sys.A;
    B2 = paraSys.sys.B2;
    B1 = paraSys.sys.B1;
    C = paraSys.sys.C;
    D2 = paraSys.sys.D2;
    D1 = paraSys.sys.D1;
    h = paraSys.sys.h;
    
%     APoly = cellfun(@(X) X.A,paraSys.poly,'UniformOutput',false);
%     B2Poly =  cellfun(@(X) X.B2,paraSys.poly,'UniformOutput',false); 
%     B1Poly =  cellfun(@(X) X.B1,paraSys.poly,'UniformOutput',false); 
%     CPoly =  cellfun(@(X) X.C,paraSys.poly,'UniformOutput',false); 
%     D2Poly =  cellfun(@(X) X.D2,paraSys.poly,'UniformOutput',false); 
%     D1Poly =  cellfun(@(X) X.D1,paraSys.poly,'UniformOutput',false);
%     
    % Simulation
    numPointsUniSpaced = paraSys.simSys.numPointsUniSpaced;
    numPointsBy2Points = paraSys.simSys.numPointsBy2Points;
    numbPointsUniSpacedSub = paraSys.simSys.numbPointsUniSpacedSub;
    onlyVertice =  paraSys.simSys.onlyVertice;
        
    % Aux
    tol = paraSys.aux.tol;
    precision = paraSys.aux.precision;
    
    % Save parameters to outPut      
    outPut.paraSys = paraSys;
    outPut.opt = opt;
     %% ------------------Discretize Sys-------------------------------------
        BChap = [B2 B1];
        DChap = [D2 D1];
        [Ad,BChapd,Cd,DChapd] = discIntSys(A,BChap,C,DChap,h);
        
        B2d = BChapd(:,1:size(B2.inf,2));
        B1d = BChapd(:,(size(B2.inf,2)+1):end);
        
        D2d = DChapd(:,1:size(D2.inf,2));
        D1d = DChapd(:,(size(D2.inf,2)+1):end);
        
        % Save parameters to outPut
        sampleSysOutPut.Ad = Ad;
        sampleSysOutPut.B2d = B2d;
        sampleSysOutPut.B1d = B1d;
        sampleSysOutPut.Cd = Cd;
        sampleSysOutPut.D2d = D2d;
        sampleSysOutPut.D1d = D1d;
        outPut.sampleSysOutPut = sampleSysOutPut;  
        
    %% ------------------Generate Polytopic Structure----------------------  
    
    APoly = genVerticesPolyFromUncerMat(Ad,tol);
    B2Poly =  genVerticesPolyFromUncerMat(B2d,tol);
    B1Poly =  genVerticesPolyFromUncerMat(B1d,tol);
    CPoly =  genVerticesPolyFromUncerMat(Cd,tol);
    D2Poly =  genVerticesPolyFromUncerMat(D2d,tol);
    D1Poly =  genVerticesPolyFromUncerMat(D1d,tol);
    
    
%     a = 1;
%     b1Len = length(B1Poly);
%     cLen = length(CPoly);
%     d2Len = length(D2Poly);
%     d1Len = length(D1Poly);

     for i=1:length(APoly)
        for j=1:length(B2Poly)
            for k=1:length(B1Poly)
                for l=1:length(CPoly)
                    for m=1:length(D2Poly)
                        for n=1:length(D1Poly)
                         indNvec = n;             
                         indMvec = (m-1)*length(D1Poly);
                         indLvec = (l-1)*length(D2Poly)*length(D1Poly);
                         indKvec = (k-1)*length(CPoly)*length(D2Poly)*length(D1Poly);
                         indJvec = (j-1)*length(B1Poly)*length(CPoly)*length(D2Poly)*length(D1Poly);
                         indIvec = (i-1)*length(B2Poly)*length(B1Poly)*length(CPoly)*length(D2Poly)*length(D1Poly);

                    
                         sysPoly{indNvec+indMvec+indLvec+indKvec+indJvec+indIvec}.A = APoly{i};
                         sysPoly{indNvec+indMvec+indLvec+indKvec+indJvec+indIvec}.B2 = B2Poly{j};
                         sysPoly{indNvec+indMvec+indLvec+indKvec+indJvec+indIvec}.B1 = B1Poly{k};
                         sysPoly{indNvec+indMvec+indLvec+indKvec+indJvec+indIvec}.C = CPoly{l};
                         sysPoly{indNvec+indMvec+indLvec+indKvec+indJvec+indIvec}.D2 = D2Poly{m};
                         sysPoly{indNvec+indMvec+indLvec+indKvec+indJvec+indIvec}.D1 = D1Poly{n};
                        end
                    end
                end
            end
        end
     end

    
    
    param.onlyVertice = onlyVertice;
    combSysPoly = genCombPoly(sysPoly, numPointsUniSpaced, numPointsBy2Points, numbPointsUniSpacedSub,param);

    combAPoly.polytopicMatrices = cellfun(@(X) X.A,combSysPoly.polytopicMatrices,'UniformOutput',false);
    combAPoly.alphaVecs = combSysPoly.alphaVecs;

    combB2Poly.polytopicMatrices = cellfun(@(X) X.B2,combSysPoly.polytopicMatrices,'UniformOutput',false);
    combB2Poly.alphaVecs = combSysPoly.alphaVecs;

    combB1Poly.polytopicMatrices = cellfun(@(X) X.B1,combSysPoly.polytopicMatrices,'UniformOutput',false);
    combB1Poly.alphaVecs = combSysPoly.alphaVecs;

    combCPoly.polytopicMatrices = cellfun(@(X) X.C,combSysPoly.polytopicMatrices,'UniformOutput',false);
    combCPoly.alphaVecs = combSysPoly.alphaVecs;

    combD2Poly.polytopicMatrices = cellfun(@(X) X.D2,combSysPoly.polytopicMatrices,'UniformOutput',false);
    combD2Poly.alphaVecs = combSysPoly.alphaVecs;

    combD1Poly.polytopicMatrices = cellfun(@(X) X.D1,combSysPoly.polytopicMatrices,'UniformOutput',false);
    combD1Poly.alphaVecs = combSysPoly.alphaVecs;

    combPoly.A = combAPoly;
    combPoly.B2 = combB2Poly;
    combPoly.B1 = combB1Poly;
    combPoly.C = combCPoly;
    combPoly.D2 = combD2Poly;
    combPoly.D1 = combD1Poly;

     % Save parameters to outPut    
    comb.Poly = combPoly;

    outPut.comb = comb;
   
    
 
        
    %% ------------------Call main Part----------------------------------------
    
    APoly = cellfun(@(X) X.A,sysPoly,'UniformOutput',false);
    B2Poly =  cellfun(@(X) X.B2,sysPoly,'UniformOutput',false); 
    B1Poly =  cellfun(@(X) X.B1,sysPoly,'UniformOutput',false); 
    CPoly =  cellfun(@(X) X.C,sysPoly,'UniformOutput',false); 
    D2Poly =  cellfun(@(X) X.D2,sysPoly,'UniformOutput',false); 
    D1Poly =  cellfun(@(X) X.D1,sysPoly,'UniformOutput',false);
    [outPutMainPart] =  discMainPart(opt,Ad,B2d,B1d,Cd,D2d,D1d,tol,precision,combPoly,folderExample,APoly,B2Poly,B1Poly,CPoly,D2Poly,D1Poly,h);
    
    %   -------------- Mid-------
% Poly Mid
    outPutMainPart.midOutPut.poly.syn.norm = outPutMainPart.midOutPut.poly.syn.norm/sqrt(h);
  
% Int Mid
    outPutMainPart.midOutPut.int.syn.norm = outPutMainPart.midOutPut.int.syn.norm/sqrt(h);
   
%   -------------- Int-------
    outPutMainPart.intOutput.syn.norm = outPutMainPart.intOutput.syn.norm/sqrt(h);
   
%   -------------- Poly-------
    outPutMainPart.polyOutput.syn.norm = outPutMainPart.polyOutput.syn.norm/sqrt(h);
   
    
    outPut.outPutMainPart = outPutMainPart;
    
    
    
    
%    ------------------ Generate polytopic system based on the interval sampled continuous system---
    APolyCont = genVerticesPolyFromUncerMat(A,tol);
    B2PolyCont =  genVerticesPolyFromUncerMat(B2,tol);
    B1PolyCont =  genVerticesPolyFromUncerMat(B1,tol);
    CPolyCont =  genVerticesPolyFromUncerMat(C,tol);
    D2PolyCont =  genVerticesPolyFromUncerMat(D2,tol);
    D1PolyCont =  genVerticesPolyFromUncerMat(D1,tol);
    
    


     for i=1:length(APolyCont)
        for j=1:length(B2PolyCont)
            for k=1:length(B1PolyCont)
                for l=1:length(CPolyCont)
                    for m=1:length(D2PolyCont)
                        for n=1:length(D1PolyCont)
                         indNvec = n;             
                         indMvec = (m-1)*length(D1PolyCont);
                         indLvec = (l-1)*length(D2PolyCont)*length(D1PolyCont);
                         indKvec = (k-1)*length(CPolyCont)*length(D2PolyCont)*length(D1PolyCont);
                         indJvec = (j-1)*length(B1PolyCont)*length(CPolyCont)*length(D2PolyCont)*length(D1PolyCont);
                         indIvec = (i-1)*length(B2PolyCont)*length(B1PolyCont)*length(CPolyCont)*length(D2PolyCont)*length(D1PolyCont);

                    
                         sysPolyCont{indNvec+indMvec+indLvec+indKvec+indJvec+indIvec}.A = APolyCont{i};
                         sysPolyCont{indNvec+indMvec+indLvec+indKvec+indJvec+indIvec}.B2 = B2PolyCont{j};
                         sysPolyCont{indNvec+indMvec+indLvec+indKvec+indJvec+indIvec}.B1 = B1PolyCont{k};
                         sysPolyCont{indNvec+indMvec+indLvec+indKvec+indJvec+indIvec}.C = CPolyCont{l};
                         sysPolyCont{indNvec+indMvec+indLvec+indKvec+indJvec+indIvec}.D2 = D2PolyCont{m};
                         sysPolyCont{indNvec+indMvec+indLvec+indKvec+indJvec+indIvec}.D1 = D1PolyCont{n};
                        end
                    end
                end
            end
        end
     end

    
    
    param.onlyVertice = onlyVertice;
    combSysPolyCont = genCombPoly(sysPolyCont, numPointsUniSpaced, numPointsBy2Points, numbPointsUniSpacedSub,param);

    combAPolyCont.polytopicMatrices = cellfun(@(X) X.A,combSysPolyCont.polytopicMatrices,'UniformOutput',false);
    combAPolyCont.alphaVecs = combSysPolyCont.alphaVecs;

    combB2PolyCont.polytopicMatrices = cellfun(@(X) X.B2,combSysPolyCont.polytopicMatrices,'UniformOutput',false);
    combB2PolyCont.alphaVecs = combSysPolyCont.alphaVecs;

    combB1PolyCont.polytopicMatrices = cellfun(@(X) X.B1,combSysPolyCont.polytopicMatrices,'UniformOutput',false);
    combB1PolyCont.alphaVecs = combSysPolyCont.alphaVecs;

    combCPolyCont.polytopicMatrices = cellfun(@(X) X.C,combSysPolyCont.polytopicMatrices,'UniformOutput',false);
    combCPolyCont.alphaVecs = combSysPolyCont.alphaVecs;

    combD2PolyCont.polytopicMatrices = cellfun(@(X) X.D2,combSysPolyCont.polytopicMatrices,'UniformOutput',false);
    combD2PolyCont.alphaVecs = combSysPolyCont.alphaVecs;

    combD1PolyCont.polytopicMatrices = cellfun(@(X) X.D1,combSysPolyCont.polytopicMatrices,'UniformOutput',false);
    combD1PolyCont.alphaVecs = combSysPolyCont.alphaVecs;

    combPolyCont.A = combAPolyCont;
    combPolyCont.B2 = combB2PolyCont;
    combPolyCont.B1 = combB1PolyCont;
    combPolyCont.C = combCPolyCont;
    combPolyCont.D2 = combD2PolyCont;
    combPolyCont.D1 = combD1PolyCont;

     % Save parameters to outPut    
    combCont.Poly = combPolyCont;

    outPut.combCont = combCont;
    
    
    
    imageName = 'KRiccatiSampledSimulation';
    axisVector = [0 15 -0.5 1.5];
    [outPutSimSampRiccati] = simulatesSampledInput(combPolyCont,h, outPutMainPart.midOutPut.riccati.syn.K,[folderExample '/' imageName],axisVector);
    
    imageName = 'KIntSampledSimulation';
    axisVector = [0 10 -0.5 1.5];
    [outPutSimSampInt] = simulatesSampledInput(combPolyCont,h, outPutMainPart.intOutPut.syn.K,[folderExample '/' imageName],axisVector);
    
    imageName = 'KPolySampledSimulation';
    axisVector = [0 10 -0.5 1.5];
    [outPutSimSampPoly] = simulatesSampledInput(combPolyCont,h, outPutMainPart.polyOutput.syn.K,[folderExample '/' imageName],axisVector);

    
    sampledSimulation.riccati = outPutSimSampRiccati;
    sampledSimulation.outPutSimSampInt = outPutSimSampInt;
    sampledSimulation.outPutSimSampPoly = outPutSimSampPoly;
    
    outPut.sampledSimulation = sampledSimulation;
end