function chamaTudo(folderExample,opt)
    
    %% ------------------Carrega os parâmetros das vars--------------------
    load([folderExample '/data']);

    
    %% ------------------Avalia a Estabilidade-----------------------------
    raioEspectralInfe = verificaEstSisHib(ACal.inf, KCal, h);
    raioEspectralSupe = verificaEstSisHib(ACal.sup, KCal, h);
    
    %% ------------------Avalia norma Matlab ------------------------------
    saidaNormMatInf = normaSistemaContinuo(A.inf+B.inf*K,E.inf,C.inf+D.inf*K,zeros(size(C.inf,1),size(E.inf,2)));
    saidaNormMatSup = normaSistemaContinuo(A.sup+B.sup*K,E.sup,C.sup+D.sup*K,zeros(size(C.sup,1),size(E.sup,2)));
    
%     saidaNormMatInf = normaSistemaContinuo(ACal.inf,ECal.inf,CCal.inf,zeros(size(CCal.inf,1),size(ECal.inf,2)));
%     saidaNormMatSup = normaSistemaContinuo(ACal.sup,ECal.sup,CCal.sup,zeros(size(CCal.inf,1),size(ECal.inf,2)));
%     
   %% ------------------Chama Avaliação Estabilidade com HInf garantida Politopico---
    
     saidaAnaHInfPoliLMILabFinal = estHInfAnaPolyLMILab(sysCalPoly,KCal,h,delta,tol);
     saidaAnaHInfPoli = estHInfAnaPoly(sysCalPoly,KCal,h,delta,tol);
     
    %% ------------------Chama Avaliação Estabilidade com HInf garantida---
    saidaEstHInfLMILab = valEstHInfLMILabInt(ACal,ECal,CCal,KCal,h,delta,tol);
     saidaEstHInf = valEstHInfInt(ACal,ECal,CCal,KCal,h,delta,tol);
 
    
    %% ------------------Simulação ----------------------------------------
    
    param.onlyVertice = onlyVertice;
    combSysPolyCont = genCombPoly(sysPolyCont, numPointsUniSpaced, numPointsBy2Points, numbPointsUniSpacedSub,param);
    
    combAPolyCont.polytopicMatrices = cellfun(@(X) X.A,combSysPolyCont.polytopicMatrices,'UniformOutput',false);
    combAPolyCont.alphaVecs = combSysPolyCont.alphaVecs;

    combBPolyCont.polytopicMatrices = cellfun(@(X) X.B2,combSysPolyCont.polytopicMatrices,'UniformOutput',false);
    combBPolyCont.alphaVecs = combSysPolyCont.alphaVecs;

    combEPolyCont.polytopicMatrices = cellfun(@(X) X.B1,combSysPolyCont.polytopicMatrices,'UniformOutput',false);
    combEPolyCont.alphaVecs = combSysPolyCont.alphaVecs;

    combCPolyCont.polytopicMatrices = cellfun(@(X) X.C,combSysPolyCont.polytopicMatrices,'UniformOutput',false);
    combCPolyCont.alphaVecs = combSysPolyCont.alphaVecs;

    combDPolyCont.polytopicMatrices = cellfun(@(X) X.D2,combSysPolyCont.polytopicMatrices,'UniformOutput',false);
    combDPolyCont.alphaVecs = combSysPolyCont.alphaVecs;

    combD1PolyCont.polytopicMatrices = cellfun(@(X) X.D1,combSysPolyCont.polytopicMatrices,'UniformOutput',false);
    combD1PolyCont.alphaVecs = combSysPolyCont.alphaVecs;
    
    
    combPolyCont.A = combAPolyCont;
    combPolyCont.B = combBPolyCont;
    combPolyCont.E = combEPolyCont;
    combPolyCont.C = combCPolyCont;
    combPolyCont.D = combDPolyCont;
    combPolyCont.D1 = combD1PolyCont;
    
  
    
    imageName = 'KIntHibrido';
    axisVector = [0 6 -1 1.5];
    isInt = 1;
    [outPutSimHibridoInt] = simulatesSampledInput(combPolyCont,h, K,[folderExample '/' imageName],axisVector,delta,tol,isInt);
    isInt = 0;
    [outPutSimHibridoIntPoly] = simulatesSampledInput(combPolyCont,h, K,[folderExample '/' imageName],axisVector,delta,tol,isInt);
    
    
    
  
    
    save([folderExample '/dataProcessado']);
       
end