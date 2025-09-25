function saida = chamaTudo(paraSys, folderExample,opt)
    tempoRefencia = datetime('now');
    %% ------------------Carrega os parâmetros das vars--------------------
        
    poly = paraSys.poly;
    ACal = paraSys.sys.ACal;
    ECal = paraSys.sys.ECal;
    CCal = paraSys.sys.CCal;
    KCal = paraSys.sys.KCal;
    
    A =  paraSys.sys.A;
    B = paraSys.sys.B;
    E = paraSys.sys.E;
    C = paraSys.sys.C;
    D = paraSys.sys.D;
    K = paraSys.sys.K;
    
    h = paraSys.sys.h;
    delta = paraSys.sys.delta;
    
    % Aux
    tol = paraSys.aux.tol;
    
    % Simulation
    numPointsUniSpaced = paraSys.simSys.numPointsUniSpaced;
    numPointsBy2Points = paraSys.simSys.numPointsBy2Points;
    numbPointsUniSpacedSub = paraSys.simSys.numbPointsUniSpacedSub;
    onlyVertice =  paraSys.simSys.onlyVertice;
    
    %%    ------------------ Generate polytopic system based on the interval sampled continuous system---
    ACalPolyCont = genVerticesPolyFromUncerMat(ACal,tol);
    ECalPolyCont =  genVerticesPolyFromUncerMat(ECal,tol);
    CCalPolyCont =  genVerticesPolyFromUncerMat(CCal,tol);

     for i=1:length(ACalPolyCont)
        for j=1:length(ECalPolyCont)
            for k=1:length(CCalPolyCont)
%                 for l=1:length(CPolyCont)
%                     for m=1:length(DPolyCont)
%                         for n=1:length(D1PolyCont)
                         indNvec = k;             
                         indMvec = (j-1)*length(CCalPolyCont);
                         indLvec = (i-1)*length(ECalPolyCont)*length(CCalPolyCont);
%                          indKvec = (k-1)*length(CPolyCont)*length(DPolyCont)*length(D1PolyCont);
%                          indJvec = (j-1)*length(EPolyCont)*length(CPolyCont)*length(DPolyCont)*length(D1PolyCont);
%                          indIvec = (i-1)*length(BPolyCont)*length(EPolyCont)*length(CPolyCont)*length(DPolyCont)*length(D1PolyCont);

                    
                         sysPoly{indNvec+indMvec+indLvec}.ACal = ACalPolyCont{i};
                         sysPoly{indNvec+indMvec+indLvec}.ECal = ECalPolyCont{j};
                         sysPoly{indNvec+indMvec+indLvec}.CCal = CCalPolyCont{k};
%                          sysPoly{indNvec+indMvec+indLvec+indKvec+indJvec+indIvec}.C = CPolyCont{l};
%                          sysPoly{indNvec+indMvec+indLvec+indKvec+indJvec+indIvec}.D = DPolyCont{m};
%                          sysPoly{indNvec+indMvec+indLvec+indKvec+indJvec+indIvec}.D1 = D1PolyCont{n};
%                         end
%                     end
%                 end
            end
        end
     end

    
    
    param.onlyVertice = onlyVertice;
    combSysPolyCont = genCombPoly(sysPoly, numPointsUniSpaced, numPointsBy2Points, numbPointsUniSpacedSub,param);

    combAPolyCont.polytopicMatrices = cellfun(@(X) X.ACal,combSysPolyCont.polytopicMatrices,'UniformOutput',false);
    combAPolyCont.alphaVecs = combSysPolyCont.alphaVecs;

    combEPolyCont.polytopicMatrices = cellfun(@(X) X.ECal,combSysPolyCont.polytopicMatrices,'UniformOutput',false);
    combEPolyCont.alphaVecs = combSysPolyCont.alphaVecs;

    combCPolyCont.polytopicMatrices = cellfun(@(X) X.CCal,combSysPolyCont.polytopicMatrices,'UniformOutput',false);
    combCPolyCont.alphaVecs = combSysPolyCont.alphaVecs;

    combPolyCont.A = combAPolyCont;
    combPolyCont.E = combEPolyCont;
    combPolyCont.C = combCPolyCont;
  
    
    APoly = cellfun(@(X) X.ACal,sysPoly,'UniformOutput',false);
    EPoly =  cellfun(@(X) X.ECal,sysPoly,'UniformOutput',false); 
    CPoly =  cellfun(@(X) X.CCal,sysPoly,'UniformOutput',false); 
    
%     poly.APoly = APoly;
%     poly.EPoly = EPoly;
%     poly.CPoly = CPoly;
    
     % Save parameters to outPut    
    saida.combCont = combPolyCont;
    saida.poly = poly;
    
    tempoAgora = datetime('now');
    disp('Gerar sistema politopico continuo, tempo de execucao');tempoAgora-tempoRefencia
    tempoRefencia = tempoAgora;
    
     %% ------------------Chama Avaliação Estabilidade com HInf garantida---
                   
    saidaNormMatInf = normaSistemaContinuo(A.inf+B.inf*K,E.inf,C.inf+D.inf*K,zeros(size(C.inf,1),size(E.inf,2)));
     saidaAnaHInf = estHInfAnaPolyLMILab(poly,KCal,h,delta,tol);
   
    sedumi.saidaAnaHInf = saidaAnaHInf;
    
    
    
    
    saida.LMILab = sedumi;
    
    tempoAgora = datetime('now');
    disp('Analise LMI Lab Hibrido, tempo de execucao');tempoAgora-tempoRefencia
    tempoRefencia = tempoAgora;
     
    %% ------------------Chama Avaliação Estabilidade com HInf garantida---
%     saidaEstHInfSintase = estHInfSint(A,B,E,C,D,h,delta,tol);
%     KCal = [eye(nx) zeros(nx,nu); saidaEstHInfSintase.K zeros(nu)];
%     raioEspectralInf = verificaEstSisHib(ACal.inf, KCal, h);
%     raioEspectralSup = verificaEstSisHib(ACal.sup, KCal, h);
%     saidaNormMatInf = normaSistemaContinuo(A.inf+B.inf*saidaEstHInfSintase.K,E.inf,C.inf+D.inf*saidaEstHInfSintase.K,zeros(size(C.inf,1),size(E.inf,2)));
%     saidaNormMatSup = normaSistemaContinuo(A.sup+B.sup*saidaEstHInfSintase.K,E.sup,C.sup+D.sup*saidaEstHInfSintase.K,zeros(size(C.sup,1),size(E.sup,2)));
%     saidaEstHInfLMILab = valEstHInfLMILabInt(ACal,ECal,CCal,KCal,h,delta,tol);
%     saidaEstHInf = valEstHInfInt(ACal,ECal,CCal,KCal,h,delta,tol);
%     save([folderExample '/saida']);
    
    % Salva para a saída
%     yalmip.saidaEstHInfSintase = saidaEstHInfSintase;
%     yalmip.raioEspectralInf = raioEspectralInf;
%     yalmip.raioEspectralSup = raioEspectralSup;
%     yalmip.saidaNormMatInf = saidaNormMatInf;
%     yalmip.saidaNormMatSup = saidaNormMatSup;
%     yalmip.saidaEstHInfLMILab = saidaEstHInfLMILab;
%     yalmip.saidaEstHInf = saidaEstHInf;
%     saida.yalmip = yalmip;
%     
%     tempoAgora = datetime('now');
%     disp('Sintase Yalmip Lab Hibrido, tempo de execucao');tempoAgora-tempoRefencia
%     tempoRefencia = tempoAgora;
   
    
    %% Simulacao ganho híbrido
%     imageName = 'KIntHibrido';
%     axisVector = [0 10 -0.5 1.5];
%     [outPutSimHibridoInt] = simulatesSampledInput(combPolyCont,h, saidaEstHInfSintaseLMILab.K,[folderExample '/' imageName],axisVector,delta/2,tol);
%     
%     saida.outPutSimHibridoInt = outPutSimHibridoInt;
%    
%     tempoAgora = datetime('now');
%     disp('Simulacao temporal, tempo de execucao');tempoAgora-tempoRefencia
%     tempoRefencia = tempoAgora;
    
    save([folderExample '/saida']);
  
    
    
    
  
end