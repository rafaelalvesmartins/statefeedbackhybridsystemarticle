function resultado = chamaTudoNomCont(paraSys, folderExample,opt)
    addpath('../HInf - Análise - Intervalar/funcoes');
    addpath('../Hinf - Análise/funcoes');
    
    tempoRefencia = datetime('now');
    %% ------------------Carrega os parâmetros das vars--------------------
    % Sistema
    A = paraSys.sys.A;
    B = paraSys.sys.B;
    E = paraSys.sys.E;
    C = paraSys.sys.C;
    D = paraSys.sys.D;
    sysPolyCont = paraSys.sys.sysPolyCont;
    
    
    %     Tamanho das matrizes
    nx = length(A.inf);
    nu = size(B.inf,2);
        
    ACal = paraSys.sys.ACal;
    ECal = paraSys.sys.ECal;
    CCal = paraSys.sys.CCal;
    h = paraSys.sys.h;
    delta = paraSys.sys.delta;
    
    % Aux
    tol = paraSys.aux.tol;
    
    
    % Simulation
    numPointsUniSpaced = paraSys.simSys.numPointsUniSpaced;
    numPointsBy2Points = paraSys.simSys.numPointsBy2Points;
    numbPointsUniSpacedSub = paraSys.simSys.numbPointsUniSpacedSub;
    onlyVertice =  paraSys.simSys.onlyVertice;
    
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

     save([folderExample '/saida']);
    combPolyCont.A = combAPolyCont;
    combPolyCont.B = combBPolyCont;
    combPolyCont.E = combEPolyCont;
    combPolyCont.C = combCPolyCont;
    combPolyCont.D = combDPolyCont;
    combPolyCont.D1 = combD1PolyCont;

     % Save parameters to outPut    
    saida.combCont = combPolyCont;
    
    tempoAgora = datetime('now');
    disp('Gerar sistema politopico continuo, tempo de execucao');tempoAgora-tempoRefencia
    tempoRefencia = tempoAgora;
    
   
    
     %% ------------------Chama Avaliação Estabilidade com HInf garantida---
    saidaEstHInfSintaseLMILab = estHInfSintLMILab(A,B,E,C,D,h,delta,tol);
    KCal = [eye(nx) zeros(nx,nu); saidaEstHInfSintaseLMILab.K zeros(nu)];
    raioEspectralInfLMILab = verificaEstSisHib(ACal.inf, KCal, h);
    raioEspectralSupLMILab = verificaEstSisHib(ACal.sup, KCal, h);
    saidaNormMatInfLMILab = normaSistemaContinuo(A.inf+B.inf*saidaEstHInfSintaseLMILab.K,E.inf,C.inf+D.inf*saidaEstHInfSintaseLMILab.K,zeros(size(C.inf,1),size(E.inf,2)));
    saidaNormMatSupLMILab = normaSistemaContinuo(A.sup+B.sup*saidaEstHInfSintaseLMILab.K,E.sup,C.sup+D.sup*saidaEstHInfSintaseLMILab.K,zeros(size(C.sup,1),size(E.sup,2)));
    saidaEstHInfLMILab = valEstHInfLMILabInt(ACal,ECal,CCal,KCal,h,delta,tol);
    saidaEstHInf = valEstHInfInt(ACal,ECal,CCal,KCal,h,delta,tol);
    save([folderExample '/saida']);

    LMILab.saidaEstHInfSintaseLMILab = saidaEstHInfSintaseLMILab;
    LMILab.raioEspectralInf = raioEspectralInfLMILab;
    LMILab.raioEspectralSup = raioEspectralSupLMILab;
    LMILab.saidaNormMatInf = saidaNormMatInfLMILab;
    LMILab.saidaNormMatSup = saidaNormMatSupLMILab;
    LMILab.saidaEstHInfLMILab = saidaEstHInfLMILab;
    LMILab.saidaEstHInf = saidaEstHInf;

    saida.LMILab = LMILab;

    tempoAgora = datetime('now');
    disp('Sintase LMI Lab Hibrido, tempo de execucao');tempoAgora-tempoRefencia
    tempoRefencia = tempoAgora;
    
    %% ------------------Comparação 1: Controlador Nominal Sampled-Data---
    A_nominal = paraSys.sys.A_nominal;
    B_nominal = paraSys.sys.B_nominal;
    E_nominal = paraSys.sys.E_nominal;
    C_nominal = paraSys.sys.C_nominal;
    D_nominal = paraSys.sys.D_nominal;

    % Converte para formato intervalar trivial (para usar funções existentes)
    A_nominal_int = infsup(A_nominal, A_nominal);
    B_nominal_int = infsup(B_nominal, B_nominal);
    E_nominal_int = infsup(E_nominal, E_nominal);
    C_nominal_int = infsup(C_nominal, C_nominal);
    D_nominal_int = infsup(D_nominal, D_nominal);

    saidaEstHInfNominal = estHInfSintLMILab(A_nominal_int,B_nominal_int,E_nominal_int,C_nominal_int,D_nominal_int,h,delta,tol);

    % Testa controlador nominal no sistema incerto
    KCal_nominal = [eye(nx) zeros(nx,nu); saidaEstHInfNominal.K zeros(nu)];
    saidaEstHInf_nominal_incerto = valEstHInfInt(ACal,ECal,CCal,KCal_nominal,h,delta,tol);

    nominal_comparison.saidaEstHInfNominal = saidaEstHInfNominal;
    nominal_comparison.saidaEstHInf_nominal_incerto = saidaEstHInf_nominal_incerto;
    saida.nominal_comparison = nominal_comparison;

    tempoAgora = datetime('now');
    disp('Controlador Nominal, tempo de execucao');tempoAgora-tempoRefencia
    tempoRefencia = tempoAgora;

    %% ------------------Comparação 2: Controlador Contínuo Robusto Discretizado---
    % Projeta controlador contínuo robusto (usando sistema intervalar)
    sys_cont_robust = paraSys.sys.sys_continuous_robust;

    % Síntese contínua (adaptada para caso contínuo)
    try
        saidaEstHInfContinuous = estHInfSintLMILab(sys_cont_robust.A,sys_cont_robust.B,sys_cont_robust.E,sys_cont_robust.C,sys_cont_robust.D,inf,delta,tol);
    catch
        % Se a função não aceitar h=inf, use um valor muito grande
        saidaEstHInfContinuous = estHInfSintLMILab(sys_cont_robust.A,sys_cont_robust.B,sys_cont_robust.E,sys_cont_robust.C,sys_cont_robust.D,1000,delta,tol);
    end

    % Discretização do controlador (implementação ZOH implícita via sampled-data)
    KCal_discretized = [eye(nx) zeros(nx,nu); saidaEstHInfContinuous.K zeros(nu)];
    saidaEstHInf_discretized = valEstHInfInt(ACal,ECal,CCal,KCal_discretized,h,delta,tol);

    continuous_comparison.saidaEstHInfContinuous = saidaEstHInfContinuous;
    continuous_comparison.saidaEstHInf_discretized = saidaEstHInf_discretized;
    saida.continuous_comparison = continuous_comparison;

    tempoAgora = datetime('now');
    disp('Controlador Continuo Discretizado, tempo de execucao');tempoAgora-tempoRefencia
    tempoRefencia = tempoAgora;
     
    %% ------------------Chama Avaliação Estabilidade com HInf garantida Politópica---
    % Laço para converter estrutura de matriz
    n = length(sysPolyCont);
    APoly = cell(n,1);
    BPoly = cell(n,1);
    EPoly = cell(n,1);
    CPoly = cell(n,1);
    DPoly = cell(n,1);
    for c = 1:n
        APoly{c} = sysPolyCont{c}.A;
        BPoly{c} = sysPolyCont{c}.B2;
        EPoly{c} = sysPolyCont{c}.B1;
        CPoly{c} = sysPolyCont{c}.C;
        DPoly{c} = sysPolyCont{c}.D2;
    end
    
    poly.APoly = APoly;
    poly.BPoly = BPoly;
    poly.EPoly = EPoly;
    poly.CPoly = CPoly;
    poly.DPoly = DPoly;
    
   %% ------------------Chama Avaliação Estabilidade com HInf garantida Politópica---
    
    saidaEstHInfSintaseLMILabPoly = estHInfSintPolyLMILab(poly,h,delta,tol);
    KCal = [eye(nx) zeros(nx,nu); saidaEstHInfSintaseLMILabPoly.K zeros(nu)];
    raioEspectralInf = verificaEstSisHib(ACal.inf, KCal, h);
    raioEspectralSup = verificaEstSisHib(ACal.sup, KCal, h);
    saidaNormMatInf = normaSistemaContinuo(A.inf+B.inf*saidaEstHInfSintaseLMILabPoly.K,E.inf,C.inf+D.inf*saidaEstHInfSintaseLMILabPoly.K,zeros(size(C.inf,1),size(E.inf,2)));
    saidaNormMatSup = normaSistemaContinuo(A.sup+B.sup*saidaEstHInfSintaseLMILabPoly.K,E.sup,C.sup+D.sup*saidaEstHInfSintaseLMILabPoly.K,zeros(size(C.sup,1),size(E.sup,2)));
    saidaEstHInfLMILab = valEstHInfLMILabInt(ACal,ECal,CCal,KCal,h,delta,tol);
    saidaEstHInf = valEstHInfInt(ACal,ECal,CCal,KCal,h,delta,tol);
    save([folderExample '/saida']);
    
    % Salva para a saída
    politopo.saidaEstHInfSintaseLMILabPoly = saidaEstHInfSintaseLMILabPoly;
    politopo.raioEspectralInf = raioEspectralInf;
    politopo.raioEspectralSup = raioEspectralSup;
    politopo.saidaNormMatInf = saidaNormMatInf;
    politopo.saidaNormMatSup = saidaNormMatSup;
    politopo.saidaEstHInfLMILab = saidaEstHInfLMILab;
    politopo.saidaEstHInf = saidaEstHInf;
    saida.politopo = politopo;
    
    tempoAgora = datetime('now');
    disp('Sintase LMIL Politopo, tempo de execucao');tempoAgora-tempoRefencia
    tempoRefencia = tempoAgora;
    
    
    %%  ------------------Chama Avaliação Estabilidade com HInf garantida Politópica Sedumu---
                                    
    saidaEstHInfSintasePoly = estHInfSintPoly(poly,h,delta,tol);
    KCal = [eye(nx) zeros(nx,nu); saidaEstHInfSintasePoly.K zeros(nu)];
    raioEspectralInf = verificaEstSisHib(ACal.inf, KCal, h);
    raioEspectralSup = verificaEstSisHib(ACal.sup, KCal, h);
    saidaNormMatInf = normaSistemaContinuo(A.inf+B.inf*saidaEstHInfSintasePoly.K,E.inf,C.inf+D.inf*saidaEstHInfSintasePoly.K,zeros(size(C.inf,1),size(E.inf,2)));
    saidaNormMatSup = normaSistemaContinuo(A.sup+B.sup*saidaEstHInfSintasePoly.K,E.sup,C.sup+D.sup*saidaEstHInfSintasePoly.K,zeros(size(C.sup,1),size(E.sup,2)));
    saidaEstHInfLMILab = valEstHInfLMILabInt(ACal,ECal,CCal,KCal,h,delta,tol);
    saidaEstHInf = valEstHInfInt(ACal,ECal,CCal,KCal,h,delta,tol);
    save([folderExample '/saida']);
    
    % Salva para a saída
    politopoSedumi.saidaEstHInfSintaseLMILabPoly = saidaEstHInfSintasePoly;
    politopoSedumi.raioEspectralInf = raioEspectralInf;
    politopoSedumi.raioEspectralSup = raioEspectralSup;
    politopoSedumi.saidaNormMatInf = saidaNormMatInf;
    politopoSedumi.saidaNormMatSup = saidaNormMatSup;
    politopoSedumi.saidaEstHInfLMILab = saidaEstHInfLMILab;
    politopoSedumi.saidaEstHInf = saidaEstHInf;
    saida.politopoSedumu = politopoSedumi;
    
    tempoAgora = datetime('now');
    disp('Sintase LMIL Politopo Sedumu, tempo de execucao');tempoAgora-tempoRefencia
    tempoRefencia = tempoAgora;

    %% ------------------Análise Comparativa de Estabilidade---
    % Calcula taxas de estabilidade para os três métodos
    resultado_estabilidade = analisaEstabilidadeComparativa(combPolyCont, ...
        saidaEstHInfSintaseLMILab.K, saidaEstHInfNominal.K, saidaEstHInfContinuous.K, h);
    
    saida.resultado_estabilidade = resultado_estabilidade;
   
    
    %% Simulacao ganho híbrido
    imageName = 'KIntHibrido';
    axisVector = [0 10 -0.5 1.5];
    flagIsPoly = false;
    [outPutSimHibridoInt] = simulatesSampledInput(combPolyCont,h, saidaEstHInfSintaseLMILab.K,[folderExample '/' imageName],axisVector,delta,tol,flagIsPoly);
    
    saida.outPutSimHibridoInt = outPutSimHibridoInt;
   
    tempoAgora = datetime('now');
    disp('Simulacao temporal, tempo de execucao');tempoAgora-tempoRefencia
    tempoRefencia = tempoAgora;
    
    
     %% Simulacao ganho híbrido Politopo
    imageName = 'KIntHibridoPoly';
    axisVector = [0 10 -0.5 1.5];
    flagIsPoly = true;
    [outPutSimHibridoIntPoly] = simulatesSampledInput(combPolyCont,h, saidaEstHInfSintaseLMILabPoly.K,[folderExample '/' imageName],axisVector,delta,tol,flagIsPoly);
    
    saida.outPutSimHibridoIntPoly = outPutSimHibridoIntPoly;
   
    tempoAgora = datetime('now');
    disp('Simulacao temporal Politopo, tempo de execucao');tempoAgora-tempoRefencia
    tempoRefencia = tempoAgora;

    %% Simulação comparativa dos três métodos
    imageName = 'ComparacaoMetodos';
    axisVector = [0 10 -0.5 1.5];

    % Simular os três controladores
    flagIsPoly = false;
    [outPutSim_Nominal] = simulatesSampledInput(combPolyCont,h, saidaEstHInfNominal.K,[folderExample '/' imageName '_Nominal'],axisVector,delta,tol,flagIsPoly);
    [outPutSim_Discretized] = simulatesSampledInput(combPolyCont,h, saidaEstHInfContinuous.K,[folderExample '/' imageName '_Discretized'],axisVector,delta,tol,flagIsPoly);

    saida.outPutSim_Nominal = outPutSim_Nominal;
    saida.outPutSim_Discretized = outPutSim_Discretized;

    tempoAgora = datetime('now');
    disp('Simulacao comparativa, tempo de execucao');tempoAgora-tempoRefencia
    tempoRefencia = tempoAgora;
    
    save([folderExample '/saida']);
end

%% Função auxiliar para análise de estabilidade
function resultado = analisaEstabilidadeComparativa(combPolyCont, K_proposto, K_nominal, K_discretizado, h)
    % Analisa estabilidade dos três controladores
    
    nRealizacoes = length(combPolyCont.A.polytopicMatrices);
    estavel_proposto = zeros(nRealizacoes, 1);
    estavel_nominal = zeros(nRealizacoes, 1);
    estavel_discretizado = zeros(nRealizacoes, 1);
    
    custos_proposto = zeros(nRealizacoes, 1);
    custos_nominal = zeros(nRealizacoes, 1);
    custos_discretizado = zeros(nRealizacoes, 1);
    
    for i = 1:nRealizacoes
        A_real = combPolyCont.A.polytopicMatrices{i};
        B_real = combPolyCont.B.polytopicMatrices{i};
        C_real = combPolyCont.C.polytopicMatrices{i};
        D_real = combPolyCont.D.polytopicMatrices{i};
        E_real = combPolyCont.E.polytopicMatrices{i};
        
        % Sistema em malha fechada para cada controlador
        Acl_proposto = A_real + B_real * K_proposto;
        Acl_nominal = A_real + B_real * K_nominal;
        Acl_discretizado = A_real + B_real * K_discretizado;
        
        % Verifica estabilidade
        estavel_proposto(i) = max(real(eig(Acl_proposto))) < -1e-6;
        estavel_nominal(i) = max(real(eig(Acl_nominal))) < -1e-6;
        estavel_discretizado(i) = max(real(eig(Acl_discretizado))) < -1e-6;
        
        % Calcula custos H-infinito (se estável)
        if estavel_proposto(i)
            try
                Ccl_proposto = C_real + D_real * K_proposto;
                custos_proposto(i) = normaSistemaContinuo(Acl_proposto, E_real, Ccl_proposto, zeros(size(C_real,1),size(E_real,2)));
            catch
                custos_proposto(i) = inf;
            end
        else
            custos_proposto(i) = inf;
        end
        
        if estavel_nominal(i)
            try
                Ccl_nominal = C_real + D_real * K_nominal;
                custos_nominal(i) = normaSistemaContinuo(Acl_nominal, E_real, Ccl_nominal, zeros(size(C_real,1),size(E_real,2)));
            catch
                custos_nominal(i) = inf;
            end
        else
            custos_nominal(i) = inf;
        end
        
        if estavel_discretizado(i)
            try
                Ccl_discretizado = C_real + D_real * K_discretizado;
                custos_discretizado(i) = normaSistemaContinuo(Acl_discretizado, E_real, Ccl_discretizado, zeros(size(C_real,1),size(E_real,2)));
            catch
                custos_discretizado(i) = inf;
            end
        else
            custos_discretizado(i) = inf;
        end
    end
    
    resultado.taxa_estabilidade_proposto = sum(estavel_proposto) / nRealizacoes * 100;
    resultado.taxa_estabilidade_nominal = sum(estavel_nominal) / nRealizacoes * 100;
    resultado.taxa_estabilidade_discretizado = sum(estavel_discretizado) / nRealizacoes * 100;
    
    % Custos apenas para sistemas estáveis
    custos_proposto_finitos = custos_proposto(isfinite(custos_proposto));
    custos_nominal_finitos = custos_nominal(isfinite(custos_nominal));
    custos_discretizado_finitos = custos_discretizado(isfinite(custos_discretizado));
    
    if ~isempty(custos_proposto_finitos)
        resultado.custo_pior_caso_proposto = max(custos_proposto_finitos);
        resultado.custo_medio_proposto = mean(custos_proposto_finitos);
    else
        resultado.custo_pior_caso_proposto = inf;
        resultado.custo_medio_proposto = inf;
    end
    
    if ~isempty(custos_nominal_finitos)
        resultado.custo_pior_caso_nominal = max(custos_nominal_finitos);
        resultado.custo_medio_nominal = mean(custos_nominal_finitos);
    else
        resultado.custo_pior_caso_nominal = inf;
        resultado.custo_medio_nominal = inf;
    end
    
    if ~isempty(custos_discretizado_finitos)
        resultado.custo_pior_caso_discretizado = max(custos_discretizado_finitos);
        resultado.custo_medio_discretizado = mean(custos_discretizado_finitos);
    else
        resultado.custo_pior_caso_discretizado = inf;
        resultado.custo_medio_discretizado = inf;
    end
    resultad.saida = saida;
end