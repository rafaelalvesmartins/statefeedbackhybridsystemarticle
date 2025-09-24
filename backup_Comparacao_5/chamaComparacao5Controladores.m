function saida = chamaComparacao5Controladores()
%   Comparação de 5 métodos de síntese para o Exemplo 6
%
%   Métodos comparados:
%   1. Nominal: estHInfSintLMILabPrec
%   2. Híbrido Intervalar: estHInfSintLMILab
%   3. Híbrido Politópico: estHInfSintPolyLMILab
%   4. Discreto Intervalar: ganho fixo K = [10.8005 24.0549]
%   5. Discreto Politópico: ganho fixo K = [9.7260 22.8904]

    addpath('../HInf - Análise - Intervalar/funcoes');
    addpath('../HInf - Análise - Intervalar/funcoes/Diversas');

    tempoInicial = datetime('now');
    fprintf('=== COMPARAÇÃO DE 5 CONTROLADORES - EXEMPLO 6 ===\n\n');

    %% Carrega o sistema do Exemplo 6
    fprintf('Carregando sistema do Exemplo 6...\n');
    paraSys = genSaveDataEx6Comparacao();

    % Sistema intervalar
    A = paraSys.sys.A;
    B = paraSys.sys.B;
    E = paraSys.sys.E;
    C = paraSys.sys.C;
    D = paraSys.sys.D;
    sysPolyCont = paraSys.sys.sysPolyCont;

    % Parâmetros do sistema
    nx = length(A.inf);
    nu = size(B.inf,2);
    h = paraSys.sys.h;
    delta = paraSys.sys.delta;
    tol = paraSys.aux.tol;

    % Matrizes caligráficas
    ACal = paraSys.sys.ACal;
    ECal = paraSys.sys.ECal;
    CCal = paraSys.sys.CCal;

    % Matrizes nominais (valores centrais) para síntese nominal
    A_nominal = (A.inf + A.sup) / 2;
    B_nominal = (B.inf + B.sup) / 2;
    E_nominal = (E.inf + E.sup) / 2;
    C_nominal = (C.inf + C.sup) / 2;
    D_nominal = (D.inf + D.sup) / 2;

    fprintf('Sistema carregado com sucesso!\n\n');

    %% Prepara sistema politópico para métodos 3 e 5
    fprintf('Preparando sistema politópico...\n');
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

    fprintf('Sistema politópico preparado!\n\n');

    %% 1. SÍNTESE NOMINAL
    fprintf('=== MÉTODO 1: SÍNTESE NOMINAL ===\n');
    tempoRef = datetime('now');
    try
        saidaNominal = estHInfSintLMILabPrec(A_nominal, B_nominal, E_nominal, C_nominal, D_nominal, h, delta, tol);

        % Avaliação do desempenho
        KCal_nominal = [eye(nx) zeros(nx,nu); saidaNominal.K zeros(nu)];
        raioEspInf_nominal = verificaEstSisHib(ACal.inf, KCal_nominal, h);
        raioEspSup_nominal = verificaEstSisHib(ACal.sup, KCal_nominal, h);
        normInf_nominal = normaSistemaContinuo(A.inf+B.inf*saidaNominal.K, E.inf, C.inf+D.inf*saidaNominal.K, zeros(size(C.inf,1),size(E.inf,2)));
        normSup_nominal = normaSistemaContinuo(A.sup+B.sup*saidaNominal.K, E.sup, C.sup+D.sup*saidaNominal.K, zeros(size(C.sup,1),size(E.sup,2)));

        % Estrutura de resultados
        nominal.sintese = saidaNominal;
        nominal.K = saidaNominal.K;
        nominal.raioEspectralInf = raioEspInf_nominal;
        nominal.raioEspectralSup = raioEspSup_nominal;
        nominal.normaSistemaInf = normInf_nominal;
        nominal.normaSistemaSup = normSup_nominal;
        nominal.sucesso = true;

        fprintf('Síntese nominal concluída com sucesso!\n');
        fprintf('Ganho K: [%.4f %.4f]\n', saidaNominal.K(1), saidaNominal.K(2));

    catch ME
        fprintf('ERRO na síntese nominal: %s\n', ME.message);
        nominal.sucesso = false;
        nominal.erro = ME.message;
    end

    tempo = datetime('now') - tempoRef;
    fprintf('Tempo de execução: %.2f segundos\n\n', seconds(tempo));

    %% 2. SÍNTESE HÍBRIDA INTERVALAR
    fprintf('=== MÉTODO 2: SÍNTESE HÍBRIDA INTERVALAR ===\n');
    tempoRef = datetime('now');
    try
        saidaHibInt = estHInfSintLMILab(A, B, E, C, D, h, delta, tol);

        % Avaliação do desempenho
        KCal_hibInt = [eye(nx) zeros(nx,nu); saidaHibInt.K zeros(nu)];
        raioEspInf_hibInt = verificaEstSisHib(ACal.inf, KCal_hibInt, h);
        raioEspSup_hibInt = verificaEstSisHib(ACal.sup, KCal_hibInt, h);
        normInf_hibInt = normaSistemaContinuo(A.inf+B.inf*saidaHibInt.K, E.inf, C.inf+D.inf*saidaHibInt.K, zeros(size(C.inf,1),size(E.inf,2)));
        normSup_hibInt = normaSistemaContinuo(A.sup+B.sup*saidaHibInt.K, E.sup, C.sup+D.sup*saidaHibInt.K, zeros(size(C.sup,1),size(E.sup,2)));

        % Estrutura de resultados
        hibridoIntervalar.sintese = saidaHibInt;
        hibridoIntervalar.K = saidaHibInt.K;
        hibridoIntervalar.raioEspectralInf = raioEspInf_hibInt;
        hibridoIntervalar.raioEspectralSup = raioEspSup_hibInt;
        hibridoIntervalar.normaSistemaInf = normInf_hibInt;
        hibridoIntervalar.normaSistemaSup = normSup_hibInt;
        hibridoIntervalar.sucesso = true;

        fprintf('Síntese híbrida intervalar concluída com sucesso!\n');
        fprintf('Ganho K: [%.4f %.4f]\n', saidaHibInt.K(1), saidaHibInt.K(2));

    catch ME
        fprintf('ERRO na síntese híbrida intervalar: %s\n', ME.message);
        hibridoIntervalar.sucesso = false;
        hibridoIntervalar.erro = ME.message;
    end

    tempo = datetime('now') - tempoRef;
    fprintf('Tempo de execução: %.2f segundos\n\n', seconds(tempo));

    %% 3. SÍNTESE HÍBRIDA POLITÓPICA
    fprintf('=== MÉTODO 3: SÍNTESE HÍBRIDA POLITÓPICA ===\n');
    tempoRef = datetime('now');
    try
        saidaHibPoly = estHInfSintPolyLMILab(poly, h, delta, tol);

        % Avaliação do desempenho
        KCal_hibPoly = [eye(nx) zeros(nx,nu); saidaHibPoly.K zeros(nu)];
        raioEspInf_hibPoly = verificaEstSisHib(ACal.inf, KCal_hibPoly, h);
        raioEspSup_hibPoly = verificaEstSisHib(ACal.sup, KCal_hibPoly, h);
        normInf_hibPoly = normaSistemaContinuo(A.inf+B.inf*saidaHibPoly.K, E.inf, C.inf+D.inf*saidaHibPoly.K, zeros(size(C.inf,1),size(E.inf,2)));
        normSup_hibPoly = normaSistemaContinuo(A.sup+B.sup*saidaHibPoly.K, E.sup, C.sup+D.sup*saidaHibPoly.K, zeros(size(C.sup,1),size(E.sup,2)));

        % Estrutura de resultados
        hibridoPolitopico.sintese = saidaHibPoly;
        hibridoPolitopico.K = saidaHibPoly.K;
        hibridoPolitopico.raioEspectralInf = raioEspInf_hibPoly;
        hibridoPolitopico.raioEspectralSup = raioEspSup_hibPoly;
        hibridoPolitopico.normaSistemaInf = normInf_hibPoly;
        hibridoPolitopico.normaSistemaSup = normSup_hibPoly;
        hibridoPolitopico.sucesso = true;

        fprintf('Síntese híbrida politópica concluída com sucesso!\n');
        fprintf('Ganho K: [%.4f %.4f]\n', saidaHibPoly.K(1), saidaHibPoly.K(2));

    catch ME
        fprintf('ERRO na síntese híbrida politópica: %s\n', ME.message);
        hibridoPolitopico.sucesso = false;
        hibridoPolitopico.erro = ME.message;
    end

    tempo = datetime('now') - tempoRef;
    fprintf('Tempo de execução: %.2f segundos\n\n', seconds(tempo));

    %% 4. CONTROLE DISCRETO INTERVALAR (GANHO FIXO)
    fprintf('=== MÉTODO 4: CONTROLE DISCRETO INTERVALAR (GANHO FIXO) ===\n');
    tempoRef = datetime('now');
    try
        % Ganho fixo da literatura
        K_discInt = [10.8005 24.0549];

        % Avaliação do desempenho
        KCal_discInt = [eye(nx) zeros(nx,nu); K_discInt zeros(nu)];
        raioEspInf_discInt = verificaEstSisHib(ACal.inf, KCal_discInt, h);
        raioEspSup_discInt = verificaEstSisHib(ACal.sup, KCal_discInt, h);
        normInf_discInt = normaSistemaContinuo(A.inf+B.inf*K_discInt, E.inf, C.inf+D.inf*K_discInt, zeros(size(C.inf,1),size(E.inf,2)));
        normSup_discInt = normaSistemaContinuo(A.sup+B.sup*K_discInt, E.sup, C.sup+D.sup*K_discInt, zeros(size(C.sup,1),size(E.sup,2)));

        % Estrutura de resultados
        discretoIntervalar.K = K_discInt;
        discretoIntervalar.raioEspectralInf = raioEspInf_discInt;
        discretoIntervalar.raioEspectralSup = raioEspSup_discInt;
        discretoIntervalar.normaSistemaInf = normInf_discInt;
        discretoIntervalar.normaSistemaSup = normSup_discInt;
        discretoIntervalar.sucesso = true;
        discretoIntervalar.tipo = 'Ganho fixo da literatura';

        fprintf('Avaliação do controle discreto intervalar concluída!\n');
        fprintf('Ganho K: [%.4f %.4f]\n', K_discInt(1), K_discInt(2));

    catch ME
        fprintf('ERRO na avaliação do controle discreto intervalar: %s\n', ME.message);
        discretoIntervalar.sucesso = false;
        discretoIntervalar.erro = ME.message;
    end

    tempo = datetime('now') - tempoRef;
    fprintf('Tempo de execução: %.2f segundos\n\n', seconds(tempo));

    %% 5. CONTROLE DISCRETO POLITÓPICO (GANHO FIXO)
    fprintf('=== MÉTODO 5: CONTROLE DISCRETO POLITÓPICO (GANHO FIXO) ===\n');
    tempoRef = datetime('now');
    try
        % Ganho fixo do artigo
        K_discPoly = [9.7260 22.8904];

        % Avaliação do desempenho
        KCal_discPoly = [eye(nx) zeros(nx,nu); K_discPoly zeros(nu)];
        raioEspInf_discPoly = verificaEstSisHib(ACal.inf, KCal_discPoly, h);
        raioEspSup_discPoly = verificaEstSisHib(ACal.sup, KCal_discPoly, h);
        normInf_discPoly = normaSistemaContinuo(A.inf+B.inf*K_discPoly, E.inf, C.inf+D.inf*K_discPoly, zeros(size(C.inf,1),size(E.inf,2)));
        normSup_discPoly = normaSistemaContinuo(A.sup+B.sup*K_discPoly, E.sup, C.sup+D.sup*K_discPoly, zeros(size(C.sup,1),size(E.sup,2)));

        % Estrutura de resultados
        discretoPolitopico.K = K_discPoly;
        discretoPolitopico.raioEspectralInf = raioEspInf_discPoly;
        discretoPolitopico.raioEspectralSup = raioEspSup_discPoly;
        discretoPolitopico.normaSistemaInf = normInf_discPoly;
        discretoPolitopico.normaSistemaSup = normSup_discPoly;
        discretoPolitopico.sucesso = true;
        discretoPolitopico.tipo = 'Ganho fixo do artigo';

        fprintf('Avaliação do controle discreto politópico concluída!\n');
        fprintf('Ganho K: [%.4f %.4f]\n', K_discPoly(1), K_discPoly(2));

    catch ME
        fprintf('ERRO na avaliação do controle discreto politópico: %s\n', ME.message);
        discretoPolitopico.sucesso = false;
        discretoPolitopico.erro = ME.message;
    end

    tempo = datetime('now') - tempoRef;
    fprintf('Tempo de execução: %.2f segundos\n\n', seconds(tempo));

    %% Salva todos os resultados
    fprintf('=== SALVANDO RESULTADOS ===\n');

    saida.sistema = paraSys.sys;
    saida.parametros.h = h;
    saida.parametros.delta = delta;
    saida.parametros.tol = tol;

    % Resultados dos 5 métodos
    saida.metodos.nominal = nominal;
    saida.metodos.hibridoIntervalar = hibridoIntervalar;
    saida.metodos.hibridoPolitopico = hibridoPolitopico;
    saida.metodos.discretoIntervalar = discretoIntervalar;
    saida.metodos.discretoPolitopico = discretoPolitopico;

    % Resumo comparativo
    fprintf('=== RESUMO COMPARATIVO ===\n');
    metodosNomes = {'Nominal', 'Híb.Int', 'Híb.Poly', 'Disc.Int', 'Disc.Poly'};
    metodosStruct = {nominal, hibridoIntervalar, hibridoPolitopico, discretoIntervalar, discretoPolitopico};

    fprintf('Método           | Ganho K1     | Ganho K2     | Raio Esp.(Inf) | Raio Esp.(Sup) | Sucesso\n');
    fprintf('-----------------|--------------|--------------|----------------|----------------|---------\n');

    for i = 1:length(metodosNomes)
        if metodosStruct{i}.sucesso
            fprintf('%-16s | %10.4f | %10.4f | %12.6f | %12.6f | %s\n', ...
                metodosNomes{i}, metodosStruct{i}.K(1), metodosStruct{i}.K(2), ...
                metodosStruct{i}.raioEspectralInf, metodosStruct{i}.raioEspectralSup, 'OK');
        else
            fprintf('%-16s | %10s | %10s | %12s | %12s | ERRO\n', ...
                metodosNomes{i}, 'N/A', 'N/A', 'N/A', 'N/A');
        end
    end

    % Tempo total
    tempoTotal = datetime('now') - tempoInicial;
    fprintf('\nTempo total de execução: %.2f segundos\n', seconds(tempoTotal));
    saida.tempoTotal = seconds(tempoTotal);

    fprintf('\n=== COMPARAÇÃO FINALIZADA ===\n');
    fprintf('Todos os resultados salvos na estrutura "saida"\n');

end