function saida = genSaveDataComparacao()
%   Sistema Mass-Spring-Damper para comparação de controladores
%   Baseado no sistema que já funciona na pasta SampledInput
%
%   Sistema incerto com:
%   - Incertezas: ±5% no amortecedor, ±2% na mola
%   - h = 0.02 (período de amostragem)
%   - Sistema nominal com valores centrais

    % Adicionar caminhos necessários
    addpath('../HInf - Análise - Intervalar/funcoes');
    addpath('../HInf - Análise - Intervalar/funcoes/Diversas');
    addpath('../funcoes');

    intvalinit('DisplayInfSup');

    % Parâmetros nominais do sistema massa-mola-amortecedor
    m1 = 0.5;
    m1 = infsup(m1, m1);

    m2 = 1;
    m2 = infsup(m2, m2);

    % Amortecedor com incerteza ±5%
    b = 0.2;
    b = infsup(b*0.95, b*1.05);

    k1 = 12;
    k1 = infsup(k1, k1);

    % Mola com incerteza ±2%
    k2 = 7;
    k2 = infsup(k2*0.98, k2*1.02);

    % Matriz A do sistema
    A = [0              0       1       0;
         0              0       0       1;
         (-k2-k1)/m1    k2/m1   -b/m1   b/m1;
         k2/m2          -k2/m2  b/m2    -b/m2];

    % Matriz B (entrada de controle)
    B = [0 0 0 1/m2]';

    % Matriz C (saída)
    C = [0 10 0 0;
         0  0 0 1;
         0  0 0 0];
    C = infsup(C, C);

    % Matriz D (feedthrough da entrada de controle)
    D = [0 0 1]';
    D = infsup(D, D);

    % Matriz E (entrada de perturbação)
    E = [0 0 1/m1 0]';

    % Tamanho das matrizes
    nx = length(A.inf);
    nu = size(B.inf, 2);
    nw = size(E.inf, 2);

    % Matrizes caligráficas
    ACal = [A B; zeros(nu, nx) zeros(nu, nu)];
    ECal = [E; zeros(nu, nw)];
    CCal = [C D];

    % Gerar politopo manualmente
    b_vals = [b.inf b.sup];
    k2_vals = [k2.inf k2.sup];

    idx = 0;
    for i = 1:length(b_vals)
        for j = 1:length(k2_vals)
            idx = idx + 1;

            sysPolyCont{idx}.A = [0              0       1       0;
                                  0              0       0       1;
                                  (-k2_vals(j)-k1.sup)/m1.sup    k2_vals(j)/m1.sup   -b_vals(i)/m1.sup   b_vals(i)/m1.sup;
                                  k2_vals(j)/m2.sup          -k2_vals(j)/m2.sup  b_vals(i)/m2.sup    -b_vals(i)/m2.sup];

            sysPolyCont{idx}.B2 = [0 0 0 1/m2.sup]';
            sysPolyCont{idx}.B1 = [0 0 1/m1.sup 0]';
            sysPolyCont{idx}.C = [0 10 0 0;
                                  0  0 0 1;
                                  0  0 0 0];
            sysPolyCont{idx}.D2 = [0 0 1]';
            sysPolyCont{idx}.D1 = [0 0 0]';
        end
    end

    % Parâmetros do sistema
    h = 0.02;  % Período de amostragem
    qtdDiv = 5;
    delta = h/qtdDiv;

    % Parâmetros para matrizes politópicas
    numPointsUniSpaced = 5;
    numPointsBy2Points = 2;
    numbPointsUniSpacedSub = 2;
    onlyVertice = 1;

    % Tolerância
    tol = 1e-5;

    % Sistema nominal (valores centrais)
    A_nominal = [0              0       1       0;
                 0              0       0       1;
                 (-7-12)/0.5    7/0.5   -0.2/0.5   0.2/0.5;
                 7/1            -7/1    0.2/1      -0.2/1];

    B_nominal = [0 0 0 1/1]';
    E_nominal = [0 0 1/0.5 0]';
    C_nominal = [0 10 0 0;
                 0  0 0 1;
                 0  0 0 0];
    D_nominal = [0 0 1]';

    % Matrizes caligráficas nominais
    ACal_nominal = [A_nominal B_nominal; zeros(nu, nx) zeros(nu, nu)];
    ECal_nominal = [E_nominal; zeros(nu, nw)];
    CCal_nominal = [C_nominal D_nominal];

    % Criar pasta se não existir
    if(exist(mfilename) ~= 7)
        mkdir(mfilename);
    end

    % Sistema intervalar
    sys.A = A;
    sys.B = B;
    sys.E = E;
    sys.C = C;
    sys.D = D;
    sys.ACal = ACal;
    sys.ECal = ECal;
    sys.CCal = CCal;
    sys.h = h;
    sys.delta = delta;
    sys.sysPolyCont = sysPolyCont;

    % Sistema nominal
    sysNominal.A = A_nominal;
    sysNominal.B = B_nominal;
    sysNominal.E = E_nominal;
    sysNominal.C = C_nominal;
    sysNominal.D = D_nominal;
    sysNominal.ACal = ACal_nominal;
    sysNominal.ECal = ECal_nominal;
    sysNominal.CCal = CCal_nominal;
    sysNominal.h = h;
    sysNominal.delta = delta;

    % Parâmetros das incertezas (para referência)
    parametros.m1 = m1;
    parametros.m2 = m2;
    parametros.b = b;
    parametros.k1 = k1;
    parametros.k2 = k2;
    parametros.h = h;
    parametros.delta = delta;
    parametros.tol = tol;

    % Estrutura de saída
    saida.sys = sys;
    saida.sysNominal = sysNominal;
    saida.parametros = parametros;

    % Parâmetros auxiliares
    saida.aux.tol = tol;

    % Simulação
    sim.numPointsUniSpaced = numPointsUniSpaced;
    sim.numPointsBy2Points = numPointsBy2Points;
    sim.numbPointsUniSpacedSub = numbPointsUniSpacedSub;
    sim.onlyVertice = onlyVertice;
    saida.simSys = sim;

    fprintf('Sistema mass-spring-damper carregado:\n');
    fprintf('  - Estados: %d\n', nx);
    fprintf('  - Entradas de controle: %d\n', nu);
    fprintf('  - Entradas de perturbação: %d\n', nw);
    fprintf('  - Saídas: %d\n', size(C.inf, 1));
    fprintf('  - Vértices do politopo: %d\n', length(sysPolyCont));
    fprintf('  - Período de amostragem: %.3f s\n', h);
    fprintf('  - Incertezas: b ± %.1f%%, k2 ± %.1f%%\n', 5, 2);

end