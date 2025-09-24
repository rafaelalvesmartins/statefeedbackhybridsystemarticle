function main_massa_mola_completo()
%   MAIN_MASSA_MOLA_COMPLETO - Síntese completa de controladores H∞
%
%   Sistema massa-mola-amortecedor baseado no Exemplo 12 (genSaveDataEx12)
%   Seguindo a lógica de G.W. Gabriel, J.C. Geromel & K.M. Grigoriadis
%
%   ETAPA 1: Síntese Nominal (sem incertezas)
%   ETAPA 2: Síntese Robusta (com incertezas)
%
%   Métodos utilizados:
%   - Abordagem Intervalar: estHInfSintLMILab
%   - Abordagem Politópica: estHInfSintPolyLMILab
%   - Abordagem Nominal: estHInfSintLMILabPrec
%
%   Autor: Baseado em genSaveDataEx12
%   Data: 2025

    clear; clc;
    fprintf('=== SÍNTESE COMPLETA DE CONTROLADORES H∞ - SISTEMA MASSA-MOLA ===\n\n');

    % Adicionar caminhos necessários
    addpath('funcoes');
    addpath('functions');
    addpath('functions/Optimize');
    addpath('HInf - Análise - Intervalar/funcoes');
    addpath('HInf - Análise - Intervalar/funcoes/Diversas');

    % Inicializar INTLAB
    intvalinit('DisplayInfSup');

    %% =================================================================
    %% DEFINIÇÃO DO SISTEMA (baseado em genSaveDataEx12)
    %% =================================================================

    % Parâmetros físicos seguindo genSaveDataEx12
    m1 = 0.5;  % massa 1 [kg]
    m2 = 1.0;  % massa 2 [kg]
    k1 = 12;   % constante da mola 1 [N/m]
    k2 = 7;    % constante da mola 2 [N/m]
    b = 0.2;   % coeficiente de amortecimento [N.s/m]

    % Parâmetros de síntese
    h = 0.02;       % Período de amostragem (como no exemplo)
    qtdDiv = 5;
    delta = h/qtdDiv;
    tol = 1e-5;     % Tolerância

    fprintf('Parâmetros do sistema (baseado em genSaveDataEx12):\n');
    fprintf('  Massa 1 (m1): %.1f kg\n', m1);
    fprintf('  Massa 2 (m2): %.1f kg\n', m2);
    fprintf('  Mola 1 (k1): %.0f N/m\n', k1);
    fprintf('  Mola 2 (k2): %.0f N/m\n', k2);
    fprintf('  Amortecimento (b): %.1f N.s/m\n', b);
    fprintf('  Período de amostragem (h): %.3f s\n', h);
    fprintf('  Delta: %.4f\n', delta);
    fprintf('  Tolerância: %.0e\n\n', tol);

    %% =================================================================
    %% ETAPA 1: SÍNTESE NOMINAL (SEM INCERTEZAS)
    %% =================================================================

    fprintf('=====================================\n');
    fprintf('ETAPA 1: SÍNTESE NOMINAL (SEM INCERTEZAS)\n');
    fprintf('=====================================\n\n');

    % Matrizes intervalares nominais (sem incerteza)
    m1_nom = infsup(m1, m1);
    m2_nom = infsup(m2, m2);
    k1_nom = infsup(k1, k1);
    k2_nom = infsup(k2, k2);
    b_nom = infsup(b, b);

    % Sistema nominal seguindo estrutura do genSaveDataEx12
    A_nom = [0              0       1       0;
             0              0       0       1;
             (-k2_nom-k1_nom)/m1_nom    k2_nom/m1_nom   -b_nom/m1_nom   b_nom/m1_nom;
             k2_nom/m2_nom          -k2_nom/m2_nom  b_nom/m2_nom    -b_nom/m2_nom];

    B_nom = [0; 0; 0; 1/m2_nom];  % Entrada de controle na massa 2

    E_nom = [0; 0; 1/m1_nom; 0];  % Entrada de perturbação na massa 1

    C_nom = [0 10 0 0;            % Saída 1: 10x posição massa 1
             0  0 0 1;            % Saída 2: velocidade massa 2
             0  0 0 0];           % Saída 3: zero

    C_nom = infsup(C_nom, C_nom);

    D_nom = [0; 0; 1];            % Feedthrough
    D_nom = infsup(D_nom, D_nom);

    fprintf('Matrizes nominais calculadas (estrutura do genSaveDataEx12).\n\n');

    %% MÉTODO 1: Abordagem Intervalar Nominal
    fprintf('--- MÉTODO 1: Abordagem Intervalar (Nominal) ---\n');
    try
        tic;
        resultado_int_nom = estHInfSintLMILab(A_nom, B_nom, E_nom, C_nom, D_nom, h, delta, tol);
        tempo_int_nom = toc;

        if isfield(resultado_int_nom, 'factivel') && resultado_int_nom.factivel
            K_int_nom = resultado_int_nom.K;
            gamma_int_nom = resultado_int_nom.gamma;
            sucesso_int_nom = true;

            fprintf('✓ Síntese bem-sucedida!\n');
            fprintf('  Ganho K: [%.6f %.6f %.6f %.6f]\n', K_int_nom(1), K_int_nom(2), K_int_nom(3), K_int_nom(4));
            fprintf('  Norma γ: %.6f\n', gamma_int_nom);
            fprintf('  Tempo: %.4f s\n\n', tempo_int_nom);
        else
            sucesso_int_nom = false;
            fprintf('✗ Síntese falhou!\n\n');
        end

    catch ME
        sucesso_int_nom = false;
        fprintf('✗ Erro: %s\n\n', ME.message);
    end

    %% MÉTODO 2: Abordagem Politópica Nominal
    fprintf('--- MÉTODO 2: Abordagem Politópica (Nominal) ---\n');
    try
        % Crear politopo nominal com um único vértice
        poly_nom.APoly = {[0              0       1       0;
                           0              0       0       1;
                           (-k2-k1)/m1    k2/m1   -b/m1   b/m1;
                           k2/m2          -k2/m2  b/m2    -b/m2]};

        poly_nom.BPoly = {[0; 0; 0; 1/m2]};
        poly_nom.EPoly = {[0; 0; 1/m1; 0]};
        poly_nom.CPoly = {[0 10 0 0; 0 0 0 1; 0 0 0 0]};
        poly_nom.DPoly = {[0; 0; 1]};

        tic;
        resultado_poly_nom = estHInfSintPolyLMILab(poly_nom, h, delta, tol);
        tempo_poly_nom = toc;

        if isfield(resultado_poly_nom, 'factivel') && resultado_poly_nom.factivel
            K_poly_nom = resultado_poly_nom.K;
            gamma_poly_nom = resultado_poly_nom.gamma;
            sucesso_poly_nom = true;

            fprintf('✓ Síntese bem-sucedida!\n');
            fprintf('  Ganho K: [%.6f %.6f %.6f %.6f]\n', K_poly_nom(1), K_poly_nom(2), K_poly_nom(3), K_poly_nom(4));
            fprintf('  Norma γ: %.6f\n', gamma_poly_nom);
            fprintf('  Tempo: %.4f s\n\n', tempo_poly_nom);
        else
            sucesso_poly_nom = false;
            fprintf('✗ Síntese falhou!\n\n');
        end

    catch ME
        sucesso_poly_nom = false;
        fprintf('✗ Erro: %s\n\n', ME.message);
    end

    %% MÉTODO 3: Abordagem Nominal Clássica
    fprintf('--- MÉTODO 3: Abordagem Nominal (estHInfSintLMILabPrec) ---\n');
    try
        % Valores centrais para abordagem nominal
        A_nom_central = [0              0       1       0;
                         0              0       0       1;
                         (-k2-k1)/m1    k2/m1   -b/m1   b/m1;
                         k2/m2          -k2/m2  b/m2    -b/m2];

        B_nom_central = [0; 0; 0; 1/m2];
        E_nom_central = [0; 0; 1/m1; 0];
        C_nom_central = [0 10 0 0; 0 0 0 1; 0 0 0 0];
        D_nom_central = [0; 0; 1];

        tic;
        resultado_nom = estHInfSintLMILabPrec(A_nom_central, B_nom_central, E_nom_central, C_nom_central, D_nom_central, h, delta, tol);
        tempo_nom = toc;

        if isfield(resultado_nom, 'factivel') && resultado_nom.factivel
            K_nom = resultado_nom.K;
            gamma_nom = resultado_nom.gamma;
            sucesso_nom = true;

            fprintf('✓ Síntese bem-sucedida!\n');
            fprintf('  Ganho K: [%.6f %.6f %.6f %.6f]\n', K_nom(1), K_nom(2), K_nom(3), K_nom(4));
            fprintf('  Norma γ: %.6f\n', gamma_nom);
            fprintf('  Tempo: %.4f s\n\n', tempo_nom);
        else
            sucesso_nom = false;
            fprintf('✗ Síntese falhou!\n\n');
        end

    catch ME
        sucesso_nom = false;
        fprintf('✗ Erro: %s\n\n', ME.message);
    end

    %% =================================================================
    %% ETAPA 2: SÍNTESE ROBUSTA (COM INCERTEZAS)
    %% =================================================================

    fprintf('=====================================\n');
    fprintf('ETAPA 2: SÍNTESE ROBUSTA (COM INCERTEZAS)\n');
    fprintf('=====================================\n\n');

    % Definir incertezas seguindo genSaveDataEx12
    m1_unc = infsup(m1, m1);                    % sem incerteza em m1
    m2_unc = infsup(m2, m2);                    % sem incerteza em m2
    b_unc = infsup(b*0.95, b*1.05);            % b com ±5% de incerteza
    k1_unc = infsup(k1, k1);                   % sem incerteza em k1
    k2_unc = infsup(k2*0.98, k2*1.02);        % k2 com ±2% de incerteza

    fprintf('Parâmetros com incertezas (seguindo genSaveDataEx12):\n');
    fprintf('  Massa 1 (m1): [%.2f, %.2f] kg (sem incerteza)\n', m1_unc.inf, m1_unc.sup);
    fprintf('  Massa 2 (m2): [%.2f, %.2f] kg (sem incerteza)\n', m2_unc.inf, m2_unc.sup);
    fprintf('  Amortecimento (b): [%.3f, %.3f] N.s/m (±5%%)\n', b_unc.inf, b_unc.sup);
    fprintf('  Mola 1 (k1): [%.1f, %.1f] N/m (sem incerteza)\n', k1_unc.inf, k1_unc.sup);
    fprintf('  Mola 2 (k2): [%.2f, %.2f] N/m (±2%%)\n\n', k2_unc.inf, k2_unc.sup);

    % Sistema com incertezas
    A_unc = [0                          0              1              0;
             0                          0              0              1;
             (-k2_unc-k1_unc)/m1_unc   k2_unc/m1_unc  -b_unc/m1_unc  b_unc/m1_unc;
             k2_unc/m2_unc             -k2_unc/m2_unc  b_unc/m2_unc  -b_unc/m2_unc];

    B_unc = [0; 0; 0; 1/m2_unc];
    E_unc = [0; 0; 1/m1_unc; 0];
    C_unc = infsup(C_nom_central, C_nom_central);  % sem incerteza
    D_unc = infsup(D_nom_central, D_nom_central);  % sem incerteza

    fprintf('Matrizes intervalares com incertezas calculadas.\n\n');

    %% MÉTODO 4: Abordagem Intervalar Robusta
    fprintf('--- MÉTODO 4: Abordagem Intervalar (Robusta) ---\n');
    try
        tic;
        resultado_int_rob = estHInfSintLMILab(A_unc, B_unc, E_unc, C_unc, D_unc, h, delta, tol);
        tempo_int_rob = toc;

        if isfield(resultado_int_rob, 'factivel') && resultado_int_rob.factivel
            K_int_rob = resultado_int_rob.K;
            gamma_int_rob = resultado_int_rob.gamma;
            sucesso_int_rob = true;

            fprintf('✓ Síntese bem-sucedida!\n');
            fprintf('  Ganho K: [%.6f %.6f %.6f %.6f]\n', K_int_rob(1), K_int_rob(2), K_int_rob(3), K_int_rob(4));
            fprintf('  Norma γ: %.6f\n', gamma_int_rob);
            fprintf('  Tempo: %.4f s\n\n', tempo_int_rob);
        else
            sucesso_int_rob = false;
            fprintf('✗ Síntese falhou!\n\n');
        end

    catch ME
        sucesso_int_rob = false;
        fprintf('✗ Erro: %s\n\n', ME.message);
    end

    %% MÉTODO 5: Abordagem Politópica Robusta
    fprintf('--- MÉTODO 5: Abordagem Politópica (Robusta) ---\n');
    try
        % Gerar vértices do politopo seguindo genSaveDataEx12
        b_vals = [b_unc.inf, b_unc.sup];
        k2_vals = [k2_unc.inf, k2_unc.sup];

        % Politopo com vértices das incertezas
        idx = 0;
        for i = 1:length(b_vals)
            for j = 1:length(k2_vals)
                idx = idx + 1;

                poly_rob.APoly{idx} = [0                                  0                 1                  0;
                                       0                                  0                 0                  1;
                                       (-k2_vals(j)-k1)/m1               k2_vals(j)/m1    -b_vals(i)/m1      b_vals(i)/m1;
                                       k2_vals(j)/m2                     -k2_vals(j)/m2    b_vals(i)/m2      -b_vals(i)/m2];

                poly_rob.BPoly{idx} = [0; 0; 0; 1/m2];    % B2
                poly_rob.EPoly{idx} = [0; 0; 1/m1; 0];    % B1
                poly_rob.CPoly{idx} = [0 10 0 0; 0 0 0 1; 0 0 0 0];
                poly_rob.DPoly{idx} = [0; 0; 1];          % D2
            end
        end

        fprintf('Politopo gerado com %d vértices (seguindo genSaveDataEx12).\n', idx);

        tic;
        resultado_poly_rob = estHInfSintPolyLMILab(poly_rob, h, delta, tol);
        tempo_poly_rob = toc;

        if isfield(resultado_poly_rob, 'factivel') && resultado_poly_rob.factivel
            K_poly_rob = resultado_poly_rob.K;
            gamma_poly_rob = resultado_poly_rob.gamma;
            sucesso_poly_rob = true;

            fprintf('✓ Síntese bem-sucedida!\n');
            fprintf('  Ganho K: [%.6f %.6f %.6f %.6f]\n', K_poly_rob(1), K_poly_rob(2), K_poly_rob(3), K_poly_rob(4));
            fprintf('  Norma γ: %.6f\n', gamma_poly_rob);
            fprintf('  Tempo: %.4f s\n\n', tempo_poly_rob);
        else
            sucesso_poly_rob = false;
            fprintf('✗ Síntese falhou!\n\n');
        end

    catch ME
        sucesso_poly_rob = false;
        fprintf('✗ Erro: %s\n\n', ME.message);
    end

    %% MÉTODO 6: Abordagem Nominal com Sistema Central
    fprintf('--- MÉTODO 6: Abordagem Nominal (Sistema Central das Incertezas) ---\n');
    try
        % Sistema central das incertezas
        b_central = (b_unc.inf + b_unc.sup) / 2;
        k2_central = (k2_unc.inf + k2_unc.sup) / 2;

        A_central = [0                            0                   1                   0;
                     0                            0                   0                   1;
                     (-k2_central-k1)/m1         k2_central/m1       -b_central/m1       b_central/m1;
                     k2_central/m2               -k2_central/m2       b_central/m2       -b_central/m2];

        B_central = [0; 0; 0; 1/m2];
        E_central = [0; 0; 1/m1; 0];

        tic;
        resultado_central = estHInfSintLMILabPrec(A_central, B_central, E_central, C_nom_central, D_nom_central, h, delta, tol);
        tempo_central = toc;

        if isfield(resultado_central, 'factivel') && resultado_central.factivel
            K_central = resultado_central.K;
            gamma_central = resultado_central.gamma;
            sucesso_central = true;

            fprintf('✓ Síntese bem-sucedida!\n');
            fprintf('  Ganho K: [%.6f %.6f %.6f %.6f]\n', K_central(1), K_central(2), K_central(3), K_central(4));
            fprintf('  Norma γ: %.6f\n', gamma_central);
            fprintf('  Tempo: %.4f s\n\n', tempo_central);
        else
            sucesso_central = false;
            fprintf('✗ Síntese falhou!\n\n');
        end

    catch ME
        sucesso_central = false;
        fprintf('✗ Erro: %s\n\n', ME.message);
    end

    %% =================================================================
    %% ETAPA 3: SÍNTESE CONTÍNUA (COM INCERTEZAS)
    %% =================================================================

    fprintf('=====================================\n');
    fprintf('ETAPA 3: SÍNTESE CONTÍNUA (COM INCERTEZAS)\n');
    fprintf('=====================================\n');
    fprintf('Usando funções: synHInfKContIntLMILab e synHInfKContPoly\n\n');

    % Preparar matrizes e parâmetros para síntese contínua
    fprintf('Preparando matrizes para síntese contínua...\n');

    % Matrizes para síntese contínua (usando sistema com incertezas)
    A_cont_unc = A_unc;
    B2_cont_unc = B_unc;          % Entrada de controle (B2)
    B1_cont_unc = E_unc;          % Entrada de perturbação (B1)
    C_cont_unc = C_unc;
    D2_cont_unc = D_unc;          % Feedthrough do controle (D2)

    % Criar matriz D1 (feedthrough da perturbação) com dimensões corretas
    num_saidas = size(C_cont_unc.inf, 1);        % número de saídas
    num_perturbacoes = size(B1_cont_unc.inf, 2); % número de perturbações
    D1_cont_unc = infsup(zeros(num_saidas, num_perturbacoes), zeros(num_saidas, num_perturbacoes));

    % Inicializar estrutura de parâmetros corretamente
    param_cont = struct();
    param_cont.tol = tol;
    param_cont.maxNormOfK = 1e8;  % Campo necessário para synHInfKContPoly

    % Inicializar variáveis de sucesso e resultados
    sucesso_cont_int_rob = false;
    sucesso_cont_poly_rob = false;
    K_cont_int_rob = [];
    gamma_cont_int_rob = NaN;
    tempo_cont_int_rob = 0;
    K_cont_poly_rob = [];
    gamma_cont_poly_rob = NaN;
    tempo_cont_poly_rob = 0;

    % Variáveis para métodos alternativos contínuos (versões não-LMILab)
    sucesso_cont_int_alt = false;
    sucesso_cont_poly_alt = false;
    K_cont_int_alt = [];
    gamma_cont_int_alt = NaN;
    tempo_cont_int_alt = 0;
    K_cont_poly_alt = [];
    gamma_cont_poly_alt = NaN;
    tempo_cont_poly_alt = 0;

    fprintf('Dimensões das matrizes contínuas:\n');
    fprintf('  A: %dx%d, B2: %dx%d, B1: %dx%d\n', size(A_cont_unc.inf,1), size(A_cont_unc.inf,2), size(B2_cont_unc.inf,1), size(B2_cont_unc.inf,2), size(B1_cont_unc.inf,1), size(B1_cont_unc.inf,2));
    fprintf('  C: %dx%d, D2: %dx%d, D1: %dx%d\n\n', size(C_cont_unc.inf,1), size(C_cont_unc.inf,2), size(D2_cont_unc.inf,1), size(D2_cont_unc.inf,2), size(D1_cont_unc.inf,1), size(D1_cont_unc.inf,2));

    %% MÉTODO 7: Abordagem Intervalar Contínua (Robusta)
    fprintf('--- MÉTODO 7: Abordagem Intervalar Contínua (Robusta) ---\n');
    try

        tic;
        resultado_cont_int_rob = synHInfKContIntLMILab(A_cont_unc, B2_cont_unc, B1_cont_unc, C_cont_unc, D2_cont_unc, D1_cont_unc, param_cont);
        tempo_cont_int_rob = toc;

        if isfield(resultado_cont_int_rob, 'feas') && resultado_cont_int_rob.feas == 1
            K_cont_int_rob = resultado_cont_int_rob.K;
            gamma_cont_int_rob = resultado_cont_int_rob.mu;
            sucesso_cont_int_rob = true;

            fprintf('✓ Síntese bem-sucedida!\n');
            fprintf('  Ganho K: [%.6f %.6f %.6f %.6f]\n', K_cont_int_rob(1), K_cont_int_rob(2), K_cont_int_rob(3), K_cont_int_rob(4));
            fprintf('  Norma γ: %.6f\n', gamma_cont_int_rob);
            fprintf('  Tempo: %.4f s\n\n', tempo_cont_int_rob);
        else
            sucesso_cont_int_rob = false;
            fprintf('✗ Síntese falhou!\n\n');
        end

    catch ME
        sucesso_cont_int_rob = false;
        fprintf('✗ Erro: %s\n\n', ME.message);
    end

    %% MÉTODO 8: Abordagem Politópica Contínua (Robusta)
    fprintf('--- MÉTODO 8: Abordagem Politópica Contínua (Robusta) ---\n');
    try
        % Preparar politopo para síntese contínua
        fprintf('Preparando politopo contínuo...\n');

        % Usar o mesmo politopo robusto, mas adicionar D1Poly (feedthrough da perturbação)
        poly_cont_rob = poly_rob;  % Copiar estrutura existente

        % Adicionar D1Poly (feedthrough da perturbação) para cada vértice
        % Dimensões: num_saidas x num_perturbacoes
        d1_matriz = zeros(num_saidas, num_perturbacoes);
        for idx_vert = 1:length(poly_cont_rob.APoly)
            poly_cont_rob.D1Poly{idx_vert} = d1_matriz;  % D1 = zeros com dimensões corretas
        end

        fprintf('Politopo contínuo preparado com %d vértices.\n', length(poly_cont_rob.APoly));

        tic;
        resultado_cont_poly_rob = synHInfKContPoly(poly_cont_rob.APoly, poly_cont_rob.BPoly, poly_cont_rob.EPoly, poly_cont_rob.CPoly, poly_cont_rob.DPoly, poly_cont_rob.D1Poly, param_cont);
        tempo_cont_poly_rob = toc;

        if isfield(resultado_cont_poly_rob, 'feas') && resultado_cont_poly_rob.feas == 1
            K_cont_poly_rob = resultado_cont_poly_rob.K;
            gamma_cont_poly_rob = resultado_cont_poly_rob.norm;
            sucesso_cont_poly_rob = true;

            fprintf('✓ Síntese bem-sucedida!\n');
            fprintf('  Ganho K: [%.6f %.6f %.6f %.6f]\n', K_cont_poly_rob(1), K_cont_poly_rob(2), K_cont_poly_rob(3), K_cont_poly_rob(4));
            fprintf('  Norma γ: %.6f\n', gamma_cont_poly_rob);
            fprintf('  Tempo: %.4f s\n\n', tempo_cont_poly_rob);
        else
            sucesso_cont_poly_rob = false;
            fprintf('✗ Síntese falhou!\n\n');
        end

    catch ME
        sucesso_cont_poly_rob = false;
        fprintf('✗ Erro: %s\n\n', ME.message);
    end

    %% MÉTODO 9: Abordagem Intervalar Contínua Alternativa (synHInfKContInt)
    fprintf('--- MÉTODO 9: Abordagem Intervalar Contínua Alternativa (synHInfKContInt) ---\n');
    try
        tic;
        resultado_cont_int_alt = synHInfKContInt(A_cont_unc, B2_cont_unc, B1_cont_unc, C_cont_unc, D2_cont_unc, D1_cont_unc, param_cont);
        tempo_cont_int_alt = toc;

        if isfield(resultado_cont_int_alt, 'feas') && resultado_cont_int_alt.feas == 1
            K_cont_int_alt = resultado_cont_int_alt.K;
            gamma_cont_int_alt = resultado_cont_int_alt.norm;
            sucesso_cont_int_alt = true;

            fprintf('✓ Síntese bem-sucedida!\n');
            fprintf('  Ganho K: [%.6f %.6f %.6f %.6f]\n', K_cont_int_alt(1), K_cont_int_alt(2), K_cont_int_alt(3), K_cont_int_alt(4));
            fprintf('  Norma γ: %.6f\n', gamma_cont_int_alt);
            fprintf('  Tempo: %.4f s\n\n', tempo_cont_int_alt);
        else
            sucesso_cont_int_alt = false;
            fprintf('✗ Síntese falhou!\n\n');
        end

    catch ME
        sucesso_cont_int_alt = false;
        fprintf('✗ Erro: %s\n\n', ME.message);
    end

    %% MÉTODO 10: Abordagem Politópica Contínua Alternativa (synHInfKContPoly)
    fprintf('--- MÉTODO 10: Abordagem Politópica Contínua Alternativa (synHInfKContPoly) ---\n');
    try
        tic;
        % Usar o mesmo politopo contínuo preparado anteriormente
        resultado_cont_poly_alt = synHInfKContPoly(poly_cont_rob.APoly, poly_cont_rob.BPoly, poly_cont_rob.EPoly, poly_cont_rob.CPoly, poly_cont_rob.DPoly, poly_cont_rob.D1Poly, param_cont);
        tempo_cont_poly_alt = toc;

        if isfield(resultado_cont_poly_alt, 'feas') && resultado_cont_poly_alt.feas == 1
            K_cont_poly_alt = resultado_cont_poly_alt.K;
            gamma_cont_poly_alt = resultado_cont_poly_alt.norm;
            sucesso_cont_poly_alt = true;

            fprintf('✓ Síntese bem-sucedida!\n');
            fprintf('  Ganho K: [%.6f %.6f %.6f %.6f]\n', K_cont_poly_alt(1), K_cont_poly_alt(2), K_cont_poly_alt(3), K_cont_poly_alt(4));
            fprintf('  Norma γ: %.6f\n', gamma_cont_poly_alt);
            fprintf('  Tempo: %.4f s\n\n', tempo_cont_poly_alt);
        else
            sucesso_cont_poly_alt = false;
            fprintf('✗ Síntese falhou!\n\n');
        end

    catch ME
        sucesso_cont_poly_alt = false;
        fprintf('✗ Erro: %s\n\n', ME.message);
    end

    %% =================================================================
    %% RESUMO DOS RESULTADOS
    %% =================================================================

    fprintf('=====================================\n');
    fprintf('RESUMO FINAL DOS RESULTADOS\n');
    fprintf('=====================================\n\n');

    % Tabela resumo
    fprintf('MÉTODO                           | STATUS | γ²        | TEMPO (s) | GANHO K (1x4)\n');
    fprintf('--------------------------------|--------|-----------|-----------|----------------------------------------\n');

    % Resultados nominais
    if sucesso_int_nom
        fprintf('1. Intervalar (Nominal)          | %-6s | %9.6f | %8.4f | [%.4f %.4f %.4f %.4f]\n', ...
            'OK', gamma_int_nom^2, tempo_int_nom, K_int_nom(1), K_int_nom(2), K_int_nom(3), K_int_nom(4));
    else
        fprintf('1. Intervalar (Nominal)          | %-6s | %9s | %8s | %s\n', ...
            'FALHOU', 'N/A', 'N/A', 'N/A');
    end

    if sucesso_poly_nom
        fprintf('2. Politópica (Nominal)          | %-6s | %9.6f | %8.4f | [%.4f %.4f %.4f %.4f]\n', ...
            'OK', gamma_poly_nom^2, tempo_poly_nom, K_poly_nom(1), K_poly_nom(2), K_poly_nom(3), K_poly_nom(4));
    else
        fprintf('2. Politópica (Nominal)          | %-6s | %9s | %8s | %s\n', ...
            'FALHOU', 'N/A', 'N/A', 'N/A');
    end

    if sucesso_nom
        fprintf('3. Nominal Clássica              | %-6s | %9.6f | %8.4f | [%.4f %.4f %.4f %.4f]\n', ...
            'OK', gamma_nom^2, tempo_nom, K_nom(1), K_nom(2), K_nom(3), K_nom(4));
    else
        fprintf('3. Nominal Clássica              | %-6s | %9s | %8s | %s\n', ...
            'FALHOU', 'N/A', 'N/A', 'N/A');
    end

    fprintf('--------------------------------|--------|-----------|-----------|----------------------------------------\n');

    % Resultados robustos
    if sucesso_int_rob
        fprintf('4. Intervalar (Robusta)          | %-6s | %9.6f | %8.4f | [%.4f %.4f %.4f %.4f]\n', ...
            'OK', gamma_int_rob^2, tempo_int_rob, K_int_rob(1), K_int_rob(2), K_int_rob(3), K_int_rob(4));
    else
        fprintf('4. Intervalar (Robusta)          | %-6s | %9s | %8s | %s\n', ...
            'FALHOU', 'N/A', 'N/A', 'N/A');
    end

    if sucesso_poly_rob
        fprintf('5. Politópica (Robusta)          | %-6s | %9.6f | %8.4f | [%.4f %.4f %.4f %.4f]\n', ...
            'OK', gamma_poly_rob^2, tempo_poly_rob, K_poly_rob(1), K_poly_rob(2), K_poly_rob(3), K_poly_rob(4));
    else
        fprintf('5. Politópica (Robusta)          | %-6s | %9s | %8s | %s\n', ...
            'FALHOU', 'N/A', 'N/A', 'N/A');
    end

    if sucesso_central
        fprintf('6. Nominal (Sistema Central)     | %-6s | %9.6f | %8.4f | [%.4f %.4f %.4f %.4f]\n', ...
            'OK', gamma_central^2, tempo_central, K_central(1), K_central(2), K_central(3), K_central(4));
    else
        fprintf('6. Nominal (Sistema Central)     | %-6s | %9s | %8s | %s\n', ...
            'FALHOU', 'N/A', 'N/A', 'N/A');
    end

    fprintf('--------------------------------|--------|-----------|-----------|----------------------------------------\n');

    % Resultados contínuos
    if sucesso_cont_int_rob && ~isempty(K_cont_int_rob) && length(K_cont_int_rob) >= 4
        fprintf('7. Intervalar (Contínua)         | %-6s | %9.6f | %8.4f | [%.4f %.4f %.4f %.4f]\n', ...
            'OK', gamma_cont_int_rob^2, tempo_cont_int_rob, K_cont_int_rob(1), K_cont_int_rob(2), K_cont_int_rob(3), K_cont_int_rob(4));
    else
        fprintf('7. Intervalar (Contínua)         | %-6s | %9s | %8s | %s\n', ...
            'FALHOU', 'N/A', 'N/A', 'N/A');
    end

    if sucesso_cont_poly_rob && ~isempty(K_cont_poly_rob) && length(K_cont_poly_rob) >= 4
        fprintf('8. Politópica (Contínua)         | %-6s | %9.6f | %8.4f | [%.4f %.4f %.4f %.4f]\n', ...
            'OK', gamma_cont_poly_rob^2, tempo_cont_poly_rob, K_cont_poly_rob(1), K_cont_poly_rob(2), K_cont_poly_rob(3), K_cont_poly_rob(4));
    else
        fprintf('8. Politópica (Contínua)         | %-6s | %9s | %8s | %s\n', ...
            'FALHOU', 'N/A', 'N/A', 'N/A');
    end

    if sucesso_cont_int_alt && ~isempty(K_cont_int_alt) && length(K_cont_int_alt) >= 4
        fprintf('9. Intervalar (Cont-Alt)         | %-6s | %9.6f | %8.4f | [%.4f %.4f %.4f %.4f]\n', ...
            'OK', gamma_cont_int_alt^2, tempo_cont_int_alt, K_cont_int_alt(1), K_cont_int_alt(2), K_cont_int_alt(3), K_cont_int_alt(4));
    else
        fprintf('9. Intervalar (Cont-Alt)         | %-6s | %9s | %8s | %s\n', ...
            'FALHOU', 'N/A', 'N/A', 'N/A');
    end

    if sucesso_cont_poly_alt && ~isempty(K_cont_poly_alt) && length(K_cont_poly_alt) >= 4
        fprintf('10. Politópica (Cont-Alt)        | %-6s | %9.6f | %8.4f | [%.4f %.4f %.4f %.4f]\n', ...
            'OK', gamma_cont_poly_alt^2, tempo_cont_poly_alt, K_cont_poly_alt(1), K_cont_poly_alt(2), K_cont_poly_alt(3), K_cont_poly_alt(4));
    else
        fprintf('10. Politópica (Cont-Alt)        | %-6s | %9s | %8s | %s\n', ...
            'FALHOU', 'N/A', 'N/A', 'N/A');
    end

    %% =================================================================
    %% SALVAMENTO DOS RESULTADOS
    %% =================================================================

    fprintf('\n=====================================\n');
    fprintf('SALVANDO RESULTADOS\n');
    fprintf('=====================================\n\n');

    % Estruturar dados para salvamento
    resultados.parametros.h = h;
    resultados.parametros.delta = delta;
    resultados.parametros.tol = tol;
    resultados.parametros.m1 = m1;
    resultados.parametros.m2 = m2;
    resultados.parametros.k1 = k1;
    resultados.parametros.k2 = k2;
    resultados.parametros.b = b;

    % Sistema baseado em genSaveDataEx12
    resultados.sistema.baseadoEm = 'genSaveDataEx12';
    resultados.sistema.estrutura = 'massa-mola-amortecedor 4x4';
    resultados.sistema.controle = 'entrada única na massa 2';
    resultados.sistema.perturbacao = 'entrada única na massa 1';

    % Ganhos obtidos
    if sucesso_int_nom
        resultados.ganhos.nominal.intervalar.K = K_int_nom;
        resultados.ganhos.nominal.intervalar.gamma = gamma_int_nom;
        resultados.ganhos.nominal.intervalar.tempo = tempo_int_nom;
    end

    if sucesso_poly_nom
        resultados.ganhos.nominal.politopica.K = K_poly_nom;
        resultados.ganhos.nominal.politopica.gamma = gamma_poly_nom;
        resultados.ganhos.nominal.politopica.tempo = tempo_poly_nom;
    end

    if sucesso_nom
        resultados.ganhos.nominal.classica.K = K_nom;
        resultados.ganhos.nominal.classica.gamma = gamma_nom;
        resultados.ganhos.nominal.classica.tempo = tempo_nom;
    end

    if sucesso_int_rob
        resultados.ganhos.robusto.intervalar.K = K_int_rob;
        resultados.ganhos.robusto.intervalar.gamma = gamma_int_rob;
        resultados.ganhos.robusto.intervalar.tempo = tempo_int_rob;
    end

    if sucesso_poly_rob
        resultados.ganhos.robusto.politopica.K = K_poly_rob;
        resultados.ganhos.robusto.politopica.gamma = gamma_poly_rob;
        resultados.ganhos.robusto.politopica.tempo = tempo_poly_rob;
    end

    if sucesso_central
        resultados.ganhos.robusto.central.K = K_central;
        resultados.ganhos.robusto.central.gamma = gamma_central;
        resultados.ganhos.robusto.central.tempo = tempo_central;
    end

    % Resultados das sínteses contínuas
    if sucesso_cont_int_rob && ~isempty(K_cont_int_rob)
        resultados.ganhos.continuo.intervalar.K = K_cont_int_rob;
        resultados.ganhos.continuo.intervalar.gamma = gamma_cont_int_rob;
        resultados.ganhos.continuo.intervalar.tempo = tempo_cont_int_rob;
    end

    if sucesso_cont_poly_rob && ~isempty(K_cont_poly_rob)
        resultados.ganhos.continuo.politopica.K = K_cont_poly_rob;
        resultados.ganhos.continuo.politopica.gamma = gamma_cont_poly_rob;
        resultados.ganhos.continuo.politopica.tempo = tempo_cont_poly_rob;
    end

    % Resultados das sínteses contínuas alternativas
    if sucesso_cont_int_alt && ~isempty(K_cont_int_alt)
        resultados.ganhos.continuo_alternativo.intervalar.K = K_cont_int_alt;
        resultados.ganhos.continuo_alternativo.intervalar.gamma = gamma_cont_int_alt;
        resultados.ganhos.continuo_alternativo.intervalar.tempo = tempo_cont_int_alt;
    end

    if sucesso_cont_poly_alt && ~isempty(K_cont_poly_alt)
        resultados.ganhos.continuo_alternativo.politopica.K = K_cont_poly_alt;
        resultados.ganhos.continuo_alternativo.politopica.gamma = gamma_cont_poly_alt;
        resultados.ganhos.continuo_alternativo.politopica.tempo = tempo_cont_poly_alt;
    end

    % Salvar workspace
    save('resultados_completos_massa_mola.mat', 'resultados');

    fprintf('Resultados salvos em: resultados_completos_massa_mola.mat\n\n');
    fprintf('=== SÍNTESE COMPLETA FINALIZADA COM SUCESSO ===\n');
    fprintf('Sistema baseado na estrutura do genSaveDataEx12\n');
    fprintf('10 métodos de síntese implementados:\n');
    fprintf('  - Métodos 1-3: Sínteses nominais\n');
    fprintf('  - Métodos 4-6: Sínteses robustas (amostrado)\n');
    fprintf('  - Métodos 7-8: Sínteses robustas (contínuo - LMILab)\n');
    fprintf('  - Métodos 9-10: Sínteses robustas (contínuo - alternativo)\n\n');

end