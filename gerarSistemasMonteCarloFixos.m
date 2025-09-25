function [sistemas_monte_carlo] = gerarSistemasMonteCarloFixos(parametros, numSimulacoes, seed)
%   GERARSISTEMASMONTECARLOFIXOS - Gera sistemas aleatórios fixos para Monte Carlo
%
%   Esta função gera uma vez só os sistemas aleatórios que serão usados
%   para testar TODOS os controladores, garantindo consistência na comparação.
%
%   Parâmetros:
%   - parametros: estrutura com parâmetros do sistema (m1, m2, k1, k2, b)
%   - numSimulacoes: número de sistemas aleatórios a gerar
%   - seed: semente para reprodutibilidade (opcional)
%
%   Saída:
%   - sistemas_monte_carlo: cell array com os sistemas gerados

    if nargin < 3
        seed = 42;  % Semente padrão para reprodutibilidade
    end

    if nargin < 2
        numSimulacoes = 100;
    end

    fprintf('      [Gerador Monte Carlo] Gerando %d sistemas aleatórios fixos...\n', numSimulacoes);
    fprintf('      Usando semente: %d (para reprodutibilidade)\n', seed);

    % Definir semente para reprodutibilidade
    try
        % Primeiro resetar para o gerador padrão
        rng('default');
        % Depois definir a semente
        rng(seed, 'twister');
    catch ME
        % Se ainda falhar, usar método legado
        warning('Usando método legado para semente aleatória: %s', ME.message);
        rand('state', seed);
        randn('state', seed);
    end

    % Extrair parâmetros do sistema
    m1 = parametros.m1;
    m2 = parametros.m2;
    k1 = parametros.k1;
    k2_nominal = parametros.k2;
    b_nominal = parametros.b;

    % Definir faixas de incerteza
    k2_min = k2_nominal * 0.98;  % -2%
    k2_max = k2_nominal * 1.02;  % +2%
    b_min = b_nominal * 0.95;    % -5%
    b_max = b_nominal * 1.05;    % +5%

    fprintf('      Incertezas: k2=[%.3f, %.3f], b=[%.4f, %.4f]\n', k2_min, k2_max, b_min, b_max);

    % Gerar parâmetros aleatórios UMA ÚNICA VEZ
    k2_valores = k2_min + (k2_max - k2_min) * rand(numSimulacoes, 1);
    b_valores = b_min + (b_max - b_min) * rand(numSimulacoes, 1);

    % Inicializar estrutura de saída
    sistemas_monte_carlo = cell(numSimulacoes, 1);

    % Gerar sistemas para cada combinação de parâmetros
    for i = 1:numSimulacoes
        k2_i = k2_valores(i);
        b_i = b_valores(i);

        % Construir sistema para estes parâmetros
        A_i = [0                      0                   1              0;
               0                      0                   0              1;
               (-k2_i-k1)/m1         k2_i/m1            -b_i/m1         b_i/m1;
               k2_i/m2               -k2_i/m2            b_i/m2         -b_i/m2];

        B_i = [0; 0; 0; 1/m2];           % Entrada de controle
        E_i = [0; 0; 1/m1; 0];           % Entrada de perturbação
        C_i = [0 10 0 0; 0 0 0 1; 0 0 0 0];  % Saídas
        D_i = [0; 0; 1];                 % Feedthrough
        D1_i = [0; 0; 0];                % Feedthrough perturbação

        % Armazenar sistema
        sistemas_monte_carlo{i}.A = A_i;
        sistemas_monte_carlo{i}.B = B_i;
        sistemas_monte_carlo{i}.E = E_i;
        sistemas_monte_carlo{i}.C = C_i;
        sistemas_monte_carlo{i}.D = D_i;
        sistemas_monte_carlo{i}.D1 = D1_i;
        sistemas_monte_carlo{i}.k2 = k2_i;
        sistemas_monte_carlo{i}.b = b_i;
        sistemas_monte_carlo{i}.indice = i;
    end

    % Estatísticas dos parâmetros gerados
    k2_mean = mean(k2_valores);
    k2_std = std(k2_valores);
    b_mean = mean(b_valores);
    b_std = std(b_valores);

    fprintf('      ✓ Sistemas gerados com sucesso!\n');
    fprintf('      Estatísticas: k2_mean=%.3f±%.3f, b_mean=%.4f±%.4f\n', ...
        k2_mean, k2_std, b_mean, b_std);

    % Verificar cobertura das faixas
    k2_coverage = (max(k2_valores) - min(k2_valores)) / (k2_max - k2_min) * 100;
    b_coverage = (max(b_valores) - min(b_valores)) / (b_max - b_min) * 100;

    fprintf('      Cobertura das faixas: k2=%.1f%%, b=%.1f%%\n', k2_coverage, b_coverage);

end