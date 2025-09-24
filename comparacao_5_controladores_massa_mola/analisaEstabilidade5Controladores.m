function resultado = analisaEstabilidade5Controladores(controladores, numRealizacoes)
%   Análise de estabilidade Monte Carlo para os 5 controladores
%
%   Entrada:
%   - controladores: estrutura com os 5 controladores (saída de comparacao5Controladores)
%   - numRealizacoes: número de realizações aleatórias (padrão: 1000)
%
%   Saída:
%   - resultado: estrutura com taxas de estabilidade e normas H-infinity

    if nargin < 2
        numRealizacoes = 1000;
    end

    fprintf('=== ANÁLISE DE ESTABILIDADE MONTE CARLO ===\n');
    fprintf('Número de realizações: %d\n\n', numRealizacoes);

    % Nomes dos controladores
    nomes = {'Hibrido_Int', 'Hibrido_Poly', 'Amostrado_Int', 'Amostrado_Poly', 'Nominal'};
    numControladores = length(nomes);

    % Adicionar caminhos necessários
    addpath('../HInf - Análise - Intervalar/funcoes');
    addpath('../HInf - Análise - Intervalar/funcoes/Diversas');
    addpath('../funcoes');

    % Carregar sistema de referência
    sistema = genSaveDataComparacao();
    h = sistema.parametros.h;

    % Parâmetros das incertezas
    m1 = 0.5;
    m2 = 1.0;
    k1 = 12.0;

    % Intervalos de incerteza
    b_min = 0.2 * 0.95;    % -5%
    b_max = 0.2 * 1.05;    % +5%
    k2_min = 7.0 * 0.98;   % -2%
    k2_max = 7.0 * 1.02;   % +2%

    % Inicializar variáveis de resultado
    taxaEstabilidade = zeros(1, numControladores);
    normaPiorCaso = zeros(1, numControladores);
    normaMedia = zeros(1, numControladores);
    desviosPadrao = zeros(1, numControladores);
    contadorEstavel = zeros(1, numControladores);
    normasHinf = cell(1, numControladores);

    for i = 1:numControladores
        normasHinf{i} = [];
    end

    %% Simulação Monte Carlo
    fprintf('Iniciando simulação Monte Carlo...\n');
    for realizacao = 1:numRealizacoes
        if mod(realizacao, 100) == 0
            fprintf('Realizacao %d/%d...\n', realizacao, numRealizacoes);
        end

        % Gerar parâmetros aleatórios
        b_real = b_min + (b_max - b_min) * rand();
        k2_real = k2_min + (k2_max - k2_min) * rand();

        % Matriz A real para esta realização
        A_real = [0              0       1       0;
                  0              0       0       1;
                  (-k2_real-k1)/m1    k2_real/m1   -b_real/m1   b_real/m1;
                  k2_real/m2          -k2_real/m2  b_real/m2    -b_real/m2];

        B_real = [0 0 0 1/m2]';
        E_real = [0 0 1/m1 0]';
        C_real = [0 10 0 0; 0 0 0 1; 0 0 0 0];

        % Testar cada controlador
        for ctrl = 1:numControladores
            nome = nomes{ctrl};

            if ~controladores.(nome).factivel
                continue; % Pular se a síntese falhou
            end

            K_ctrl = controladores.(nome).K;

            try
                % Verificar estabilidade
                A_cl = A_real + B_real * K_ctrl;
                eigenvals = eig(A_cl);
                estavel = all(real(eigenvals) < -1e-6); % Margem de estabilidade

                if estavel
                    contadorEstavel(ctrl) = contadorEstavel(ctrl) + 1;

                    % Calcular norma H-infinity
                    normaHinf = calculaNormaHinfContinuo(A_cl, E_real, C_real);

                    if ~isnan(normaHinf) && normaHinf < inf
                        normasHinf{ctrl} = [normasHinf{ctrl}, normaHinf];
                    end
                end

            catch
                % Se houver erro, considera como instável
                continue;
            end
        end
    end

    %% Calcular estatísticas
    fprintf('\nCalculando estatísticas finais...\n');

    for ctrl = 1:numControladores
        nome = nomes{ctrl};

        if ~controladores.(nome).factivel
            taxaEstabilidade(ctrl) = 0;
            normaPiorCaso(ctrl) = inf;
            normaMedia(ctrl) = inf;
            desviosPadrao(ctrl) = 0;
            continue;
        end

        % Taxa de estabilidade
        taxaEstabilidade(ctrl) = (contadorEstavel(ctrl) / numRealizacoes) * 100;

        % Estatísticas das normas
        if ~isempty(normasHinf{ctrl})
            normaPiorCaso(ctrl) = max(normasHinf{ctrl});
            normaMedia(ctrl) = mean(normasHinf{ctrl});
            desviosPadrao(ctrl) = std(normasHinf{ctrl});
        else
            normaPiorCaso(ctrl) = inf;
            normaMedia(ctrl) = inf;
            desviosPadrao(ctrl) = 0;
        end
    end

    %% Estruturar resultados
    resultado.taxaEstabilidade = taxaEstabilidade;
    resultado.normaPiorCaso = normaPiorCaso;
    resultado.normaMedia = normaMedia;
    resultado.desviosPadrao = desviosPadrao;
    resultado.nomes = {'Híbrido Int', 'Híbrido Poly', 'Amostrado Int', 'Amostrado Poly', 'Nominal'};
    resultado.numRealizacoes = numRealizacoes;
    resultado.contadorEstavel = contadorEstavel;
    resultado.normasCompletas = normasHinf;

    %% Exibir resultados
    fprintf('\n=== RESULTADOS DA ANÁLISE DE ESTABILIDADE ===\n');
    fprintf('Controlador           | Taxa Estab. | Pior Caso  | Média      | Desvio     | Realizações\n');
    fprintf('----------------------|-------------|------------|------------|------------|-----------\n');

    for ctrl = 1:numControladores
        nome = nomes{ctrl};
        nomeDisplay = resultado.nomes{ctrl};

        if controladores.(nome).factivel
            if isinf(normaPiorCaso(ctrl))
                piorStr = '∞';
                mediaStr = '∞';
                desvioStr = 'N/A';
            else
                piorStr = sprintf('%.3f', normaPiorCaso(ctrl));
                mediaStr = sprintf('%.3f', normaMedia(ctrl));
                desvioStr = sprintf('%.3f', desviosPadrao(ctrl));
            end

            fprintf('%-20s | %9.1f%% | %10s | %10s | %10s | %d/%d\n', ...
                nomeDisplay, taxaEstabilidade(ctrl), piorStr, mediaStr, ...
                desvioStr, contadorEstavel(ctrl), numRealizacoes);
        else
            fprintf('%-20s | %9s | %10s | %10s | %10s | SÍNTESE FALHOU\n', ...
                nomeDisplay, 'N/A', 'N/A', 'N/A', 'N/A');
        end
    end

    fprintf('\nAnálise de estabilidade finalizada!\n');

end

%% Função auxiliar para cálculo da norma H-infinity
function normaHinf = calculaNormaHinfContinuo(A, E, C)
%   Calcula norma H-infinity do sistema contínuo usando método existente
%   A: matriz de estado em malha fechada
%   E: matriz de entrada de perturbação
%   C: matriz de saída

    try
        % Verificar estabilidade
        eigenvals = eig(A);
        if any(real(eigenvals) >= -1e-6)
            normaHinf = inf;
            return;
        end

        % Usar função existente normaSistemaContinuo se disponível
        if exist('normaSistemaContinuo', 'file')
            D1 = zeros(size(C, 1), size(E, 2)); % Sem feedthrough da perturbação
            normaHinf = normaSistemaContinuo(A, E, C, D1);
        else
            % Método alternativo usando Control Systems Toolbox
            sys = ss(A, E, C, zeros(size(C, 1), size(E, 2)));
            normaHinf = norm(sys, inf);
        end

        % Verificar se o resultado é válido
        if isnan(normaHinf) || normaHinf <= 0
            normaHinf = inf;
        end

    catch
        normaHinf = inf;
    end
end