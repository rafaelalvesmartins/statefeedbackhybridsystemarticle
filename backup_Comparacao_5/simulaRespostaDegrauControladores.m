function resultado = simulaRespostaDegrauControladores(ganhosControladores, numRealizacoes, tempoSimulacao)
%   Simulação da resposta ao degrau para os 5 controladores
%   com análise Monte Carlo no Sistema do Exemplo 6
%
%   Entrada:
%   - ganhosControladores: struct com os 5 ganhos dos controladores
%   - numRealizacoes: número de realizações aleatórias (padrão: 100)
%   - tempoSimulacao: tempo total de simulação em segundos (padrão: 10)
%
%   Saída:
%   - resultado: struct com respostas, normas L2 e gráficos

    if nargin < 2
        numRealizacoes = 100;
    end
    if nargin < 3
        tempoSimulacao = 10;
    end

    fprintf('=== SIMULAÇÃO RESPOSTA AO DEGRAU - 5 CONTROLADORES ===\n');
    fprintf('Sistema: Exemplo 6\n');
    fprintf('Número de realizações: %d\n', numRealizacoes);
    fprintf('Tempo de simulação: %.1f segundos\n\n', tempoSimulacao);

    %% Parâmetros do Sistema do Exemplo 6
    % Definição das incertezas intervalares
    A_intervals = [infsup(0.9, 1.1), infsup(-1.1, -0.9);
                   infsup(-0.5, 0.5), infsup(3.8, 4.2)];

    Bw = [1.0; -1.0];  % Entrada de perturbação
    Bu = [1.0; -1.0];  % Entrada de controle
    C = [eye(2); zeros(2,2)];  % Matriz de saída
    h = 0.1;  % Período de amostragem

    % Dimensões do sistema
    nx = 2;
    nu = 1;

    %% Definição dos controladores
    nomeControladores = {'Nominal', 'Discreto Int', 'Discreto Poly', 'Hibrido Int', 'Hibrido Poly'};
    numControladores = length(nomeControladores);

    % Extração dos ganhos
    ganhos = extrairGanhos(ganhosControladores, nomeControladores);

    %% Configuração da simulação
    dt = 0.01;  % Passo de integração
    t = 0:dt:tempoSimulacao;
    N = length(t);

    % Entrada degrau unitário por 1 segundo
    entrada = zeros(size(t));
    entrada(t <= 1) = 1;

    %% Inicialização das variáveis
    respostasCompletas = cell(numControladores, numRealizacoes);
    normasL2 = cell(numControladores, 1);
    contadorEstavel = zeros(numControladores, 1);

    for ctrl = 1:numControladores
        normasL2{ctrl} = [];
    end

    %% Monte Carlo - Simulação para cada realização
    fprintf('Iniciando simulações Monte Carlo...\n');

    for realizacao = 1:numRealizacoes
        if mod(realizacao, 10) == 0
            fprintf('Processando realização %d/%d...\n', realizacao, numRealizacoes);
        end

        % Gerar realização aleatória do sistema
        A_real = geraMatrizAleatoria(A_intervals);

        % Testar cada controlador
        for ctrl = 1:numControladores
            if isempty(ganhos{ctrl})
                respostasCompletas{ctrl, realizacao} = nan(N, nx);
                continue;
            end

            K_ctrl = ganhos{ctrl};

            try
                % Sistema em malha fechada contínuo
                Acl = A_real + Bu * K_ctrl;
                Bcl = Bu;  % Entrada de referência
                Ccl = C(1:2, :);  % Apenas os dois primeiros estados (saídas principais)
                Dcl = zeros(2, 1);

                % Verificar estabilidade
                eigenvals = eig(Acl);
                if any(real(eigenvals) >= 0)
                    respostasCompletas{ctrl, realizacao} = nan(N, 2);
                    continue;
                end

                % Criar sistema LTI
                sys_cl = ss(Acl, Bcl, Ccl, Dcl);

                % Simular resposta
                [y_sim, t_sim] = lsim(sys_cl, entrada, t);

                % Armazenar resposta (apenas as 2 primeiras saídas)
                respostasCompletas{ctrl, realizacao} = y_sim;
                contadorEstavel(ctrl) = contadorEstavel(ctrl) + 1;

                % Calcular norma L2 da resposta
                normaL2_resposta = sqrt(trapz(t, sum(y_sim.^2, 2)) * dt);
                normasL2{ctrl} = [normasL2{ctrl}, normaL2_resposta];

            catch ME
                % Se houver erro na simulação
                respostasCompletas{ctrl, realizacao} = nan(N, 2);
                fprintf('Erro na simulação - Controlador %d, Realização %d: %s\n', ctrl, realizacao, ME.message);
            end
        end
    end

    %% Cálculo de estatísticas
    fprintf('\nCalculando estatísticas das respostas...\n');

    respostasMedia = cell(numControladores, 1);
    respostasStd = cell(numControladores, 1);
    normasL2_stats = struct();

    for ctrl = 1:numControladores
        if contadorEstavel(ctrl) > 0
            % Coletar respostas válidas
            respostas_validas = [];
            for real = 1:numRealizacoes
                resp = respostasCompletas{ctrl, real};
                if ~any(isnan(resp(:)))
                    if isempty(respostas_validas)
                        respostas_validas = resp;
                    else
                        respostas_validas = cat(3, respostas_validas, resp);
                    end
                end
            end

            if ~isempty(respostas_validas)
                % Calcular média e desvio padrão
                respostasMedia{ctrl} = mean(respostas_validas, 3);
                respostasStd{ctrl} = std(respostas_validas, 0, 3);

                % Estatísticas das normas L2
                normasL2_stats.(sprintf('ctrl%d', ctrl)).media = mean(normasL2{ctrl});
                normasL2_stats.(sprintf('ctrl%d', ctrl)).std = std(normasL2{ctrl});
                normasL2_stats.(sprintf('ctrl%d', ctrl)).max = max(normasL2{ctrl});
                normasL2_stats.(sprintf('ctrl%d', ctrl)).min = min(normasL2{ctrl});
            else
                respostasMedia{ctrl} = nan(N, 2);
                respostasStd{ctrl} = zeros(N, 2);
            end
        else
            respostasMedia{ctrl} = nan(N, 2);
            respostasStd{ctrl} = zeros(N, 2);
        end
    end

    %% Geração dos gráficos
    fprintf('Gerando gráficos...\n');
    figHandle = geraGraficosComparacao(t, entrada, respostasMedia, respostasStd, nomeControladores, contadorEstavel, numRealizacoes);

    %% Estruturação dos resultados
    resultado.tempo = t;
    resultado.entrada = entrada;
    resultado.respostasMedia = respostasMedia;
    resultado.respostasStd = respostasStd;
    resultado.nomeControladores = nomeControladores;
    resultado.contadorEstavel = contadorEstavel;
    resultado.percentualEstavel = (contadorEstavel / numRealizacoes) * 100;
    resultado.normasL2 = normasL2;
    resultado.normasL2_stats = normasL2_stats;
    resultado.figHandle = figHandle;
    resultado.numRealizacoes = numRealizacoes;

    %% Exibição dos resultados
    fprintf('\n=== RESULTADOS DA SIMULAÇÃO ===\n');
    fprintf('Controlador           | %% Estável | Norma L2 (Média) | Norma L2 (Std) | Realizações Estáveis\n');
    fprintf('----------------------|----------|------------------|----------------|--------------------\n');

    for ctrl = 1:numControladores
        if contadorEstavel(ctrl) > 0
            fprintf('%-20s | %7.1f%% | %15.4f | %13.4f | %d/%d\n', ...
                nomeControladores{ctrl}, resultado.percentualEstavel(ctrl), ...
                normasL2_stats.(sprintf('ctrl%d', ctrl)).media, ...
                normasL2_stats.(sprintf('ctrl%d', ctrl)).std, ...
                contadorEstavel(ctrl), numRealizacoes);
        else
            fprintf('%-20s | %7.1f%% | %15s | %13s | %d/%d\n', ...
                nomeControladores{ctrl}, 0.0, 'N/A', 'N/A', 0, numRealizacoes);
        end
    end

    fprintf('\nSimulação finalizada!\n');
    fprintf('Gráfico salvo como figura MATLAB.\n');

end

%% Funções auxiliares

function ganhos = extrairGanhos(ganhosControladores, nomeControladores)
    numControladores = length(nomeControladores);
    ganhos = cell(1, numControladores);

    if isstruct(ganhosControladores)
        if isfield(ganhosControladores, 'nominal') && ganhosControladores.nominal.sucesso
            ganhos{1} = ganhosControladores.nominal.K;
        else
            ganhos{1} = [];
        end

        if isfield(ganhosControladores, 'discretoIntervalar') && ganhosControladores.discretoIntervalar.sucesso
            ganhos{2} = ganhosControladores.discretoIntervalar.K;
        else
            ganhos{2} = [];
        end

        if isfield(ganhosControladores, 'discretoPolitopico') && ganhosControladores.discretoPolitopico.sucesso
            ganhos{3} = ganhosControladores.discretoPolitopico.K;
        else
            ganhos{3} = [];
        end

        if isfield(ganhosControladores, 'hibridoIntervalar') && ganhosControladores.hibridoIntervalar.sucesso
            ganhos{4} = ganhosControladores.hibridoIntervalar.K;
        else
            ganhos{4} = [];
        end

        if isfield(ganhosControladores, 'hibridoPolitopico') && ganhosControladores.hibridoPolitopico.sucesso
            ganhos{5} = ganhosControladores.hibridoPolitopico.K;
        else
            ganhos{5} = [];
        end
    end
end

function A_real = geraMatrizAleatoria(A_intervals)
    [m, n] = size(A_intervals);
    A_real = zeros(m, n);

    for i = 1:m
        for j = 1:n
            inf_val = A_intervals(i,j).inf;
            sup_val = A_intervals(i,j).sup;
            A_real(i,j) = inf_val + (sup_val - inf_val) * rand();
        end
    end
end

function figHandle = geraGraficosComparacao(t, entrada, respostasMedia, respostasStd, nomeControladores, contadorEstavel, numRealizacoes)
    figHandle = figure('Position', [100, 100, 1400, 1000]);

    cores = {'b', 'r', 'g', 'm', 'c'};  % Cores para cada controlador

    for ctrl = 1:length(nomeControladores)
        % Subplot para cada controlador
        subplot(3, 2, ctrl);

        if contadorEstavel(ctrl) > 0 && ~any(isnan(respostasMedia{ctrl}(:)))
            % Plotar entrada (referência)
            yyaxis right;
            plot(t, entrada, 'k--', 'LineWidth', 1, 'DisplayName', 'Referência');
            ylabel('Entrada');
            ylim([-0.1, 1.2]);

            yyaxis left;
            % Plotar resposta média ± desvio padrão para cada estado
            for estado = 1:size(respostasMedia{ctrl}, 2)
                media = respostasMedia{ctrl}(:, estado);
                desvio = respostasStd{ctrl}(:, estado);

                % Área do desvio padrão
                fill([t, fliplr(t)], [media + desvio; flipud(media - desvio)]', ...
                     cores{ctrl}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
                     'HandleVisibility', 'off');
                hold on;

                % Linha da média
                plot(t, media, 'Color', cores{ctrl}, 'LineWidth', 2, ...
                     'DisplayName', sprintf('x_%d', estado));
            end

            ylabel('Estados');
            xlabel('Tempo (s)');
            title(sprintf('%s\n(%d/%d realizações estáveis = %.1f%%)', ...
                  nomeControladores{ctrl}, contadorEstavel(ctrl), numRealizacoes, ...
                  (contadorEstavel(ctrl)/numRealizacoes)*100));
            grid on;
            legend('Location', 'best');

        else
            % Controlador instável ou sem dados
            text(0.5, 0.5, sprintf('%s\nInstável ou sem dados\n(%d/%d realizações)', ...
                 nomeControladores{ctrl}, contadorEstavel(ctrl), numRealizacoes), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 12, 'Color', 'red');
            xlim([0, max(t)]);
            ylim([0, 1]);
            title(nomeControladores{ctrl});
            grid on;
        end
    end

    % Subplot adicional com comparação das normas L2
    subplot(3, 2, 6);
    normasL2_medias = zeros(1, length(nomeControladores));
    normasL2_stds = zeros(1, length(nomeControladores));

    for ctrl = 1:length(nomeControladores)
        if contadorEstavel(ctrl) > 0
            normasL2_medias(ctrl) = mean(respostasMedia{ctrl}(:));
            normasL2_stds(ctrl) = std(respostasMedia{ctrl}(:));
        else
            normasL2_medias(ctrl) = 0;
            normasL2_stds(ctrl) = 0;
        end
    end

    bar(normasL2_medias);
    hold on;
    errorbar(1:length(nomeControladores), normasL2_medias, normasL2_stds, 'k.');

    set(gca, 'XTickLabel', nomeControladores);
    xtickangle(45);
    ylabel('Norma L2 Média');
    title('Comparação Normas L2');
    grid on;

    sgtitle('Comparação de Respostas ao Degrau - 5 Controladores', 'FontSize', 16, 'FontWeight', 'bold');

    % Salvar figura
    saveas(figHandle, 'comparacao_resposta_degrau_controladores.fig');
    saveas(figHandle, 'comparacao_resposta_degrau_controladores.png');
end