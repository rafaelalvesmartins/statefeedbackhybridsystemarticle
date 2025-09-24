function resultado = simulacaoTemporal5Controladores(controladores, numRealizacoes, tempoSimulacao)
%   Simulação temporal com resposta ao degrau para os 5 controladores
%
%   Entrada:
%   - controladores: estrutura com os 5 controladores
%   - numRealizacoes: número de realizações aleatórias (padrão: 50)
%   - tempoSimulacao: tempo total de simulação em segundos (padrão: 15)
%
%   Saída:
%   - resultado: estrutura com métricas temporais e gráficos

    if nargin < 2
        numRealizacoes = 50;
    end
    if nargin < 3
        tempoSimulacao = 15;
    end

    fprintf('=== SIMULAÇÃO TEMPORAL - 5 CONTROLADORES ===\n');
    fprintf('Número de realizações: %d\n', numRealizacoes);
    fprintf('Tempo de simulação: %.1f segundos\n\n', tempoSimulacao);

    % Nomes dos controladores
    nomes = {'Hibrido_Int', 'Hibrido_Poly', 'Amostrado_Int', 'Amostrado_Poly', 'Nominal'};
    nomesDisplay = {'Híbrido Int', 'Híbrido Poly', 'Amostrado Int', 'Amostrado Poly', 'Nominal'};
    numControladores = length(nomes);

    % Adicionar caminhos necessários
    addpath('../HInf - Análise - Intervalar/funcoes');
    addpath('../HInf - Análise - Intervalar/funcoes/Diversas');
    addpath('../funcoes');

    % Carregar sistema
    sistema = genSaveDataComparacao();
    h = sistema.parametros.h;

    % Parâmetros das incertezas
    m1 = 0.5; m2 = 1.0; k1 = 12.0;
    b_min = 0.2 * 0.95; b_max = 0.2 * 1.05;
    k2_min = 7.0 * 0.98; k2_max = 7.0 * 1.02;

    % Configuração da simulação
    dt = 0.01; % Passo de integração
    t = 0:dt:tempoSimulacao;
    N = length(t);

    % Entrada degrau unitário por 1 segundo
    entrada = zeros(size(t));
    entrada(t <= 1) = 1;

    % Inicializar variáveis
    respostasCompletas = cell(numControladores, numRealizacoes);
    normasL2 = cell(numControladores, 1);
    temposAcomodacao = cell(numControladores, 1);
    overshoots = cell(numControladores, 1);
    contadorEstavel = zeros(numControladores, 1);

    for i = 1:numControladores
        normasL2{i} = [];
        temposAcomodacao{i} = [];
        overshoots{i} = [];
    end

    %% Simulação Monte Carlo
    fprintf('Iniciando simulações...\n');
    for realizacao = 1:numRealizacoes
        if mod(realizacao, 10) == 0
            fprintf('Realizacao %d/%d...\n', realizacao, numRealizacoes);
        end

        % Gerar sistema aleatório
        b_real = b_min + (b_max - b_min) * rand();
        k2_real = k2_min + (k2_max - k2_min) * rand();

        A_real = [0              0       1       0;
                  0              0       0       1;
                  (-k2_real-k1)/m1    k2_real/m1   -b_real/m1   b_real/m1;
                  k2_real/m2          -k2_real/m2  b_real/m2    -b_real/m2];

        B_real = [0 0 0 1/m2]';
        C_real = [0 10 0 0; 0 0 0 1]; % Apenas primeiros 2 estados como saída

        % Testar cada controlador
        for ctrl = 1:numControladores
            nome = nomes{ctrl};

            if ~controladores.(nome).factivel
                respostasCompletas{ctrl, realizacao} = nan(N, 2);
                continue;
            end

            K_ctrl = controladores.(nome).K;

            try
                % Sistema em malha fechada
                A_cl = A_real + B_real * K_ctrl;
                eigenvals = eig(A_cl);

                % Verificar estabilidade
                if any(real(eigenvals) >= -1e-6)
                    respostasCompletas{ctrl, realizacao} = nan(N, 2);
                    continue;
                end

                % Simular resposta
                sys_cl = ss(A_cl, B_real, C_real, zeros(size(C_real, 1), size(B_real, 2)));
                [y_sim, ~] = lsim(sys_cl, entrada, t);

                % Armazenar resposta
                respostasCompletas{ctrl, realizacao} = y_sim;
                contadorEstavel(ctrl) = contadorEstavel(ctrl) + 1;

                % Calcular métricas
                normaL2 = sqrt(trapz(t, sum(y_sim.^2, 2)) * dt);
                normasL2{ctrl} = [normasL2{ctrl}, normaL2];

                % Tempo de acomodação (2% da resposta final)
                tempoAcomod = calcularTempoAcomodacao(t, y_sim(:, 1), 0.02);
                temposAcomodacao{ctrl} = [temposAcomodacao{ctrl}, tempoAcomod];

                % Overshoot máximo
                overshoot = calcularOvershoot(y_sim(:, 1));
                overshoots{ctrl} = [overshoots{ctrl}, overshoot];

            catch ME
                respostasCompletas{ctrl, realizacao} = nan(N, 2);
                fprintf('Erro na simulação - Controlador %s, Realizacao %d: %s\n', ...
                    nome, realizacao, ME.message);
            end
        end
    end

    %% Calcular estatísticas
    fprintf('\nCalculando estatísticas das respostas...\n');

    respostasMedia = cell(numControladores, 1);
    respostasStd = cell(numControladores, 1);
    metricas = struct();

    for ctrl = 1:numControladores
        nome = nomes{ctrl};

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
                respostasMedia{ctrl} = mean(respostas_validas, 3);
                respostasStd{ctrl} = std(respostas_validas, 0, 3);
            else
                respostasMedia{ctrl} = nan(N, 2);
                respostasStd{ctrl} = zeros(N, 2);
            end

            % Métricas
            metricas.(nome).normaL2_media = mean(normasL2{ctrl});
            metricas.(nome).normaL2_std = std(normasL2{ctrl});
            metricas.(nome).tempoAcomodacao_media = mean(temposAcomodacao{ctrl});
            metricas.(nome).tempoAcomodacao_std = std(temposAcomodacao{ctrl});
            metricas.(nome).overshoot_media = mean(overshoots{ctrl});
            metricas.(nome).overshoot_std = std(overshoots{ctrl});

        else
            respostasMedia{ctrl} = nan(N, 2);
            respostasStd{ctrl} = zeros(N, 2);

            metricas.(nome).normaL2_media = inf;
            metricas.(nome).normaL2_std = 0;
            metricas.(nome).tempoAcomodacao_media = inf;
            metricas.(nome).tempoAcomodacao_std = 0;
            metricas.(nome).overshoot_media = inf;
            metricas.(nome).overshoot_std = 0;
        end
    end

    %% Gerar gráficos
    fprintf('Gerando gráficos...\n');
    figHandle = gerarGraficosComparativos(t, entrada, respostasMedia, respostasStd, ...
        nomesDisplay, contadorEstavel, numRealizacoes, metricas, nomes);

    %% Estruturar resultados
    resultado.tempo = t;
    resultado.entrada = entrada;
    resultado.respostasMedia = respostasMedia;
    resultado.respostasStd = respostasStd;
    resultado.nomesControladores = nomesDisplay;
    resultado.contadorEstavel = contadorEstavel;
    resultado.percentualEstavel = (contadorEstavel / numRealizacoes) * 100;
    resultado.metricas = metricas;
    resultado.figHandle = figHandle;
    resultado.numRealizacoes = numRealizacoes;

    %% Exibir resultados
    fprintf('\n=== RESULTADOS DA SIMULAÇÃO TEMPORAL ===\n');
    fprintf('Controlador           | %% Estável | Norma L2   | Tempo Acom. | Overshoot  | Realizações\n');
    fprintf('----------------------|----------|------------|-------------|------------|-----------\n');

    for ctrl = 1:numControladores
        nome = nomes{ctrl};
        nomeDisplay = nomesDisplay{ctrl};

        if contadorEstavel(ctrl) > 0
            fprintf('%-20s | %7.1f%% | %8.3f | %9.3f | %8.1f%% | %d/%d\n', ...
                nomeDisplay, resultado.percentualEstavel(ctrl), ...
                metricas.(nome).normaL2_media, ...
                metricas.(nome).tempoAcomodacao_media, ...
                metricas.(nome).overshoot_media * 100, ...
                contadorEstavel(ctrl), numRealizacoes);
        else
            fprintf('%-20s | %7.1f%% | %8s | %9s | %8s | %d/%d\n', ...
                nomeDisplay, 0.0, 'N/A', 'N/A', 'N/A', 0, numRealizacoes);
        end
    end

    fprintf('\nSimulação temporal finalizada!\n');

end

%% Funções auxiliares

function tempoAcomod = calcularTempoAcomodacao(t, resposta, tolerancia)
    % Calcula tempo de acomodação (resposta dentro de ±tolerancia do valor final)

    try
        valorFinal = resposta(end);
        limiteInf = valorFinal * (1 - tolerancia);
        limiteSup = valorFinal * (1 + tolerancia);

        % Encontrar último ponto fora da banda de tolerância
        foraLimites = (resposta < limiteInf) | (resposta > limiteSup);
        indiceUltimo = find(foraLimites, 1, 'last');

        if isempty(indiceUltimo)
            tempoAcomod = 0;
        else
            tempoAcomod = t(indiceUltimo);
        end
    catch
        tempoAcomod = inf;
    end
end

function overshoot = calcularOvershoot(resposta)
    % Calcula overshoot máximo em relação ao valor final

    try
        valorFinal = resposta(end);
        if valorFinal ~= 0
            overshoot = (max(resposta) - valorFinal) / abs(valorFinal);
            overshoot = max(0, overshoot); % Apenas overshoot positivo
        else
            overshoot = 0;
        end
    catch
        overshoot = inf;
    end
end

function figHandle = gerarGraficosComparativos(t, entrada, respostasMedia, respostasStd, ...
    nomesDisplay, contadorEstavel, numRealizacoes, metricas, nomes)

    figHandle = figure('Position', [100, 100, 1600, 1200]);
    cores = {'b', 'r', 'g', 'm', 'c'};

    % Subplots para cada controlador
    for ctrl = 1:length(nomesDisplay)
        subplot(3, 2, ctrl);

        if contadorEstavel(ctrl) > 0 && ~any(isnan(respostasMedia{ctrl}(:)))
            % Plotar entrada
            yyaxis right;
            plot(t, entrada, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Referência');
            ylabel('Entrada');
            ylim([-0.2, 1.3]);

            yyaxis left;
            % Plotar resposta com área de incerteza
            for estado = 1:min(2, size(respostasMedia{ctrl}, 2))
                media = respostasMedia{ctrl}(:, estado);
                desvio = respostasStd{ctrl}(:, estado);

                % Área do desvio padrão
                fill([t, fliplr(t)], [media + desvio; flipud(media - desvio)]', ...
                     cores{ctrl}, 'FaceAlpha', 0.25, 'EdgeColor', 'none', ...
                     'HandleVisibility', 'off');
                hold on;

                % Linha da média
                plot(t, media, 'Color', cores{ctrl}, 'LineWidth', 2.5, ...
                     'DisplayName', sprintf('x_%d', estado));
            end

            ylabel('Estados');
            xlabel('Tempo (s)');

            % Título com métricas
            titulo = sprintf('%s (%.1f%% estável)\nL2: %.3f | Acom: %.2fs | OS: %.1f%%', ...
                nomesDisplay{ctrl}, (contadorEstavel(ctrl)/numRealizacoes)*100, ...
                metricas.(nomes{ctrl}).normaL2_media, ...
                metricas.(nomes{ctrl}).tempoAcomodacao_media, ...
                metricas.(nomes{ctrl}).overshoot_media * 100);
            title(titulo, 'FontSize', 10);

            grid on;
            legend('Location', 'best');
            xlim([0, max(t)]);

        else
            % Controlador instável ou sem dados
            text(0.5, 0.5, sprintf('%s\nInstável ou Síntese Falhou\n(%d/%d realizações)', ...
                 nomesDisplay{ctrl}, contadorEstavel(ctrl), numRealizacoes), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontSize', 12, 'Color', 'red');
            xlim([0, max(t)]);
            ylim([0, 1]);
            title(nomesDisplay{ctrl});
            grid on;
        end
    end

    % Subplot comparativo de métricas
    subplot(3, 2, 6);
    metricasMatrix = zeros(length(nomesDisplay), 3);

    for ctrl = 1:length(nomesDisplay)
        if contadorEstavel(ctrl) > 0
            metricasMatrix(ctrl, 1) = metricas.(nomes{ctrl}).normaL2_media;
            metricasMatrix(ctrl, 2) = metricas.(nomes{ctrl}).tempoAcomodacao_media;
            metricasMatrix(ctrl, 3) = metricas.(nomes{ctrl}).overshoot_media * 100;
        else
            metricasMatrix(ctrl, :) = [0, 0, 0];
        end
    end

    bar(metricasMatrix);
    set(gca, 'XTickLabel', nomesDisplay);
    xtickangle(45);
    ylabel('Valores das Métricas');
    legend({'Norma L2', 'Tempo Acomod. (s)', 'Overshoot (%)'}, 'Location', 'best');
    title('Comparação de Métricas Temporais');
    grid on;

    sgtitle('Simulação Temporal - Comparação de 5 Controladores', ...
        'FontSize', 16, 'FontWeight', 'bold');

    % Salvar figura
    saveas(figHandle, 'comparacao_5_controladores_massa_mola/simulacao_temporal_5_controladores.fig');
    saveas(figHandle, 'comparacao_5_controladores_massa_mola/simulacao_temporal_5_controladores.png');
end