function scriptComparacaoCompleta()
%   Script principal para comparação completa dos 5 controladores
%   Sistema mass-spring-damper com incertezas
%
%   Executa sequencialmente:
%   1. Carregamento do sistema
%   2. Síntese dos 5 controladores
%   3. Análise de estabilidade Monte Carlo
%   4. Simulação temporal com resposta ao degrau
%   5. Geração de tabelas LaTeX
%   6. Salvamento de todos os resultados

    fprintf('========================================================\n');
    fprintf('       COMPARAÇÃO COMPLETA DE 5 CONTROLADORES\n');
    fprintf('         Sistema Mass-Spring-Damper Incerto\n');
    fprintf('========================================================\n\n');

    % Adicionar caminhos necessários
    addpath('../HInf - Análise - Intervalar/funcoes');
    addpath('../HInf - Análise - Intervalar/funcoes/Diversas');
    addpath('../funcoes');
    addpath('..');

    % Mudar para o diretório da comparação
    % cd('comparacao_5_controladores_massa_mola');

    tempoInicial = datetime('now');

    %% ETAPA 1: Carregamento e verificação do sistema
    fprintf('ETAPA 1: Carregando sistema mass-spring-damper...\n');
    try
        sistema = genSaveDataComparacao();
        fprintf('✓ Sistema carregado com sucesso!\n');
        fprintf('  - Estados: %d | Entradas: %d | Saídas: %d\n', ...
            length(sistema.sys.A.inf), size(sistema.sys.B.inf, 2), size(sistema.sys.C.inf, 1));
        fprintf('  - Período de amostragem: %.3f s\n', sistema.parametros.h);
        fprintf('  - Incertezas: amortecedor ±5%%, mola ±2%%\n\n');
    catch ME
        fprintf('✗ ERRO no carregamento: %s\n', ME.message);
        return;
    end

    %% ETAPA 2: Síntese dos controladores
    fprintf('ETAPA 2: Executando síntese dos 5 controladores...\n');
    try
        controladores = comparacao5Controladores();
        fprintf('✓ Síntese de controladores concluída!\n\n');

        % Verificar quais sínteses foram bem-sucedidas
        nomes = fieldnames(controladores);
        sucessos = 0;
        for i = 1:length(nomes)
            if strcmp(nomes{i}, 'sistema') || strcmp(nomes{i}, 'nomeControladores')
                continue;
            end
            if controladores.(nomes{i}).factivel
                sucessos = sucessos + 1;
            end
        end
        fprintf('  Sínteses bem-sucedidas: %d/5\n\n', sucessos);

    catch ME
        fprintf('✗ ERRO na síntese: %s\n', ME.message);
        return;
    end

    %% ETAPA 3: Análise de estabilidade Monte Carlo
    fprintf('ETAPA 3: Executando análise de estabilidade (1000 realizações)...\n');
    try
        estabilidade = analisaEstabilidade5Controladores(controladores, 1000);
        fprintf('✓ Análise de estabilidade concluída!\n\n');
    catch ME
        fprintf('✗ ERRO na análise de estabilidade: %s\n', ME.message);
        return;
    end

    %% ETAPA 4: Simulação temporal
    fprintf('ETAPA 4: Executando simulação temporal (50 realizações)...\n');
    try
        simulacao = simulacaoTemporal5Controladores(controladores, 50, 15);
        fprintf('✓ Simulação temporal concluída!\n\n');
    catch ME
        fprintf('✗ ERRO na simulação temporal: %s\n', ME.message);
        return;
    end

    %% ETAPA 5: Geração de tabelas LaTeX
    fprintf('ETAPA 5: Gerando tabelas LaTeX para o artigo...\n');
    try
        % Compilar dados para as tabelas
        dadosTabelas = compilarDadosTabelas(controladores, estabilidade, simulacao);

        % Gerar tabelas
        tabelaComplexidade = gerarTabelaComplexidade(dadosTabelas);
        tabelaPerformance = gerarTabelaPerformance(dadosTabelas);

        % Salvar tabelas
        salvarTabelasLaTeX(tabelaComplexidade, tabelaPerformance);

        fprintf('✓ Tabelas LaTeX geradas e salvas!\n\n');
    catch ME
        fprintf('✗ ERRO na geração de tabelas: %s\n', ME.message);
        return;
    end

    %% ETAPA 6: Salvamento completo dos resultados
    fprintf('ETAPA 6: Salvando resultados completos...\n');
    try
        % Estrutura final com todos os resultados
        resultadosCompletos.sistema = sistema;
        resultadosCompletos.controladores = controladores;
        resultadosCompletos.estabilidade = estabilidade;
        resultadosCompletos.simulacao = simulacao;
        resultadosCompletos.tabelas.complexidade = tabelaComplexidade;
        resultadosCompletos.tabelas.performance = tabelaPerformance;
        resultadosCompletos.timestamp = datestr(now);

        % Salvar em arquivo MAT
        save('resultados_comparacao_completa.mat', 'resultadosCompletos', '-v7.3');

        % Gerar gráficos adicionais de alta qualidade
        gerarGraficosArtigo(resultadosCompletos);

        fprintf('✓ Resultados salvos em "resultados_comparacao_completa.mat"\n');
        fprintf('✓ Gráficos de alta qualidade gerados\n\n');
    catch ME
        fprintf('✗ ERRO no salvamento: %s\n', ME.message);
        return;
    end

    %% Resumo final
    tempoTotal = datetime('now') - tempoInicial;
    fprintf('========================================================\n');
    fprintf('                    RESUMO FINAL\n');
    fprintf('========================================================\n');
    fprintf('Tempo total de execução: %.1f segundos\n\n', seconds(tempoTotal));

    fprintf('Arquivos gerados:\n');
    fprintf('  📊 resultados_comparacao_completa.mat\n');
    fprintf('  📄 tabela_complexidade_massa_mola.tex\n');
    fprintf('  📄 tabela_performance_massa_mola.tex\n');
    fprintf('  📈 simulacao_temporal_5_controladores.fig/.png\n');
    fprintf('  📈 comparacao_metricas_artigo.fig/.png\n\n');

    % Exibir melhor controlador por critério
    fprintf('🏆 RANKING DOS CONTROLADORES:\n');
    exibirRanking(estabilidade, simulacao);

    fprintf('\n========================================================\n');
    fprintf('           COMPARAÇÃO FINALIZADA COM SUCESSO!\n');
    fprintf('========================================================\n');

    cd('..'); % Voltar ao diretório principal
end

%% Função para compilar dados das tabelas
function dados = compilarDadosTabelas(controladores, estabilidade, simulacao)

    nomes = {'Hibrido_Int', 'Hibrido_Poly', 'Amostrado_Int', 'Amostrado_Poly', 'Nominal'};
    nomesDisplay = {'Híbrido Int', 'Híbrido Poly', 'Amostrado Int', 'Amostrado Poly', 'Nominal'};

    dados.nomes = nomes;
    dados.nomesDisplay = nomesDisplay;
    dados.numControladores = length(nomes);

    % Dados de complexidade
    dados.ganhos = cell(1, dados.numControladores);
    dados.variaveis = zeros(1, dados.numControladores);
    dados.restricoes = zeros(1, dados.numControladores);
    dados.tempoSintese = zeros(1, dados.numControladores);
    dados.normaGarantida = zeros(1, dados.numControladores);

    % Dados de performance
    dados.taxaEstabilidade = estabilidade.taxaEstabilidade;
    dados.normaPiorCaso = estabilidade.normaPiorCaso;
    dados.normaMedia = estabilidade.normaMedia;
    dados.desviosPadrao = estabilidade.desviosPadrao;

    for i = 1:dados.numControladores
        nome = nomes{i};
        if controladores.(nome).factivel
            dados.ganhos{i} = controladores.(nome).K;
            dados.variaveis(i) = controladores.(nome).variaveis;
            dados.restricoes(i) = controladores.(nome).restricoes;
            dados.tempoSintese(i) = controladores.(nome).tempo;
            dados.normaGarantida(i) = controladores.(nome).mu;
        else
            dados.ganhos{i} = [];
            dados.variaveis(i) = 0;
            dados.restricoes(i) = 0;
            dados.tempoSintese(i) = NaN;
            dados.normaGarantida(i) = NaN;
        end
    end
end

%% Função para gerar tabela de complexidade
function tabela = gerarTabelaComplexidade(dados)

    tabela = sprintf('\\begin{table}[htbp]\n');
    tabela = [tabela, sprintf('\\centering\n')];
    tabela = [tabela, sprintf('\\caption{Complexidade Computacional - Sistema Mass-Spring-Damper}\n')];
    tabela = [tabela, sprintf('\\label{tab:complexidade_massa_mola}\n')];
    tabela = [tabela, sprintf('\\begin{tabular}{lccccc}\n')];
    tabela = [tabela, sprintf('\\hline\n')];
    tabela = [tabela, sprintf('Método & Ganho $K$ & Variáveis & Restrições & Tempo (s) & $\\gamma$ \\\\\n')];
    tabela = [tabela, sprintf('\\hline\n')];

    for i = 1:dados.numControladores
        nome = dados.nomesDisplay{i};

        % Ganho
        if ~isempty(dados.ganhos{i})
            if length(dados.ganhos{i}) <= 4
                ganhoStr = sprintf('$[%s]$', strjoin(arrayfun(@(x) sprintf('%.2f', x), dados.ganhos{i}, 'UniformOutput', false), ', '));
            else
                ganhoStr = sprintf('$[%.2f, \\ldots]$', dados.ganhos{i}(1));
            end
        else
            ganhoStr = 'N/A';
        end

        % Variáveis e restrições
        if dados.variaveis(i) == 0
            varStr = 'N/A';
            restStr = 'N/A';
        else
            varStr = sprintf('%d', dados.variaveis(i));
            restStr = sprintf('%d', dados.restricoes(i));
        end

        % Tempo
        if isnan(dados.tempoSintese(i))
            tempoStr = 'N/A';
        else
            tempoStr = sprintf('%.3f', dados.tempoSintese(i));
        end

        % Norma garantida
        if isnan(dados.normaGarantida(i))
            normaStr = 'N/A';
        else
            normaStr = sprintf('%.3f', dados.normaGarantida(i));
        end

        tabela = [tabela, sprintf('%s & %s & %s & %s & %s & %s \\\\\n', ...
                  nome, ganhoStr, varStr, restStr, tempoStr, normaStr)];
    end

    tabela = [tabela, sprintf('\\hline\n')];
    tabela = [tabela, sprintf('\\end{tabular}\n')];
    tabela = [tabela, sprintf('\\end{table}\n')];
end

%% Função para gerar tabela de performance
function tabela = gerarTabelaPerformance(dados)

    tabela = sprintf('\\begin{table}[htbp]\n');
    tabela = [tabela, sprintf('\\centering\n')];
    tabela = [tabela, sprintf('\\caption{Análise de Performance Monte Carlo (1000 realizações)}\n')];
    tabela = [tabela, sprintf('\\label{tab:performance_massa_mola}\n')];
    tabela = [tabela, sprintf('\\begin{tabular}{lcccc}\n')];
    tabela = [tabela, sprintf('\\hline\n')];
    tabela = [tabela, sprintf('Método & Taxa Estab. (\\%%) & Pior Caso & Norma Média & Desvio Padrão \\\\\n')];
    tabela = [tabela, sprintf('\\hline\n')];

    for i = 1:dados.numControladores
        nome = dados.nomesDisplay{i};

        % Taxa de estabilidade
        taxaStr = sprintf('%.1f', dados.taxaEstabilidade(i));

        % Pior caso
        if isinf(dados.normaPiorCaso(i)) || dados.normaPiorCaso(i) > 1000
            piorCasoStr = '$\\infty$';
        else
            piorCasoStr = sprintf('%.3f', dados.normaPiorCaso(i));
        end

        % Norma média
        if isinf(dados.normaMedia(i)) || dados.normaMedia(i) > 1000
            mediaStr = '$\\infty$';
        else
            mediaStr = sprintf('%.3f', dados.normaMedia(i));
        end

        % Desvio padrão
        if dados.desviosPadrao(i) == 0 || isnan(dados.desviosPadrao(i))
            desvioStr = 'N/A';
        else
            desvioStr = sprintf('%.3f', dados.desviosPadrao(i));
        end

        tabela = [tabela, sprintf('%s & %s & %s & %s & %s \\\\\n', ...
                  nome, taxaStr, piorCasoStr, mediaStr, desvioStr)];
    end

    tabela = [tabela, sprintf('\\hline\n')];
    tabela = [tabela, sprintf('\\end{tabular}\n')];
    tabela = [tabela, sprintf('\\end{table}\n')];
end

%% Função para salvar tabelas LaTeX
function salvarTabelasLaTeX(tabelaComplexidade, tabelaPerformance)

    % Salvar tabela de complexidade
    fileID = fopen('tabela_complexidade_massa_mola.tex', 'w');
    fprintf(fileID, '%s', tabelaComplexidade);
    fclose(fileID);

    % Salvar tabela de performance
    fileID = fopen('tabela_performance_massa_mola.tex', 'w');
    fprintf(fileID, '%s', tabelaPerformance);
    fclose(fileID);
end

%% Função para gerar gráficos adicionais para o artigo
function gerarGraficosArtigo(resultados)

    % Gráfico comparativo de métricas
    figure('Position', [100, 100, 1200, 800]);

    nomes = resultados.estabilidade.nomes;
    taxas = resultados.estabilidade.taxaEstabilidade;
    normas = resultados.estabilidade.normaMedia;
    normas(isinf(normas)) = 100; % Limitar para visualização

    % Subplot 1: Taxa de estabilidade
    subplot(2, 2, 1);
    bar(taxas, 'FaceColor', [0.3, 0.6, 0.9]);
    set(gca, 'XTickLabel', nomes);
    xtickangle(45);
    ylabel('Taxa de Estabilidade (%)');
    title('Taxa de Estabilidade por Controlador');
    grid on;

    % Subplot 2: Norma média H-infinity
    subplot(2, 2, 2);
    bar(normas, 'FaceColor', [0.9, 0.4, 0.3]);
    set(gca, 'XTickLabel', nomes);
    xtickangle(45);
    ylabel('Norma H∞ Média');
    title('Norma H∞ Média por Controlador');
    grid on;

    % Subplot 3: Tempo de síntese
    tempos = zeros(1, length(nomes));
    for i = 1:length(nomes)
        nome = resultados.controladores.nomeControladores{i};
        nome = strrep(nome, ' ', '_');
        nome = strrep(nome, 'í', 'i');
        nome = strrep(nome, 'ú', 'u');
        if resultados.controladores.(nome).factivel
            tempos(i) = resultados.controladores.(nome).tempo;
        end
    end

    subplot(2, 2, 3);
    bar(tempos, 'FaceColor', [0.4, 0.8, 0.4]);
    set(gca, 'XTickLabel', nomes);
    xtickangle(45);
    ylabel('Tempo de Síntese (s)');
    title('Tempo de Síntese por Controlador');
    grid on;

    % Subplot 4: Comparação geral (radar plot simplificado)
    subplot(2, 2, 4);
    metricas_norm = [taxas/100; (100-normas)/100; (1-tempos/max(tempos))]'; % Normalizar para [0,1]
    plot(metricas_norm', 'LineWidth', 2, 'Marker', 'o');
    legend(nomes, 'Location', 'best');
    set(gca, 'XTickLabel', {'Estabilidade', 'Performance', 'Rapidez'});
    ylabel('Métrica Normalizada');
    title('Comparação Geral (Normalizada)');
    grid on;

    sgtitle('Análise Comparativa dos 5 Controladores', 'FontSize', 16, 'FontWeight', 'bold');

    % Salvar
    saveas(gcf, 'comparacao_metricas_artigo.fig');
    saveas(gcf, 'comparacao_metricas_artigo.png');
end

%% Função para exibir ranking
function exibirRanking(estabilidade, simulacao)

    nomes = estabilidade.nomes;
    [~, rankEstab] = sort(estabilidade.taxaEstabilidade, 'descend');
    [~, rankNorma] = sort(estabilidade.normaMedia, 'ascend');

    fprintf('  Por Estabilidade: ');
    for i = 1:min(3, length(rankEstab))
        fprintf('%s (%.1f%%) ', nomes{rankEstab(i)}, estabilidade.taxaEstabilidade(rankEstab(i)));
    end
    fprintf('\n');

    fprintf('  Por Performance:  ');
    for i = 1:min(3, length(rankNorma))
        if ~isinf(estabilidade.normaMedia(rankNorma(i)))
            fprintf('%s (%.3f) ', nomes{rankNorma(i)}, estabilidade.normaMedia(rankNorma(i)));
        end
    end
    fprintf('\n');
end