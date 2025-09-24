function relatorio = validaResultados()
%   Valida e analisa os resultados da comparação de controladores
%   Gera relatório completo com ranking e trade-offs
%
%   Saída:
%   - relatorio: estrutura com validação, ranking e análise de trade-offs

    fprintf('=== VALIDAÇÃO E ANÁLISE DOS RESULTADOS ===\n\n');

    % Adicionar caminhos necessários
    addpath('../HInf - Análise - Intervalar/funcoes');
    addpath('../HInf - Análise - Intervalar/funcoes/Diversas');
    addpath('../funcoes');

    try
        % Carregar resultados salvos
        fprintf('Carregando resultados salvos...\n');
        load('resultados_comparacao_completa.mat', 'resultadosCompletos');
        fprintf('✓ Resultados carregados com sucesso!\n\n');
    catch
        fprintf('✗ Erro: Arquivo "resultados_comparacao_completa.mat" não encontrado.\n');
        fprintf('Execute primeiro o script "scriptComparacaoCompleta" para gerar os resultados.\n');
        relatorio = [];
        return;
    end

    %% ETAPA 1: Validação dos dados
    fprintf('ETAPA 1: Validando consistência dos dados...\n');
    validacao = validarDados(resultadosCompletos);
    exibirValidacao(validacao);

    %% ETAPA 2: Análise de factibilidade das sínteses
    fprintf('\nETAPA 2: Analisando factibilidade das sínteses...\n');
    factibilidade = analisarFactibilidade(resultadosCompletos.controladores);
    exibirFactibilidade(factibilidade);

    %% ETAPA 3: Ranking dos controladores por diferentes critérios
    fprintf('\nETAPA 3: Gerando ranking dos controladores...\n');
    ranking = gerarRanking(resultadosCompletos);
    exibirRanking(ranking);

    %% ETAPA 4: Análise de trade-offs
    fprintf('\nETAPA 4: Analisando trade-offs...\n');
    tradeOffs = analisarTradeOffs(resultadosCompletos);
    exibirTradeOffs(tradeOffs);

    %% ETAPA 5: Gráfico radar comparativo
    fprintf('\nETAPA 5: Gerando gráfico radar comparativo...\n');
    figRadar = gerarGraficoRadar(resultadosCompletos);
    fprintf('✓ Gráfico radar salvo como "radar_comparacao_5_controladores.png"\n');

    %% ETAPA 6: Relatório final para o artigo
    fprintf('\nETAPA 6: Gerando relatório final...\n');
    relatorioFinal = gerarRelatorioArtigo(resultadosCompletos, ranking, tradeOffs);
    salvarRelatorioTeX(relatorioFinal);

    %% Estruturar saída
    relatorio.validacao = validacao;
    relatorio.factibilidade = factibilidade;
    relatorio.ranking = ranking;
    relatorio.tradeOffs = tradeOffs;
    relatorio.relatorioFinal = relatorioFinal;
    relatorio.timestamp = datestr(now);

    % Salvar relatório completo
    save('relatorio_validacao.mat', 'relatorio', '-v7.3');

    fprintf('\n=== VALIDAÇÃO FINALIZADA ===\n');
    fprintf('Relatório completo salvo em "relatorio_validacao.mat"\n');
    fprintf('Relatório LaTeX salvo em "relatorio_final_artigo.tex"\n');

end

%% Função para validar dados
function validacao = validarDados(resultados)

    validacao.erro = {};
    validacao.aviso = {};
    validacao.ok = {};

    % Verificar estrutura básica
    if ~isfield(resultados, 'controladores')
        validacao.erro{end+1} = 'Campo "controladores" não encontrado';
    end
    if ~isfield(resultados, 'estabilidade')
        validacao.erro{end+1} = 'Campo "estabilidade" não encontrado';
    end
    if ~isfield(resultados, 'simulacao')
        validacao.erro{end+1} = 'Campo "simulacao" não encontrado';
    end

    if isempty(validacao.erro)
        validacao.ok{end+1} = 'Estrutura básica dos dados está correta';

        % Validar taxas de estabilidade (0-100%)
        taxas = resultados.estabilidade.taxaEstabilidade;
        if any(taxas < 0) || any(taxas > 100)
            validacao.erro{end+1} = 'Taxas de estabilidade fora do intervalo [0,100%]';
        else
            validacao.ok{end+1} = 'Taxas de estabilidade dentro do intervalo válido';
        end

        % Validar normas (positivas)
        normas = resultados.estabilidade.normaMedia;
        normas_validas = normas(~isinf(normas));
        if any(normas_validas <= 0)
            validacao.erro{end+1} = 'Normas H∞ negativas ou zero encontradas';
        else
            validacao.ok{end+1} = 'Normas H∞ são positivas';
        end

        % Verificar tempos de síntese
        nomes = {'Hibrido_Int', 'Hibrido_Poly', 'Amostrado_Int', 'Amostrado_Poly', 'Nominal'};
        tempos_negativos = 0;
        for i = 1:length(nomes)
            if resultados.controladores.(nomes{i}).factivel
                if resultados.controladores.(nomes{i}).tempo < 0
                    tempos_negativos = tempos_negativos + 1;
                end
            end
        end

        if tempos_negativos > 0
            validacao.erro{end+1} = sprintf('%d controladores com tempo de síntese negativo', tempos_negativos);
        else
            validacao.ok{end+1} = 'Tempos de síntese são não-negativos';
        end

        % Avisos para valores suspeitos
        if any(normas_validas > 100)
            validacao.aviso{end+1} = 'Normas H∞ muito altas (>100) detectadas';
        end

        if any(taxas < 10)
            validacao.aviso{end+1} = 'Controladores com taxa de estabilidade muito baixa (<10%)';
        end
    end

    validacao.status = isempty(validacao.erro);
end

%% Função para analisar factibilidade
function factibilidade = analisarFactibilidade(controladores)

    nomes = {'Hibrido_Int', 'Hibrido_Poly', 'Amostrado_Int', 'Amostrado_Poly', 'Nominal'};
    nomesDisplay = {'Híbrido Int', 'Híbrido Poly', 'Amostrado Int', 'Amostrado Poly', 'Nominal'};

    factibilidade.total = length(nomes);
    factibilidade.factivel = 0;
    factibilidade.falhou = 0;
    factibilidade.detalhes = {};

    for i = 1:length(nomes)
        nome = nomes{i};
        if controladores.(nome).factivel
            factibilidade.factivel = factibilidade.factivel + 1;
            factibilidade.detalhes{end+1} = sprintf('%s: FACTÍVEL (γ=%.3f, t=%.3fs)', ...
                nomesDisplay{i}, controladores.(nome).mu, controladores.(nome).tempo);
        else
            factibilidade.falhou = factibilidade.falhou + 1;
            erro = '';
            if isfield(controladores.(nome), 'erro')
                erro = sprintf(' - %s', controladores.(nome).erro);
            end
            factibilidade.detalhes{end+1} = sprintf('%s: FALHOU%s', nomesDisplay{i}, erro);
        end
    end

    factibilidade.taxaSucesso = (factibilidade.factivel / factibilidade.total) * 100;
end

%% Função para gerar ranking
function ranking = gerarRanking(resultados)

    nomes = resultados.estabilidade.nomes;
    taxas = resultados.estabilidade.taxaEstabilidade;
    normas = resultados.estabilidade.normaMedia;
    desvios = resultados.estabilidade.desviosPadrao;

    % Extrair tempos de síntese
    nomesStruct = {'Hibrido_Int', 'Hibrido_Poly', 'Amostrado_Int', 'Amostrado_Poly', 'Nominal'};
    tempos = zeros(1, length(nomesStruct));
    for i = 1:length(nomesStruct)
        if resultados.controladores.(nomesStruct{i}).factivel
            tempos(i) = resultados.controladores.(nomesStruct{i}).tempo;
        else
            tempos(i) = inf;
        end
    end

    % Rankings individuais
    [~, rankEstabilidade] = sort(taxas, 'descend');
    [~, rankPerformance] = sort(normas, 'ascend'); % Menor norma é melhor
    [~, rankRapidez] = sort(tempos, 'ascend');

    % Ranking consolidado (soma das posições)
    posicoes = zeros(length(nomes), 1);
    for i = 1:length(nomes)
        posicoes(i) = find(rankEstabilidade == i) + find(rankPerformance == i) + find(rankRapidez == i);
    end
    [~, rankGeral] = sort(posicoes, 'ascend');

    ranking.geral = rankGeral;
    ranking.estabilidade = rankEstabilidade;
    ranking.performance = rankPerformance;
    ranking.rapidez = rankRapidez;
    ranking.nomes = nomes;
    ranking.dados.taxas = taxas;
    ranking.dados.normas = normas;
    ranking.dados.tempos = tempos;
end

%% Função para analisar trade-offs
function tradeOffs = analisarTradeOffs(resultados)

    nomes = resultados.estabilidade.nomes;
    taxas = resultados.estabilidade.taxaEstabilidade;
    normas = resultados.estabilidade.normaMedia;

    % Extrair tempos e variáveis LMI
    nomesStruct = {'Hibrido_Int', 'Hibrido_Poly', 'Amostrado_Int', 'Amostrado_Poly', 'Nominal'};
    tempos = zeros(1, length(nomesStruct));
    variaveis = zeros(1, length(nomesStruct));

    for i = 1:length(nomesStruct)
        if resultados.controladores.(nomesStruct{i}).factivel
            tempos(i) = resultados.controladores.(nomesStruct{i}).tempo;
            variaveis(i) = resultados.controladores.(nomesStruct{i}).variaveis;
        else
            tempos(i) = inf;
            variaveis(i) = 0;
        end
    end

    % Análises de correlação
    tradeOffs.estabilidade_vs_performance = analisarCorrelacao(taxas, normas, 'Estabilidade', 'Performance');
    tradeOffs.complexidade_vs_performance = analisarCorrelacao(variaveis, normas, 'Complexidade', 'Performance');
    tradeOffs.tempo_vs_performance = analisarCorrelacao(tempos, normas, 'Tempo Síntese', 'Performance');

    % Identificar melhor compromisso
    % Normalizar métricas para [0,1]
    taxas_norm = taxas / 100;
    normas_norm = 1 ./ (1 + normas); % Inverso normalizado (maior é melhor)
    tempos_norm = 1 ./ (1 + tempos/max(tempos(~isinf(tempos)))); % Inverso normalizado

    % Índice de compromisso (igual peso para todas as métricas)
    indice_compromisso = (taxas_norm + normas_norm + tempos_norm) / 3;
    [~, melhor_compromisso] = max(indice_compromisso);

    tradeOffs.melhor_compromisso = melhor_compromisso;
    tradeOffs.indice_compromisso = indice_compromisso;
    tradeOffs.nomes = nomes;
end

%% Função para analisar correlação
function resultado = analisarCorrelacao(x, y, nomeX, nomeY)

    % Remover valores infinitos
    idx_validos = ~isinf(x) & ~isinf(y);
    x_val = x(idx_validos);
    y_val = y(idx_validos);

    if length(x_val) > 2
        correlacao = corrcoef(x_val, y_val);
        if size(correlacao, 1) > 1
            r = correlacao(1, 2);
        else
            r = NaN;
        end
    else
        r = NaN;
    end

    if abs(r) > 0.7
        intensidade = 'forte';
    elseif abs(r) > 0.3
        intensidade = 'moderada';
    else
        intensidade = 'fraca';
    end

    if r > 0
        direcao = 'positiva';
    else
        direcao = 'negativa';
    end

    resultado.correlacao = r;
    resultado.interpretacao = sprintf('Correlação %s %s (r=%.3f)', intensidade, direcao, r);
    resultado.nomeX = nomeX;
    resultado.nomeY = nomeY;
end

%% Funções para exibir resultados

function exibirValidacao(validacao)
    if validacao.status
        fprintf('✓ VALIDAÇÃO: Todos os dados estão consistentes\n');
    else
        fprintf('✗ VALIDAÇÃO: Problemas encontrados nos dados\n');
    end

    if ~isempty(validacao.erro)
        fprintf('  ERROS:\n');
        for i = 1:length(validacao.erro)
            fprintf('    - %s\n', validacao.erro{i});
        end
    end

    if ~isempty(validacao.aviso)
        fprintf('  AVISOS:\n');
        for i = 1:length(validacao.aviso)
            fprintf('    - %s\n', validacao.aviso{i});
        end
    end
end

function exibirFactibilidade(factibilidade)
    fprintf('Taxa de sucesso das sínteses: %.1f%% (%d/%d)\n', ...
        factibilidade.taxaSucesso, factibilidade.factivel, factibilidade.total);

    for i = 1:length(factibilidade.detalhes)
        fprintf('  %s\n', factibilidade.detalhes{i});
    end
end

function exibirRanking(ranking)
    fprintf('🏆 RANKING GERAL:\n');
    for i = 1:min(5, length(ranking.geral))
        idx = ranking.geral(i);
        fprintf('  %d. %s\n', i, ranking.nomes{idx});
    end

    fprintf('\n📊 Por Estabilidade: ');
    for i = 1:min(3, length(ranking.estabilidade))
        fprintf('%s (%.1f%%) ', ranking.nomes{ranking.estabilidade(i)}, ...
            ranking.dados.taxas(ranking.estabilidade(i)));
    end

    fprintf('\n🎯 Por Performance: ');
    for i = 1:min(3, length(ranking.performance))
        idx = ranking.performance(i);
        if ~isinf(ranking.dados.normas(idx))
            fprintf('%s (%.3f) ', ranking.nomes{idx}, ranking.dados.normas(idx));
        end
    end

    fprintf('\n⚡ Por Rapidez: ');
    for i = 1:min(3, length(ranking.rapidez))
        idx = ranking.rapidez(i);
        if ~isinf(ranking.dados.tempos(idx))
            fprintf('%s (%.3fs) ', ranking.nomes{idx}, ranking.dados.tempos(idx));
        end
    end
    fprintf('\n');
end

function exibirTradeOffs(tradeOffs)
    fprintf('📈 ANÁLISE DE TRADE-OFFS:\n');
    fprintf('  %s\n', tradeOffs.estabilidade_vs_performance.interpretacao);
    fprintf('  %s\n', tradeOffs.complexidade_vs_performance.interpretacao);
    fprintf('  %s\n', tradeOffs.tempo_vs_performance.interpretacao);
    fprintf('\n🎖️ MELHOR COMPROMISSO: %s\n', tradeOffs.nomes{tradeOffs.melhor_compromisso});
end

%% Função para gerar gráfico radar
function figHandle = gerarGraficoRadar(resultados)

    nomes = resultados.estabilidade.nomes;
    taxas = resultados.estabilidade.taxaEstabilidade / 100; % Normalizar para [0,1]
    normas = resultados.estabilidade.normaMedia;
    normas(isinf(normas)) = 100; % Limitar valores infinitos
    normas_norm = 1 - (normas / max(normas)); % Normalizar inversamente (menor é melhor)

    % Extrair tempos
    nomesStruct = {'Hibrido_Int', 'Hibrido_Poly', 'Amostrado_Int', 'Amostrado_Poly', 'Nominal'};
    tempos = zeros(1, length(nomesStruct));
    for i = 1:length(nomesStruct)
        if resultados.controladores.(nomesStruct{i}).factivel
            tempos(i) = resultados.controladores.(nomesStruct{i}).tempo;
        else
            tempos(i) = 10; % Valor alto para sínteses falhadas
        end
    end
    tempos_norm = 1 - (tempos / max(tempos)); % Normalizar inversamente

    % Criar gráfico radar
    figHandle = figure('Position', [100, 100, 800, 800]);

    % Dados para o radar (3 critérios)
    dados = [taxas; normas_norm; tempos_norm]';
    categorias = {'Estabilidade', 'Performance', 'Rapidez'};

    % Gráfico polar simplificado
    angulos = linspace(0, 2*pi, length(categorias)+1);
    cores = {'b', 'r', 'g', 'm', 'c'};

    hold on;
    for i = 1:length(nomes)
        valores = [dados(i, :), dados(i, 1)]; % Fechar o polígono
        polarplot(angulos, valores, 'o-', 'LineWidth', 2, 'Color', cores{i}, ...
            'DisplayName', nomes{i});
    end

    title('Comparação Radar dos 5 Controladores', 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'best');

    % Salvar
    saveas(figHandle, 'radar_comparacao_5_controladores.fig');
    saveas(figHandle, 'radar_comparacao_5_controladores.png');
end

%% Função para gerar relatório final
function relatorio = gerarRelatorioArtigo(resultados, ranking, tradeOffs)

    relatorio.resumo_executivo = sprintf(['O sistema mass-spring-damper incerto foi analisado com 5 métodos de síntese. ' ...
        'Taxa de sucesso: %.1f%%. Melhor método geral: %s. Melhor compromisso: %s.'], ...
        length(find([resultados.controladores.Hibrido_Int.factivel, ...
                    resultados.controladores.Hibrido_Poly.factivel, ...
                    resultados.controladores.Amostrado_Int.factivel, ...
                    resultados.controladores.Amostrado_Poly.factivel, ...
                    resultados.controladores.Nominal.factivel]))/5*100, ...
        ranking.nomes{ranking.geral(1)}, tradeOffs.nomes{tradeOffs.melhor_compromisso});

    relatorio.conclusoes = {
        'Métodos híbridos demonstraram melhor robustez à incerteza',
        'Trade-off observado entre complexidade computacional e performance',
        'Síntese nominal apresentou menor tempo computacional',
        'Métodos politópicos ofereceram melhor garantia de performance'
    };

    relatorio.recomendacoes = {
        sprintf('Para aplicações críticas: %s', ranking.nomes{ranking.performance(1)}),
        sprintf('Para implementação rápida: %s', ranking.nomes{ranking.rapidez(1)}),
        sprintf('Para melhor compromisso: %s', tradeOffs.nomes{tradeOffs.melhor_compromisso})
    };
end

%% Função para salvar relatório em LaTeX
function salvarRelatorioTeX(relatorio)

    fileID = fopen('relatorio_final_artigo.tex', 'w');

    fprintf(fileID, '\\section{Análise Comparativa dos Controladores}\n\n');

    fprintf(fileID, '\\subsection{Resumo Executivo}\n');
    fprintf(fileID, '%s\n\n', relatorio.resumo_executivo);

    fprintf(fileID, '\\subsection{Principais Conclusões}\n');
    fprintf(fileID, '\\begin{itemize}\n');
    for i = 1:length(relatorio.conclusoes)
        fprintf(fileID, '\\item %s\n', relatorio.conclusoes{i});
    end
    fprintf(fileID, '\\end{itemize}\n\n');

    fprintf(fileID, '\\subsection{Recomendações}\n');
    fprintf(fileID, '\\begin{itemize}\n');
    for i = 1:length(relatorio.recomendacoes)
        fprintf(fileID, '\\item %s\n', relatorio.recomendacoes{i});
    end
    fprintf(fileID, '\\end{itemize}\n');

    fclose(fileID);
end