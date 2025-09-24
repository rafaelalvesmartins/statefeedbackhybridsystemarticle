function resultado = comparacaoCompletaControladores()
%   Comparação completa dos 5 métodos de síntese para o Exemplo 6
%
%   Executa sequencialmente:
%   1. Carregamento do sistema
%   2. Síntese dos controladores
%   3. Teste de estabilidade Monte Carlo
%   4. Simulação de resposta ao degrau
%   5. Geração de tabelas LaTeX
%
%   Saída:
%   - resultado: struct com todos os resultados e tabelas LaTeX

    fprintf('=== COMPARAÇÃO COMPLETA DE CONTROLADORES - EXEMPLO 6 ===\n\n');

    % Adicionar caminhos necessários
    addpath('../HInf - Análise - Intervalar/funcoes');
    addpath('../HInf - Análise - Intervalar/funcoes/Diversas');
    addpath('../funcoes')

    tempoInicial = datetime('now');

    %% 1. Carregamento do sistema
    fprintf('ETAPA 1: Carregando sistema do Exemplo 6...\n');
    try
        sistema = genSaveDataEx6Comparacao();
        fprintf('✓ Sistema carregado com sucesso!\n\n');
    catch ME
        fprintf('✗ ERRO no carregamento do sistema: %s\n', ME.message);
        resultado.erro = ME.message;
        return;
    end

    %% 2. Síntese dos controladores
    fprintf('ETAPA 2: Executando síntese dos 5 controladores...\n');
    try
        controladores = chamaComparacao5Controladores();
        fprintf('✓ Síntese de controladores concluída!\n\n');

        % Exibir resumo dos ganhos
        fprintf('Ganhos obtidos:\n');
        metodosNomes = fieldnames(controladores.metodos);
        for i = 1:length(metodosNomes)
            metodo = controladores.metodos.(metodosNomes{i});
            if metodo.sucesso
                fprintf('  %s: K = [%.4f, %.4f]\n', metodosNomes{i}, metodo.K(1), metodo.K(2));
            else
                fprintf('  %s: FALHOU\n', metodosNomes{i});
            end
        end
        fprintf('\n');

    catch ME
        fprintf('✗ ERRO na síntese de controladores: %s\n', ME.message);
        resultado.erro = ME.message;
        return;
    end

    %% 3. Teste de estabilidade Monte Carlo
    fprintf('ETAPA 3: Executando análise de estabilidade Monte Carlo...\n');
    try
        estabilidade = avaliaMonteCarloControladores(controladores.metodos, 1000);
        fprintf('✓ Análise de estabilidade concluída!\n\n');
    catch ME
        fprintf('✗ ERRO na análise de estabilidade: %s\n', ME.message);
        resultado.erro = ME.message;
        return;
    end

    %% 4. Simulação de resposta ao degrau
    fprintf('ETAPA 4: Executando simulação de resposta ao degrau...\n');
    try
        simulacao = simulaRespostaDegrauControladores(controladores.metodos, 100, 10);
        fprintf('✓ Simulação de resposta ao degrau concluída!\n\n');
    catch ME
        fprintf('✗ ERRO na simulação de resposta: %s\n', ME.message);
        resultado.erro = ME.message;
        return;
    end

    %% 5. Geração de tabelas LaTeX
    fprintf('ETAPA 5: Gerando tabelas LaTeX...\n');
    try
        % Compilar dados para as tabelas
        dadosCompletos = compilarDadosTabelas(controladores, estabilidade, simulacao, sistema);

        % Gerar tabelas LaTeX
        tabelaComplexidade = gerarTabelaComplexidade(dadosCompletos);
        tabelaPerformance = gerarTabelaPerformance(dadosCompletos);

        % Salvar tabelas em arquivos
        salvarTabelasTeX(tabelaComplexidade, tabelaPerformance);

        fprintf('✓ Tabelas LaTeX geradas e salvas!\n\n');
    catch ME
        fprintf('✗ ERRO na geração de tabelas: %s\n', ME.message);
        resultado.erro = ME.message;
        return;
    end

    %% Compilação dos resultados finais
    tempoTotal = datetime('now') - tempoInicial;

    resultado.sistema = sistema;
    resultado.controladores = controladores;
    resultado.estabilidade = estabilidade;
    resultado.simulacao = simulacao;
    resultado.tabelas.complexidade = tabelaComplexidade;
    resultado.tabelas.performance = tabelaPerformance;
    resultado.tempoTotal = seconds(tempoTotal);

    %% Exibição do resumo final
    fprintf('=== RESUMO FINAL ===\n');
    fprintf('Tempo total de execução: %.1f segundos\n', seconds(tempoTotal));
    fprintf('\nTabelas LaTeX geradas:\n');
    fprintf('  - tabela_complexidade_ex6.tex\n');
    fprintf('  - tabela_performance_ex6.tex\n');

    fprintf('\n=== COMPARAÇÃO COMPLETA FINALIZADA COM SUCESSO ===\n');

end

%% Função para compilar dados das tabelas
function dados = compilarDadosTabelas(controladores, estabilidade, simulacao, sistema)

    % Nomes dos controladores na ordem correta
    nomes = {'Nominal', 'Discreto Int', 'Discreto Poly', 'Híbrido Int', 'Híbrido Poly'};
    metodosStruct = {'nominal', 'discretoIntervalar', 'discretoPolitopico', 'hibridoIntervalar', 'hibridoPolitopico'};

    dados.nomes = nomes;
    dados.numControladores = length(nomes);

    % Extrair dados de complexidade
    dados.variaveisLMI = [2, 0, 0, 8, 16];  % Estimativas baseadas no tipo de síntese
    dados.tempoSintese = zeros(1, dados.numControladores);
    dados.normaGarantida = zeros(1, dados.numControladores);
    dados.ganhos = cell(1, dados.numControladores);

    for i = 1:dados.numControladores
        metodo = controladores.metodos.(metodosStruct{i});
        if metodo.sucesso
            if isfield(metodo, 'sintese') && isfield(metodo.sintese, 'tempoExecucao')
                dados.tempoSintese(i) = metodo.sintese.tempoExecucao;
            else
                dados.tempoSintese(i) = NaN;
            end

            if isfield(metodo, 'sintese') && isfield(metodo.sintese, 'gamma')
                dados.normaGarantida(i) = metodo.sintese.gamma;
            else
                dados.normaGarantida(i) = NaN;
            end

            dados.ganhos{i} = metodo.K;
        else
            dados.tempoSintese(i) = NaN;
            dados.normaGarantida(i) = NaN;
            dados.ganhos{i} = [];
        end
    end

    % Extrair dados de performance
    dados.taxaEstabilidade = estabilidade.percentualEstavel;
    dados.piorCaso = estabilidade.normaPiorCaso;
    dados.normaMedia = estabilidade.normaMedia;

    % Calcular desvio padrão das normas L2 da simulação
    dados.desvioNormaL2 = zeros(1, dados.numControladores);
    for i = 1:dados.numControladores
        if isfield(simulacao.normasL2_stats, sprintf('ctrl%d', i))
            dados.desvioNormaL2(i) = simulacao.normasL2_stats.(sprintf('ctrl%d', i)).std;
        else
            dados.desvioNormaL2(i) = NaN;
        end
    end

end

%% Função para gerar tabela de complexidade
function tabela = gerarTabelaComplexidade(dados)

    tabela = sprintf('\\begin{table}[htbp]\n');
    tabela = [tabela, sprintf('\\centering\n')];
    tabela = [tabela, sprintf('\\caption{Complexidade Computacional dos Métodos de Síntese - Exemplo 6}\n')];
    tabela = [tabela, sprintf('\\label{tab:complexidade_ex6}\n')];
    tabela = [tabela, sprintf('\\begin{tabular}{lcccc}\n')];
    tabela = [tabela, sprintf('\\hline\n')];
    tabela = [tabela, sprintf('Método & Ganho $K$ & Variáveis LMI & Tempo (s) & $\\gamma$ Garantido \\\\\n')];
    tabela = [tabela, sprintf('\\hline\n')];

    for i = 1:dados.numControladores
        nome = dados.nomes{i};

        % Ganho
        if ~isempty(dados.ganhos{i})
            ganhoStr = sprintf('[%.3f, %.3f]', dados.ganhos{i}(1), dados.ganhos{i}(2));
        else
            ganhoStr = 'N/A';
        end

        % Variáveis LMI
        varLMI = dados.variaveisLMI(i);
        if varLMI == 0
            varLMIStr = 'N/A';
        else
            varLMIStr = sprintf('%d', varLMI);
        end

        % Tempo de síntese
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

        tabela = [tabela, sprintf('%s & %s & %s & %s & %s \\\\\n', ...
                  nome, ganhoStr, varLMIStr, tempoStr, normaStr)];
    end

    tabela = [tabela, sprintf('\\hline\n')];
    tabela = [tabela, sprintf('\\end{tabular}\n')];
    tabela = [tabela, sprintf('\\end{table}\n')];

end

%% Função para gerar tabela de performance
function tabela = gerarTabelaPerformance(dados)

    tabela = sprintf('\\begin{table}[htbp]\n');
    tabela = [tabela, sprintf('\\centering\n')];
    tabela = [tabela, sprintf('\\caption{Análise de Performance Monte Carlo - Exemplo 6 (1000 realizações)}\n')];
    tabela = [tabela, sprintf('\\label{tab:performance_ex6}\n')];
    tabela = [tabela, sprintf('\\begin{tabular}{lcccc}\n')];
    tabela = [tabela, sprintf('\\hline\n')];
    tabela = [tabela, sprintf('Método & Taxa Estab. (\\%%) & Pior Caso & Norma Média & Desvio Padrão \\\\\n')];
    tabela = [tabela, sprintf('\\hline\n')];

    for i = 1:dados.numControladores
        nome = dados.nomes{i};

        % Taxa de estabilidade
        taxaStr = sprintf('%.1f', dados.taxaEstabilidade(i));

        % Pior caso
        if isinf(dados.piorCaso(i))
            piorCasoStr = '$\\infty$';
        else
            piorCasoStr = sprintf('%.3f', dados.piorCaso(i));
        end

        % Norma média
        if isinf(dados.normaMedia(i))
            mediaStr = '$\\infty$';
        else
            mediaStr = sprintf('%.3f', dados.normaMedia(i));
        end

        % Desvio padrão
        if isnan(dados.desvioNormaL2(i))
            desvioStr = 'N/A';
        else
            desvioStr = sprintf('%.3f', dados.desvioNormaL2(i));
        end

        tabela = [tabela, sprintf('%s & %s & %s & %s & %s \\\\\n', ...
                  nome, taxaStr, piorCasoStr, mediaStr, desvioStr)];
    end

    tabela = [tabela, sprintf('\\hline\n')];
    tabela = [tabela, sprintf('\\end{tabular}\n')];
    tabela = [tabela, sprintf('\\end{table}\n')];

end

%% Função para salvar tabelas em arquivos .tex
function salvarTabelasTeX(tabelaComplexidade, tabelaPerformance)

    % Salvar tabela de complexidade
    fileID = fopen('tabela_complexidade_ex6.tex', 'w');
    if fileID == -1
        error('Não foi possível criar o arquivo tabela_complexidade_ex6.tex');
    end
    fprintf(fileID, '%s', tabelaComplexidade);
    fclose(fileID);

    % Salvar tabela de performance
    fileID = fopen('tabela_performance_ex6.tex', 'w');
    if fileID == -1
        error('Não foi possível criar o arquivo tabela_performance_ex6.tex');
    end
    fprintf(fileID, '%s', tabelaPerformance);
    fclose(fileID);

    fprintf('Arquivos LaTeX salvos:\n');
    fprintf('  - tabela_complexidade_ex6.tex\n');
    fprintf('  - tabela_performance_ex6.tex\n');

end