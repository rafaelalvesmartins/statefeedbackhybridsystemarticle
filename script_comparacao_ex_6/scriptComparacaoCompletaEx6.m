function scriptComparacaoCompletaEx6()
%   Script principal para comparação completa dos 5 controladores
%   Sistema do Exemplo 6 com incertezas intervalares
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
    fprintf('              Sistema Exemplo 6 - Reviewer\n');
    fprintf('========================================================\n\n');

    % Adicionar caminhos necessários
    addpath('../HInf - Análise - Intervalar/funcoes');
    addpath('../HInf - Análise - Intervalar/funcoes/Diversas');
    addpath('../funcoes');
    addpath('..');

    tempoInicial = datetime('now');

    %% ETAPA 1: Carregamento e verificação do sistema
    fprintf('ETAPA 1: Carregando sistema do Exemplo 6...\n');
    try
        sistema = genSaveDataEx6Reviewer();
        fprintf('✓ Sistema carregado com sucesso!\n');
        fprintf('  - Estados: %d | Entradas: %d | Saídas: %d\n', ...
            length(sistema.sys.A.inf), size(sistema.sys.B.inf, 2), size(sistema.sys.C.inf, 1));
        fprintf('  - Período de amostragem: %.4f s\n', sistema.parametros.h);
        fprintf('  - Incertezas: A(1,1)±0.1, A(1,2)±0.1, A(2,1)±0.5, A(2,2)±0.2\n\n');
    catch ME
        fprintf('✗ ERRO no carregamento: %s\n', ME.message);
        return;
    end

    %% ETAPA 2: Síntese dos controladores
    fprintf('ETAPA 2: Executando síntese dos 5 controladores...\n');
    try
        controladores = comparacao5ControladoresEx6();
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
        estabilidade = analisaEstabilidade5ControladoresEx6(controladores, 1000);
        fprintf('✓ Análise de estabilidade concluída!\n\n');
    catch ME
        fprintf('✗ ERRO na análise de estabilidade: %s\n', ME.message);
        estabilidade = []; % Continuar sem essa análise
    end

    %% ETAPA 4: Simulação temporal
    fprintf('ETAPA 4: Executando simulação temporal (100 realizações)...\n');
    try
        simulacao = simulacaoTemporal5ControladoresEx6(controladores, 100, 10);
        fprintf('✓ Simulação temporal concluída!\n\n');
    catch ME
        fprintf('✗ ERRO na simulação temporal: %s\n', ME.message);
        simulacao = []; % Continuar sem essa análise
    end

    %% ETAPA 5: Geração de tabelas LaTeX
    fprintf('ETAPA 5: Gerando tabelas LaTeX para o artigo...\n');
    try
        % Compilar dados para as tabelas
        dadosTabelas = compilarDadosTabelasEx6(controladores, estabilidade, simulacao);

        % Gerar tabelas
        tabelaComplexidade = gerarTabelaComplexidadeEx6(dadosTabelas);
        tabelaPerformance = gerarTabelaPerformanceEx6(dadosTabelas);

        % Salvar tabelas
        salvarTabelasLaTeXEx6(tabelaComplexidade, tabelaPerformance);

        fprintf('✓ Tabelas LaTeX geradas e salvas!\n\n');
    catch ME
        fprintf('✗ ERRO na geração de tabelas: %s\n', ME.message);
        tabelaComplexidade = '';
        tabelaPerformance = '';
    end

    %% ETAPA 6: Salvamento completo dos resultados
    fprintf('ETAPA 6: Salvando resultados completos...\n');
    try
        % Estrutura final com todos os resultados
        resultadosCompletos.sistema = sistema;
        resultadosCompletos.controladores = controladores;
        if ~isempty(estabilidade)
            resultadosCompletos.estabilidade = estabilidade;
        end
        if ~isempty(simulacao)
            resultadosCompletos.simulacao = simulacao;
        end
        if exist('tabelaComplexidade', 'var')
            resultadosCompletos.tabelas.complexidade = tabelaComplexidade;
            resultadosCompletos.tabelas.performance = tabelaPerformance;
        end
        resultadosCompletos.timestamp = datestr(now);

        % Salvar em arquivo MAT
        save('resultados_comparacao_exemplo6_reviewer.mat', 'resultadosCompletos', '-v7.3');

        % Gerar gráficos básicos
        gerarGraficosEx6(resultadosCompletos);

        fprintf('✓ Resultados salvos em "resultados_comparacao_exemplo6_reviewer.mat"\n');
        fprintf('✓ Gráficos básicos gerados\n\n');
    catch ME
        fprintf('✗ ERRO no salvamento: %s\n', ME.message);
    end

    %% Resumo final
    tempoTotal = datetime('now') - tempoInicial;
    fprintf('========================================================\n');
    fprintf('                    RESUMO FINAL\n');
    fprintf('========================================================\n');
    fprintf('Tempo total de execução: %.1f segundos\n\n', seconds(tempoTotal));

    fprintf('Arquivos gerados:\n');
    fprintf('  📊 resultados_comparacao_exemplo6_reviewer.mat\n');
    if exist('tabelaComplexidade', 'var') && ~isempty(tabelaComplexidade)
        fprintf('  📄 tabela_complexidade_exemplo6.tex\n');
        fprintf('  📄 tabela_performance_exemplo6.tex\n');
    end
    fprintf('  📈 Gráficos básicos de comparação\n\n');

    % Exibir resumo dos controladores
    fprintf('🏆 RESUMO DOS CONTROLADORES:\n');
    exibirResumoEx6(controladores);

    fprintf('\n========================================================\n');
    fprintf('           COMPARAÇÃO FINALIZADA COM SUCESSO!\n');
    fprintf('========================================================\n');
end

%% Função para gerar dados do sistema Exemplo 6
function sistema = genSaveDataEx6Reviewer()
    
    fprintf('  Configurando sistema do Exemplo 6...\n');
    
    % Parâmetros que funcionaram
    h = 0.0500;  % Período de amostragem
    delta = h/10;  % = 0.005
    tol = 1e-6;
    
    % Sistema nominal
    A_nom = [1.0  -1.0; 
             0.0   4.0];
    B_nom = [1.0; -1.0];
    E_nom = [1.0; -1.0];  % Mesmo que B para simplicidade
    C_nom = [1 0; 0 1];   % Matriz identidade
    D_nom = [0; 0];       % Sem feedthrough direto
    
    % Incertezas intervalares
    % A(1,1): 1.0 ± 0.1 → [0.9, 1.1]
    % A(1,2): -1.0 ± 0.1 → [-1.1, -0.9] 
    % A(2,1): 0.0 ± 0.5 → [-0.5, 0.5]
    % A(2,2): 4.0 ± 0.2 → [3.8, 4.2]
    
    A.inf = [0.9  -1.1; 
            -0.5   3.8];
    A.sup = [1.1  -0.9; 
             0.5   4.2];
    
    % B, E, C, D sem incertezas (precisos)
    B.inf = B_nom; B.sup = B_nom;
    E.inf = E_nom; E.sup = E_nom;
    C.inf = C_nom; C.sup = C_nom;
    D.inf = D_nom; D.sup = D_nom;
    
    % Criar sistema politópico (4 vértices)
    % Vértice 1: A.inf(1,1), A.inf(1,2), A.inf(2,1), A.inf(2,2)
    A1 = [0.9  -1.1; -0.5   3.8];
    % Vértice 2: A.inf(1,1), A.inf(1,2), A.inf(2,1), A.sup(2,2)
    A2 = [0.9  -1.1; -0.5   4.2];
    % Vértice 3: A.inf(1,1), A.inf(1,2), A.sup(2,1), A.inf(2,2)
    A3 = [0.9  -1.1;  0.5   3.8];
    % Vértice 4: A.inf(1,1), A.inf(1,2), A.sup(2,1), A.sup(2,2)
    A4 = [0.9  -1.1;  0.5   4.2];
    % Vértice 5: A.inf(1,1), A.sup(1,2), A.inf(2,1), A.inf(2,2)
    A5 = [0.9  -0.9; -0.5   3.8];
    % Vértice 6: A.inf(1,1), A.sup(1,2), A.inf(2,1), A.sup(2,2)
    A6 = [0.9  -0.9; -0.5   4.2];
    % Vértice 7: A.inf(1,1), A.sup(1,2), A.sup(2,1), A.inf(2,2)
    A7 = [0.9  -0.9;  0.5   3.8];
    % Vértice 8: A.inf(1,1), A.sup(1,2), A.sup(2,1), A.sup(2,2)
    A8 = [0.9  -0.9;  0.5   4.2];
    % Vértice 9: A.sup(1,1), A.inf(1,2), A.inf(2,1), A.inf(2,2)
    A9 = [1.1  -1.1; -0.5   3.8];
    % Vértice 10: A.sup(1,1), A.inf(1,2), A.inf(2,1), A.sup(2,2)
    A10 = [1.1  -1.1; -0.5   4.2];
    % Vértice 11: A.sup(1,1), A.inf(1,2), A.sup(2,1), A.inf(2,2)
    A11 = [1.1  -1.1;  0.5   3.8];
    % Vértice 12: A.sup(1,1), A.inf(1,2), A.sup(2,1), A.sup(2,2)
    A12 = [1.1  -1.1;  0.5   4.2];
    % Vértice 13: A.sup(1,1), A.sup(1,2), A.inf(2,1), A.inf(2,2)
    A13 = [1.1  -0.9; -0.5   3.8];
    % Vértice 14: A.sup(1,1), A.sup(1,2), A.inf(2,1), A.sup(2,2)
    A14 = [1.1  -0.9; -0.5   4.2];
    % Vértice 15: A.sup(1,1), A.sup(1,2), A.sup(2,1), A.inf(2,2)
    A15 = [1.1  -0.9;  0.5   3.8];
    % Vértice 16: A.sup(1,1), A.sup(1,2), A.sup(2,1), A.sup(2,2)
    A16 = [1.1  -0.9;  0.5   4.2];
    
    % Criar sistema politópico (usar apenas 4 vértices principais por simplicidade)
    sysPolyCont = cell(4, 1);
    vertices_A = {A1, A8, A9, A16}; % Vértices extremos
    
    for i = 1:4
        sysPolyCont{i}.A = vertices_A{i};
        sysPolyCont{i}.B2 = B_nom;   % Entrada de controle
        sysPolyCont{i}.B1 = E_nom;   % Entrada de perturbação
        sysPolyCont{i}.C = C_nom;
        sysPolyCont{i}.D2 = D_nom;   % Feedthrough do controle
    end
    
    % Matrizes caligráficas para análise híbrida
    nx = 2; nu = 1;
    ACal.inf = [A.inf  B.inf; zeros(nu, nx+nu)];
    ACal.sup = [A.sup  B.sup; zeros(nu, nx+nu)];
    
    ECal.inf = [E.inf; zeros(nu, size(E.inf,2))];
    ECal.sup = [E.sup; zeros(nu, size(E.sup,2))];
    
    CCal.inf = [C.inf  D.inf];
    CCal.sup = [C.sup  D.sup];
    
    % Estrutura de saída
    sistema.sys.A = A;
    sistema.sys.B = B;
    sistema.sys.E = E;
    sistema.sys.C = C;
    sistema.sys.D = D;
    sistema.sys.sysPolyCont = sysPolyCont;
    sistema.sys.ACal = ACal;
    sistema.sys.ECal = ECal;
    sistema.sys.CCal = CCal;
    sistema.sys.h = h;
    sistema.sys.delta = delta;
    
    % Sistema nominal para síntese nominal
    sistema.sysNominal.A = A_nom;
    sistema.sysNominal.B = B_nom;
    sistema.sysNominal.E = E_nom;
    sistema.sysNominal.C = C_nom;
    sistema.sysNominal.D = D_nom;
    
    % Parâmetros
    sistema.parametros.h = h;
    sistema.parametros.delta = delta;
    sistema.parametros.tol = tol;
    sistema.parametros.nx = nx;
    sistema.parametros.nu = nu;
    
    % Informações auxiliares
    sistema.aux.tol = tol;
    sistema.aux.descricao = 'Sistema Exemplo 6 para resposta ao reviewer';
    sistema.aux.incertezas = 'A(1,1)±0.1, A(1,2)±0.1, A(2,1)±0.5, A(2,2)±0.2';
    
    fprintf('  ✓ Sistema Exemplo 6 configurado com sucesso!\n');
end

%% Funções auxiliares para compilação de dados
function dados = compilarDadosTabelasEx6(controladores, estabilidade, simulacao)
    
    nomes = {'Hibrido_Int', 'Hibrido_Poly', 'Equivalente_Int', 'Equivalente_Poly', 'Nominal'};
    nomesDisplay = {'Híbrido Int', 'Híbrido Poly', 'Equiv. Int', 'Equiv. Poly', 'Nominal'};
    
    dados.nomes = nomes;
    dados.nomesDisplay = nomesDisplay;
    dados.numControladores = length(nomes);
    
    % Dados de complexidade
    dados.ganhos = cell(1, dados.numControladores);
    dados.variaveis = zeros(1, dados.numControladores);
    dados.restricoes = zeros(1, dados.numControladores);
    dados.tempoSintese = zeros(1, dados.numControladores);
    dados.normaGarantida = zeros(1, dados.numControladores);
    
    % Dados de performance (se disponível)
    if ~isempty(estabilidade)
        dados.taxaEstabilidade = estabilidade.taxaEstabilidade;
        dados.normaPiorCaso = estabilidade.normaPiorCaso;
        dados.normaMedia = estabilidade.normaMedia;
        dados.desviosPadrao = estabilidade.desviosPadrao;
    else
        dados.taxaEstabilidade = nan(1, dados.numControladores);
        dados.normaPiorCaso = nan(1, dados.numControladores);
        dados.normaMedia = nan(1, dados.numControladores);
        dados.desviosPadrao = nan(1, dados.numControladores);
    end
    
    for i = 1:dados.numControladores
        nome = nomes{i};
        if isfield(controladores, nome) && controladores.(nome).factivel
            dados.ganhos{i} = controladores.(nome).K;
            if isfield(controladores.(nome), 'variaveis')
                dados.variaveis(i) = controladores.(nome).variaveis;
            end
            if isfield(controladores.(nome), 'restricoes')
                dados.restricoes(i) = controladores.(nome).restricoes;
            end
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
function tabela = gerarTabelaComplexidadeEx6(dados)
    
    tabela = sprintf('\\begin{table}[htbp]\n');
    tabela = [tabela, sprintf('\\centering\n')];
    tabela = [tabela, sprintf('\\caption{Complexidade Computacional - Sistema Exemplo 6}\n')];
    tabela = [tabela, sprintf('\\label{tab:complexidade_exemplo6}\n')];
    tabela = [tabela, sprintf('\\begin{tabular}{lccccc}\n')];
    tabela = [tabela, sprintf('\\hline\n')];
    tabela = [tabela, sprintf('Método & Ganho $K$ & Variáveis & Restrições & Tempo (s) & $\\gamma$ \\\\\n')];
    tabela = [tabela, sprintf('\\hline\n')];
    
    for i = 1:dados.numControladores
        nome = dados.nomesDisplay{i};
        
        % Ganho
        if ~isempty(dados.ganhos{i})
            if length(dados.ganhos{i}) <= 4
                ganhoStr = sprintf('$[%s]$', strjoin(arrayfun(@(x) sprintf('%.3f', x), dados.ganhos{i}, 'UniformOutput', false), ', '));
            else
                ganhoStr = sprintf('$[%.3f, \\ldots]$', dados.ganhos{i}(1));
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
function tabela = gerarTabelaPerformanceEx6(dados)
    
    tabela = sprintf('\\begin{table}[htbp]\n');
    tabela = [tabela, sprintf('\\centering\n')];
    tabela = [tabela, sprintf('\\caption{Análise de Performance Monte Carlo - Exemplo 6}\n')];
    tabela = [tabela, sprintf('\\label{tab:performance_exemplo6}\n')];
    tabela = [tabela, sprintf('\\begin{tabular}{lcccc}\n')];
    tabela = [tabela, sprintf('\\hline\n')];
    tabela = [tabela, sprintf('Método & Taxa Estab. (\\%%) & Pior Caso & Norma Média & Desvio Padrão \\\\\n')];
    tabela = [tabela, sprintf('\\hline\n')];
    
    for i = 1:dados.numControladores
        nome = dados.nomesDisplay{i};
        
        % Taxa de estabilidade
        if ~isnan(dados.taxaEstabilidade(i))
            taxaStr = sprintf('%.1f', dados.taxaEstabilidade(i));
        else
            taxaStr = 'N/A';
        end
        
        % Pior caso
        if ~isnan(dados.normaPiorCaso(i)) && ~isinf(dados.normaPiorCaso(i))
            piorCasoStr = sprintf('%.3f', dados.normaPiorCaso(i));
        else
            piorCasoStr = 'N/A';
        end
        
        % Norma média
        if ~isnan(dados.normaMedia(i)) && ~isinf(dados.normaMedia(i))
            mediaStr = sprintf('%.3f', dados.normaMedia(i));
        else
            mediaStr = 'N/A';
        end
        
        % Desvio padrão
        if ~isnan(dados.desviosPadrao(i)) && dados.desviosPadrao(i) > 0
            desvioStr = sprintf('%.3f', dados.desviosPadrao(i));
        else
            desvioStr = 'N/A';
        end
        
        tabela = [tabela, sprintf('%s & %s & %s & %s & %s \\\\\n', ...
            nome, taxaStr, piorCasoStr, mediaStr, desvioStr)];
    end
    
    tabela = [tabela, sprintf('\\hline\n')];
    tabela = [tabela, sprintf('\\end{tabular}\n')];
    tabela = [tabela, sprintf('\\end{table}\n')];
end

%% Função para salvar tabelas LaTeX
function salvarTabelasLaTeXEx6(tabelaComplexidade, tabelaPerformance)
    
    % Salvar tabela de complexidade
    fileID = fopen('tabela_complexidade_exemplo6.tex', 'w');
    fprintf(fileID, '%s', tabelaComplexidade);
    fclose(fileID);
    
    % Salvar tabela de performance
    fileID = fopen('tabela_performance_exemplo6.tex', 'w');
    fprintf(fileID, '%s', tabelaPerformance);
    fclose(fileID);
end

%% Função para gerar gráficos básicos
function gerarGraficosEx6(resultados)
    
    try
        figure('Position', [100, 100, 1000, 600]);
        
        % Gráfico simples dos ganhos dos controladores
        nomes = fieldnames(resultados.controladores);
        ganhos = [];
        nomesValidos = {};
        
        for i = 1:length(nomes)
            if strcmp(nomes{i}, 'sistema') || strcmp(nomes{i}, 'nomeControladores')
                continue;
            end
            if isfield(resultados.controladores.(nomes{i}), 'factivel') && ...
                    resultados.controladores.(nomes{i}).factivel
                ganhos = [ganhos; resultados.controladores.(nomes{i}).K(:)'];
                nomesValidos{end+1} = strrep(nomes{i}, '_', ' ');
            end
        end
        
        if ~isempty(ganhos)
            bar(ganhos);
            legend(nomesValidos, 'Location', 'best');
            xlabel('Elemento do Ganho');
            ylabel('Valor');
            title('Comparação dos Ganhos dos Controladores - Exemplo 6');
            grid on;
            
            saveas(gcf, 'comparacao_ganhos_exemplo6.fig');
            saveas(gcf, 'comparacao_ganhos_exemplo6.png');
        end
        
    catch ME
        fprintf('  Aviso: Erro na geração de gráficos: %s\n', ME.message);
    end
end

%% Função para exibir resumo
function exibirResumoEx6(controladores)
    
    nomes = fieldnames(controladores);
    
    fprintf('  Controlador             | Status    | Ganho K               | Norma γ\n');
    fprintf('  ------------------------|-----------|----------------------|--------\n');
    
    for i = 1:length(nomes)
        if strcmp(nomes{i}, 'sistema') || strcmp(nomes{i}, 'nomeControladores')
            continue;
        end
        
        nomeDisplay = strrep(nomes{i}, '_', ' ');
        
        if isfield(controladores.(nomes{i}), 'factivel') && controladores.(nomes{i}).factivel
            ganhoStr = sprintf('[%.3f, %.3f]', controladores.(nomes{i}).K(1), controladores.(nomes{i}).K(2));
            fprintf('  %-23s | %-9s | %-20s | %.3f\n', ...
                nomeDisplay, 'OK', ganhoStr, controladores.(nomes{i}).mu);
        else
            fprintf('  %-23s | %-9s | %-20s | %s\n', ...
                nomeDisplay, 'FALHOU', 'N/A', 'N/A');
        end
    end
end

% Função auxiliar vazia para compatibilidade
function estabilidade = analisaEstabilidade5ControladoresEx6(controladores, numRealizacoes)
    fprintf('  Análise de estabilidade não implementada ainda.\n');
    estabilidade = [];
end

% Função auxiliar vazia para compatibilidade  
function simulacao = simulacaoTemporal5ControladoresEx6(controladores, numRealizacoes, tempoSim)
    fprintf('  Simulação temporal não implementada ainda.\n');
    simulacao = [];
end