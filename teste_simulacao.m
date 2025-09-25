function teste_simulacao()
%   TESTE_SIMULACAO - Função para testar as funções de simulação
%
%   Esta função cria dados de teste e verifica se as funções estão funcionando

    fprintf('=== TESTE DAS FUNÇÕES DE SIMULAÇÃO ===\n\n');

    % Tentar carregar dados reais primeiro
    try
        fprintf('Tentando carregar dados reais...\n');
        load('resultados_completos_massa_mola.mat', 'resultados');
        fprintf('✓ Dados reais carregados!\n');

        % Usar os dados reais se disponíveis
        use_real_data = true;

    catch
        fprintf('⚠ Dados reais não encontrados. Criando dados de teste...\n');
        use_real_data = false;
    end

    if use_real_data
        % TESTE COM DADOS REAIS
        fprintf('\n--- EXECUTANDO main_simulacao_completa ---\n');
        try
            main_simulacao_completa();  % Função não retorna saída
            fprintf('✓ Script principal executado com sucesso!\n');

            % Tentar carregar resultados salvos
            fprintf('\n=== VERIFICANDO RESULTADOS SALVOS ===\n');
            try
                if exist('resultados_resposta_revisor.mat', 'file')
                    load('resultados_resposta_revisor.mat', 'saida_revisor');
                    fprintf('✓ Resultados encontrados para %d controladores\n', length(saida_revisor.nomes_controladores));
                else
                    fprintf('⚠ Arquivo de resultados não encontrado\n');
                end
            catch
                fprintf('⚠ Erro ao carregar resultados\n');
            end

        catch ME
            fprintf('✗ ERRO no script principal: %s\n', ME.message);
            if ~isempty(ME.stack)
                fprintf('   Linha: %d, Arquivo: %s\n', ME.stack(1).line, ME.stack(1).name);
            end
        end

    else
        % TESTE COM DADOS SINTÉTICOS
        fprintf('\n--- CRIANDO DADOS DE TESTE ---\n');

        % Criar um sistema simples 2x2 para teste
        nx = 2; nu = 1; nw = 1;

        % Criar politopo de teste com 2 vértices
        combPoly = struct();

        % Vértice 1
        combPoly.A.alphaVecs{1} = 1;
        combPoly.A.polytopicMatrices{1} = [0 1; -2 -3];
        combPoly.B.alphaVecs{1} = 1;
        combPoly.B.polytopicMatrices{1} = [0; 1];
        combPoly.E.alphaVecs{1} = 1;
        combPoly.E.polytopicMatrices{1} = [1; 0];
        combPoly.C.alphaVecs{1} = 1;
        combPoly.C.polytopicMatrices{1} = [1 0];
        combPoly.D.alphaVecs{1} = 1;
        combPoly.D.polytopicMatrices{1} = [0];
        combPoly.D1.alphaVecs{1} = 1;
        combPoly.D1.polytopicMatrices{1} = [0];

        % Vértice 2 (ligeiramente diferente)
        combPoly.A.alphaVecs{2} = 2;
        combPoly.A.polytopicMatrices{2} = [0 1; -2.2 -3.1];
        combPoly.B.alphaVecs{2} = 2;
        combPoly.B.polytopicMatrices{2} = [0; 1];
        combPoly.E.alphaVecs{2} = 2;
        combPoly.E.polytopicMatrices{2} = [1; 0];
        combPoly.C.alphaVecs{2} = 2;
        combPoly.C.polytopicMatrices{2} = [1 0];
        combPoly.D.alphaVecs{2} = 2;
        combPoly.D.polytopicMatrices{2} = [0];
        combPoly.D1.alphaVecs{2} = 2;
        combPoly.D1.polytopicMatrices{2} = [0];

        % Parâmetros de teste
        h = 0.1;
        K = [2 1];  % Ganho estabilizante simples
        imageName = 'teste';
        axisVector = [0 10 -1 1];
        delta = 0.01;
        tol = 1e-6;

        fprintf('Sistema de teste: nx=%d, nu=%d, nw=%d, %d vértices\n', nx, nu, nw, length(combPoly.A.alphaVecs));

        % TESTE 1: Versão simplificada
        fprintf('\n--- TESTANDO VERSÃO SIMPLIFICADA ---\n');
        try
            tic;
            output1 = simulatesSampledInputSimplified(combPoly, h, K, imageName, axisVector, delta, tol, false);
            tempo1 = toc;

            fprintf('✓ Versão simplificada OK (%.2f s)\n', tempo1);
            if isfield(output1, 'gamma') && isfield(output1.gamma, 'maxGamma')
                fprintf('  γ máximo: %.6f\n', output1.gamma.maxGamma);
                fprintf('  Taxa sucesso: %.1f%%\n', output1.lmi_success_rate*100);
            end

        catch ME
            fprintf('✗ ERRO na versão simplificada: %s\n', ME.message);
        end

        % TESTE 2: Versão completa
        fprintf('\n--- TESTANDO VERSÃO COMPLETA ---\n');
        try
            tic;
            output2 = simulatesSampledInputCompleto(combPoly, h, K, imageName, axisVector, delta, tol, false);
            tempo2 = toc;

            fprintf('✓ Versão completa OK (%.2f s)\n', tempo2);
            if isfield(output2, 'maxCost') && isfinite(output2.maxCost)
                fprintf('  Custo máximo: %.6f\n', output2.maxCost);
                fprintf('  Custo médio: %.6f\n', output2.meanCost);
            end
            if isfield(output2, 'gamma') && isfield(output2.gamma, 'maxGamma')
                fprintf('  γ máximo: %.6f\n', output2.gamma.maxGamma);
                fprintf('  Taxa sucesso LMI: %.1f%%\n', output2.lmi_success_rate*100);
                fprintf('  Taxa sucesso sim: %.1f%%\n', output2.sim_success_rate*100);
            end

        catch ME
            fprintf('✗ ERRO na versão completa: %s\n', ME.message);
        end
    end

    fprintf('\n=== TESTE FINALIZADO ===\n');
end