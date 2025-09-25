%% TESTE RÁPIDO DA SIMULAÇÃO ADAPTADA DO SIMULINK
clc; clear; close all;

fprintf('=== TESTE DA SIMULAÇÃO ADAPTADA ===\n\n');

% Setup paths
if exist('setup_paths.m', 'file')
    setup_paths;
else
    fprintf('⚠ setup_paths.m não encontrado, continuando...\n');
end

try
    % Carregar dados de teste
    if exist('resultados_completos_massa_mola.mat', 'file')
        fprintf('Carregando dados reais...\n');
        load('resultados_completos_massa_mola.mat');
        fprintf('✓ Dados carregados!\n\n');
    else
        error('Arquivo de dados não encontrado!');
    end

    % Verificar se temos controladores
    if exist('resultados', 'var') && isfield(resultados, 'ganhos')
        fprintf('Testando controlador dos resultados salvos...\n');

        % Extrair controlador nominal clássico
        ganhos = resultados.ganhos;
        if isfield(ganhos, 'nominal') && isfield(ganhos.nominal, 'classica') && isfield(ganhos.nominal.classica, 'K')
            K = ganhos.nominal.classica.K;
            gamma_teorico = ganhos.nominal.classica.gamma;
            fprintf('  Controlador nominal clássico encontrado (γ=%.6f)\n', gamma_teorico);
        elseif isfield(ganhos, 'nominal') && isfield(ganhos.nominal, 'intervalar') && isfield(ganhos.nominal.intervalar, 'K')
            K = ganhos.nominal.intervalar.K;
            gamma_teorico = ganhos.nominal.intervalar.gamma;
            fprintf('  Controlador nominal intervalar encontrado (γ=%.6f)\n', gamma_teorico);
        else
            error('Controlador nominal K não encontrado');
        end

        % Parâmetros do sistema
        if isfield(resultados, 'parametros')
            params = resultados.parametros;
            h = params.h;
            delta = params.delta;
            tol = 1e-6;
        else
            h = 0.02;
            delta = 0.004;
            tol = 1e-6;
        end

        fprintf('  K = [%.6f %.6f %.6f %.6f]\n', K);
        fprintf('  h = %.4f, delta = %.6f\n', h, delta);

        % Preparar politopo de teste simples (2 vértices para teste rápido)
        fprintf('\nPreparando sistema politópico de teste...\n');

        % Sistema massa-mola: 2 massas conectadas por molas
        m1 = 0.5; m2 = 1.0; k1 = 12; k2 = 7; b = 0.20;

        % Vértices do politopo (variação em k2 e b)
        vertices = [
            6.86, 0.19;  % k2_min, b_min
            7.14, 0.21;  % k2_max, b_max
        ];

        combPoly.A.alphaVecs = {[1 0], [0 1]}; % Pesos dos vértices
        combPoly.A.polytopicMatrices = {};
        combPoly.B.polytopicMatrices = {};
        combPoly.E.polytopicMatrices = {};
        combPoly.C.polytopicMatrices = {};
        combPoly.D.polytopicMatrices = {};

        for v = 1:size(vertices, 1)
            k2_v = vertices(v, 1);
            b_v = vertices(v, 2);

            % Matrizes do sistema massa-mola
            A_v = [0 1 0 0;
                   -k1/m1 -b_v/m1 k2_v/m1 0;
                   0 0 0 1;
                   k2_v/m2 0 -k2_v/m2 0];

            B_v = [0; 1/m1; 0; 0];  % Entrada de controle
            E_v = [0; 0.1/m1; 0; 0.1/m2];  % Entrada de perturbação
            C_v = [1 0 0 0; 0 0 1 0];  % Saídas: posições das massas
            D_v = [0; 0];  % Sem feedthrough

            combPoly.A.polytopicMatrices{v} = A_v;
            combPoly.B.polytopicMatrices{v} = B_v;
            combPoly.E.polytopicMatrices{v} = E_v;
            combPoly.C.polytopicMatrices{v} = C_v;
            combPoly.D.polytopicMatrices{v} = D_v;

            fprintf('  Vértice %d: k2=%.3f, b=%.4f\n', v, k2_v, b_v);
        end

        fprintf('✓ Politopo preparado: %d vértices\n\n', length(combPoly.A.alphaVecs));

        % Testar simulação corrigida
        fprintf('--- TESTE COM SIMULAÇÃO CORRIGIDA ---\n');
        tic;
        try
            output = simulatesSampledInputFixed(combPoly, h, K, 'teste', [-2 2 -1 1], delta, tol, false);
            tempo_exec = toc;

            fprintf('✓ Simulação concluída: %.2f s\n', tempo_exec);
            fprintf('  Custo máximo: %.6f\n', output.maxCost);
            fprintf('  Custo médio: %.6f ± %.6f\n', output.meanCost, output.stdCost);

            if isfield(output, 'gamma')
                fprintf('  γ médio (LMI): %.6f\n', output.gamma.meanGamma);
                fprintf('  γ máximo (LMI): %.6f\n', output.gamma.maxGamma);
            end

            fprintf('  Taxa sucesso LMI: %.1f%%\n', output.lmi_success_rate * 100);
            fprintf('  Taxa sucesso simulação: %.1f%%\n', output.sim_success_rate * 100);
            fprintf('  Método: %s\n', output.metodo);

            if output.sim_success_rate > 0
                fprintf('✓ SIMULAÇÃO BEM-SUCEDIDA!\n');
            else
                fprintf('⚠ Simulação com problemas de estabilidade.\n');
            end

        catch ME
            tempo_exec = toc;
            fprintf('✗ ERRO no teste: %s\n', ME.message);
            fprintf('   Tempo até erro: %.2f s\n', tempo_exec);
        end

    else
        fprintf('⚠ Controladores não encontrados nos dados.\n');
    end

catch ME
    fprintf('✗ ERRO GERAL: %s\n', ME.message);
end

fprintf('\n=== TESTE FINALIZADO ===\n');