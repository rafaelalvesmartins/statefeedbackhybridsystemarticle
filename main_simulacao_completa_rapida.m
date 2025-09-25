function main_simulacao_completa_rapida()
%   MAIN_SIMULACAO_COMPLETA_RAPIDA - Versão mais rápida com menos simulações
%
%   Esta é uma versão otimizada para execução mais rápida, usando menos
%   simulações Monte Carlo (30 em vez de 100) e com melhor tratamento de timeouts.

    clear; clc; close all;

    fprintf('==================================================================================\n');
    fprintf('                      SIMULAÇÃO MONTE CARLO RÁPIDA\n');
    fprintf('     Versão otimizada com 30 simulações por controlador\n');
    fprintf('==================================================================================\n\n');

    % Carregar resultados das sínteses
    fprintf('Carregando controladores sintetizados...\n');
    try
        load('resultados_completos_massa_mola.mat', 'resultados');
        fprintf('✓ Controladores carregados com sucesso!\n\n');
    catch ME
        error('ERRO: Execute primeiro main_massa_mola_completo.m para gerar os controladores');
    end

    % Extrair controladores
    K_nominal_sampled = resultados.ganhos.nominal.classica.K;
    K_robust_interval = resultados.ganhos.robusto.intervalar.K;
    K_robust_polytopic = resultados.ganhos.robusto.politopica.K;
    K_continuous_sampled = resultados.ganhos.continuo.intervalar.K;

    % Parâmetros do sistema
    m1 = resultados.parametros.m1;
    m2 = resultados.parametros.m2;
    k1 = resultados.parametros.k1;
    k2 = resultados.parametros.k2;
    b = resultados.parametros.b;
    h = resultados.parametros.h;
    delta = resultados.parametros.delta;
    tol = 1e-6;

    % Configuração rápida
    numSimulacoesMC = 30;  % Reduzido para execução mais rápida
    fprintf('Configuração Monte Carlo RÁPIDA: %d simulações por controlador\n\n');

    % Preparar parâmetros
    parametros_sistema.m1 = m1;
    parametros_sistema.m2 = m2;
    parametros_sistema.k1 = k1;
    parametros_sistema.k2 = k2;
    parametros_sistema.b = b;

    % Gerar sistemas Monte Carlo fixos
    fprintf('Gerando %d sistemas Monte Carlo...\n', numSimulacoesMC);
    seed_monte_carlo = 54321;
    sistemas_monte_carlo_fixos = gerarSistemasMonteCarloFixos(parametros_sistema, numSimulacoesMC, seed_monte_carlo);

    % Controladores para análise
    nomes_controladores = {
        'H∞ Nominal Sampled';
        'H∞ Robust Interval';
        'H∞ Robust Polytopic';
        'H∞ Continuous Sampled'
    };

    controladores_K = {K_nominal_sampled; K_robust_interval; K_robust_polytopic; K_continuous_sampled};
    resultados_simulacao = cell(4, 1);
    axisVector = [0 15 -0.5 1.5];

    fprintf('\n==================================================================================\n');
    fprintf('                    EXECUTANDO SIMULAÇÕES RÁPIDAS\n');
    fprintf('==================================================================================\n\n');

    % Loop pelos controladores
    for ctrl = 1:4
        fprintf('--- SIMULAÇÃO %d: %s ---\n', ctrl, nomes_controladores{ctrl});

        try
            imageName = sprintf('SimRapida_%d_%s', ctrl, strrep(nomes_controladores{ctrl}, ' ', '_'));
            flagIsPoly = (ctrl == 3);  % Apenas o politópico usa análise politópica

            tic;
            fprintf('    Testando controlador em %d sistemas...\n', numSimulacoesMC);

            output_sim = simulatesSampledInputMonteCarlo(sistemas_monte_carlo_fixos, h, controladores_K{ctrl}, ...
                imageName, axisVector, delta, tol, flagIsPoly);

            tempo_sim = toc;
            resultados_simulacao{ctrl} = output_sim;

            fprintf('  ✓ Simulação concluída: %.2f s\n', tempo_sim);

            if isfield(output_sim, 'maxCost') && isfinite(output_sim.maxCost)
                fprintf('    Custo máximo: %.6f\n', output_sim.maxCost);
                fprintf('    Custo médio: %.6f ± %.6f\n', output_sim.meanCost, output_sim.stdCost);
            end

            if isfield(output_sim, 'gamma') && isfield(output_sim.gamma, 'meanGamma')
                fprintf('    γ médio: %.6f, máx: %.6f\n', output_sim.gamma.meanGamma, output_sim.gamma.maxGamma);
            end

            fprintf('    Taxa sucesso: %.1f%% simulações\n', output_sim.sim_success_rate * 100);

        catch ME
            fprintf('  ✗ ERRO: %s\n', ME.message);
            resultados_simulacao{ctrl} = [];
        end
        fprintf('\n');
    end

    % Resultados resumidos
    fprintf('==================================================================================\n');
    fprintf('                         RESULTADOS RESUMIDOS\n');
    fprintf('==================================================================================\n\n');

    fprintf('%-25s | %-12s | %-12s | %-12s\n', 'CONTROLADOR', 'Custo Máx', 'Custo Médio', 'Taxa Sucesso');
    fprintf('%s\n', repmat('-', 70, 1));

    for ctrl = 1:4
        if ~isempty(resultados_simulacao{ctrl})
            sim_result = resultados_simulacao{ctrl};

            if isfield(sim_result, 'maxCost')
                fprintf('%-25s | %10.6f | %10.6f | %9.1f%%\n', ...
                    nomes_controladores{ctrl}, sim_result.maxCost, sim_result.meanCost, ...
                    sim_result.sim_success_rate * 100);
            else
                fprintf('%-25s | %10s | %10s | %9.1f%%\n', ...
                    nomes_controladores{ctrl}, 'N/A', 'N/A', sim_result.sim_success_rate * 100);
            end
        else
            fprintf('%-25s | %10s | %10s | %9s\n', nomes_controladores{ctrl}, 'FALHOU', 'FALHOU', '0.0%');
        end
    end

    % Análise comparativa rápida
    fprintf('\n📊 ANÁLISE COMPARATIVA RÁPIDA:\n');
    custos_validos = [];
    nomes_validos = {};

    for ctrl = 1:4
        if ~isempty(resultados_simulacao{ctrl}) && isfield(resultados_simulacao{ctrl}, 'maxCost') && ...
           isfinite(resultados_simulacao{ctrl}.maxCost)
            custos_validos(end+1) = resultados_simulacao{ctrl}.maxCost;
            nomes_validos{end+1} = nomes_controladores{ctrl};
        end
    end

    if ~isempty(custos_validos)
        [melhor_custo, melhor_idx] = min(custos_validos);
        fprintf('   🏆 MELHOR MÉTODO: %s (custo=%.6f)\n', nomes_validos{melhor_idx}, melhor_custo);

        % Mostrar ranking
        [custos_ordenados, indices_ordenados] = sort(custos_validos);
        fprintf('   📈 RANKING (melhor → pior):\n');
        for i = 1:length(custos_ordenados)
            fprintf('      %d. %s: %.6f\n', i, nomes_validos{indices_ordenados(i)}, custos_ordenados(i));
        end
    else
        fprintf('   ⚠ Nenhum resultado válido obtido\n');
    end

    % Salvar resultados rápidos
    saida_rapida.resultados_simulacao = resultados_simulacao;
    saida_rapida.nomes_controladores = nomes_controladores;
    saida_rapida.controladores_K = controladores_K;
    saida_rapida.sistemas_monte_carlo = sistemas_monte_carlo_fixos;
    saida_rapida.parametros.numSimulacoesMC = numSimulacoesMC;
    saida_rapida.parametros.seed_monte_carlo = seed_monte_carlo;
    saida_rapida.metodo = 'Monte Carlo Rápido';
    saida_rapida.timestamp = datestr(now);

    try
        save('resultados_monte_carlo_rapido.mat', 'saida_rapida');
        fprintf('\n✓ Resultados salvos em: resultados_monte_carlo_rapido.mat\n');
    catch
        fprintf('\n⚠ Erro ao salvar resultados\n');
    end

    fprintf('\n==================================================================================\n');
    fprintf('                    SIMULAÇÃO RÁPIDA FINALIZADA\n');
    fprintf('                            %s\n', datestr(now));
    fprintf('==================================================================================\n');

end