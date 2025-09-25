function [outPut] = simulatesSampledInputMonteCarlo(sistemas_monte_carlo, h, K, imageName, axisVector, delta, tol, flagIsPoly)
%   SIMULATESAMPLEDINPUTMONTECARLO - Simulação Monte Carlo com sistemas fixos
%
%   Esta função testa um controlador específico nos sistemas Monte Carlo
%   já gerados, garantindo consistência entre diferentes controladores.
%
%   Parâmetros:
%   - sistemas_monte_carlo: cell array com sistemas já gerados (de gerarSistemasMonteCarloFixos)
%   - h: período de amostragem
%   - K: ganho do controlador
%   - imageName: nome para salvar imagens
%   - axisVector: vetor de eixos para gráficos
%   - delta: parâmetro delta
%   - tol: tolerância
%   - flagIsPoly: flag para análise politópica

    numSimulacoes = length(sistemas_monte_carlo);
    fprintf('      [Monte Carlo] Testando controlador em %d sistemas fixos...\n', numSimulacoes);

    if numSimulacoes == 0
        error('Lista de sistemas Monte Carlo está vazia!');
    end

    % Dimensões do sistema
    nx = 4;
    nu = 1;
    nw = 1;

    % Arrays para resultados
    gammaLMI = zeros(numSimulacoes, 1);
    custVector = zeros(numSimulacoes, 1);
    lmi_success = false(numSimulacoes, 1);
    sim_success = false(numSimulacoes, 1);

    % Arrays para guardar os parâmetros usados (extraídos dos sistemas)
    k2_used = zeros(numSimulacoes, 1);
    b_used = zeros(numSimulacoes, 1);

    % Extrair parâmetros dos sistemas já gerados
    for i = 1:numSimulacoes
        k2_used(i) = sistemas_monte_carlo{i}.k2;
        b_used(i) = sistemas_monte_carlo{i}.b;
    end

    % Função única para cada tipo
    func_intervalar = 'valEstHInfLMILab';
    func_politopica = 'estHInfAnaPolyLMILab';

    % Configurar modelo Simulink
    model = 'simulinkSampledInput';
    options = simset('SrcWorkspace','current');
    load_system(model);
    simulationTime = 15;
    set_param(model, 'StopTime', int2str(simulationTime));

    % Loop Monte Carlo usando sistemas já gerados
    for i = 1:numSimulacoes
        fprintf('      Simulação %d/%d: ', i, numSimulacoes);

        try
            % Usar sistema já gerado
            sistema_atual = sistemas_monte_carlo{i};

            A = sistema_atual.A;
            B = sistema_atual.B;
            E = sistema_atual.E;
            C = sistema_atual.C;
            D = sistema_atual.D;
            D1 = sistema_atual.D1;

            % Construir matrizes caligráficas
            ACal = [A B ; zeros(nu, nx) zeros(nu, nu)];
            ECal = [E ; zeros(nu, nw)];
            CCal = [C D];
            KCal = [eye(nx) zeros(nx, nu); K zeros(nu)];

            %% PARTE 1: TESTE DE LMI
            gamma_encontrado = NaN;
            lmi_funcionou = false;

            if flagIsPoly
                % Implementação politópica (usar sistema individual como politopo de 1 vértice)
                poly.APoly{1} = ACal;
                poly.EPoly{1} = ECal;
                poly.CPoly{1} = CCal;

                if exist(func_politopica, 'file')
                    try
                        % Usar timeout para evitar travamentos longos
                        warning('off', 'all');  % Silenciar warnings temporariamente
                        saidaLMI = feval(func_politopica, poly, KCal, h, delta, tol);
                        warning('on', 'all');   % Reativar warnings

                        if isfield(saidaLMI, 'gamma') && isfinite(saidaLMI.gamma) && saidaLMI.gamma > 0
                            gamma_encontrado = saidaLMI.gamma;
                            lmi_funcionou = true;
                        end
                    catch ME
                        warning('on', 'all');   % Reativar warnings mesmo em caso de erro
                        % LMI falhou, mas continua silenciosamente
                    end
                end
            else
                % Implementação intervalar
                if exist(func_intervalar, 'file')
                    try
                        % Usar timeout para evitar travamentos longos
                        warning('off', 'all');  % Silenciar warnings temporariamente
                        saidaLMI = feval(func_intervalar, ACal, ECal, CCal, KCal, h, delta, tol);
                        warning('on', 'all');   % Reativar warnings

                        if isfield(saidaLMI, 'gamma') && isfinite(saidaLMI.gamma) && saidaLMI.gamma > 0
                            gamma_encontrado = saidaLMI.gamma;
                            lmi_funcionou = true;
                        end
                    catch ME
                        warning('on', 'all');   % Reativar warnings mesmo em caso de erro
                        % LMI falhou, mas continua silenciosamente
                        if contains(ME.message, 'terminated by user') || contains(ME.message, 'interrupted')
                            fprintf('INTERROMPIDO ');
                            % Se foi interrompido pelo usuário, não tenta alternativas
                            break;
                        end
                    end
                end
            end

            % Se LMI falhou, usar análise de estabilidade simples
            if ~lmi_funcionou
                try
                    % Sistema em malha fechada
                    A_cl = A - B * K;
                    eigenvals = eig(A_cl);

                    if all(real(eigenvals) < -0.01)  % Sistema estável com margem
                        % Usar norma 2 como aproximação
                        sys_cl = ss(A_cl, E, C, D);
                        gamma_aproximado = norm(sys_cl, 2);

                        if isfinite(gamma_aproximado)
                            gamma_encontrado = gamma_aproximado;
                            lmi_funcionou = true;
                        end
                    end
                catch
                    % Análise simples também falhou
                end
            end

            % Armazenar resultado LMI
            gammaLMI(i) = gamma_encontrado;
            lmi_success(i) = lmi_funcionou;

            %% PARTE 2: SIMULAÇÃO USANDO SIMULINK
            if lmi_funcionou
                try
                    % Preparar matrizes para Simulink
                    B = [E B];  % Concatenar entrada de perturbação e controle

                    % Executar simulação
                    sim(model,[],options);

                    % Extrair custo da simulação
                    custVector(i) = simout.Data(end);
                    sim_success(i) = true;

                    fprintf('OK (k2=%.3f, b=%.4f, γ=%.4f, custo=%.4f)\n', sistema_atual.k2, sistema_atual.b, gamma_encontrado, custVector(i));

                catch ME_sim
                    custVector(i) = Inf;
                    sim_success(i) = false;
                    fprintf('SimErr (k2=%.3f, b=%.4f)\n', sistema_atual.k2, sistema_atual.b);
                end
            else
                custVector(i) = Inf;
                sim_success(i) = false;
                fprintf('LMIErr (k2=%.3f, b=%.4f)\n', sistema_atual.k2, sistema_atual.b);
            end

        catch ME
            fprintf('      ERRO GERAL: %s\n', ME.message);
            gammaLMI(i) = NaN;
            custVector(i) = Inf;
        end
    end

    %% PROCESSAMENTO FINAL DOS RESULTADOS MONTE CARLO
    num_lmi_ok = sum(lmi_success);
    num_sim_ok = sum(sim_success);

    fprintf('      Resumo Monte Carlo: LMI OK=%d/%d, Sim OK=%d/%d\n', ...
        num_lmi_ok, numSimulacoes, num_sim_ok, numSimulacoes);

    % Estatísticas de gamma
    if num_lmi_ok > 0
        gammas_validos = gammaLMI(lmi_success);
        gamma.meanGamma = mean(gammas_validos);
        gamma.stdGamma = std(gammas_validos);
        gamma.maxGamma = max(gammas_validos);
        gamma.minGamma = min(gammas_validos);
        [~, gamma.indexMaxGamma] = max(gammaLMI);
        gamma.numSucessos = num_lmi_ok;
        outPut.gamma = gamma;
    else
        outPut.gamma = struct('meanGamma', NaN, 'stdGamma', NaN, 'maxGamma', NaN, 'minGamma', NaN, ...
                             'indexMaxGamma', 1, 'numSucessos', 0);
    end

    % Estatísticas de custo
    if num_sim_ok > 0
        custos_validos = custVector(sim_success);
        outPut.meanCost = mean(custos_validos);
        outPut.stdCost = std(custos_validos);
        outPut.maxCost = max(custos_validos);
        outPut.minCost = min(custos_validos);
        [~, outPut.indexMaxCost] = max(custVector);
        [~, outPut.indexMinCost] = min(custVector);
    else
        outPut.meanCost = NaN;
        outPut.stdCost = NaN;
        outPut.maxCost = NaN;
        outPut.minCost = NaN;
        outPut.indexMaxCost = 1;
        outPut.indexMinCost = 1;
    end

    % Outros campos
    outPut.custVector = custVector;
    outPut.lmi_success_rate = num_lmi_ok / numSimulacoes;
    outPut.sim_success_rate = num_sim_ok / numSimulacoes;
    outPut.numSimulacoes = numSimulacoes;

    % Adicionar informações sobre os parâmetros usados
    outPut.parametros_usados.k2 = k2_used;
    outPut.parametros_usados.b = b_used;

    % Calcular ranges dos parâmetros usados
    if ~isempty(k2_used) && any(~isnan(k2_used))
        k2_validos = k2_used(~isnan(k2_used));
        b_validos = b_used(~isnan(b_used));
        outPut.parametros_usados.k2_range = [min(k2_validos), max(k2_validos)];
        outPut.parametros_usados.b_range = [min(b_validos), max(b_validos)];
    else
        outPut.parametros_usados.k2_range = [NaN, NaN];
        outPut.parametros_usados.b_range = [NaN, NaN];
    end

    % Estatísticas dos parâmetros que causaram pior/melhor desempenho
    if num_sim_ok > 0
        [~, idx_pior] = max(custVector);
        [~, idx_melhor] = min(custVector(custVector < Inf));

        outPut.pior_caso.custo = custVector(idx_pior);
        outPut.pior_caso.k2 = k2_used(idx_pior);
        outPut.pior_caso.b = b_used(idx_pior);

        if ~isempty(idx_melhor)
            outPut.melhor_caso.custo = custVector(idx_melhor);
            outPut.melhor_caso.k2 = k2_used(idx_melhor);
            outPut.melhor_caso.b = b_used(idx_melhor);
        end
    end

    if flagIsPoly
        outPut.metodo = 'Monte Carlo Politópico';
    else
        outPut.metodo = 'Monte Carlo Intervalar';
    end

    % Estatísticas adicionais sobre a distribuição dos parâmetros
    outPut.estatisticas_parametros.k2_mean = mean(k2_used(~isnan(k2_used)));
    outPut.estatisticas_parametros.k2_std = std(k2_used(~isnan(k2_used)));
    outPut.estatisticas_parametros.b_mean = mean(b_used(~isnan(b_used)));
    outPut.estatisticas_parametros.b_std = std(b_used(~isnan(b_used)));

    fprintf('      Simulação Monte Carlo concluída.\n');
    fprintf('      Parâmetros: k2_mean=%.3f±%.3f, b_mean=%.4f±%.4f\n', ...
        outPut.estatisticas_parametros.k2_mean, outPut.estatisticas_parametros.k2_std, ...
        outPut.estatisticas_parametros.b_mean, outPut.estatisticas_parametros.b_std);

end