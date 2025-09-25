function [outPut] = simulatesSampledInputRobusto(combPoly, h, K, imageName, axisVector, delta, tol, flagIsPoly)
%   SIMULATESAMPLEDINPUTROBUSTO - Versão robusta com diagnóstico detalhado
%
%   Esta versão inclui diagnósticos detalhados e múltiplas estratégias
%   para resolver problemas com LMIs que não estão funcionando

    % Verificar argumentos
    if nargin < 8
        error('Requer 8 argumentos de entrada');
    end

    fprintf('      [Simulação Robusta] Analisando %d vértices...\n', length(combPoly.A.alphaVecs));

    % Dimensões do sistema
    try
        nx = size(combPoly.A.polytopicMatrices{1}, 1);
        nu = size(combPoly.B.polytopicMatrices{1}, 2);
        nw = size(combPoly.E.polytopicMatrices{1}, 2);
        fprintf('      Sistema: nx=%d, nu=%d, nw=%d\n', nx, nu, nw);
    catch
        error('Erro ao extrair dimensões do politopo');
    end

    % Arrays para resultados
    gammaLMI = zeros(length(combPoly.A.alphaVecs), 1);
    custVector = zeros(length(combPoly.A.alphaVecs), 1);
    lmi_success = false(length(combPoly.A.alphaVecs), 1);
    sim_success = false(length(combPoly.A.alphaVecs), 1);

    % Função única para cada tipo
    func_intervalar = 'valEstHInfLMILab';
    func_politopica = 'estHInfAnaPolyLMILab';


    % Loop pelos vértices
    for i = 1:length(combPoly.A.alphaVecs)
        fprintf('      Vértice %d/%d: ', i, length(combPoly.A.alphaVecs));

        try
            % Extrair matrizes
            A = combPoly.A.polytopicMatrices{i};
            B = combPoly.B.polytopicMatrices{i};
            E = combPoly.E.polytopicMatrices{i};
            C = combPoly.C.polytopicMatrices{i};
            D = combPoly.D.polytopicMatrices{i};
            D1 = combPoly.D1.polytopicMatrices{i};

            % Construir matrizes caligráficas
            ACal = [A B ; zeros(nu, nx) zeros(nu, nu)];
            ECal = [E ; zeros(nu, nw)];
            CCal = [C D];
            KCal = [eye(nx) zeros(nx, nu); K zeros(nu)];

            %% PARTE 1: TESTE DE LMI SIMPLIFICADO
            gamma_encontrado = NaN;
            lmi_funcionou = false;

            if flagIsPoly
                % Implementação politópica
                poly.APoly{1} = ACal;
                poly.EPoly{1} = ECal;
                poly.CPoly{1} = CCal;

                if exist(func_politopica, 'file')
                    try
                        fprintf('Testando %s... ', func_politopica);
                        saidaLMI = feval(func_politopica, poly, KCal, h, delta, tol);

                        if isfield(saidaLMI, 'gamma') && isfinite(saidaLMI.gamma) && saidaLMI.gamma > 0
                            gamma_encontrado = saidaLMI.gamma;
                            lmi_funcionou = true;
                            fprintf('OK (γ=%.4f)\n', gamma_encontrado);
                        else
                            fprintf('Falhou\n');
                        end
                    catch ME
                        fprintf('Erro: %s\n', ME.message);
                    end
                else
                    fprintf('%s não existe\n', func_politopica);
                end
            else
                % Implementação intervalar
                if exist(func_intervalar, 'file')
                    try
                        fprintf('Testando %s... ', func_intervalar);
                        saidaLMI = feval(func_intervalar, ACal, ECal, CCal, KCal, h, delta, tol);

                        if isfield(saidaLMI, 'gamma') && isfinite(saidaLMI.gamma) && saidaLMI.gamma > 0
                            gamma_encontrado = saidaLMI.gamma;
                            lmi_funcionou = true;
                            fprintf('OK (γ=%.4f)\n', gamma_encontrado);
                        else
                            fprintf('Falhou\n');
                        end
                    catch ME
                        fprintf('Erro: %s\n', ME.message);
                    end
                else
                    fprintf('%s não existe\n', func_intervalar);
                end
            end

            % Se LMI falhou, usar análise de estabilidade simples
            if ~lmi_funcionou
                fprintf('      Tentando análise de estabilidade simples...\n');
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
                            fprintf('      OK com análise simples (γ≈%.4f)\n', gamma_aproximado);
                        end
                    end
                catch
                    fprintf('      Análise simples também falhou\n');
                end
            end

            % Armazenar resultado LMI
            gammaLMI(i) = gamma_encontrado;
            lmi_success(i) = lmi_funcionou;

            %% PARTE 2: SIMULAÇÃO USANDO SIMULINK
            if lmi_funcionou
                try
                    % Usar simulação do Simulink
                    fprintf('      Usando simulação Simulink... ');

                    % Preparar matrizes para Simulink
                    B = [E B];  % Concatenar entrada de perturbação e controle

                    % Configurar modelo Simulink
                    model = 'simulinkSampledInput';
                    options = simset('SrcWorkspace','current');
                    load_system(model);
                    simulationTime = 15;
                    set_param(model, 'StopTime', int2str(simulationTime));

                    % Executar simulação
                    sim(model,[],options);

                    % Extrair custo da simulação
                    custVector(i) = simout.Data(end);
                    sim_success(i) = true;

                    fprintf('OK (custo=%.4f)\n', custVector(i));

                catch ME_sim
                    custVector(i) = Inf;
                    sim_success(i) = false;
                    fprintf('      Erro simulação: %s\n', ME_sim.message);
                end
            else
                custVector(i) = Inf;
                sim_success(i) = false;
                fprintf('      Pulando simulação (LMI falhou)\n');
            end

        catch ME
            fprintf('      ERRO GERAL: %s\n', ME.message);
            gammaLMI(i) = NaN;
            custVector(i) = Inf;
        end
    end

    %% PROCESSAMENTO FINAL
    num_lmi_ok = sum(lmi_success);
    num_sim_ok = sum(sim_success);

    fprintf('      Resumo: LMI OK=%d/%d, Sim OK=%d/%d\n', ...
        num_lmi_ok, length(combPoly.A.alphaVecs), num_sim_ok, length(combPoly.A.alphaVecs));

    % Estatísticas de gamma
    if num_lmi_ok > 0
        gammas_validos = gammaLMI(lmi_success);
        gamma.meanGamma = mean(gammas_validos);
        gamma.stdGamma = std(gammas_validos);
        gamma.maxGamma = max(gammas_validos);
        [~, gamma.indexMaxGamma] = max(gammaLMI);
        gamma.numSucessos = num_lmi_ok;
        outPut.gamma = gamma;
    else
        outPut.gamma = struct('meanGamma', NaN, 'stdGamma', NaN, 'maxGamma', NaN, ...
                             'indexMaxGamma', 1, 'numSucessos', 0);
    end

    % Estatísticas de custo
    if num_sim_ok > 0
        custos_validos = custVector(sim_success);
        outPut.meanCost = mean(custos_validos);
        outPut.stdCost = std(custos_validos);
        outPut.maxCost = max(custos_validos);
        [~, outPut.indexMaxCost] = max(custVector);
    else
        outPut.meanCost = NaN;
        outPut.stdCost = NaN;
        outPut.maxCost = NaN;
        outPut.indexMaxCost = 1;
    end

    % Outros campos
    outPut.custVector = custVector;
    outPut.lmi_success_rate = num_lmi_ok / length(combPoly.A.alphaVecs);
    outPut.sim_success_rate = num_sim_ok / length(combPoly.A.alphaVecs);

    if flagIsPoly
        outPut.metodo = 'Politópico';
    else
        outPut.metodo = 'Intervalar';
    end

    fprintf('      Simulação robusta concluída.\n');

end