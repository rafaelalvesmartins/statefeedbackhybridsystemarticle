function [outPut] = simulatesSampledInputNumerico(combPoly, h, K, imageName, axisVector, delta, tol, flagIsPoly)
%   SIMULATESAMPLEDINPUTNUMERICO - Simulação puramente numérica (sem Simulink)
%
%   Esta função implementa simulação amostrada usando apenas integração numérica,
%   eliminando dependências do Simulink e resolvendo problemas de estabilidade.
%
%   Parâmetros:
%   - combPoly: estrutura do politopo
%   - h: período de amostragem
%   - K: ganho do controlador
%   - imageName: nome da imagem (não usado)
%   - axisVector: limites dos eixos (não usado)
%   - delta: parâmetro de discretização
%   - tol: tolerância
%   - flagIsPoly: flag para tipo de análise

    % Verificar argumentos
    if nargin < 8
        error('simulatesSampledInputNumerico: Requer 8 argumentos de entrada.');
    end

    % Verificar estrutura
    if ~isstruct(combPoly) || ~isfield(combPoly, 'A') || ~isfield(combPoly.A, 'alphaVecs')
        error('combPoly deve ser uma estrutura com campo A.alphaVecs');
    end

    fprintf('      [Simulação Numérica Pura] Processando %d vértices...\n', length(combPoly.A.alphaVecs));

    % Dimensões do sistema
    try
        nx = size(combPoly.A.polytopicMatrices{1}, 1);
        if isfield(combPoly, 'B') && isfield(combPoly.B, 'polytopicMatrices')
            nu = size(combPoly.B.polytopicMatrices{1}, 2);
        else
            nu = 1; % Assumir entrada escalar
        end
        if isfield(combPoly, 'E') && isfield(combPoly.E, 'polytopicMatrices')
            nw = size(combPoly.E.polytopicMatrices{1}, 2);
        else
            nw = 1; % Assumir perturbação escalar
        end
    catch
        error('Erro ao extrair dimensões do politopo');
    end

    fprintf('      Sistema: nx=%d, nu=%d, nw=%d\n', nx, nu, nw);

    % Parâmetros de simulação
    simulationTime = 30;  % Tempo total (conforme Simulink)
    dt = 0.01;           % Passo de integração
    t = 0:dt:simulationTime;
    numSteps = length(t);

    % Arrays para resultados
    gammaLMI = zeros(length(combPoly.A.alphaVecs), 1);
    custVector = zeros(length(combPoly.A.alphaVecs), 1);
    lmi_success = false(length(combPoly.A.alphaVecs), 1);
    states = cell(length(combPoly.A.alphaVecs), 1);

    % Loop pelos vértices
    for i = 1:length(combPoly.A.alphaVecs)
        fprintf('      Processando vértice %d/%d...\n', i, length(combPoly.A.alphaVecs));

        try
            % Extrair matrizes do vértice i
            A = combPoly.A.polytopicMatrices{i};
            B = combPoly.B.polytopicMatrices{i};
            if isfield(combPoly, 'E') && isfield(combPoly.E, 'polytopicMatrices')
                E = combPoly.E.polytopicMatrices{i};
            else
                E = B; % Usar B como E se não existir
            end
            C = combPoly.C.polytopicMatrices{i};
            D = combPoly.D.polytopicMatrices{i};

            %% ANÁLISE LMI
            % Construir matrizes caligráficas
            ACal = [A B ; zeros(nu, nx) zeros(nu, nu)];
            ECal = [E ; zeros(nu, nw)];
            CCal = [C D];
            KCal = [eye(nx) zeros(nx, nu); K zeros(nu)];

            % Executar LMI
            lmi_sucesso = false;
            try
                if flagIsPoly
                    % Análise politópica
                    poly.APoly{1} = ACal;
                    poly.EPoly{1} = ECal;
                    poly.CPoly{1} = CCal;

                    try
                        saidaLMI = estHInfAnaPolyLMILab(poly, KCal, h, delta, tol);
                        lmi_sucesso = true;
                    catch
                        try
                            saidaLMI = estHInfAnaPolyLMILabOriginal(poly, KCal, h, delta, tol);
                            lmi_sucesso = true;
                        catch
                            saidaLMI.gamma = NaN;
                        end
                    end
                else
                    % Análise intervalar
                    try
                        saidaLMI = valEstHInfLMILab(ACal, ECal, CCal, KCal, h, delta, tol);
                        lmi_sucesso = true;
                    catch
                        try
                            saidaLMI = valEstHInfLMILabSemInt(ACal, ECal, CCal, KCal, h, delta, tol);
                            lmi_sucesso = true;
                        catch
                            try
                                saidaLMI = valEstHInfLMILabInt(ACal, ECal, CCal, KCal, h, delta, tol);
                                lmi_sucesso = true;
                            catch
                                saidaLMI.gamma = NaN;
                            end
                        end
                    end
                end
            catch
                saidaLMI.gamma = NaN;
            end

            % Armazenar γ da LMI
            if lmi_sucesso && isfield(saidaLMI, 'gamma') && isfinite(saidaLMI.gamma)
                gammaLMI(i) = saidaLMI.gamma;
                lmi_success(i) = true;
                fprintf('        LMI OK (γ=%.4f)\n', saidaLMI.gamma);
            else
                gammaLMI(i) = NaN;
                lmi_success(i) = false;
                fprintf('        LMI FALHOU\n');
            end

            %% SIMULAÇÃO NUMÉRICA CORRIGIDA
            try
                fprintf('        Iniciando simulação numérica...\n');

                % Estados e sinais
                x = zeros(nx, numSteps);
                x(:, 1) = zeros(nx, 1);  % Estado inicial zero

                % Sinal de referência: degrau unitário em t=1s
                r_signal = zeros(numSteps, 1);
                r_signal(t >= 1) = 1;

                % Perturbação: impulso em t=2s
                w_signal = zeros(numSteps, nw);
                impulse_indices = find(t >= 2 & t < 2.1); % Impulso curto
                if ~isempty(impulse_indices)
                    w_signal(impulse_indices, :) = 0.1; % Amplitude moderada
                end

                % Variáveis de controle
                u_control = zeros(nu, numSteps);
                y = zeros(size(C, 1), numSteps);

                % Variáveis de amostragem
                last_sample_time = -h; % Forçar primeira amostragem

                % Simulação principal
                for k = 1:numSteps-1
                    current_time = t(k);

                    % Verificar se é momento de amostragem
                    if (current_time - last_sample_time) >= h - dt/2
                        % Atualizar controle baseado no estado atual
                        u_control(:, k) = -K * x(:, k); % Realimentação de estado
                        last_sample_time = current_time;
                    else
                        % Zero-order hold: manter controle anterior
                        u_control(:, k) = u_control(:, max(1, k-1));
                    end

                    % Entrada total: controle + referência (se aplicável)
                    % Assumindo que a referência entra via perturbação ou diretamente
                    u_total = u_control(:, k);
                    w_total = w_signal(k, :)' + r_signal(k) * 0.1; % Referência como perturbação

                    % Dinâmica do sistema: dx/dt = Ax + Bu + Ew
                    x_dot = A * x(:, k) + B * u_total + E * w_total;

                    % Integração (Euler)
                    x(:, k+1) = x(:, k) + dt * x_dot;

                    % Saída: y = Cx + Du + Dw (simplificado)
                    y(:, k) = C * x(:, k) + D * u_total;
                end

                % Última saída
                y(:, end) = C * x(:, end) + D * u_control(:, end);

                % Verificação de estabilidade
                max_state = max(max(abs(x)));
                if max_state > 1e3
                    fprintf('        ⚠ Sistema instável (max|x|=%.2e)\n', max_state);
                    custVector(i) = Inf;
                else
                    % Custo quadrático: integral da norma da saída
                    cost_integrand = sum(y.^2, 1);
                    custVector(i) = trapz(t, cost_integrand);
                    fprintf('        ✓ Simulação estável, custo=%.6f\n', custVector(i));
                end

                % Armazenar estados
                states{i} = x';

            catch ME_sim
                fprintf('        ✗ Erro na simulação: %s\n', ME_sim.message);
                custVector(i) = Inf;
                states{i} = NaN(numSteps, nx);
            end

        catch ME
            fprintf('        Vértice %d: ERRO GERAL - %s\n', i, ME.message);
            gammaLMI(i) = NaN;
            custVector(i) = Inf;
            lmi_success(i) = false;
            states{i} = NaN(numSteps, nx);
        end
    end

    %% PROCESSAMENTO DOS RESULTADOS

    % Filtrar resultados válidos
    custos_validos = custVector(isfinite(custVector));
    gammas_validos = gammaLMI(lmi_success & isfinite(gammaLMI));

    numSucessosLMI = sum(lmi_success);
    numSucessosSim = sum(isfinite(custVector));

    fprintf('      LMIs bem-sucedidas: %d/%d\n', numSucessosLMI, length(combPoly.A.alphaVecs));
    fprintf('      Simulações válidas: %d/%d\n', numSucessosSim, length(combPoly.A.alphaVecs));

    %% RESULTADOS FINAIS
    if numSucessosSim > 0
        % Estatísticas de custo (simulação temporal)
        meanCost = mean(custos_validos);
        stdCost = std(custos_validos);
        [maxCost, ~] = max(custos_validos);
        [~, indexMaxCost] = max(custVector);

        outPut.custVector = custVector;
        outPut.meanCost = meanCost;
        outPut.stdCost = stdCost;
        outPut.maxCost = maxCost;
        outPut.indexMaxCost = indexMaxCost;

        fprintf('      Custo médio: %.6f, máximo: %.6f\n', meanCost, maxCost);

    else
        % Nenhuma simulação válida
        outPut.custVector = custVector;
        outPut.meanCost = NaN;
        outPut.stdCost = NaN;
        outPut.maxCost = NaN;
        outPut.indexMaxCost = 1;
    end

    if numSucessosLMI > 0
        % Estatísticas de γ (LMI)
        meanGamma = mean(gammas_validos);
        stdGamma = std(gammas_validos);
        [maxGamma, ~] = max(gammas_validos);
        [~, indexMaxGamma] = max(gammaLMI);

        gamma.meanGamma = meanGamma;
        gamma.stdGamma = stdGamma;
        gamma.maxGamma = maxGamma;
        gamma.indexMaxGamma = indexMaxGamma;
        gamma.numSucessos = numSucessosLMI;

        outPut.gamma = gamma;

        fprintf('      γ médio: %.6f, máximo: %.6f\n', meanGamma, maxGamma);

    else
        % Nenhuma LMI válida
        outPut.gamma = struct('meanGamma', NaN, 'stdGamma', NaN, 'maxGamma', NaN, ...
                             'indexMaxGamma', 1, 'numSucessos', 0);
    end

    % Informações adicionais
    outPut.time = t';
    outPut.states = states;
    outPut.lmi_success_rate = numSucessosLMI / length(combPoly.A.alphaVecs);
    outPut.sim_success_rate = numSucessosSim / length(combPoly.A.alphaVecs);

    if flagIsPoly
        outPut.metodo = 'Politópico (Numérico)';
    else
        outPut.metodo = 'Intervalar (Numérico)';
    end

    fprintf('      Simulação numérica concluída.\n');

end