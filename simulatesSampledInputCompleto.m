function [outPut] = simulatesSampledInputCompleto(combPoly, h, K, imageName, axisVector, delta, tol, flagIsPoly)
%   SIMULATESAMPLEDINPUTCOMPLETO - Versão completa com simulação temporal
%
%   Esta versão combina análise LMI com simulação temporal usando integração
%   numérica (sem depender do Simulink), conforme exemplos do projeto.
%
%   Baseada nos exemplos: imageName = 'KIntHibrido'; flagIsPoly = false;
%
%   Uso: output = simulatesSampledInputCompleto(combPoly, h, K, imageName, axisVector, delta, tol, flagIsPoly)
%
%   Parâmetros:
%   - combPoly: estrutura do politopo com campos A, B, E, C, D, D1
%   - h: período de amostragem
%   - K: ganho do controlador (vetor linha)
%   - imageName: nome da imagem (string)
%   - axisVector: limites dos eixos [xmin xmax ymin ymax]
%   - delta: parâmetro de discretização
%   - tol: tolerância
%   - flagIsPoly: true para politópico, false para intervalar

    % Verificar número de argumentos
    if nargin < 8
        error('simulatesSampledInputCompleto: Requer 8 argumentos de entrada.\nUso: simulatesSampledInputCompleto(combPoly, h, K, imageName, axisVector, delta, tol, flagIsPoly)');
    end

    % Verificar se combPoly tem a estrutura esperada
    if ~isstruct(combPoly) || ~isfield(combPoly, 'A') || ~isfield(combPoly.A, 'alphaVecs')
        error('combPoly deve ser uma estrutura com campo A.alphaVecs');
    end

    fprintf('      [Simulação Completa] Analisando %d vértices com LMI + Temporal...\n', length(combPoly.A.alphaVecs));

    % Dimensões do sistema
    try
        nx = size(combPoly.A.polytopicMatrices{1}, 1);
        nu = size(combPoly.B.polytopicMatrices{1}, 2);
        nw = size(combPoly.E.polytopicMatrices{1}, 2);
    catch
        error('Erro ao extrair dimensões do politopo');
    end

    fprintf('      Sistema: nx=%d, nu=%d, nw=%d\n', nx, nu, nw);

    % Parâmetros de simulação temporal
    simulationTime = 15;  % segundos (conforme original)
    dt = 0.01;           % passo de integração
    t = 0:dt:simulationTime;
    numSteps = length(t);

    % Arrays para armazenar resultados
    gammaLMI = zeros(length(combPoly.A.alphaVecs), 1);
    custVector = zeros(length(combPoly.A.alphaVecs), 1);
    lmi_success = false(length(combPoly.A.alphaVecs), 1);
    states = cell(length(combPoly.A.alphaVecs), 1);

    % Sinal de entrada (degrau unitário em t=1s)
    u_ref = ones(numSteps, 1);
    u_ref(t < 1) = 0;

    % Perturbação (impulso em t=2s)
    w_dist = zeros(numSteps, nw);
    impulse_idx = find(t >= 2, 1);
    if ~isempty(impulse_idx) && impulse_idx <= numSteps
        w_dist(impulse_idx, :) = 1;  % impulso unitário
    end

    % Loop pelos vértices do politopo
    for i = 1:length(combPoly.A.alphaVecs)
        fprintf('      Processando vértice %d/%d...\n', i, length(combPoly.A.alphaVecs));

        try
            % Extrair matrizes do vértice i
            A = combPoly.A.polytopicMatrices{i};
            B = combPoly.B.polytopicMatrices{i};
            E = combPoly.E.polytopicMatrices{i};
            C = combPoly.C.polytopicMatrices{i};
            D = combPoly.D.polytopicMatrices{i};

            %% PARTE 1: ANÁLISE LMI (conforme metodologia híbrida)
            % Construir matrizes caligráficas
            ACal = [A B ; zeros(nu, nx) zeros(nu, nu)];
            ECal = [E ; zeros(nu, nw)];
            CCal = [C D];
            KCal = [eye(nx) zeros(nx, nu); K zeros(nu)];

            % Executar LMI conforme flag
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
                            % Falha na LMI politópica
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
            else
                gammaLMI(i) = NaN;
                lmi_success(i) = false;
            end

            %% PARTE 2: SIMULAÇÃO TEMPORAL (usando integração numérica)
            try
                % Sistema em malha fechada amostrado
                A_cl = A - B * K;

                % Verificar estabilidade
                eigenvals = eig(A_cl);
                if all(real(eigenvals) < -1e-6)  % Sistema estável

                    % Estado inicial
                    x = zeros(nx, numSteps);
                    x(:, 1) = zeros(nx, 1);  % condição inicial zero

                    % Saída
                    y = zeros(size(C, 1), numSteps);

                    % Controlador amostrado (ZOH)
                    u_control = zeros(nu, numSteps);

                    % Simulação com amostragem
                    for k = 1:numSteps-1
                        % Controle amostrado (atualizado a cada h segundos)
                        if mod(t(k), h) < dt || k == 1
                            % Atualizar controle (referência - realimentação)
                            u_control(:, k) = K * x(:, k);
                        else
                            % Manter controle anterior (ZOH)
                            u_control(:, k) = u_control(:, k-1);
                        end

                        % Entrada total
                        u_total = u_ref(k) - u_control(:, k);

                        % Integração (Euler)
                        x_dot = A * x(:, k) + B * u_total + E * w_dist(k, :)';
                        x(:, k+1) = x(:, k) + dt * x_dot;

                        % Saída
                        y(:, k) = C * x(:, k) + D * u_total;
                    end

                    % Última saída
                    y(:, end) = C * x(:, end) + D * (u_ref(end) - u_control(:, end));

                    % Calcular custo quadrático (norma H2 aproximada)
                    cost_integrand = sum(y.^2, 1);  % soma sobre todas as saídas
                    custVector(i) = trapz(t, cost_integrand);  % integral temporal

                    % Armazenar estados para posterior análise
                    states{i} = x';

                else
                    % Sistema instável
                    custVector(i) = Inf;
                    states{i} = NaN(numSteps, nx);
                end

            catch ME_sim
                fprintf('        Erro na simulação temporal: %s\n', ME_sim.message);
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
        [maxCost, idx_max] = max(custos_validos);
        [~, indexMaxCost] = max(custVector);  % Index no array original

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
    outPut.time = t;
    outPut.states = states;
    outPut.lmi_success_rate = numSucessosLMI / length(combPoly.A.alphaVecs);
    outPut.sim_success_rate = numSucessosSim / length(combPoly.A.alphaVecs);

    if flagIsPoly
        outPut.metodo = 'Politópico';
    else
        outPut.metodo = 'Intervalar';
    end

    fprintf('      Simulação completa concluída.\n');

end