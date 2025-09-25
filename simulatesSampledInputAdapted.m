function [outPut] = simulatesSampledInputAdapted(combPoly, h, K, imageName, axisVector, delta, tol, flagIsPoly)
%   SIMULATESAMPLEDINPUTADAPTED - Adaptação da simulação Simulink original
%
%   Esta função adapta a simulação original do Simulink (simulinkSampledInput)
%   para trabalhar com a estrutura de dados atual, mantendo a fidelidade
%   ao modelo testado.
%
%   Parâmetros:
%   - combPoly: estrutura do politopo
%   - h: período de amostragem
%   - K: ganho do controlador
%   - imageName: nome da imagem
%   - axisVector: limites dos eixos
%   - delta: parâmetro de discretização
%   - tol: tolerância
%   - flagIsPoly: flag para tipo de análise

    % Verificar argumentos
    if nargin < 8
        error('simulatesSampledInputAdapted: Requer 8 argumentos de entrada.');
    end

    % Verificar estrutura
    if ~isstruct(combPoly) || ~isfield(combPoly, 'A') || ~isfield(combPoly.A, 'alphaVecs')
        error('combPoly deve ser uma estrutura com campo A.alphaVecs');
    end

    fprintf('      [Simulação Adaptada Simulink] Processando %d vértices...\n', length(combPoly.A.alphaVecs));

    % Verificar se o modelo Simulink está disponível
    model = 'simulinkSampledInput';
    try
        load_system(model);
        simulink_available = true;
        fprintf('      Modelo Simulink carregado com sucesso.\n');
    catch
        simulink_available = false;
        fprintf('      ⚠ Modelo Simulink não disponível. Usando simulação numérica.\n');
    end

    % Dimensões do sistema
    try
        nx = size(combPoly.A.polytopicMatrices{1}, 1);
        if isfield(combPoly, 'B') && isfield(combPoly.B, 'polytopicMatrices')
            % Estrutura atual (B unificado)
            nu = size(combPoly.B.polytopicMatrices{1}, 2);
            current_structure = true;
        elseif isfield(combPoly, 'B2') && isfield(combPoly.B2, 'polytopicMatrices')
            % Estrutura original do Simulink (B1, B2 separados)
            nu = size(combPoly.B2.polytopicMatrices{1}, 2);
            current_structure = false;
        else
            error('Estrutura de entrada B não reconhecida');
        end

        if isfield(combPoly, 'E') && isfield(combPoly.E, 'polytopicMatrices')
            nw = size(combPoly.E.polytopicMatrices{1}, 2);
        else
            nw = nu; % Assumir mesmo tamanho se E não existir
        end
    catch
        error('Erro ao extrair dimensões do politopo');
    end

    fprintf('      Sistema: nx=%d, nu=%d, nw=%d\n', nx, nu, nw);

    % Configurar simulação
    if simulink_available
        % Parâmetros do Simulink (conforme original)
        simulationTime = 30;  % Tempo original do Simulink
        options = simset('SrcWorkspace', 'current');
        set_param(model, 'StopTime', int2str(simulationTime));
    else
        % Parâmetros para simulação numérica
        simulationTime = 30;  % Manter consistência
        dt = 0.01;
        t = 0:dt:simulationTime;
        numSteps = length(t);
    end

    % Arrays para resultados
    gammaLMI = zeros(length(combPoly.A.alphaVecs), 1);
    custVector = zeros(length(combPoly.A.alphaVecs), 1);
    lmi_success = false(length(combPoly.A.alphaVecs), 1);
    states = cell(length(combPoly.A.alphaVecs), 1);
    time = [];
    lenTime = 0;

    % Loop pelos vértices
    for i = 1:length(combPoly.A.alphaVecs)
        fprintf('      Processando vértice %d/%d...\n', i, length(combPoly.A.alphaVecs));

        try
            % Extrair matrizes do vértice i
            A = combPoly.A.polytopicMatrices{i};

            % Adaptar estrutura de dados
            if current_structure
                % Estrutura atual: B, E, C, D unificados
                B_control = combPoly.B.polytopicMatrices{i};
                if isfield(combPoly, 'E') && isfield(combPoly.E, 'polytopicMatrices')
                    B_dist = combPoly.E.polytopicMatrices{i};
                else
                    B_dist = B_control; % Usar B como E se não existir
                end
                C = combPoly.C.polytopicMatrices{i};
                D_control = combPoly.D.polytopicMatrices{i};

                % Converter para formato original do Simulink (CORREÇÃO CRÍTICA)
                B1 = B_dist;   % Entrada de perturbação (w)
                B2 = B_control; % Entrada de controle (u)
                D1 = zeros(size(C, 1), size(B_dist, 2));  % Sem feedthrough da perturbação
                D2 = D_control; % Feedthrough do controle

                % Para análise LMI, usar matrizes de controle
                B = B_control;
                E = B_dist;
                D = D_control;

            else
                % Estrutura original do Simulink (CORREÇÃO DO BUG CRÍTICO)
                B1 = combPoly.B1.polytopicMatrices{i}; % CORREÇÃO: era B2 antes
                B2 = combPoly.B2.polytopicMatrices{i};
                C = combPoly.C.polytopicMatrices{i};
                D1 = combPoly.D1.polytopicMatrices{i};
                D2 = combPoly.D2.polytopicMatrices{i};

                % Para análise LMI, criar matrizes unificadas
                B = B2;  % Entrada de controle
                E = B1;  % Entrada de perturbação
                D = D2;  % Feedthrough do controle
            end

            %% ANÁLISE LMI (usando estrutura unificada)
            % Construir matrizes caligráficas
            ACal = [A B ; zeros(nu, nx) zeros(nu, nu)];
            ECal = [E ; zeros(nu, size(E, 2))];
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
            else
                gammaLMI(i) = NaN;
                lmi_success(i) = false;
            end

            %% SIMULAÇÃO (Simulink ou numérica)
            try
                if simulink_available
                    % Usar Simulink (conforme implementação original)
                    % Definir variáveis no workspace para o Simulink
                    assignin('base', 'A', A);
                    assignin('base', 'B1', B1);
                    assignin('base', 'B2', B2);
                    assignin('base', 'B', [B1 B2]); % Matriz combinada para o modelo
                    assignin('base', 'C', C);
                    assignin('base', 'D1', D1);
                    assignin('base', 'D2', D2);
                    assignin('base', 'K', K);
                    assignin('base', 'h', h);

                    % Executar simulação
                    sim(model, [], options);

                    % Obter resultados (assumindo que o modelo produz 'simout')
                    simout = evalin('base', 'simout');
                    custVector(i) = simout.Data(end);

                    % Obter tempo e estados se disponíveis
                    if i == 1
                        time = simout.Time;
                        lenTime = length(time);
                    end

                    % Obter estados se disponíveis
                    try
                        simoutStates = evalin('base', 'simoutStates');
                        states{i} = normalizeSizeVec(simoutStates.Data, lenTime);
                    catch
                        states{i} = NaN(lenTime, nx);
                    end

                else
                    % Simulação numérica melhorada (baseada no modelo Simulink)
                    fprintf('        Usando simulação numérica...\n');

                    % Estados iniciais
                    x = zeros(nx, numSteps);
                    x(:, 1) = zeros(nx, 1);  % Condição inicial zero

                    % Sinais de entrada (conforme modelo Simulink original)
                    % Referência: degrau unitário em t=1s
                    r_signal = ones(numSteps, 1);
                    r_signal(t < 1) = 0;

                    % Perturbação: impulso em t=2s
                    w_signal = zeros(numSteps, size(B1, 2));
                    impulse_idx = find(t >= 2, 1);
                    if ~isempty(impulse_idx)
                        w_signal(impulse_idx, :) = 0.1; % Impulso menor para estabilidade
                    end

                    % Saída e controle
                    y = zeros(size(C, 1), numSteps);
                    u_control = zeros(size(B2, 2), numSteps);
                    e_signal = zeros(size(B2, 2), numSteps); % Erro

                    % Simulação com sistema amostrado corretamente
                    sample_counter = 0;
                    last_sample_time = 0;

                    for k = 1:numSteps-1
                        % Verificar se é momento de amostragem
                        current_time = t(k);
                        if (current_time - last_sample_time) >= (h - dt/2) || k == 1
                            % Momento de amostragem - atualizar controle
                            e_signal(:, k) = r_signal(k) - (C * x(:, k))'; % Erro de rastreamento
                            u_control(:, k) = K * x(:, k); % Lei de controle por realimentação de estado
                            last_sample_time = current_time;
                            sample_counter = sample_counter + 1;
                        else
                            % Entre amostras - manter controle anterior (ZOH)
                            e_signal(:, k) = e_signal(:, k-1);
                            u_control(:, k) = u_control(:, k-1);
                        end

                        % Entrada total do sistema: referência + controle realimentado
                        % Assumindo que a estrutura é: x' = Ax + B2*u + B1*w, y = Cx + D2*u + D1*w
                        % onde u é o sinal de controle total
                        u_total = -u_control(:, k); % Sinal de controle (negativo por realimentação)

                        % Dinâmica do sistema contínuo
                        x_dot = A * x(:, k) + B2 * u_total + B1 * w_signal(k, :)';

                        % Integração numérica (Euler)
                        x(:, k+1) = x(:, k) + dt * x_dot;

                        % Calcular saída
                        y(:, k) = C * x(:, k) + D2 * u_total + D1 * w_signal(k, :)';
                    end

                    % Última saída
                    u_total_final = -u_control(:, end);
                    y(:, end) = C * x(:, end) + D2 * u_total_final + D1 * w_signal(end, :)';

                    % Verificar estabilidade a posteriori
                    final_states = x(:, end);
                    if any(abs(final_states) > 1e3)
                        fprintf('        ⚠ Sistema pode estar instável (estados finais grandes)\n');
                        custVector(i) = Inf;
                    else
                        % Calcular custo quadrático (integral da norma da saída)
                        cost_integrand = sum(y.^2, 1);  % Soma das saídas ao quadrado
                        custVector(i) = trapz(t, cost_integrand);  % Integral no tempo
                        fprintf('        ✓ Simulação estável, custo = %.6f\n', custVector(i));
                    end

                    % Armazenar estados
                    states{i} = x';
                    if i == 1
                        time = t';
                        lenTime = length(time);
                    end
                end

            catch ME_sim
                fprintf('        Erro na simulação: %s\n', ME_sim.message);
                custVector(i) = Inf;
                if lenTime > 0
                    states{i} = NaN(lenTime, nx);
                else
                    states{i} = NaN(1000, nx); % Tamanho padrão
                end
            end

        catch ME
            fprintf('        Vértice %d: ERRO GERAL - %s\n', i, ME.message);
            gammaLMI(i) = NaN;
            custVector(i) = Inf;
            lmi_success(i) = false;
            if lenTime > 0
                states{i} = NaN(lenTime, nx);
            else
                states{i} = NaN(1000, nx); % Tamanho padrão
            end
        end
    end

    % Fechar modelo Simulink se foi aberto
    if simulink_available
        try
            close_system(model, 0);
        catch
            % Ignorar erros ao fechar
        end
    end

    %% PROCESSAMENTO DOS RESULTADOS (conforme implementação original)

    % Filtrar resultados válidos
    custos_validos = custVector(isfinite(custVector));
    gammas_validos = gammaLMI(lmi_success & isfinite(gammaLMI));

    numSucessosLMI = sum(lmi_success);
    numSucessosSim = sum(isfinite(custVector));

    fprintf('      LMIs bem-sucedidas: %d/%d\n', numSucessosLMI, length(combPoly.A.alphaVecs));
    fprintf('      Simulações válidas: %d/%d\n', numSucessosSim, length(combPoly.A.alphaVecs));

    % Calcular estatísticas (conforme original)
    if numSucessosSim > 0
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
        outPut.custVector = custVector;
        outPut.meanCost = NaN;
        outPut.stdCost = NaN;
        outPut.maxCost = NaN;
        outPut.indexMaxCost = 1;
    end

    if numSucessosLMI > 0
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
        outPut.gamma = struct('meanGamma', NaN, 'stdGamma', NaN, 'maxGamma', NaN, ...
                             'indexMaxGamma', 1, 'numSucessos', 0);
    end

    % Informações adicionais (conforme original)
    outPut.time = time;
    outPut.states = states;
    outPut.lmi_success_rate = numSucessosLMI / length(combPoly.A.alphaVecs);
    outPut.sim_success_rate = numSucessosSim / length(combPoly.A.alphaVecs);

    if simulink_available
        if flagIsPoly
            outPut.metodo = 'Politópico (Simulink)';
        else
            outPut.metodo = 'Intervalar (Simulink)';
        end
    else
        if flagIsPoly
            outPut.metodo = 'Politópico (Numérico)';
        else
            outPut.metodo = 'Intervalar (Numérico)';
        end
    end

    fprintf('      Simulação adaptada concluída.\n');

end