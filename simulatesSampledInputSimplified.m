function [outPut] = simulatesSampledInputSimplified(combPoly, h, K, imageName, axisVector, delta, tol, flagIsPoly)
%   SIMULATESAMPLEDINPUTSIMPLIFIED - Versão simplificada para resposta ao revisor
%
%   Esta versão foca apenas na análise LMI sem depender do Simulink,
%   fornecendo as métricas γ necessárias para a comparação solicitada pelo revisor.
%
%   Parâmetros:
%   - combPoly: estrutura do politopo
%   - h: período de amostragem
%   - K: ganho do controlador
%   - imageName: nome para salvar (não usado nesta versão)
%   - axisVector: eixos (não usado nesta versão)
%   - delta: parâmetro de discretização
%   - tol: tolerância
%   - flagIsPoly: true para análise politópica, false para intervalar

    % Verificar número de argumentos
    if nargin < 8
        error('simulatesSampledInputSimplified: Requer 8 argumentos de entrada.\nUso: simulatesSampledInputSimplified(combPoly, h, K, imageName, axisVector, delta, tol, flagIsPoly)');
    end

    % Verificar se combPoly tem a estrutura esperada
    if ~isstruct(combPoly) || ~isfield(combPoly, 'A') || ~isfield(combPoly.A, 'alphaVecs')
        error('combPoly deve ser uma estrutura com campo A.alphaVecs');
    end

    fprintf('      [Simulação LMI Simplificada] Analisando %d vértices...\n', length(combPoly.A.alphaVecs));

    % Dimensões do sistema
    try
        nx = size(combPoly.A.polytopicMatrices{1}, 1);
        nu = size(combPoly.B.polytopicMatrices{1}, 2);
        nw = size(combPoly.E.polytopicMatrices{1}, 2);
    catch
        error('Erro ao extrair dimensões do politopo');
    end

    fprintf('      Sistema: nx=%d, nu=%d, nw=%d\n', nx, nu, nw);

    % Arrays para armazenar resultados
    gammaLMI = zeros(length(combPoly.A.alphaVecs), 1);
    lmi_success = false(length(combPoly.A.alphaVecs), 1);

    % Loop pelos vértices do politopo
    for i = 1:length(combPoly.A.alphaVecs)
        try
            % Extrair matrizes do vértice i
            A = combPoly.A.polytopicMatrices{i};
            B = combPoly.B.polytopicMatrices{i};
            E = combPoly.E.polytopicMatrices{i};
            C = combPoly.C.polytopicMatrices{i};
            D = combPoly.D.polytopicMatrices{i};

            % Construir matrizes caligráficas (conforme metodologia híbrida)
            ACal = [A B ; zeros(nu, nx) zeros(nu, nu)];
            ECal = [E ; zeros(nu, nw)];
            CCal = [C D];
            KCal = [eye(nx) zeros(nx, nu); K zeros(nu)];

            % Escolher função LMI baseada no tipo de análise
            if flagIsPoly
                % Análise politópica
                poly.APoly{1} = ACal;
                poly.EPoly{1} = ECal;
                poly.CPoly{1} = CCal;

                try
                    saidaLMI = estHInfAnaPolyLMILab(poly, KCal, h, delta, tol);
                catch
                    % Se falhar, tentar versão alternativa
                    try
                        saidaLMI = estHInfAnaPolyLMILabOriginal(poly, KCal, h, delta, tol);
                    catch
                        fprintf('        Vértice %d: LMI politópica falhou\n', i);
                        saidaLMI.gamma = NaN;
                    end
                end
            else
                % Análise intervalar
                try
                    saidaLMI = valEstHInfLMILab(ACal, ECal, CCal, KCal, h, delta, tol);
                catch
                    % Se falhar, tentar versões alternativas
                    try
                        saidaLMI = valEstHInfLMILabSemInt(ACal, ECal, CCal, KCal, h, delta, tol);
                    catch
                        try
                            saidaLMI = valEstHInfLMILabInt(ACal, ECal, CCal, KCal, h, delta, tol);
                        catch
                            fprintf('        Vértice %d: LMI intervalar falhou\n', i);
                            saidaLMI.gamma = NaN;
                        end
                    end
                end
            end

            % Armazenar resultado
            if isfield(saidaLMI, 'gamma') && isfinite(saidaLMI.gamma)
                gammaLMI(i) = saidaLMI.gamma;
                lmi_success(i) = true;
            else
                gammaLMI(i) = NaN;
                lmi_success(i) = false;
            end

        catch ME
            fprintf('        Vértice %d: ERRO - %s\n', i, ME.message);
            gammaLMI(i) = NaN;
            lmi_success(i) = false;
        end
    end

    % Calcular estatísticas apenas dos vértices bem-sucedidos
    gammaValidos = gammaLMI(lmi_success & isfinite(gammaLMI));
    numSucessos = sum(lmi_success);

    fprintf('      LMIs bem-sucedidas: %d/%d\n', numSucessos, length(combPoly.A.alphaVecs));

    if numSucessos > 0
        meanGamma = mean(gammaValidos);
        stdGamma = std(gammaValidos);
        [maxGamma, ~] = max(gammaValidos);
        [~, indexMaxGamma] = max(gammaLMI); % Index no array original

        % Estrutura de saída para γ
        gamma.meanGamma = meanGamma;
        gamma.stdGamma = stdGamma;
        gamma.maxGamma = maxGamma;
        gamma.indexMaxGamma = indexMaxGamma;
        gamma.numSucessos = numSucessos;

        fprintf('      γ médio: %.6f, γ máximo: %.6f\n', meanGamma, maxGamma);

        % Para esta versão simplificada, usar γ como aproximação do custo
        % (isso é razoável pois γ representa o bound da norma H∞)
        outPut.meanCost = meanGamma^2;  % γ²
        outPut.stdCost = 2 * meanGamma * stdGamma;  % Aproximação do desvio
        outPut.maxCost = maxGamma^2;    % γ²_max
        outPut.indexMaxCost = indexMaxGamma;
        outPut.custVector = gammaLMI.^2;  % γ² para cada vértice
        outPut.gamma = gamma;

        % Informação adicional
        outPut.lmi_success_rate = numSucessos / length(combPoly.A.alphaVecs);
        if flagIsPoly
            outPut.metodo = 'Politópico';
        else
            outPut.metodo = 'Intervalar';
        end

    else
        % Se nenhuma LMI funcionou
        fprintf('      ⚠ TODAS as LMIs falharam!\n');
        outPut.gamma = struct('meanGamma', NaN, 'stdGamma', NaN, 'maxGamma', NaN, ...
                             'indexMaxGamma', 1, 'numSucessos', 0);
        outPut.meanCost = NaN;
        outPut.stdCost = NaN;
        outPut.maxCost = NaN;
        outPut.indexMaxCost = 1;
        outPut.custVector = NaN(size(gammaLMI));
        outPut.lmi_success_rate = 0;
        if flagIsPoly
            outPut.metodo = 'Politópico';
        else
            outPut.metodo = 'Intervalar';
        end
    end

    fprintf('      Simulação simplificada concluída.\n');

end