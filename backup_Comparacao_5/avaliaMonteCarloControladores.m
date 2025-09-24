function resultado = avaliaMonteCarloControladores(ganhosControladores, numRealizacoes)
%   Avaliação Monte Carlo de 5 controladores no Sistema do Exemplo 6
%
%   Entrada:
%   - ganhosControladores: struct com os 5 ganhos dos controladores
%   - numRealizacoes: número de realizações aleatórias (padrão: 1000)
%
%   Saída:
%   - resultado: struct com normaPiorCaso, normaMedia e nomeControladores

    if nargin < 2
        numRealizacoes = 1000;
    end

    fprintf('=== AVALIAÇÃO MONTE CARLO DE CONTROLADORES ===\n');
    fprintf('Sistema: Exemplo 6\n');
    fprintf('Número de realizações: %d\n\n', numRealizacoes);

    %% Parâmetros do Sistema do Exemplo 6
    % Definição das incertezas intervalares
    A_intervals = [infsup(0.9, 1.1), infsup(-1.1, -0.9);
                   infsup(-0.5, 0.5), infsup(3.8, 4.2)];

    Bw = [1.0; -1.0];
    Bu = [1.0; -1.0];
    C = [eye(2); zeros(2,2)];
    Dw = [0; 0; 0; 1];
    Du = [0; 0; 1; 0];
    h = 0.1;
    delta = h/5;
    tol = 1e-5;

    % Dimensões do sistema
    nx = 2;
    nu = 1;
    nw = 1;

    %% Definição dos controladores
    nomeControladores = {'Nominal', 'Discreto Int', 'Discreto Poly', 'Hibrido Int', 'Hibrido Poly'};
    numControladores = length(nomeControladores);

    % Extração dos ganhos
    if isstruct(ganhosControladores)
        ganhos = cell(1, numControladores);
        if isfield(ganhosControladores, 'nominal') && ganhosControladores.nominal.sucesso
            ganhos{1} = ganhosControladores.nominal.K;
        else
            ganhos{1} = [];
        end

        if isfield(ganhosControladores, 'discretoIntervalar') && ganhosControladores.discretoIntervalar.sucesso
            ganhos{2} = ganhosControladores.discretoIntervalar.K;
        else
            ganhos{2} = [];
        end

        if isfield(ganhosControladores, 'discretoPolitopico') && ganhosControladores.discretoPolitopico.sucesso
            ganhos{3} = ganhosControladores.discretoPolitopico.K;
        else
            ganhos{3} = [];
        end

        if isfield(ganhosControladores, 'hibridoIntervalar') && ganhosControladores.hibridoIntervalar.sucesso
            ganhos{4} = ganhosControladores.hibridoIntervalar.K;
        else
            ganhos{4} = [];
        end

        if isfield(ganhosControladores, 'hibridoPolitopico') && ganhosControladores.hibridoPolitopico.sucesso
            ganhos{5} = ganhosControladores.hibridoPolitopico.K;
        else
            ganhos{5} = [];
        end
    else
        error('ganhosControladores deve ser uma estrutura com os campos: nominal, discretoIntervalar, discretoPolitopico, hibridoIntervalar, hibridoPolitopico');
    end

    %% Inicialização das variáveis de resultado
    normasHinf = cell(1, numControladores);
    contadorEstavel = zeros(1, numControladores);

    for i = 1:numControladores
        normasHinf{i} = [];
    end

    %% Monte Carlo - Geração de realizações aleatórias
    fprintf('Iniciando simulação Monte Carlo...\n');

    for realizacao = 1:numRealizacoes
        if mod(realizacao, 100) == 0
            fprintf('Processando realização %d/%d...\n', realizacao, numRealizacoes);
        end

        % Gerar realização aleatória do sistema
        A_real = geraMatrizAleatoria(A_intervals);

        % Matrizes caligráficas para esta realização
        ACal_real = [A_real, Bu; zeros(nu, nx), zeros(nu, nu)];
        ECal_real = [Bw; zeros(nu, nw)];
        CCal_real = [C, Du];

        % Testar cada controlador
        for ctrl = 1:numControladores
            if isempty(ganhos{ctrl})
                continue; % Pula se o controlador não está disponível
            end

            K_ctrl = ganhos{ctrl};

            % Matriz de ganho caligráfica
            KCal_ctrl = [eye(nx), zeros(nx, nu); K_ctrl, zeros(nu)];

            try
                % Verificar estabilidade
                estavel = verificaEstabilidade(ACal_real, KCal_ctrl, h);

                if estavel
                    contadorEstavel(ctrl) = contadorEstavel(ctrl) + 1;

                    % Calcular norma H-infinity
                    normaHinf = calculaNormaHinf(A_real, Bu, Bw, C, Du, Dw, K_ctrl);

                    if ~isnan(normaHinf) && normaHinf < inf
                        normasHinf{ctrl} = [normasHinf{ctrl}, normaHinf];
                    end
                end
            catch
                % Se houver erro, considera como instável
                continue;
            end
        end
    end

    %% Cálculo dos resultados
    fprintf('\nCalculando estatísticas finais...\n');

    normaPiorCaso = zeros(1, numControladores);
    normaMedia = zeros(1, numControladores);

    for ctrl = 1:numControladores
        if ~isempty(normasHinf{ctrl})
            normaPiorCaso(ctrl) = max(normasHinf{ctrl});
            normaMedia(ctrl) = mean(normasHinf{ctrl});
        else
            normaPiorCaso(ctrl) = inf;
            normaMedia(ctrl) = inf;
        end
    end

    %% Estruturação dos resultados
    resultado.normaPiorCaso = normaPiorCaso;
    resultado.normaMedia = normaMedia;
    resultado.nomeControladores = nomeControladores;
    resultado.numRealizacoes = numRealizacoes;
    resultado.contadorEstavel = contadorEstavel;
    resultado.percentualEstavel = (contadorEstavel / numRealizacoes) * 100;
    resultado.normasCompletas = normasHinf;

    %% Exibição dos resultados
    fprintf('\n=== RESULTADOS DA AVALIAÇÃO MONTE CARLO ===\n');
    fprintf('Controlador           | %% Estável | Norma Média | Pior Caso  | Realizações Estáveis\n');
    fprintf('----------------------|----------|-------------|------------|-----------------\n');

    for ctrl = 1:numControladores
        if ~isempty(ganhos{ctrl})
            fprintf('%-20s | %7.1f%% | %10.4f | %9.4f | %d/%d\n', ...
                nomeControladores{ctrl}, resultado.percentualEstavel(ctrl), ...
                normaMedia(ctrl), normaPiorCaso(ctrl), ...
                contadorEstavel(ctrl), numRealizacoes);
        else
            fprintf('%-20s | %7s | %10s | %9s | N/A\n', ...
                nomeControladores{ctrl}, 'N/A', 'N/A', 'N/A');
        end
    end

    fprintf('\nAvaliação Monte Carlo finalizada!\n');

end

%% Funções auxiliares

function A_real = geraMatrizAleatoria(A_intervals)
    % Gera uma realização aleatória da matriz A dentro dos intervalos
    [m, n] = size(A_intervals);
    A_real = zeros(m, n);

    for i = 1:m
        for j = 1:n
            inf_val = A_intervals(i,j).inf;
            sup_val = A_intervals(i,j).sup;
            A_real(i,j) = inf_val + (sup_val - inf_val) * rand();
        end
    end
end

function estavel = verificaEstabilidade(ACal, KCal, h)
    % Verifica estabilidade do sistema híbrido amostrado
    try
        % Matriz do sistema em malha fechada
        Acl = ACal + [zeros(size(ACal,1), size(KCal,2)-size(ACal,2)), zeros(size(ACal,1), size(KCal,1)-size(ACal,1))] * KCal;

        % Matriz de transição discreta aproximada
        Ad = expm(Acl * h);

        % Verificar se todos os autovalores estão dentro do círculo unitário
        eigenvals = eig(Ad);
        raioEspectral = max(abs(eigenvals));

        estavel = (raioEspectral < 1.0);
    catch
        estavel = false;
    end
end

function normaHinf = calculaNormaHinf(A, Bu, Bw, C, Du, Dw, K)
    % Calcula a norma H-infinity do sistema em malha fechada
    try
        % Sistema em malha fechada
        Acl = A + Bu * K;
        Bcl = Bw;
        Ccl = C + Du * K;
        Dcl = Dw;

        % Verificar estabilidade
        eigenvals = eig(Acl);
        if any(real(eigenvals) >= 0)
            normaHinf = inf;
            return;
        end

        % Criar sistema LTI e calcular norma H-infinity
        sys = ss(Acl, Bcl, Ccl, Dcl);
        normaHinf = norm(sys, inf);

    catch
        normaHinf = inf;
    end
end