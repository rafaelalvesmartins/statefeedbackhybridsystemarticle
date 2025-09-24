function saida = comparacao5ControladoresEx6()
%   Comparação de 5 métodos de síntese para o Exemplo 6
%   Sistema com incertezas intervalares/politópicas
%
%   Métodos implementados:
%   1. Híbrido Intervalar - estHInfSintLMILab
%   2. Híbrido Politópico - estHInfSintPolyLMILab
%   3. Amostrado Intervalar - estHInfSintSampledInt
%   4. Amostrado Politópico - estHInfSintSampledPoly
%   5. Nominal - estHInfSintLMILabPrec
%
%   Saída:
%   - saida: estrutura com resultados dos 5 controladores

    fprintf('=== SÍNTESE DE 5 CONTROLADORES - EXEMPLO 6 ===\n\n');

    % Adicionar caminhos necessários
    addpath('../HInf - Análise - Intervalar/funcoes');
    addpath('../HInf - Análise - Intervalar/funcoes/Diversas');
    addpath('../funcoes');
    addpath('../functions');
    addpath('../functions/Optimize/');
    addpath('../functions/INTLAB/');
    addpath('../functions/INTLAB/exp');
    addpath('../functions/Varied/');
    
    % Carregar sistema
    fprintf('Carregando sistema do Exemplo 6...\n');
    sistema = genSaveDataEx6Comparacao();

    % Extrair parâmetros do sistema
    A = sistema.sys.A;
    B = sistema.sys.B;
    E = sistema.sys.E;
    C = sistema.sys.C;
    D = sistema.sys.D;
    sysPolyCont = sistema.sys.sysPolyCont;
    h = sistema.parametros.h;
    delta = sistema.parametros.delta;
    tol = sistema.parametros.tol;

    % Sistema nominal
    A_nom = sistema.sysNominal.A;
    B_nom = sistema.sysNominal.B;
    E_nom = sistema.sysNominal.E;
    C_nom = sistema.sysNominal.C;
    D_nom = sistema.sysNominal.D;

    % Dimensões
    nx = size(A.inf, 1);
    nu = size(B.inf, 2);

    % Inicializar estruturas de resultado
    nomes = {'Hibrido_Int', 'Hibrido_Poly', 'Amostrado_Int', 'Amostrado_Poly', 'Nominal'};
    numControladores = length(nomes);

    for i = 1:numControladores
        saida.(nomes{i}) = struct();
    end

    %% 1. HÍBRIDO INTERVALAR
    fprintf('\n=== MÉTODO 1: HÍBRIDO INTERVALAR ===\n');
    tempoInicio = tic;
    try
        resultado1 = estHInfSintLMILab(A, B, E, C, D, h, delta, tol);
        tempoSintese = toc(tempoInicio);

        saida.Hibrido_Int.K = resultado1.K;
        saida.Hibrido_Int.mu = resultado1.gamma;
        saida.Hibrido_Int.tempo = tempoSintese;
        if isfield(resultado1, 'numVariaveis')
            saida.Hibrido_Int.variaveis = resultado1.numVariaveis;
        else
            saida.Hibrido_Int.variaveis = NaN;
        end
        if isfield(resultado1, 'numRestricoes')
            saida.Hibrido_Int.restricoes = resultado1.numRestricoes;
        else
            saida.Hibrido_Int.restricoes = NaN;
        end
        saida.Hibrido_Int.factivel = true;

        fprintf('✓ Síntese bem-sucedida!\n');
        fprintf('  Ganho K: [%s]\n', mat2str(resultado1.K, 4));
        fprintf('  Norma garantida: %.4f\n', resultado1.gamma);
        fprintf('  Tempo: %.3f s\n', tempoSintese);

    catch ME
        saida.Hibrido_Int.factivel = false;
        saida.Hibrido_Int.erro = ME.message;
        saida.Hibrido_Int.K = [];
        saida.Hibrido_Int.mu = inf;
        saida.Hibrido_Int.tempo = toc(tempoInicio);
        saida.Hibrido_Int.variaveis = 0;
        saida.Hibrido_Int.restricoes = 0;
        fprintf('✗ ERRO: %s\n', ME.message);
    end

    %% 2. HÍBRIDO POLITÓPICO
    fprintf('\n=== MÉTODO 2: HÍBRIDO POLITÓPICO ===\n');
    tempoInicio = tic;
    try
        % Preparar sistema politópico
        n = length(sysPolyCont);
        poly.APoly = cell(n, 1);
        poly.BPoly = cell(n, 1);
        poly.EPoly = cell(n, 1);
        poly.CPoly = cell(n, 1);
        poly.DPoly = cell(n, 1);

        for c = 1:n
            poly.APoly{c} = sysPolyCont{c}.A;
            poly.BPoly{c} = sysPolyCont{c}.B2;
            poly.EPoly{c} = sysPolyCont{c}.B1;
            poly.CPoly{c} = sysPolyCont{c}.C;
            poly.DPoly{c} = sysPolyCont{c}.D2;
        end

        resultado2 = estHInfSintPolyLMILab(poly, h, delta, tol);
        tempoSintese = toc(tempoInicio);

        saida.Hibrido_Poly.K = resultado2.K;
        saida.Hibrido_Poly.mu = resultado2.gamma;
        saida.Hibrido_Poly.tempo = tempoSintese;
        if isfield(resultado2, 'numVariaveis')
            saida.Hibrido_Poly.variaveis = resultado2.numVariaveis;
        else
            saida.Hibrido_Poly.variaveis = NaN;
        end
        if isfield(resultado2, 'numRestricoes')
            saida.Hibrido_Poly.restricoes = resultado2.numRestricoes;
        else
            saida.Hibrido_Poly.restricoes = NaN;
        end
        saida.Hibrido_Poly.factivel = true;

        fprintf('✓ Síntese bem-sucedida!\n');
        fprintf('  Ganho K: [%s]\n', mat2str(resultado2.K, 4));
        fprintf('  Norma garantida: %.4f\n', resultado2.gamma);
        fprintf('  Tempo: %.3f s\n', tempoSintese);

    catch ME
        saida.Hibrido_Poly.factivel = false;
        saida.Hibrido_Poly.erro = ME.message;
        saida.Hibrido_Poly.K = [];
        saida.Hibrido_Poly.mu = inf;
        saida.Hibrido_Poly.tempo = toc(tempoInicio);
        saida.Hibrido_Poly.variaveis = 0;
        saida.Hibrido_Poly.restricoes = 0;
        fprintf('✗ ERRO: %s\n', ME.message);
    end

    %% 3. AMOSTRADO INTERVALAR
    fprintf('\n=== MÉTODO 3: AMOSTRADO INTERVALAR ===\n');
    tempoInicio = tic;
    try
        % Usar função de síntese amostrada intervalar
        resultado3 = estHInfSintSampledInt(A, B, E, C, D, h, tol);
        tempoSintese = toc(tempoInicio);

        saida.Amostrado_Int.K = resultado3.K;
        saida.Amostrado_Int.mu = resultado3.gamma;
        saida.Amostrado_Int.tempo = tempoSintese;
        saida.Amostrado_Int.variaveis = resultado3.numVariaveis;
        saida.Amostrado_Int.restricoes = resultado3.numRestricoes;
        saida.Amostrado_Int.factivel = true;

        fprintf('✓ Síntese bem-sucedida!\n');
        fprintf('  Ganho K: [%s]\n', mat2str(resultado3.K, 4));
        fprintf('  Norma garantida: %.4f\n', resultado3.gamma);
        fprintf('  Tempo: %.3f s\n', tempoSintese);

    catch ME
        saida.Amostrado_Int.factivel = false;
        saida.Amostrado_Int.erro = ME.message;
        saida.Amostrado_Int.K = [];
        saida.Amostrado_Int.mu = inf;
        saida.Amostrado_Int.tempo = toc(tempoInicio);
        saida.Amostrado_Int.variaveis = 0;
        saida.Amostrado_Int.restricoes = 0;
        fprintf('✗ ERRO: %s\n', ME.message);
    end

    %% 4. AMOSTRADO POLITÓPICO
    fprintf('\n=== MÉTODO 4: AMOSTRADO POLITÓPICO ===\n');
    tempoInicio = tic;
    try
        % Usar função de síntese amostrada politópica
        resultado4 = estHInfSintSampledPoly(poly, h, tol);
        tempoSintese = toc(tempoInicio);

        saida.Amostrado_Poly.K = resultado4.K;
        saida.Amostrado_Poly.mu = resultado4.gamma;
        saida.Amostrado_Poly.tempo = tempoSintese;
        saida.Amostrado_Poly.variaveis = resultado4.numVariaveis;
        saida.Amostrado_Poly.restricoes = resultado4.numRestricoes;
        saida.Amostrado_Poly.factivel = true;

        fprintf('✓ Síntese bem-sucedida!\n');
        fprintf('  Ganho K: [%s]\n', mat2str(resultado4.K, 4));
        fprintf('  Norma garantida: %.4f\n', resultado4.gamma);
        fprintf('  Tempo: %.3f s\n', tempoSintese);

    catch ME
        saida.Amostrado_Poly.factivel = false;
        saida.Amostrado_Poly.erro = ME.message;
        saida.Amostrado_Poly.K = [];
        saida.Amostrado_Poly.mu = inf;
        saida.Amostrado_Poly.tempo = toc(tempoInicio);
        saida.Amostrado_Poly.variaveis = 0;
        saida.Amostrado_Poly.restricoes = 0;
        fprintf('✗ ERRO: %s\n', ME.message);
    end

    %% 5. NOMINAL
    fprintf('\n=== MÉTODO 5: NOMINAL ===\n');
    tempoInicio = tic;
    try
        resultado5 = estHInfSintLMILabPrec(A_nom, B_nom, E_nom, C_nom, D_nom, h, delta, tol);
        tempoSintese = toc(tempoInicio);

        saida.Nominal.K = resultado5.K;
        saida.Nominal.mu = resultado5.gamma;
        saida.Nominal.tempo = tempoSintese;
        if isfield(resultado5, 'numVariaveis')
            saida.Nominal.variaveis = resultado5.numVariaveis;
        else
            saida.Nominal.variaveis = NaN;
        end
        if isfield(resultado5, 'numRestricoes')
            saida.Nominal.restricoes = resultado5.numRestricoes;
        else
            saida.Nominal.restricoes = NaN;
        end
        saida.Nominal.factivel = true;

        fprintf('✓ Síntese bem-sucedida!\n');
        fprintf('  Ganho K: [%s]\n', mat2str(resultado5.K, 4));
        fprintf('  Norma garantida: %.4f\n', resultado5.gamma);
        fprintf('  Tempo: %.3f s\n', tempoSintese);

    catch ME
        saida.Nominal.factivel = false;
        saida.Nominal.erro = ME.message;
        saida.Nominal.K = [];
        saida.Nominal.mu = inf;
        saida.Nominal.tempo = toc(tempoInicio);
        saida.Nominal.variaveis = 0;
        saida.Nominal.restricoes = 0;
        fprintf('✗ ERRO: %s\n', ME.message);
    end

    %% Resumo final
    fprintf('\n=== RESUMO DA SÍNTESE ===\n');
    fprintf('Método               | Status    | Ganho K                    | Norma γ | Tempo (s)\n');
    fprintf('---------------------|-----------|----------------------------|---------|----------\n');

    for i = 1:numControladores
        nome = nomes{i};
        nomeDisplay = strrep(nome, '_', ' ');
        if saida.(nome).factivel
            if length(saida.(nome).K) <= 4
                ganhoStr = mat2str(saida.(nome).K, 3);
            else
                ganhoStr = sprintf('[%.3f ...%dx1]', saida.(nome).K(1), length(saida.(nome).K));
            end
            fprintf('%-20s | %-9s | %-26s | %7.3f | %8.3f\n', ...
                nomeDisplay, 'OK', ganhoStr, ...
                saida.(nome).mu, saida.(nome).tempo);
        else
            fprintf('%-20s | %-9s | %-26s | %7s | %8s\n', ...
                nomeDisplay, 'FALHOU', 'N/A', 'N/A', 'N/A');
        end
    end

    % Salvar sistema para referência
    saida.sistema = sistema;
    saida.nomeControladores = nomes;

    fprintf('\n=== SÍNTESE FINALIZADA ===\n');

end

%% Funções auxiliares para síntese amostrada

function resultado = estHInfSintSampledInt(A, B, E, C, D, h, tol)
    % Síntese H-infinity para sistema amostrado intervalar
    % Usa função discIntSys para discretização adequada

    fprintf('  Executando síntese amostrada intervalar...\n');

    try
        % Discretizar sistema intervalar usando discIntSys
        % Construir matriz caligráfica seguindo convenção B2/B1 e D2/D1
        % B2 = entrada de controle, B1 = entrada de perturbação
        BChap = [B E];  % [B2 B1] = [B E] - controle primeiro, depois perturbação
        DChap = [D zeros(size(C.inf,1), size(E.inf,2))];  % [D2 D1] = [D 0] - feedthrough do controle primeiro

        % Discretizar usando discIntSys (funciona diretamente com intervalos)
        [Ad, BChapd, Cd, DChapd] = discIntSys(A, BChap, C, DChap, h);

        % Separar matrizes conforme convenção B2/B1
        B2d = BChapd(:, 1:size(B.inf,2));                    % Entrada de controle (B2)
        B1d = BChapd(:, (size(B.inf,2)+1):end);             % Entrada de perturbação (B1)

        % Separar matrizes feedthrough conforme convenção D2/D1
        D2d = DChapd(:, 1:size(D.inf,2));                    % Feedthrough do controle (D2)
        D1d = DChapd(:, (size(D.inf,2)+1):end);             % Feedthrough da perturbação (D1)

        % Renomear para compatibilidade com resto do código
        Bd = B2d;  % Entrada de controle
        Ed = B1d;  % Entrada de perturbação
        Dd = D2d;  % Feedthrough do controle

        % Resolver problema LMI discreto intervalar usando synHInfKDiscIntWG
        nx = size(A.inf, 1);
        nu = size(B.inf, 2);
        nw = size(Ed.inf, 2);
        ny = size(Cd.inf, 1);

        % Parâmetros para síntese
        param.tol = tol;

        % Chamar síntese H-infinity discreta intervalar
        sintese = synHInfKDiscIntWG(Ad, Bd, Ed, Cd, Dd, D1d, param);

        if sintese.feas == 1
            % Síntese bem-sucedida
            resultado.K = sintese.K;
            resultado.gamma = sintese.norm;
            resultado.numVariaveis = sintese.var;
            resultado.numRestricoes = sintese.line;
            resultado.factivel = true;
            resultado.metodo = 'Amostrado Intervalar';
            resultado.tempoSintese = sintese.cpusec;
            resultado.delta = sintese.delta;
        else
            % Síntese falhou
            resultado.K = [];
            resultado.gamma = inf;
            resultado.numVariaveis = 0;
            resultado.numRestricoes = 0;
            resultado.factivel = false;
            resultado.metodo = 'Amostrado Intervalar';
            resultado.erro = 'Síntese LMI não factível';
        end

        fprintf('    ✓ Síntese amostrada intervalar concluída\n');

    catch ME
        fprintf('    ✗ Erro na síntese amostrada intervalar: %s\n', ME.message);

        % Resultado de falha
        resultado.K = [];
        resultado.gamma = inf;
        resultado.numVariaveis = 0;
        resultado.numRestricoes = 0;
        resultado.factivel = false;
        resultado.erro = ME.message;
        resultado.metodo = 'Amostrado Intervalar';
    end
end

function resultado = estHInfSintSampledPoly(poly, h, tol)
    % Síntese H-infinity para sistema amostrado politópico
    % Usa função discIntSys para discretização adequada de cada vértice

    fprintf('  Executando síntese amostrada politópica...\n');

    try
        n = length(poly.APoly);
        AdPoly = cell(n, 1);
        BdPoly = cell(n, 1);
        EdPoly = cell(n, 1);
        CdPoly = cell(n, 1);
        DdPoly = cell(n, 1);

        % Discretizar cada vértice do politopo usando discIntSys
        for i = 1:n
            A = poly.APoly{i};
            B = poly.BPoly{i};
            E = poly.EPoly{i};
            C = poly.CPoly{i};
            D = poly.DPoly{i};

            % Construir matriz caligráfica seguindo convenção B2/B1 e D2/D1
            % B2 = entrada de controle, B1 = entrada de perturbação
            BChap = [B E];  % [B2 B1] = [B E] - controle primeiro, depois perturbação
            DChap = [D zeros(size(C,1), size(E,2))];  % [D2 D1] = [D 0] - feedthrough do controle primeiro

            % Discretizar usando discIntSys
            [Ad, BChapd, Cd, DChapd] = discIntSys(A, BChap, C, DChap, h);

            % Armazenar matrizes discretas
            AdPoly{i} = Ad;
            CdPoly{i} = Cd;

            % Separar matrizes conforme convenção B2/B1
            B2d = BChapd(:, 1:size(B,2));                    % Entrada de controle (B2)
            B1d = BChapd(:, (size(B,2)+1):end);             % Entrada de perturbação (B1)

            % Separar matrizes feedthrough conforme convenção D2/D1
            D2d = DChapd(:, 1:size(D,2));                    % Feedthrough do controle (D2)
            D1d = DChapd(:, (size(D,2)+1):end);             % Feedthrough da perturbação (D1)

            % Armazenar seguindo convenção original do código
            BdPoly{i} = B2d;  % Entrada de controle
            EdPoly{i} = B1d;  % Entrada de perturbação
            DdPoly{i} = D2d;  % Feedthrough do controle
        end

        % Sistema politópico discreto
        polyDiscreto.APoly = AdPoly;
        polyDiscreto.BPoly = BdPoly;
        polyDiscreto.EPoly = EdPoly;
        polyDiscreto.CPoly = CdPoly;
        polyDiscreto.DPoly = DdPoly;

        % Resolver problema LMI politópico discreto usando synHInfKDiscPoly
        nx = size(poly.APoly{1}, 1);
        nu = size(poly.BPoly{1}, 2);
        nw = size(EdPoly{1}, 2);
        ny = size(CdPoly{1}, 1);

        % Parâmetros para síntese
        param.tol = tol;

        % Preparar matrizes para síntese politópica (seguindo convenção B2/B1, D2/D1)
        APoly_disc = polyDiscreto.APoly;
        B2Poly_disc = polyDiscreto.BPoly;  % Entrada de controle
        B1Poly_disc = polyDiscreto.EPoly;  % Entrada de perturbação
        CPoly_disc = polyDiscreto.CPoly;
        D2Poly_disc = polyDiscreto.DPoly;  % Feedthrough do controle

        % Criar D1Poly para feedthrough da perturbação (zeros para este sistema)
        D1Poly_disc = cell(n, 1);
        for i = 1:n
            D1Poly_disc{i} = zeros(size(CdPoly{i}, 1), size(EdPoly{i}, 2));
        end

        % Chamar síntese H-infinity discreta politópica
        sintese = synHInfKDiscPoly(APoly_disc, B2Poly_disc, B1Poly_disc, CPoly_disc, D2Poly_disc, D1Poly_disc, param);

        if sintese.feas == 1
            % Síntese bem-sucedida
            resultado.K = sintese.K;
            resultado.gamma = sintese.norm;
            resultado.numVariaveis = sintese.var;
            resultado.numRestricoes = sintese.line;
            resultado.factivel = true;
            resultado.metodo = 'Amostrado Politópico';
            resultado.tempoSintese = sintese.cpusec;
            resultado.delta = sintese.delta;
        else
            % Síntese falhou
            resultado.K = [];
            resultado.gamma = inf;
            resultado.numVariaveis = 0;
            resultado.numRestricoes = 0;
            resultado.factivel = false;
            resultado.metodo = 'Amostrado Politópico';
            resultado.erro = 'Síntese LMI não factível';
        end

        fprintf('    ✓ Síntese amostrada politópica concluída\n');

    catch ME
        fprintf('    ✗ Erro na síntese amostrada politópica: %s\n', ME.message);

        % Resultado de falha
        resultado.K = [];
        resultado.gamma = inf;
        resultado.numVariaveis = 0;
        resultado.numRestricoes = 0;
        resultado.factivel = false;
        resultado.erro = ME.message;
        resultado.metodo = 'Amostrado Politópico';
    end
end