function saida = main_comparacao_exemplo6()
%   MAIN_COMPARACAO_EXEMPLO6 - Comparação de 5 métodos H∞ para Exemplo 6
%
%   Resposta ao revisor: Comparação entre métodos híbridos propostos e
%   métodos de discretização equivalente da literatura.
%
%   Métodos comparados:
%   1. Híbrido Intervalar (estHInfSintLMILab) - PROPOSTO
%   2. Híbrido Politópico (estHInfSintPolyLMILab) - PROPOSTO
%   3. Equivalente Intervalar (discretização + síntese) - LITERATURA
%   4. Equivalente Politópico (discretização + síntese) - LITERATURA
%   5. Nominal (estHInfSintLMILabPrec) - BASELINE
%
%   Parâmetros fixos conforme especificação do revisor:
%   - h = 0.0500 (período de amostragem)
%   - delta = h/10 = 0.005 (parâmetro de discretização)
%   - tol = 1e-6
%
%   Autor: Rafael Alves, Maurício Souza
%   Data: 2025

    fprintf('=== COMPARAÇÃO EXEMPLO 6 - RESPOSTA AO REVISOR ===\n\n');
    tempoInicial = datetime('now');

    % Adicionar caminhos necessários
    addpath('funcoes');
    addpath('../funcoes');
    addpath('../HInf - Análise - Intervalar/funcoes');
    addpath('../HInf - Análise - Intervalar/funcoes/Diversas');

    %% Parâmetros fixos conforme especificação
    h = 0.0500;      % Período de amostragem
    delta = h/10;    % Parâmetro de discretização (0.005)
    tol = 1e-6;      % Tolerância

    fprintf('Parâmetros do sistema:\n');
    fprintf('  h (período de amostragem): %.4f\n', h);
    fprintf('  delta (discretização): %.6f\n', delta);
    fprintf('  tol (tolerância): %.0e\n\n', tol);

    %% Carregar sistema do Exemplo 6
    fprintf('Carregando sistema do Exemplo 6...\n');
    sistemaDados = genSaveDataEx6Reviewer();

    % Extrair matrizes intervalares
    A = sistemaDados.sys.A;
    B = sistemaDados.sys.B;
    E = sistemaDados.sys.E;
    C = sistemaDados.sys.C;
    D = sistemaDados.sys.D;
    sysPolyCont = sistemaDados.sys.sysPolyCont;

    % Dimensões do sistema
    nx = length(A.inf);
    nu = size(B.inf, 2);
    nw = size(E.inf, 2);
    ny = size(C.inf, 1);

    % Sistema nominal para método baseline
    A_nom = (A.inf + A.sup) / 2;
    B_nom = (B.inf + B.sup) / 2;
    E_nom = (E.inf + E.sup) / 2;
    C_nom = (C.inf + C.sup) / 2;
    D_nom = (D.inf + D.sup) / 2;

    fprintf('Sistema carregado com sucesso!\n');
    fprintf('  Dimensões: nx=%d, nu=%d, nw=%d, ny=%d\n\n', nx, nu, nw, ny);

    % Inicializar estruturas de resultado
    nomes = {'HibridoIntervalar', 'HibridoPolitopico', 'EquivalenteIntervalar', 'EquivalentePolitopico', 'Nominal'};
    numMetodos = length(nomes);

    for i = 1:numMetodos
        saida.(nomes{i}) = struct();
        saida.(nomes{i}).sucesso = false;
        saida.(nomes{i}).tempo = 0;
        saida.(nomes{i}).metodo = nomes{i};
    end

    %% MÉTODO 1: HÍBRIDO INTERVALAR
    fprintf('=== MÉTODO 1: HÍBRIDO INTERVALAR (PROPOSTO) ===\n');
    tempoRef = tic;
    try
        resultado1 = estHInfSintLMILab(A, B, E, C, D, h, delta, tol);
        tempoSintese = toc(tempoRef);

        if isfield(resultado1, 'factivel') && resultado1.factivel
            saida.HibridoIntervalar.K = resultado1.K;
            saida.HibridoIntervalar.gamma = resultado1.gamma;
            saida.HibridoIntervalar.tempo = tempoSintese;
            saida.HibridoIntervalar.sucesso = true;

            fprintf('✓ Síntese bem-sucedida!\n');
            fprintf('  Ganho K: [%.6f %.6f]\n', resultado1.K(1), resultado1.K(2));
            fprintf('  Norma γ garantida: %.6f\n', resultado1.gamma);
            fprintf('  Tempo: %.4f s\n\n', tempoSintese);
        else
            saida.HibridoIntervalar.erro = 'Síntese não factível';
            fprintf('✗ Síntese não factível\n\n');
        end

    catch ME
        saida.HibridoIntervalar.erro = ME.message;
        fprintf('✗ ERRO: %s\n\n', ME.message);
    end

    %% MÉTODO 2: HÍBRIDO POLITÓPICO
    fprintf('=== MÉTODO 2: HÍBRIDO POLITÓPICO (PROPOSTO) ===\n');
    tempoRef = tic;
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
        tempoSintese = toc(tempoRef);

        if isfield(resultado2, 'factivel') && resultado2.factivel
            saida.HibridoPolitopico.K = resultado2.K;
            saida.HibridoPolitopico.gamma = resultado2.gamma;
            saida.HibridoPolitopico.tempo = tempoSintese;
            saida.HibridoPolitopico.sucesso = true;

            fprintf('✓ Síntese bem-sucedida!\n');
            fprintf('  Ganho K: [%.6f %.6f]\n', resultado2.K(1), resultado2.K(2));
            fprintf('  Norma γ garantida: %.6f\n', resultado2.gamma);
            fprintf('  Tempo: %.4f s\n\n', tempoSintese);
        else
            saida.HibridoPolitopico.erro = 'Síntese não factível';
            fprintf('✗ Síntese não factível\n\n');
        end

    catch ME
        saida.HibridoPolitopico.erro = ME.message;
        fprintf('✗ ERRO: %s\n\n', ME.message);
    end

    %% MÉTODO 3: EQUIVALENTE INTERVALAR
    fprintf('=== MÉTODO 3: EQUIVALENTE INTERVALAR (LITERATURA) ===\n');
    tempoRef = tic;
    try
        resultado3 = synHInfKDiscIntWG(A, B, E, C, D, h, tol);
        tempoSintese = toc(tempoRef);

        if resultado3.factivel
            saida.EquivalenteIntervalar.K = resultado3.K;
            saida.EquivalenteIntervalar.gamma = resultado3.gamma;
            saida.EquivalenteIntervalar.tempo = tempoSintese;
            saida.EquivalenteIntervalar.sucesso = true;

            fprintf('✓ Síntese bem-sucedida!\n');
            fprintf('  Ganho K: [%.6f %.6f]\n', resultado3.K(1), resultado3.K(2));
            fprintf('  Norma γ garantida: %.6f\n', resultado3.gamma);
            fprintf('  Tempo: %.4f s\n\n', tempoSintese);
        else
            saida.EquivalenteIntervalar.erro = 'Síntese não factível';
            fprintf('✗ Síntese não factível\n\n');
        end

    catch ME
        saida.EquivalenteIntervalar.erro = ME.message;
        fprintf('✗ ERRO: %s\n\n', ME.message);
    end

    %% MÉTODO 4: EQUIVALENTE POLITÓPICO
    fprintf('=== MÉTODO 4: EQUIVALENTE POLITÓPICO (LITERATURA) ===\n');
    tempoRef = tic;
    try
        resultado4 = synHInfEquivPoly(poly, h, tol);
        tempoSintese = toc(tempoRef);

        if resultado4.factivel
            saida.EquivalentePolitopico.K = resultado4.K;
            saida.EquivalentePolitopico.gamma = resultado4.gamma;
            saida.EquivalentePolitopico.tempo = tempoSintese;
            saida.EquivalentePolitopico.sucesso = true;

            fprintf('✓ Síntese bem-sucedida!\n');
            fprintf('  Ganho K: [%.6f %.6f]\n', resultado4.K(1), resultado4.K(2));
            fprintf('  Norma γ garantida: %.6f\n', resultado4.gamma);
            fprintf('  Tempo: %.4f s\n\n', tempoSintese);
        else
            saida.EquivalentePolitopico.erro = 'Síntese não factível';
            fprintf('✗ Síntese não factível\n\n');
        end

    catch ME
        saida.EquivalentePolitopico.erro = ME.message;
        fprintf('✗ ERRO: %s\n\n', ME.message);
    end

    %% MÉTODO 5: NOMINAL (BASELINE)
    fprintf('=== MÉTODO 5: NOMINAL (BASELINE) ===\n');
    tempoRef = tic;
    try
        resultado5 = estHInfSintLMILabPrec(A_nom, B_nom, E_nom, C_nom, D_nom, h, delta, tol);
        tempoSintese = toc(tempoRef);

        if isfield(resultado5, 'factivel') && resultado5.factivel
            saida.Nominal.K = resultado5.K;
            saida.Nominal.gamma = resultado5.gamma;
            saida.Nominal.tempo = tempoSintese;
            saida.Nominal.sucesso = true;

            fprintf('✓ Síntese bem-sucedida!\n');
            fprintf('  Ganho K: [%.6f %.6f]\n', resultado5.K(1), resultado5.K(2));
            fprintf('  Norma γ garantida: %.6f\n', resultado5.gamma);
            fprintf('  Tempo: %.4f s\n\n', tempoSintese);
        else
            saida.Nominal.erro = 'Síntese não factível';
            fprintf('✗ Síntese não factível\n\n');
        end

    catch ME
        saida.Nominal.erro = ME.message;
        fprintf('✗ ERRO: %s\n\n', ME.message);
    end

    %% RESUMO COMPARATIVO
    fprintf('=== TABELA COMPARATIVA PARA O REVISOR ===\n');
    fprintf('Método                    | Status | K1        | K2        | γ         | Tempo (s) | Tipo\n');
    fprintf('--------------------------|--------|-----------|-----------|-----------|-----------|----------\n');

    tipos = {'PROPOSTO', 'PROPOSTO', 'LITERATURA', 'LITERATURA', 'BASELINE'};

    for i = 1:numMetodos
        nome = nomes{i};
        if saida.(nome).sucesso
            fprintf('%-25s | %-6s | %9.6f | %9.6f | %9.6f | %9.4f | %s\n', ...
                nome, 'OK', saida.(nome).K(1), saida.(nome).K(2), ...
                saida.(nome).gamma, saida.(nome).tempo, tipos{i});
        else
            fprintf('%-25s | %-6s | %9s | %9s | %9s | %9s | %s\n', ...
                nome, 'FALHOU', 'N/A', 'N/A', 'N/A', 'N/A', tipos{i});
        end
    end

    % Salvar parâmetros do sistema
    saida.sistema.h = h;
    saida.sistema.delta = delta;
    saida.sistema.tol = tol;
    saida.sistema.nx = nx;
    saida.sistema.nu = nu;
    saida.sistema.nw = nw;
    saida.sistema.ny = ny;

    % Tempo total
    tempoTotal = datetime('now') - tempoInicial;
    saida.tempoTotal = seconds(tempoTotal);

    fprintf('\nTempo total de execução: %.2f segundos\n', seconds(tempoTotal));
    fprintf('\n=== COMPARAÇÃO FINALIZADA ===\n');
    fprintf('Dados salvos para análise posterior e resposta ao revisor.\n');

end