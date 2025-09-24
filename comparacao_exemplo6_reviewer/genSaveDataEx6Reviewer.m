function saida = genSaveDataEx6Reviewer()
%   GENSAVEDATAEX6REVIEWER - Geração dos dados do Exemplo 6 para o revisor
%
%   Sistema do Exemplo 6 conforme especificação do artigo:
%   Sistema nominal:
%   A_nom = [1.0  -1.0; 0.0   4.0];
%   B_nom = [1.0; -1.0];
%   E_nom = B_nom;
%   C_nom = [1 0; 0 1];
%   D_nom = [0; 0];
%
%   Incertezas: ±0.1 em A(1,1), ±0.1 em A(1,2), ±0.5 em A(2,1), ±0.2 em A(2,2)
%
%   Parâmetros fixos conforme especificação do revisor:
%   - h = 0.0500 (período de amostragem)
%   - delta = h/10 = 0.005 (parâmetro de discretização)
%   - tol = 1e-6
%
%   Autor: Rafael Alves, Maurício Souza
%   Data: 2025

    fprintf('Gerando dados do Exemplo 6 para resposta ao revisor...\n');

    % Inicializar INTLAB para aritmética intervalar
    intvalinit('DisplayInfSup');

    %% Definição do sistema nominal conforme especificação
    A_nom = [1.0  -1.0;
             0.0   4.0];
    B_nom = [1.0; -1.0];
    E_nom = B_nom;  % E_nom = B_nom conforme especificação
    C_nom = [1 0; 0 1];  % Matriz identidade 2x2
    D_nom = [0; 0];

    %% Definição das incertezas
    % ±0.1 em A(1,1), ±0.1 em A(1,2), ±0.5 em A(2,1), ±0.2 em A(2,2)

    % Elementos da matriz A com incertezas
    a11 = infsup(1.0 - 0.1, 1.0 + 0.1);    % 1.0 ± 0.1
    a12 = infsup(-1.0 - 0.1, -1.0 + 0.1);  % -1.0 ± 0.1
    a21 = infsup(0.0 - 0.5, 0.0 + 0.5);    % 0.0 ± 0.5
    a22 = infsup(4.0 - 0.2, 4.0 + 0.2);    % 4.0 ± 0.2

    % Matriz A intervalar
    A = [a11  a12;
         a21  a22];

    % Matrizes B, E, C, D (sem incertezas, transformadas em intervalares)
    B = infsup(B_nom, B_nom);
    E = infsup(E_nom, E_nom);
    C = infsup(C_nom, C_nom);
    D = infsup(D_nom, D_nom);

    %% Parâmetros fixos conforme especificação do revisor
    h = 0.0500;      % Período de amostragem
    delta = h/10;    % Parâmetro de discretização = 0.005
    tol = 1e-6;      % Tolerância

    %% Dimensões do sistema
    nx = size(A_nom, 1);  % Estados: 2
    nu = size(B_nom, 2);  % Entradas de controle: 1
    nw = size(E_nom, 2);  % Entradas de perturbação: 1
    ny = size(C_nom, 1);  % Saídas: 2

    %% Matrizes caligráficas para análise híbrida
    ACal = [A B ; zeros(nu, nx) zeros(nu, nu)];
    ECal = [E ; zeros(nu, nw)];
    CCal = [C D];

    %% Geração do politopo - todos os vértices das incertezas
    % Valores extremos para cada elemento incerto
    a11_vals = [a11.inf a11.sup];
    a12_vals = [a12.inf a12.sup];
    a21_vals = [a21.inf a21.sup];
    a22_vals = [a22.inf a22.sup];

    idx = 0;
    fprintf('Gerando vértices do politopo...\n');

    for i = 1:length(a11_vals)
        for j = 1:length(a12_vals)
            for k = 1:length(a21_vals)
                for l = 1:length(a22_vals)
                    idx = idx + 1;

                    % Matriz A para este vértice
                    sysPolyCont{idx}.A = [a11_vals(i)  a12_vals(j);
                                          a21_vals(k)  a22_vals(l)];

                    % Matrizes B2, B1 (controle e perturbação)
                    sysPolyCont{idx}.B2 = B_nom;  % Bu (entrada de controle)
                    sysPolyCont{idx}.B1 = E_nom;  % Bw (entrada de perturbação)

                    % Matriz C (constante)
                    sysPolyCont{idx}.C = C_nom;

                    % Matrizes D2, D1 (feedthrough)
                    sysPolyCont{idx}.D2 = D_nom;         % Du (feedthrough controle)
                    sysPolyCont{idx}.D1 = zeros(ny, nw); % Dw (feedthrough perturbação)
                end
            end
        end
    end

    fprintf('Politopo gerado com %d vértices.\n', idx);

    %% Estrutura de saída
    % Sistema intervalar
    sys.A = A;
    sys.B = B;
    sys.E = E;
    sys.C = C;
    sys.D = D;

    % Matrizes caligráficas
    sys.ACal = ACal;
    sys.ECal = ECal;
    sys.CCal = CCal;

    % Parâmetros
    sys.h = h;
    sys.delta = delta;

    % Sistema politópico
    sys.sysPolyCont = sysPolyCont;

    % Sistema nominal
    sysNominal.A = A_nom;
    sysNominal.B = B_nom;
    sysNominal.E = E_nom;
    sysNominal.C = C_nom;
    sysNominal.D = D_nom;

    % Estrutura de saída principal
    saida.sys = sys;
    saida.sysNominal = sysNominal;

    % Parâmetros auxiliares
    saida.aux.tol = tol;
    saida.aux.nx = nx;
    saida.aux.nu = nu;
    saida.aux.nw = nw;
    saida.aux.ny = ny;
    saida.aux.numVertices = idx;

    % Informações para o revisor
    saida.reviewer.especificacao = 'Exemplo 6 - Resposta ao Revisor';
    saida.reviewer.h = h;
    saida.reviewer.delta = delta;
    saida.reviewer.tol = tol;
    saida.reviewer.incertezas = 'A(1,1)±0.1, A(1,2)±0.1, A(2,1)±0.5, A(2,2)±0.2';

    fprintf('Sistema do Exemplo 6 gerado com sucesso!\n');
    fprintf('Parâmetros: h=%.4f, delta=%.6f, tol=%.0e\n', h, delta, tol);
    fprintf('Dimensões: nx=%d, nu=%d, nw=%d, ny=%d\n', nx, nu, nw, ny);
    fprintf('Politopo: %d vértices\n\n', idx);

end