function saida = genSaveDataEx6Comparacao()
% Optimal H∞ State Feedback Sampled-data Control
% Design for Markov Jump Linear Systems
% G.W. Gabriel, J.C. Geromel & K.M. Grigoriadis
%
% Exemplo 6 - Sistema com matrizes intervalares
% Adaptado para comparação de 5 controladores (resposta ao reviewer)

intvalinit('DisplayInfSup');

% Definição das incertezas para o Exemplo 6
% A = [1±0.1, -1±0.1; 0±0.5, 4±0.2]
% Elementos da matriz A com incertezas
a11 = infsup(0.9, 1.1);     % 1 ± 0.1
a12 = infsup(-1.1, -0.9);   % -1 ± 0.1
a21 = infsup(-0.5, 0.5);    % 0 ± 0.5
a22 = infsup(3.8, 4.2);     % 4 ± 0.2

% Matriz A intervalar
A = [a11 a12;
     a21 a22];

% Matrizes B (Bw = Bu)
B = [1.0; -1.0];
B = infsup(B, B);

% Matriz E (mesmo que B para este exemplo)
E = [1.0; -1.0];
E = infsup(E, E);

% Matriz C = [I; 0] (matriz 4x2)
C = [eye(2); zeros(2,2)];
C = infsup(C, C);

% Matriz D
D = [0; 0; 1; 0];  % Du
D = infsup(D, D);

% Tamanho das matrizes
nx = length(A.inf);
nu = size(B.inf, 2);
nw = size(E.inf, 2);

% Matrizes caligráficas
ACal = [A B ; zeros(nu, nx) zeros(nu, nu)];
ECal = [E ; zeros(nu, nw)];
CCal = [C D];

% Generate polytope manually - vértices das incertezas
a11_vals = [a11.inf a11.sup];
a12_vals = [a12.inf a12.sup];
a21_vals = [a21.inf a21.sup];
a22_vals = [a22.inf a22.sup];

idx = 0;
for i = 1:length(a11_vals)
    for j = 1:length(a12_vals)
        for k = 1:length(a21_vals)
            for l = 1:length(a22_vals)
                idx = idx + 1;
                
                % Matriz A para este vértice
                sysPolyCont{idx}.A = [a11_vals(i) a12_vals(j);
                                      a21_vals(k) a22_vals(l)];
                
                % Matrizes B constantes
                sysPolyCont{idx}.B2 = [1.0; -1.0];  % Bu
                sysPolyCont{idx}.B1 = [1.0; -1.0];  % Bw
                
                % Matriz C constante
                sysPolyCont{idx}.C = [eye(2); zeros(2,2)];
                
                % Matrizes D
                sysPolyCont{idx}.D2 = [0; 0; 1; 0];  % Du
                sysPolyCont{idx}.D1 = [0; 0; 0; 1];  % Dw
            end
        end
    end
end

% Outras variáveis importantes no sistema
h = 0.0500;  % Período de amostragem que funcionou (ajustado de 0.1)
qtdDiv = 10;  % Ajustado para manter delta = h/10
delta = h/qtdDiv;  % = 0.005

%% ------------------Parâmetro para as matrizes Politópicas------------
numPointsUniSpaced = 5;
numPointsBy2Points = 2;
numbPointsUniSpacedSub = 2;

% numPointsUniSpaced = 2500;
% numPointsBy2Points = 1000;
% numbPointsUniSpacedSub = 1500;

onlyVertice = 1;

%% ------------------Aux Vars------------------------------------------
tol = 1e-6;  % Ajustado para maior precisão

%% ------------------Save File-----------------------------------------
% Folder does not exist
if(exist(mfilename) ~= 7)
    mkdir(mfilename);
end

% Sys
sys.A = A;
sys.B = B;
sys.E = E;
sys.C = C;
sys.D = D;
sys.ACal = ACal;
sys.ECal = ECal;
sys.CCal = CCal;
sys.h = h;
sys.delta = delta;
sys.sysPolyCont = sysPolyCont;
saida.sys = sys;

% Aux
saida.aux.tol = tol;

% Simulacao
sim.numPointsUniSpaced = numPointsUniSpaced;
sim.numPointsBy2Points = numPointsBy2Points;
sim.numbPointsUniSpacedSub = numbPointsUniSpacedSub;
sim.onlyVertice = onlyVertice;
saida.simSys = sim;

% Adicionais para compatibilidade com script de comparação
% Sistema nominal (valores centrais) para método nominal
saida.sysNominal.A = [1.0  -1.0; 
                      0.0   4.0];
saida.sysNominal.B = [1.0; -1.0];
saida.sysNominal.E = [1.0; -1.0];
saida.sysNominal.C = [eye(2); zeros(2,2)];
saida.sysNominal.D = [0; 0; 1; 0];

% Parâmetros adicionais
saida.parametros.h = h;
saida.parametros.delta = delta;
saida.parametros.tol = tol;
saida.parametros.nx = nx;
saida.parametros.nu = nu;
saida.parametros.qtdDiv = qtdDiv;

fprintf('Sistema Exemplo 6 gerado com sucesso!\n');
fprintf('  - h = %.4f s\n', h);
fprintf('  - delta = %.4f s\n', delta);
fprintf('  - %d vértices politópicos gerados\n', length(sysPolyCont));
fprintf('  - Incertezas: A(1,1)±0.1, A(1,2)±0.1, A(2,1)±0.5, A(2,2)±0.2\n');

end