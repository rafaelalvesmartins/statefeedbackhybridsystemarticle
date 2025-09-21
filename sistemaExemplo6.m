% Script para implementação do Sistema do Exemplo 6
% Sistema incerto com matrizes intervalares usando INTLAB
%
% Sistema:
% A = [1±0.1, -1±0.1; 0±0.5, 4±0.2]
% Bw = Bu = [1.0000; -1.0000]
% C = [I; 0] (matriz identidade 2x2 com zeros embaixo)
% Dw = [0; 0; 0; 1]
% Du = [0; 0; 1; 0]
% Período de amostragem h = 0.1

clear; clc;

% Verificar se INTLAB está disponível
if ~exist('infsup', 'file')
    error('INTLAB não está instalado ou não está no path. Instale INTLAB para usar matrizes intervalares.');
end

%% Definição das matrizes intervalares usando INTLAB
fprintf('Definindo matrizes intervalares usando INTLAB...\n');

% Matriz A com incertezas
% A = [1±0.1, -1±0.1; 0±0.5, 4±0.2]
A_interval = [infsup(0.9, 1.1), infsup(-1.1, -0.9);
              infsup(-0.5, 0.5), infsup(3.8, 4.2)];

% Matrizes B (Bw = Bu)
Bw_interval = [1.0; -1.0];
Bu_interval = [1.0; -1.0];

% Matriz C = [I; 0] (identidade 2x2 com zeros embaixo)
C_interval = [eye(2); zeros(2,2)];

% Matrizes D
Dw_interval = [0; 0; 0; 1];
Du_interval = [0; 0; 1; 0];

%% Definição das matrizes precisas (valores nominais) para comparação
fprintf('Definindo matrizes precisas (valores nominais)...\n');

% Matriz A nominal (valores centrais)
A_nominal = [1.0, -1.0;
             0.0,  4.0];

% Matrizes B nominais
Bw_nominal = [1.0; -1.0];
Bu_nominal = [1.0; -1.0];

% Matriz C nominal
C_nominal = [eye(2); zeros(2,2)];

% Matrizes D nominais
Dw_nominal = [0; 0; 0; 1];
Du_nominal = [0; 0; 1; 0];

%% Parâmetros do sistema
h = 0.1;  % Período de amostragem

%% Salvando todas as matrizes na estrutura sistemaEx6
fprintf('Salvando matrizes na estrutura sistemaEx6...\n');

% Matrizes intervalares
sistemaEx6.interval.A = A_interval;
sistemaEx6.interval.Bw = Bw_interval;
sistemaEx6.interval.Bu = Bu_interval;
sistemaEx6.interval.C = C_interval;
sistemaEx6.interval.Dw = Dw_interval;
sistemaEx6.interval.Du = Du_interval;

% Matrizes nominais
sistemaEx6.nominal.A = A_nominal;
sistemaEx6.nominal.Bw = Bw_nominal;
sistemaEx6.nominal.Bu = Bu_nominal;
sistemaEx6.nominal.C = C_nominal;
sistemaEx6.nominal.Dw = Dw_nominal;
sistemaEx6.nominal.Du = Du_nominal;

% Parâmetros
sistemaEx6.parametros.h = h;

%% Exibindo informações do sistema
fprintf('\n=== SISTEMA DO EXEMPLO 6 ===\n');
fprintf('Período de amostragem: h = %.1f\n', h);

fprintf('\n--- MATRIZES INTERVALARES ---\n');
fprintf('Matriz A (intervalar):\n');
disp(A_interval);

fprintf('Matriz Bw (intervalar):\n');
disp(Bw_interval);

fprintf('Matriz Bu (intervalar):\n');
disp(Bu_interval);

fprintf('Matriz C (intervalar):\n');
disp(C_interval);

fprintf('Matriz Dw (intervalar):\n');
disp(Dw_interval);

fprintf('Matriz Du (intervalar):\n');
disp(Du_interval);

fprintf('\n--- MATRIZES NOMINAIS ---\n');
fprintf('Matriz A (nominal):\n');
disp(A_nominal);

fprintf('Matriz Bw (nominal):\n');
disp(Bw_nominal);

fprintf('Matriz Bu (nominal):\n');
disp(Bu_nominal);

fprintf('Matriz C (nominal):\n');
disp(C_nominal);

fprintf('Matriz Dw (nominal):\n');
disp(Dw_nominal);

fprintf('Matriz Du (nominal):\n');
disp(Du_nominal);

%% Informações adicionais do sistema
fprintf('\n--- INFORMAÇÕES DO SISTEMA ---\n');
fprintf('Dimensões do sistema:\n');
fprintf('  - Estados (n): %d\n', size(A_nominal, 1));
fprintf('  - Entradas de perturbação (nw): %d\n', size(Bw_nominal, 2));
fprintf('  - Entradas de controle (nu): %d\n', size(Bu_nominal, 2));
fprintf('  - Saídas (p): %d\n', size(C_nominal, 1));

fprintf('\nSistema salvo na estrutura "sistemaEx6"\n');
fprintf('Acesse as matrizes intervalares com: sistemaEx6.interval.<matriz>\n');
fprintf('Acesse as matrizes nominais com: sistemaEx6.nominal.<matriz>\n');
fprintf('Acesse os parâmetros com: sistemaEx6.parametros.<parâmetro>\n');